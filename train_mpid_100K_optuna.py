import os
import time
import pickle
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader, Subset
import torchvision.transforms as transforms
import torch.nn.functional as F
import matplotlib.pyplot as plt
import seaborn as sns
import optuna

from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import (
    classification_report,
    multilabel_confusion_matrix,
    roc_auc_score,
    roc_curve,
    hamming_loss,
    accuracy_score,
    f1_score
)

# GPU bellek parçalanmasını (VRAM fragmentation) önlemek ve CUDA Out-of-Memory (OOM) 
# hatalarını engellemek için dinamik bellek  kalkanını
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

BASE_DIR = '/home/byildirim/Desktop/wire_images_100K/wire_images/'
CSV_PATH = os.path.join(BASE_DIR, 'labels.csv')
OUTPUT_DIR = '/home/byildirim/Desktop/MPID_100k_Optuna_Results'

os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_NAMES = [
    'Elektron (e-)',
    'Foton (γ)',
    'Muon (μ-)',
    'Piyon (π)',
    'Proton (p)'
]

NUM_EPOCHS = 100
EARLY_STOP_PATIENCE = 5
THRESHOLD = 0.5            # Sigmoid çıktılarını ikili sınıfa yuvarlama eşiği (Sieving Threshold)
NUM_OPTUNA_TRIALS = 15     # Bayesian Optimization (TPE Sampler) arama uzayı simülasyon sayısı
OPTIMIZER_NAME = "Adam"

# 1. SMART VERTEX CROP 
class SmartVertexCrop(object):
    """
    LArTPC okumaları (800x800) fiziksel yapıları gereği 'sparse' (seyrek) verilerdir, 
    yani matrisin %95'inden fazlası sıfır enerjili boş arka plandır. Ağın boş piksellerle vakit kaybetmesini 
    engellemek için, sıfırdan büyük enerji piksellerinin medyanını bularak etkileşim vertex'ini (merkezini) 
    yakalıyoruz ve tam bu noktadan 512x512'lik dinamik bir kırpma yapıyoruz. 
    Eğer vertex köşeye yakınsa, matris boyutunu korumak adına dışarı taşan kısımları Zero-Padding ile besliyoruz.
    """
    def __init__(self, crop_size=512):
        self.crop_size = crop_size

    def __call__(self, img_tensor):
        c, h, w = img_tensor.shape
        total_energy = torch.sum(img_tensor, dim=0)
        nonzero_indices = torch.nonzero(total_energy > 0)

        if len(nonzero_indices) == 0:
            center_y, center_x = h // 2, w // 2
        else:
            center_y = int(torch.median(nonzero_indices[:, 0]).item())
            center_x = int(torch.median(nonzero_indices[:, 1]).item())

        half_crop = self.crop_size // 2
        y1 = max(0, center_y - half_crop)
        y2 = min(h, center_y + half_crop)
        x1 = max(0, center_x - half_crop)
        x2 = min(w, center_x + half_crop)

        cropped_img = img_tensor[:, y1:y2, x1:x2]

        pad_l = max(0, -(center_x - half_crop))
        pad_t = max(0, -(center_y - half_crop))
        pad_r = max(0, (center_x + half_crop) - w)
        pad_b = max(0, (center_y + half_crop) - h)

        if any(p > 0 for p in [pad_l, pad_r, pad_t, pad_b]):
            cropped_img = F.pad(cropped_img, (pad_l, pad_r, pad_t, pad_b), mode='constant', value=0.0)
        return cropped_img

# 2. NORMALIZATION (Fiziksel Yük Ölçeklendirmesi)
class ChargeNormalize(object):
    """
    Simülasyondan gelen ham enerji (ADC) değerleri çok geniş bir skalaya yayılabilir. 
    Ağın gradyan kararlılığını korumak için 500.0 ADC üzerindeki ekstrem enerji piklerini kırpıyor (clamp) 
    ve tüm matrisi [0, 1] aralığına normalize ederek kararlı bir girdi uzayı hazırlıyoruz.
    """
    def __init__(self, clip_value=500.0):
        self.clip_value = clip_value

    def __call__(self, tensor):
        tensor = torch.clamp(tensor, 0, self.clip_value)
        tensor = tensor / self.clip_value
        return tensor

# 3. DATASET 
class LArTPCDataset(Dataset):
    """
    Geant4 simülasyonundan çıkan seyrek formatlı C++ binary çıktılardır. RAM darboğazını önlemek amacıyla, 
    her örnek çağrıldığında on-the-fly (anlık) olarak dosyadan tel düzlemleri (U, V, Y) okunur, yoğun matrislere açılır ve 
    ardından birleştirilerek 3 kanallı (RGB gibi) bir PyTorch Tensor yapısına dönüştürülür.
    """
    def __init__(self, csv_file, wire_dir, transform=None):
        self.df = pd.read_csv(csv_file, comment='#', sep=',')
        self.wire_dir = wire_dir
        self.transform = transform
        self.sparse_dtype = np.dtype([('wire', np.int16), ('tick', np.int16), ('val', np.float32)])

    def __len__(self):
        return len(self.df)

    def _read_wire_image(self, filepath):
        img = np.zeros((800, 800), dtype=np.float32)
        if not os.path.exists(filepath): return img
        with open(filepath, 'rb') as f:
            nw = np.fromfile(f, dtype=np.int32, count=1)[0]
            nt = np.fromfile(f, dtype=np.int32, count=1)[0]
            nz = np.fromfile(f, dtype=np.int32, count=1)[0]
            if nz > 0:
                data = np.fromfile(f, dtype=self.sparse_dtype, count=nz)
                valid = ((data['wire'] >= 0) & (data['wire'] < nw) & (data['tick'] >= 0) & (data['tick'] < nt))
                img[data['tick'][valid], data['wire'][valid]] = data['val'][valid]
        return img

    def _get_truth_labels(self, truth_path):
        """
        Bir nötrino etkileşiminde birden fazla parçacık eşzamanlı 
        olarak üretilebilir (Sınıflar mutually exclusive değildir). Dolayısıyla tekil Softmax yerine, 
        her parçacığın varlığını bağımsız birer 1 veya 0 olarak kodlayan 5 boyutlu bir hedef vektör çıkarıyoruz.
        """
        if not os.path.exists(truth_path): return np.zeros(5, dtype=np.float32)
        with open(truth_path, 'rb') as f:
            n = np.fromfile(f, dtype=np.int32, count=1)[0]
            np.fromfile(f, dtype=np.int32, count=2)
            if n == 0: return np.zeros(5, dtype=np.float32)
            np.fromfile(f, dtype=np.int16, count=n * 2)
            pdgs = np.fromfile(f, dtype=np.int32, count=n)
        updgs = set(np.abs(pdgs))
        return np.array([
            1. if 11 in updgs else 0.,  # e-
            1. if 22 in updgs else 0.,  # γ
            1. if 13 in updgs else 0.,  # μ-
            1. if (211 in updgs or 111 in updgs) else 0., # π
            1. if 2212 in updgs else 0. # p
        ], dtype=np.float32)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        imgs = [self._read_wire_image(os.path.join(self.wire_dir, row[f'file_view{i}_{v}'])) for i, v in enumerate(['U', 'V', 'Y'])]
        tensor = torch.from_numpy(np.stack(imgs, axis=0))
        truth_path = os.path.join(self.wire_dir, row['file_truth_Y'])
        labels = torch.from_numpy(self._get_truth_labels(truth_path))
        if self.transform:
            tensor = self.transform(tensor)
        return tensor, labels

# 4. MODEL MİMARİSİ (Dinamik ve Esnek CNN)
class MicroBooNE_MPID_Net(nn.Module):
    """
    MİMARİ STRATEJİ:
    1. Dinamik Katman Yapısı: Optuna'dan gelecek 'num_blocks' parametresine göre kendini esnekçe inşa eder.
    2. Neden GroupNorm?: LArTPC verisi yüksek oranda seyrek olduğundan standart BatchNorm istatistikleri çöker.
       Küçük batch boyutlarında (32) GroupNorm çok daha kararlı gradyan akışları sağlar.
    3. Neden Global Average Pooling (GAP)?: Ağın son evrişim katmanından sonra doğrudan Flatten yapmak 
       milyonlarca gereksiz parametre oluşturur ve overfitting'i tetikler. GAP kullanarak parametre 
       sayısını devasa ölçüde düşürdük, model boyutunu 1.70 MB'a çekerek CUDA OOM riskini sıfırladık.
    """
    def __init__(self, dropout_rate, num_blocks):
        super(MicroBooNE_MPID_Net, self).__init__()
        layers = []
        in_channels = 3
        channel_sizes = [32, 64, 128, 256, 256]

        for i in range(num_blocks):
            out_channels = channel_sizes[i]
            layers.append(nn.Conv2d(in_channels, out_channels, kernel_size=3, padding=1))
            num_groups = 8 if out_channels <= 64 else (16 if out_channels == 128 else 32)
            layers.append(nn.GroupNorm(num_groups, out_channels))
            layers.append(nn.ReLU())
            layers.append(nn.MaxPool2d(kernel_size=2, stride=2))
            in_channels = out_channels

        self.features = nn.Sequential(*layers)
        self.global_pool = nn.AdaptiveAvgPool2d((1, 1))
        self.classifier = nn.Sequential(
            nn.Linear(in_channels, 128),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(128, 5)
        )

    def forward(self, x):
        x = self.features(x)
        x = self.global_pool(x).flatten(1)
        return self.classifier(x)

# 5. OPTUNA OBJECTIVE FONKSİYONU (Eğitim Döngüsü)
def objective(trial, train_ds, val_ds, device, best_state):
    # Search Space Tanımlıyoruz
    lr = trial.suggest_float("lr", 1e-4, 5e-3, log=True)
    dropout_rate = trial.suggest_float("dropout_rate", 0.2, 0.5, step=0.1)
    batch_size = 32
    num_blocks = trial.suggest_int("num_blocks", 3, 4)

    print(f"\n---> Trial {trial.number} | LR: {lr:.5f} | Dropout: {dropout_rate} | Blocks: {num_blocks}")

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=True)

    model = MicroBooNE_MPID_Net(dropout_rate=dropout_rate, num_blocks=num_blocks).to(device)
    
    # Sınıflar bağımsız (non-exclusive) olduğu için her sınıfa izole bir binary loss 
    # uygulayan BCEWithLogitsLoss fonksiyonunu kullanıyoruz. (Arka planda Sigmoid entegrasyonu içerir).
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    
    # Validasyon kaybındaki sert zikzakları ve dalgalanmaları önlemek amacıyla, 
    # kayıp değeri 2 epoch boyunca yerinde sayarsa öğrenme oranını otomatik olarak yarıya indiren dinamik scheduler.
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2)

    history = {'train_loss': [], 'val_loss': [], 'train_acc': [], 'val_acc': []}
    patience_counter = 0
    trial_best_loss = float('inf')

    for epoch in range(NUM_EPOCHS):
        # ---------------- TRAINING PHASE ----------------
        model.train()
        train_loss = 0
        train_preds, train_truths = [], []

        for imgs, lbls in train_loader:
            imgs, lbls = imgs.to(device), lbls.to(device)
            optimizer.zero_grad()
            out = model(imgs)
            loss = criterion(out, lbls)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            probs = torch.sigmoid(out)
            train_preds.extend((probs > THRESHOLD).float().cpu().numpy())
            train_truths.extend(lbls.cpu().numpy())

        avg_train_loss = train_loss / len(train_loader)
        train_acc = (np.array(train_truths) == np.array(train_preds)).mean()

        # ---------------- VALIDATION PHASE ----------------
        model.eval()
        val_loss = 0
        val_preds, val_truths = [], []

        with torch.no_grad():
            for imgs, lbls in val_loader:
                imgs, lbls = imgs.to(device), lbls.to(device)
                out = model(imgs)
                loss = criterion(out, lbls)
                val_loss += loss.item()
                probs = torch.sigmoid(out)
                val_preds.extend((probs > THRESHOLD).float().cpu().numpy())
                val_truths.extend(lbls.cpu().numpy())

        avg_val_loss = val_loss / len(val_loader)
        val_acc = (np.array(val_truths) == np.array(val_preds)).mean()

        history['train_loss'].append(avg_train_loss)
        history['val_loss'].append(avg_val_loss)
        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)

        print(f"Epoch [{epoch+1}/{NUM_EPOCHS}] | Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f} | Train Acc: {train_acc:.4f} | Val Acc: {val_acc:.4f}")
        
        # Her epoch sonunda validasyon loss değerini scheduler mekanizmasına besliyoruz.
        scheduler.step(avg_val_loss)

        # OPTUNA PRUNING: Eğer bu deneme (trial), önceki başarılı denemelere kıyasla gidişat 
        # olarak umut vadetmiyorsa (MedianPruner kurallarına göre), eğitimi yarıda keserek kaynak israfını önler.
        trial.report(avg_val_loss, epoch)
        if trial.should_prune():
            print(f"Trial {trial.number} PRUNED.")
            torch.cuda.empty_cache()
            raise optuna.exceptions.TrialPruned()

        # EARLY STOPPING & GLOBAL REKOR KONTROLÜ
        if avg_val_loss < trial_best_loss:
            trial_best_loss = avg_val_loss
            patience_counter = 0

            if trial_best_loss < best_state['best_loss']:
                best_state['best_loss'] = trial_best_loss
                
                # Referans sızıntılarını (shallow copy) engellemek amacıyla 
                # rekor kırdığımız andaki saf geçmiş verilerini copy.deepcopy yardımıyla belleğe mühürlüyoruz.
                import copy
                best_state['history'] = copy.deepcopy(history)

                torch.save(model.state_dict(), os.path.join(OUTPUT_DIR, 'best_mpid_model.pth'))
                print(f"[NEW BEST MODEL] Val Loss: {avg_val_loss:.4f}")
        else:
            patience_counter += 1
            if patience_counter >= EARLY_STOP_PATIENCE:
                print(f"[EARLY STOP] Epoch {epoch+1}")
                break

    torch.cuda.empty_cache()
    return trial_best_loss

# 6. MAIN (Veri Akışı ve Dağıtımı)
def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    transform_pipeline = transforms.Compose([
        SmartVertexCrop(512),
        ChargeNormalize(clip_value=500.0)
    ])

    full_ds = LArTPCDataset(CSV_PATH, BASE_DIR, transform=transform_pipeline)
    df = full_ds.df
    groups = df['file_truth_Y'].values
    indices = np.arange(len(df))

    # DATA LEAKAGE KALKANI: Aynı etkileşim kökenine sahip olayların (dosyaların) 
    # şans eseri hem train hem test setine dağılarak ağın ezberleme yapmasını engellemek amacıyla, 
    # bölünmeleri "file_truth_Y" dosya gruplarını baz alarak GroupShuffleSplit ile kümelenmiş olarak gerçekleştiriyoruz.
    gss1 = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=42)
    trainval_idx, test_idx = next(gss1.split(indices, groups=groups))

    trainval_groups = groups[trainval_idx]
    gss2 = GroupShuffleSplit(n_splits=1, test_size=0.1765, random_state=42)
    train_idx_rel, val_idx_rel = next(gss2.split(trainval_idx, groups=trainval_groups))

    train_idx = trainval_idx[train_idx_rel]
    val_idx = trainval_idx[val_idx_rel]

    train_ds = Subset(full_ds, train_idx)
    val_ds = Subset(full_ds, val_idx)
    test_ds = Subset(full_ds, test_idx)

    print("\nDATASET SPLIT")
    print(f"Train      : {len(train_ds)}")
    print(f"Validation : {len(val_ds)}")
    print(f"Test       : {len(test_ds)}")

    best_state = {'best_loss': float('inf'), 'history': {}}

    print("\nOPTUNA STARTED")
    study = optuna.create_study(direction="minimize", pruner=optuna.pruners.MedianPruner(n_warmup_steps=5))
    study.optimize(lambda trial: objective(trial, train_ds, val_ds, device, best_state), n_trials=NUM_OPTUNA_TRIALS)

    best_params = study.best_params
    best_params['batch_size'] = 32
    best_loss = study.best_value

    print("\n" + "*" * 60)
    print("OPTUNA FINISHED")
    print("*" * 60)
    print(best_params)
    print(f"Best Loss: {best_loss:.4f}")

    # TEST DEĞERLENDİRME VE  ANALİZ AŞAMASI
    print("\nLoading best model and evaluating on TEST split...")
    best_model = MicroBooNE_MPID_Net(dropout_rate=best_params['dropout_rate'], num_blocks=best_params['num_blocks']).to(device)
    best_model.load_state_dict(torch.load(os.path.join(OUTPUT_DIR, 'best_mpid_model.pth')))
    best_model.eval()

    test_loader = DataLoader(test_ds, batch_size=best_params['batch_size'], shuffle=False, num_workers=4)

    all_p, all_t, all_logits = [], [], []
    
    # GPU asenkron çalıştığı için süreci ölçmeden önce 
    # tüm CUDA iş parçacıklarının senkronize olmasını bekliyoruz.
    torch.cuda.synchronize()
    inf_start = time.time()

    with torch.no_grad():
        for imgs, lbls in test_loader:
            imgs = imgs.to(device)
            out = best_model(imgs)
            probs = torch.sigmoid(out)
            all_logits.extend(probs.cpu().numpy())
            all_p.extend((probs > THRESHOLD).cpu().numpy())
            all_t.extend(lbls.numpy())

    torch.cuda.synchronize()
    inf_end = time.time()
    inf_time_per_sample = ((inf_end - inf_start) / len(test_ds)) * 1000

    all_p = np.array(all_p).astype(int)
    all_t = np.array(all_t).astype(int)
    all_logits = np.array(all_logits)

    sns.set_theme(style="whitegrid")

    #  LOSS + ACC CURVES ÇİZİMİ
    history = best_state['history']
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].plot(history['train_loss'], label='Train Loss', lw=2)
    axes[0].plot(history['val_loss'], label='Validation Loss', lw=2)
    axes[0].set_title(f'Loss vs Epoch (Best Trial: {best_loss:.4f})', fontweight='bold')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('BCEWithLogitsLoss')
    axes[0].legend()

    axes[1].plot(history['train_acc'], label='Train Accuracy', lw=2)
    axes[1].plot(history['val_acc'], label='Validation Accuracy', lw=2)
    axes[1].set_title('Accuracy vs Epoch (Element-wise)', fontweight='bold')
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('Accuracy')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, '01_Best_Training_Curves.png'), dpi=300)
    plt.close()

    #  ROC CURVES
    plt.figure(figsize=(10, 8))
    for i, class_name in enumerate(TARGET_NAMES):
        try:
            fpr, tpr, _ = roc_curve(all_t[:, i], all_logits[:, i])
            auc_score = roc_auc_score(all_t[:, i], all_logits[:, i])
            plt.plot(fpr, tpr, lw=2, label=f'{class_name} (AUC = {auc_score:.3f})')
        except ValueError:
            pass
    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves per Particle Class', fontsize=16, fontweight='bold')
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(OUTPUT_DIR, '02_Best_ROC_AUC_Curves.png'), dpi=300)
    plt.close()

    #  CONFUSION MATRICES
    mcm = multilabel_confusion_matrix(all_t, all_p)
    fig, axes = plt.subplots(1, 5, figsize=(25, 5))
    fig.suptitle('Independent Confusion Matrices', fontsize=20, fontweight='bold')
    for matrix, ax, class_name in zip(mcm, axes, TARGET_NAMES):
        sns.heatmap(matrix, annot=True, fmt='d', cmap='Blues', ax=ax, cbar=False, annot_kws={"size": 16})
        ax.set_title(class_name, fontsize=16, fontweight='bold')
        ax.set_xlabel('Predicted Label')
        ax.set_ylabel('True Label')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, '03_Best_Confusion_Matrices.png'), dpi=300)
    plt.close()

    #  FINAL REPORT GENERATION
    print("\nWriting final report...")
    h_loss = hamming_loss(all_t, all_p)
    subset_acc = accuracy_score(all_t, all_p)
    elementwise_acc = (all_t == all_p).mean()
    micro_f1 = f1_score(all_t, all_p, average='micro', zero_division=0)
    macro_f1 = f1_score(all_t, all_p, average='macro', zero_division=0)

    with open(os.path.join(OUTPUT_DIR, 'Nihai_Rapor_Optuna.txt'), 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write(" MICROBOONE MPID REPORT\n")
        f.write("=" * 60 + "\n\n")
        f.write("--- DATASET ---\n")
        f.write(f"Total Events : {len(full_ds)}\n")
        f.write(f"Train Set    : {len(train_ds)}\n")
        f.write(f"Validation   : {len(val_ds)}\n")
        f.write(f"Test Set     : {len(test_ds)}\n\n")
        f.write("--- BEST HYPERPARAMETERS ---\n")
        f.write(f"Conv Blocks  : {best_params['num_blocks']}\n")
        f.write(f"Dropout      : {best_params['dropout_rate']}\n")
        f.write(f"Learning Rate: {best_params['lr']:.6f}\n")
        f.write(f"Batch Size   : {best_params['batch_size']}\n")
        f.write(f"Optimizer    : {OPTIMIZER_NAME}\n\n")
        f.write("--- PERFORMANCE ---\n")
        f.write(f"Subset Accuracy : %{subset_acc*100:.2f}\n")
        f.write(f"Element-wise Accuracy : %{elementwise_acc*100:.2f}\n")
        f.write(f"Hamming Loss    : {h_loss:.4f}\n")
        f.write(f"Micro F1 Score  : {micro_f1:.4f}\n")
        f.write(f"Macro F1 Score  : {macro_f1:.4f}\n\n")
        f.write("--- INFERENCE ---\n")
        f.write(f"Inference Time per Sample: {inf_time_per_sample:.2f} ms\n")
        f.write(f"Model Size: {os.path.getsize(os.path.join(OUTPUT_DIR, 'best_mpid_model.pth'))/1e6:.2f} MB\n\n")
        f.write("--- CONFUSION MATRICES ---\n\n")
        f.write(f"{'SINIF':<15} | {'TP':<10} | {'TN':<10} | {'FP':<10} | {'FN':<10}\n")
        f.write("-" * 70 + "\n")
        for i, class_name in enumerate(TARGET_NAMES):
            tn, fp, fn, tp = mcm[i].ravel()
            f.write(f"{class_name:<15} | {tp:<10} | {tn:<10} | {fp:<10} | {fn:<10}\n")
        f.write("\n--- CLASSIFICATION REPORT ---\n\n")
        f.write(classification_report(all_t, all_p, target_names=TARGET_NAMES, zero_division=0))

    print("\nAll outputs saved successfully.")

if __name__ == '__main__':
    main()