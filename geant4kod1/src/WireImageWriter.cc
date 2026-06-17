#include "WireImageWriter.hh"
#include "DetectorConstruction.hh"
#include "G4AutoLock.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <numeric>

namespace {
    G4Mutex wireImageMutex = G4MUTEX_INITIALIZER;
    constexpr G4int kMicroBooNEReadoutTicks = 9600;
}

WireImageWriter::WireImageWriter()
    : fDetResponse(nullptr),
      fOutputDir("wire_images"),
      fNWireBins(800),
      fNTickBins(800),
      fTPC_X(DetectorConstruction::fTPC_X),
      fTPC_Y(DetectorConstruction::fTPC_Y),
      fTPC_Z(DetectorConstruction::fTPC_Z),
      fApplyRecombination(true),
      fApplyDiffusion(true),
      fApplyLifetime(true),
      fApplyNoise(true),
      fApplyConvolution(true),
      fApplyDeconvolution(true),
      fLabelsFileOpen(false),
      fImagesWritten(0)
{
    fDetResponse = new DetectorResponse();
    InitializeWirePlanes();

    fSignalResponse.resize(3);
    for (G4int p = 0; p < 3; p++) {
        fSignalResponse[p] = fDetResponse->GetSignalResponse(p);
    }

    fNetKernel.resize(3);
    for (G4int p = 0; p < 3; p++) {
        fNetKernel[p] = PrecomputeNetKernel(p);
    }
}

std::vector<G4double> WireImageWriter::PrecomputeNetKernel(G4int planeID)
{
    const std::vector<G4double>& response = fSignalResponse[planeID];
    G4int kResp = response.size();

    G4int N = 256;
    while (N < 4 * kResp) N *= 2;
    G4int nFreq = N / 2 + 1;

    G4double samplingRateHz = 1e6 / fDetResponse->GetSamplingPeriod();
    G4double df = samplingRateHz / N;

    G4double NSR = 0.02;
    G4double fCutoff = 150000.0;

    std::vector<G4double> Hre(nFreq, 0.0), Him(nFreq, 0.0);
    for (G4int k = 0; k < nFreq; k++) {
        G4double sumRe = 0.0, sumIm = 0.0;
        G4double w = 2.0 * M_PI * k / N;
        for (G4int n = 0; n < kResp; n++) {
            sumRe += response[n] * std::cos(w * n);
            sumIm -= response[n] * std::sin(w * n);
        }
        Hre[k] = sumRe;
        Him[k] = sumIm;
    }

    std::vector<G4double> Fre(nFreq, 0.0);
    for (G4int k = 0; k < nFreq; k++) {
        G4double Hpower = Hre[k] * Hre[k] + Him[k] * Him[k];
        G4double wiener = Hpower / (Hpower + NSR);

        G4double f = k * df;
        G4double gaussian = std::exp(-f * f / (2.0 * fCutoff * fCutoff));

        Fre[k] = wiener * gaussian;
    }

    std::vector<G4double> kernel(N, 0.0);
    for (G4int n = 0; n < N; n++) {
        G4double val = Fre[0];
        G4double w = 2.0 * M_PI * n / N;
        for (G4int k = 1; k < nFreq - 1; k++) {
            val += 2.0 * Fre[k] * std::cos(w * k);
        }
        if (N % 2 == 0) {
            val += Fre[nFreq - 1] * ((n % 2 == 0) ? 1.0 : -1.0);
        }
        kernel[n] = val / N;
    }

    G4double maxK = 0.0;
    for (auto v : kernel) maxK = std::max(maxK, std::abs(v));
    G4double cutoff = 0.001 * maxK;

    G4int kLen = 1;
    for (G4int i = 0; i < N / 2; i++) {
        if (std::abs(kernel[i]) > cutoff) kLen = i + 1;
    }
    kLen = std::min(kLen + 5, N / 2);

    std::vector<G4double> truncated(kLen);
    for (G4int i = 0; i < kLen; i++) {
        truncated[i] = kernel[i];
    }

    return truncated;
}

WireImageWriter::~WireImageWriter()
{
    FinalizeLabelsFile();
    delete fDetResponse;
}

void WireImageWriter::SetImageDimensions(G4int nWires, G4int nTicks)
{
    fNWireBins = nWires;
    fNTickBins = nTicks;
}

void WireImageWriter::SetApplyAllEffects(G4bool apply)
{
    fApplyRecombination = apply;
    fApplyDiffusion = apply;
    fApplyLifetime = apply;
    fApplyNoise = apply;
    fApplyConvolution = apply;
    fApplyDeconvolution = apply;
}

void WireImageWriter::InitializeWirePlanes()
{
    fWirePlanes.clear();

    WirePlane uPlane;
    uPlane.planeID = 0;
    uPlane.name = "U";
    uPlane.angle = 60.0 * M_PI / 180.0;
    uPlane.nWiresReal = DetectorConstruction::fNWires_U;
    uPlane.planeZ = DetectorConstruction::fWirePlaneU_X;
    uPlane.isCollection = false;
    fWirePlanes.push_back(uPlane);

    WirePlane vPlane;
    vPlane.planeID = 1;
    vPlane.name = "V";
    vPlane.angle = -60.0 * M_PI / 180.0;
    vPlane.nWiresReal = DetectorConstruction::fNWires_V;
    vPlane.planeZ = DetectorConstruction::fWirePlaneV_X;
    vPlane.isCollection = false;
    fWirePlanes.push_back(vPlane);

    WirePlane yPlane;
    yPlane.planeID = 2;
    yPlane.name = "Y";
    yPlane.angle = 0.0;
    yPlane.nWiresReal = DetectorConstruction::fNWires_Y;
    yPlane.planeZ = DetectorConstruction::fWirePlaneY_X;
    yPlane.isCollection = true;
    fWirePlanes.push_back(yPlane);

    G4double halfY = fTPC_Y / 2.0;
    G4double halfZ = fTPC_Z / 2.0;
    fWireCoordMin.resize(3);
    fWireCoordMax.resize(3);

    fWireCoordMin[2] = -halfZ;
    fWireCoordMax[2] = halfZ;

    G4double maxCoord = halfY * std::sin(60.0 * M_PI / 180.0)
                      + halfZ * std::cos(60.0 * M_PI / 180.0);
    fWireCoordMin[0] = -maxCoord;
    fWireCoordMax[0] = maxCoord;
    fWireCoordMin[1] = -maxCoord;
    fWireCoordMax[1] = maxCoord;
}

G4double WireImageWriter::CalculateWireCoord(const WirePlane& plane,
                                              G4double y, G4double z) const
{
    if (plane.planeID == 2) {
        return z;
    } else if (plane.planeID == 0) {
        return y * std::sin(60.0 * M_PI / 180.0) + z * std::cos(60.0 * M_PI / 180.0);
    } else {
        return -y * std::sin(60.0 * M_PI / 180.0) + z * std::cos(60.0 * M_PI / 180.0);
    }
}

G4int WireImageWriter::CalculateWireNumber(const WirePlane& plane,
                                            G4double y, G4double z) const
{
    G4double wireCoord = CalculateWireCoord(plane, y, z);
    G4int pid = plane.planeID;

    G4double range = fWireCoordMax[pid] - fWireCoordMin[pid];
    if (range <= 0.0) return -1;

    G4double frac = (wireCoord - fWireCoordMin[pid]) / range;
    G4int wireNum = static_cast<G4int>(frac * plane.nWiresReal);

    if (wireNum < 0 || wireNum >= plane.nWiresReal) return -1;
    return wireNum;
}

G4double WireImageWriter::CalculateDriftDistance(const WirePlane& plane,
                                                  G4double x) const
{
    G4double halfX = fTPC_X / 2.0;
    G4double driftDist = halfX - x;
    if (driftDist < 0.0) driftDist = 0.0;
    if (driftDist > fTPC_X) driftDist = fTPC_X;
    return driftDist;
}

void WireImageWriter::InitializeLabelsFile(const G4String& filename)
{
    G4AutoLock lock(&wireImageMutex);

    mkdir(fOutputDir.c_str(), 0755);

    std::vector<G4String> subdirs = {"nue_cc", "numu_cc", "nutau_cc", "nc", "cosmic"};
    for (const auto& sub : subdirs) {
        G4String path = fOutputDir + "/" + sub;
        mkdir(path.c_str(), 0755);
    }

    G4String fullPath = fOutputDir + "/" + filename;

    bool fileExists = false;
    {
        std::ifstream test(fullPath.c_str());
        if (test.good()) {
            test.seekg(0, std::ios::end);
            fileExists = (test.tellg() > 10);
        }
    }

    if (fileExists) {
        fLabelsFile.open(fullPath.c_str(), std::ios::app);
        if (!fLabelsFile.is_open()) {
            G4cerr << "ERROR: Cannot open labels file for append: " << fullPath << G4endl;
            return;
        }
    } else {
        fLabelsFile.open(fullPath.c_str());
        if (!fLabelsFile.is_open()) {
            G4cerr << "ERROR: Cannot open labels file: " << fullPath << G4endl;
            return;
        }
        fLabelsFile << "# LAr TPC Wire x Time Image Dataset (sparse format)\n";
        fLabelsFile << "# Wire-by-wire signal formation with full detector response\n";
        fLabelsFile << "# Signal chain: Recombination -> Drift -> Diffusion -> "
                    << "NetKernel(conv+deconv) -> Downsample -> Noise -> Threshold\n";
        fLabelsFile << "# Views: U (induction, +60), V (induction, -60), Y (collection, 0)\n";
        fLabelsFile << "# Charge sparse format: header(nW,nT,nNZ) + entries(wire_i16,tick_i16,val_f32)\n";
        fLabelsFile << "# Truth COO format (Y-plane only, multi-label per pixel):\n";
        fLabelsFile << "#   header(nEntries,nW,nT) + wire[]_i16 + tick[]_i16 + pdg[]_i32 + energy[]_f32\n";
        fLabelsFile << "#   loadable with sparse_loader.load_truth(filename)\n";
        fLabelsFile << "#\n";
        fLabelsFile << "event_id,"
                    << "event_label,"
                    << "interaction_type,"
                    << "file_view0_U,"
                    << "file_view1_V,"
                    << "file_view2_Y,"
                    << "file_truth_Y,"
                    << "n_truth_entries,"
                    << "n_hits,"
                    << "visible_energy_MeV,"
                    << "n_electrons_total,"
                    << "image_nwires,"
                    << "image_nticks"
                    << "\n";
        fLabelsFile.flush();
    }

    fLabelsFileOpen = true;
}

void WireImageWriter::FinalizeLabelsFile()
{
    G4AutoLock lock(&wireImageMutex);
    if (fLabelsFileOpen && fLabelsFile.is_open()) {
        fLabelsFile.close();
        fLabelsFileOpen = false;
    }
}

void WireImageWriter::WriteEventImages(
    G4int eventID,
    const G4String& eventLabel,
    const G4String& interactionType,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4double>& hitDeDx,
    const std::vector<G4int>& hitPDG,
    const std::vector<G4int>& hitTrackID)
{
    if (hitX.empty() || hitX.size() < 3) return;

    G4int nPixels = fNWireBins * fNTickBins;

    std::ostringstream ss;
    ss << fOutputDir << "/" << eventLabel
       << "/evt_" << std::setfill('0') << std::setw(6) << eventID;
    G4String fileBase = ss.str();

    std::ostringstream ssRel;
    ssRel << eventLabel << "/evt_" << std::setfill('0') << std::setw(6) << eventID;
    G4String relBase = ssRel.str();

    G4double totalElectrons = 0.0;
    G4double visibleEnergy = 0.0;
    for (const auto& e : hitEnergy) visibleEnergy += e;

    G4String viewFiles[3];

    for (size_t p = 0; p < fWirePlanes.size(); p++) {
        const WirePlane& plane = fWirePlanes[p];

        std::vector<float> imageData(nPixels, 0.0f);
        FillWireTimeImage(plane, hitX, hitY, hitZ, hitEnergy, hitTime, hitDeDx, imageData);

        for (const auto& v : imageData) totalElectrons += v;

        G4String viewFile = fileBase + "_view" + std::to_string(p) + "_" + plane.name + ".sparse";
        WriteSparseImage(viewFile, imageData);
        viewFiles[p] = relBase + "_view" + std::to_string(p) + "_" + plane.name + ".sparse";
    }

    std::vector<TruthEntry> truthData;
    FillWireTimeTruthSparse(fWirePlanes[2], hitX, hitY, hitZ,
                             hitEnergy, hitTime, hitPDG, truthData);
    G4String truthFile    = fileBase + "_truth.coo";
    G4String truthRelFile = relBase  + "_truth.coo";
    WriteTruthCoo(truthFile, truthData);

    {
        G4AutoLock lock(&wireImageMutex);

        if (fLabelsFileOpen && fLabelsFile.is_open()) {
            fLabelsFile << eventID << ","
                        << eventLabel << ","
                        << interactionType << ","
                        << viewFiles[0] << ","
                        << viewFiles[1] << ","
                        << viewFiles[2] << ","
                        << truthRelFile << ","
                        << static_cast<int>(truthData.size()) << ","
                        << hitX.size() << ","
                        << std::fixed << std::setprecision(2) << visibleEnergy << ","
                        << std::fixed << std::setprecision(0) << totalElectrons << ","
                        << fNWireBins << ","
                        << fNTickBins
                        << "\n";
            fLabelsFile.flush();
        }
        fImagesWritten++;
    }
}

void WireImageWriter::FillWireTimeImage(
    const WirePlane& plane,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4double>& hitDeDx,
    std::vector<float>& imageData)
{
    G4int nPhysicalTicks = kMicroBooNEReadoutTicks;

    std::map<G4int, std::vector<float>> wireWaveforms;
    DepositChargeOnWires(plane, hitX, hitY, hitZ, hitEnergy, hitTime, hitDeDx,
                         wireWaveforms, nPhysicalTicks);

    if (wireWaveforms.empty()) return;

    if (fApplyConvolution) {
        const std::vector<G4double>& netK = fNetKernel[plane.planeID];
        for (auto& [wireIdx, waveform] : wireWaveforms) {
            ConvolveWireWithResponse(waveform, netK);
        }
    }

    if (!plane.isCollection) {
        G4double signalFraction = (plane.planeID == 0) ? 0.75 : 0.80;
        for (auto& [wireIdx, waveform] : wireWaveforms) {
            for (auto& v : waveform) v *= static_cast<float>(signalFraction);
        }
    }

    DownsampleToImage(wireWaveforms, nPhysicalTicks, plane.nWiresReal, imageData);

    if (fApplyNoise) {
        AddNoise(imageData, plane.planeID);
    }

    ApplyThreshold(imageData, 4.5f);
}

void WireImageWriter::DepositChargeOnWires(
    const WirePlane& plane,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4double>& hitDeDx,
    std::map<G4int, std::vector<float>>& wireWaveforms,
    G4int nPhysicalTicks)
{
    G4double driftVelocity = fDetResponse->GetDriftVelocity();
    G4double samplingPeriod = fDetResponse->GetSamplingPeriod();
    G4double halfX = fTPC_X / 2.0;

    for (size_t i = 0; i < hitX.size(); i++) {
        G4double energy = hitEnergy[i];
        if (energy < 1e-6) continue;

        G4double dEdx = (i < hitDeDx.size() && hitDeDx[i] > 0) ? hitDeDx[i] : 2.1;
        G4double nElectrons;

        if (fApplyRecombination) {
            nElectrons = fDetResponse->EnergyToElectrons(energy, dEdx);
        } else {
            nElectrons = energy / fDetResponse->GetWionization();
        }
        if (nElectrons < 1.0) continue;

        G4int wireNum = CalculateWireNumber(plane, hitY[i], hitZ[i]);
        if (wireNum < 0) continue;

        G4double driftDist = halfX - hitX[i];
        if (driftDist < 0.0) driftDist = 0.0;
        if (driftDist > fTPC_X) driftDist = fTPC_X;

        G4double driftTime = driftDist / driftVelocity;

        G4double t0 = (i < hitTime.size() && hitTime[i] >= 0.0) ? hitTime[i] : 0.0;
        G4double arrivalTime = t0 + driftTime;
        G4int centerTick = static_cast<G4int>(arrivalTime / samplingPeriod);
        if (centerTick < 0 || centerTick >= nPhysicalTicks) continue;

        if (fApplyLifetime) {
            G4double attenuation = fDetResponse->CalculateLifetimeAttenuation(driftTime);
            nElectrons *= attenuation;
        }
        if (nElectrons < 0.5) continue;

        if (fApplyDiffusion && driftTime > 0.0) {
            G4double sigmaWires = fDetResponse->CalculateTransverseDiffusion(driftTime);
            G4double sigmaTicks = fDetResponse->CalculateLongitudinalDiffusion(driftTime);

            G4int wireRadius = static_cast<G4int>(std::ceil(3.0 * sigmaWires));
            G4int tickRadius = static_cast<G4int>(std::ceil(3.0 * sigmaTicks));
            wireRadius = std::max(0, std::min(wireRadius, 5));
            tickRadius = std::max(0, std::min(tickRadius, 5));

            if (wireRadius == 0 && tickRadius == 0) {
                auto& wf = wireWaveforms[wireNum];
                if (wf.empty()) wf.resize(nPhysicalTicks, 0.0f);
                wf[centerTick] += static_cast<float>(nElectrons);
            } else {
                G4double sigW2 = sigmaWires * sigmaWires + 0.01;
                G4double sigT2 = sigmaTicks * sigmaTicks + 0.01;
                G4double weightSum = 0.0;

                struct DiffHit { G4int dw; G4int dt; G4double weight; };
                std::vector<DiffHit> diffHits;
                diffHits.reserve((2*wireRadius+1) * (2*tickRadius+1));

                for (G4int dw = -wireRadius; dw <= wireRadius; dw++) {
                    G4int w = wireNum + dw;
                    if (w < 0 || w >= plane.nWiresReal) continue;
                    G4double wWeight = std::exp(-0.5 * dw * dw / sigW2);

                    for (G4int dt = -tickRadius; dt <= tickRadius; dt++) {
                        G4int t = centerTick + dt;
                        if (t < 0 || t >= nPhysicalTicks) continue;
                        G4double tWeight = std::exp(-0.5 * dt * dt / sigT2);
                        G4double w2d = wWeight * tWeight;
                        diffHits.push_back({dw, dt, w2d});
                        weightSum += w2d;
                    }
                }

                if (weightSum > 0.0) {
                    for (const auto& dh : diffHits) {
                        G4int w = wireNum + dh.dw;
                        G4int t = centerTick + dh.dt;
                        G4double charge = nElectrons * dh.weight / weightSum;

                        auto& wf = wireWaveforms[w];
                        if (wf.empty()) wf.resize(nPhysicalTicks, 0.0f);
                        wf[t] += static_cast<float>(charge);
                    }
                }
            }
        } else {
            auto& wf = wireWaveforms[wireNum];
            if (wf.empty()) wf.resize(nPhysicalTicks, 0.0f);
            wf[centerTick] += static_cast<float>(nElectrons);
        }
    }
}

void WireImageWriter::ConvolveWireWithResponse(
    std::vector<float>& waveform,
    const std::vector<G4double>& signalResponse)
{
    G4int N = waveform.size();
    G4int K = signalResponse.size();
    if (N < 1 || K < 1) return;

    std::vector<float> convolved(N, 0.0f);

    for (G4int n = 0; n < N; n++) {
        G4double sum = 0.0;
        for (G4int k = 0; k < K; k++) {
            G4int idx = n - k;
            if (idx >= 0 && idx < N) {
                sum += waveform[idx] * signalResponse[k];
            }
        }
        convolved[n] = static_cast<float>(sum);
    }

    waveform = std::move(convolved);
}

void WireImageWriter::DownsampleToImage(
    const std::map<G4int, std::vector<float>>& wireWaveforms,
    G4int nPhysicalTicks,
    G4int nWiresReal,
    std::vector<float>& imageData)
{
    G4double wireBinWidth = static_cast<G4double>(nWiresReal) / fNWireBins;
    G4double tickBinWidth = static_cast<G4double>(nPhysicalTicks) / fNTickBins;

    for (const auto& [wireIdx, waveform] : wireWaveforms) {
        G4int wireBin = static_cast<G4int>(wireIdx / wireBinWidth);
        if (wireBin < 0 || wireBin >= fNWireBins) continue;

        for (G4int t = 0; t < nPhysicalTicks && t < static_cast<G4int>(waveform.size()); t++) {
            if (std::abs(waveform[t]) < 0.01f) continue;

            G4int tickBin = static_cast<G4int>(t / tickBinWidth);
            if (tickBin < 0 || tickBin >= fNTickBins) continue;

            G4int pixelIndex = tickBin * fNWireBins + wireBin;
            imageData[pixelIndex] += waveform[t];
        }
    }
}

void WireImageWriter::FillWireTimePDGMask(
    const WirePlane& plane,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4int>& hitPDG,
    std::vector<int32_t>& maskData)
{
    G4int nPixels = fNWireBins * fNTickBins;
    std::vector<std::map<G4int, G4double>> pixelPDGEnergy(nPixels);

    G4double driftVelocity = fDetResponse->GetDriftVelocity();
    G4double samplingPeriod = fDetResponse->GetSamplingPeriod();
    G4double halfX = fTPC_X / 2.0;
    G4double maxDriftTicks = (fTPC_X / driftVelocity) / samplingPeriod;

    for (size_t i = 0; i < hitX.size(); i++) {
        if (hitEnergy[i] < 1e-6) continue;

        G4int wireNum = CalculateWireNumber(plane, hitY[i], hitZ[i]);
        if (wireNum < 0) continue;

        G4double driftDist = halfX - hitX[i];
        if (driftDist < 0.0 || driftDist > fTPC_X) continue;

        G4double timeTick = (driftDist / driftVelocity) / samplingPeriod;

        G4int wireBin = static_cast<G4int>(
            static_cast<G4double>(wireNum) / plane.nWiresReal * fNWireBins);
        G4int tickBin = static_cast<G4int>(
            (timeTick / maxDriftTicks) * fNTickBins);

        if (wireBin < 0 || wireBin >= fNWireBins) continue;
        if (tickBin < 0 || tickBin >= fNTickBins) continue;

        G4int pixelIndex = tickBin * fNWireBins + wireBin;
        G4int pdg = (i < hitPDG.size()) ? hitPDG[i] : 0;
        pixelPDGEnergy[pixelIndex][pdg] += hitEnergy[i];
    }

    for (G4int p = 0; p < nPixels; p++) {
        if (pixelPDGEnergy[p].empty()) { maskData[p] = 0; continue; }
        G4int best = 0;
        G4double maxE = 0.0;
        for (const auto& pair : pixelPDGEnergy[p]) {
            if (pair.second > maxE) { maxE = pair.second; best = pair.first; }
        }
        maskData[p] = best;
    }
}

void WireImageWriter::FillWireTimeTruthSparse(
    const WirePlane& plane,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4int>& hitPDG,
    std::vector<TruthEntry>& truthData)
{
    truthData.clear();

    G4double driftVelocity = fDetResponse->GetDriftVelocity();
    G4double samplingPeriod = fDetResponse->GetSamplingPeriod();
    G4double halfX = fTPC_X / 2.0;
    G4int    nPhysicalTicks = kMicroBooNEReadoutTicks;

    G4double tickBinWidth = static_cast<G4double>(nPhysicalTicks) / fNTickBins;
    G4double wireBinWidth = static_cast<G4double>(plane.nWiresReal) / fNWireBins;

    std::map<std::tuple<G4int, G4int, G4int>, G4double> truthMap;

    for (size_t i = 0; i < hitX.size(); i++) {
        if (hitEnergy[i] < 1e-6) continue;

        G4int wireNum = CalculateWireNumber(plane, hitY[i], hitZ[i]);
        if (wireNum < 0) continue;

        G4double driftDist = halfX - hitX[i];
        if (driftDist < 0.0) driftDist = 0.0;
        if (driftDist > fTPC_X) driftDist = fTPC_X;

        G4double driftTime = driftDist / driftVelocity;
        G4double t0 = (i < hitTime.size() && hitTime[i] >= 0.0) ? hitTime[i] : 0.0;
        G4double arrivalTime = t0 + driftTime;
        G4int    physicalTick = static_cast<G4int>(arrivalTime / samplingPeriod);
        if (physicalTick < 0 || physicalTick >= nPhysicalTicks) continue;

        G4int wireBin = static_cast<G4int>(wireNum / wireBinWidth);
        G4int tickBin = static_cast<G4int>(physicalTick / tickBinWidth);
        if (wireBin < 0 || wireBin >= fNWireBins)  continue;
        if (tickBin < 0 || tickBin >= fNTickBins)  continue;

        G4int pdg = (i < hitPDG.size()) ? hitPDG[i] : 0;
        truthMap[std::make_tuple(wireBin, tickBin, pdg)] += hitEnergy[i];
    }

    truthData.reserve(truthMap.size());
    for (const auto& kv : truthMap) {
        const G4int wb  = std::get<0>(kv.first);
        const G4int tb  = std::get<1>(kv.first);
        const G4int pdg = std::get<2>(kv.first);
        const G4double E = kv.second;
        if (E < 1e-4) continue;
        TruthEntry te;
        te.wire = static_cast<int16_t>(wb);
        te.tick = static_cast<int16_t>(tb);
        te.pdg  = static_cast<int32_t>(pdg);
        te.energy = static_cast<float>(E);
        truthData.push_back(te);
    }
}

void WireImageWriter::WriteTruthCoo(const G4String& filename,
                                     const std::vector<TruthEntry>& truthData)
{
    std::ofstream out(filename.c_str(), std::ios::binary);
    if (!out.is_open()) {
        G4cerr << "ERROR: Cannot write truth file: " << filename << G4endl;
        return;
    }

    int32_t nEntries = static_cast<int32_t>(truthData.size());
    int32_t nW = fNWireBins;
    int32_t nT = fNTickBins;
    out.write(reinterpret_cast<const char*>(&nEntries), sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&nW),       sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&nT),       sizeof(int32_t));

    for (const auto& tp : truthData)
        out.write(reinterpret_cast<const char*>(&tp.wire),   sizeof(int16_t));
    for (const auto& tp : truthData)
        out.write(reinterpret_cast<const char*>(&tp.tick),   sizeof(int16_t));
    for (const auto& tp : truthData)
        out.write(reinterpret_cast<const char*>(&tp.pdg),    sizeof(int32_t));
    for (const auto& tp : truthData)
        out.write(reinterpret_cast<const char*>(&tp.energy), sizeof(float));

    out.close();
}

void WireImageWriter::AddNoise(std::vector<float>& imageData, G4int planeID)
{
    G4int spectrumInterval = 16;

    for (G4int w = 0; w < fNWireBins; w++) {
        std::vector<float> wireNoise(fNTickBins, 0.0f);

        if (w % spectrumInterval == 0) {
            fDetResponse->GenerateNoise(wireNoise, planeID);
        } else {
            G4double nPhysTicksPerBin = (fTPC_X / fDetResponse->GetDriftVelocity()
                                          / fDetResponse->GetSamplingPeriod()) / fNTickBins;
            G4double scaledENC = 400.0 * std::sqrt(nPhysTicksPerBin);
            G4double enc = 400.0;
            if (planeID == 0) enc *= 1.15;
            else if (planeID == 1) enc *= 1.08;

            for (G4int t = 0; t < fNTickBins; t++) {
                wireNoise[t] = static_cast<float>(G4RandGauss::shoot(0.0, enc));
            }
        }

        for (G4int t = 0; t < fNTickBins; t++) {
            imageData[t * fNWireBins + w] += wireNoise[t];
        }
    }
}

void WireImageWriter::ApplyThreshold(std::vector<float>& imageData, float nSigma)
{
    float threshold = nSigma * 400.0f;

    for (auto& v : imageData) {
        if (std::abs(v) < threshold) {
            v = 0.0f;
        }
    }
}

void WireImageWriter::WriteSparseImage(const G4String& filename,
                                        const std::vector<float>& imageData)
{
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file.is_open()) {
        G4cerr << "ERROR: Cannot write: " << filename << G4endl;
        return;
    }

    std::vector<int16_t> wireIndices;
    std::vector<int16_t> tickIndices;
    std::vector<float> values;

    for (G4int t = 0; t < fNTickBins; t++) {
        for (G4int w = 0; w < fNWireBins; w++) {
            float val = imageData[t * fNWireBins + w];
            if (val != 0.0f) {
                wireIndices.push_back(static_cast<int16_t>(w));
                tickIndices.push_back(static_cast<int16_t>(t));
                values.push_back(val);
            }
        }
    }

    int32_t nw = fNWireBins;
    int32_t nt = fNTickBins;
    int32_t nz = static_cast<int32_t>(values.size());

    file.write(reinterpret_cast<const char*>(&nw), sizeof(int32_t));
    file.write(reinterpret_cast<const char*>(&nt), sizeof(int32_t));
    file.write(reinterpret_cast<const char*>(&nz), sizeof(int32_t));

    for (size_t i = 0; i < values.size(); i++) {
        file.write(reinterpret_cast<const char*>(&wireIndices[i]), sizeof(int16_t));
        file.write(reinterpret_cast<const char*>(&tickIndices[i]), sizeof(int16_t));
        file.write(reinterpret_cast<const char*>(&values[i]), sizeof(float));
    }

    file.close();
}

void WireImageWriter::WriteBinaryFloat(const G4String& filename,
                                        const std::vector<float>& data)
{
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file.is_open()) {
        G4cerr << "ERROR: Cannot write: " << filename << G4endl;
        return;
    }
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
    file.close();
}

void WireImageWriter::WriteBinaryInt(const G4String& filename,
                                      const std::vector<int32_t>& data)
{
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file.is_open()) {
        G4cerr << "ERROR: Cannot write: " << filename << G4endl;
        return;
    }
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(int32_t));
    file.close();
}
