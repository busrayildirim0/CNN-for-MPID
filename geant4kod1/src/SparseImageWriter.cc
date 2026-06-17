#include "SparseImageWriter.hh"
#include "DetectorConstruction.hh"
#include "G4AutoLock.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <map>

namespace {
    G4Mutex sparseWriterMutex = G4MUTEX_INITIALIZER;
    constexpr G4int kMicroBooNERawTicks = 9600;
    constexpr G4int kMicroBooNEPublicTicks = 6048;
    constexpr G4int kMicroBooNETrimFront = 2400;
    constexpr G4int kMicroBooNETrimBack = 1152;
}

SparseImageWriter::SparseImageWriter()
    : fDetResponse(nullptr),
      fOutputDir("sparse_images"),
      fNWireBins(3456),
      fNTickBins(kMicroBooNEPublicTicks),
      fTPC_X(DetectorConstruction::fTPC_X),
      fTPC_Y(DetectorConstruction::fTPC_Y),
      fTPC_Z(DetectorConstruction::fTPC_Z),
      fApplyRecombination(true),
      fApplyDiffusion(true),
      fApplyLifetime(true),
      fApplyNoise(false),
      fLabelsFileOpen(false),
      fImagesWritten(0)
{
    fDetResponse = new DetectorResponse();
    InitializeWirePlanes();
}

SparseImageWriter::~SparseImageWriter()
{
    FinalizeLabelsFile();
    delete fDetResponse;
}

void SparseImageWriter::SetImageDimensions(G4int nW, G4int nT)
{
    fNWireBins = nW;
    fNTickBins = nT;
}

void SparseImageWriter::SetApplyAllEffects(G4bool b)
{
    fApplyRecombination = b;
    fApplyDiffusion = b;
    fApplyLifetime = b;

}

void SparseImageWriter::InitializeWirePlanes()
{
    fWirePlanes.clear();

    WirePlane u;
    u.planeID = 0; u.name = "U"; u.angle = 60.0 * M_PI / 180.0;
    u.nWiresReal = 2400; u.isCollection = false;
    fWirePlanes.push_back(u);

    WirePlane v;
    v.planeID = 1; v.name = "V"; v.angle = -60.0 * M_PI / 180.0;
    v.nWiresReal = 2400; v.isCollection = false;
    fWirePlanes.push_back(v);

    WirePlane y;
    y.planeID = 2; y.name = "Y"; y.angle = 0.0;
    y.nWiresReal = 3456; y.isCollection = true;
    fWirePlanes.push_back(y);
}

G4double SparseImageWriter::CalculateWireCoord(const WirePlane& plane,
                                                G4double y, G4double z) const
{

    if (plane.planeID == 2) return z;
    else if (plane.planeID == 0)
        return y * std::sin(60.0*M_PI/180.0) + z * std::cos(60.0*M_PI/180.0);
    else
        return -y * std::sin(60.0*M_PI/180.0) + z * std::cos(60.0*M_PI/180.0);
}

void SparseImageWriter::InitializeLabelsFile(const G4String& filename)
{
    G4AutoLock lock(&sparseWriterMutex);

    mkdir(fOutputDir.c_str(), 0755);
    std::vector<G4String> subs = {"nue_cc","numu_cc","nutau_cc","nc","cosmic"};
    for (auto& s : subs) mkdir((fOutputDir + "/" + s).c_str(), 0755);

    G4String path = fOutputDir + "/" + filename;
    fLabelsFile.open(path.c_str());

    if (!fLabelsFile.is_open()) {
        G4cerr << "ERROR: Cannot open " << path << G4endl;
        return;
    }

    fLabelsFile << "# LAr TPC Sparse Wire×Time Dataset\n"
                << "# Format: Sparse COO binary (.coo files)\n"
                << "# Resolution: " << fNWireBins << "×" << fNTickBins << "\n"
                << "# Detector: Recombination + Drift + Diffusion (no noise in sparse)\n"
                << "#\n"
                << "event_id,"
                << "event_label,"
                << "interaction_type,"
                << "file_view0,"
                << "file_view1,"
                << "file_view2,"
                << "file_truth,"
                << "n_hits,"
                << "visible_energy_MeV,"
                << "n_pixels_view0,"
                << "n_pixels_view1,"
                << "n_pixels_view2,"
                << "n_truth_entries,"
                << "file_size_KB,"
                << "view0_nwires,"
                << "view1_nwires,"
                << "view2_nwires,"
                << "image_nticks"
                << "\n";
    fLabelsFile.flush();
    fLabelsFileOpen = true;
}

void SparseImageWriter::FinalizeLabelsFile()
{
    G4AutoLock lock(&sparseWriterMutex);
    if (fLabelsFileOpen && fLabelsFile.is_open()) {
        fLabelsFile.close();
        fLabelsFileOpen = false;
    }
}

void SparseImageWriter::WriteEventSparse(
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
    if (hitX.size() < 3) return;

    std::ostringstream ss;
    ss << fOutputDir << "/" << eventLabel
       << "/evt_" << std::setfill('0') << std::setw(6) << eventID;
    G4String fileBase = ss.str();

    std::ostringstream ssRel;
    ssRel << eventLabel << "/evt_" << std::setfill('0') << std::setw(6) << eventID;
    G4String relBase = ssRel.str();

    G4double visEnergy = 0;
    for (auto& e : hitEnergy) visEnergy += e;

    std::vector<SparsePixel> sparseViews[3];
    G4String viewFiles[3];
    G4int totalFileSize = 0;

    for (int p = 0; p < 3; p++) {
        FillSparseView(fWirePlanes[p], hitX, hitY, hitZ, hitEnergy, hitTime, hitDeDx, sparseViews[p]);

        G4String fname = fileBase + "_view" + std::to_string(p) + ".coo";
        viewFiles[p] = relBase + "_view" + std::to_string(p) + ".coo";

        std::ofstream out(fname.c_str(), std::ios::binary);
        if (out.is_open()) {

            int32_t nPix = sparseViews[p].size();
            int32_t nW = fWirePlanes[p].nWiresReal;
            int32_t nT = fNTickBins;
            out.write(reinterpret_cast<char*>(&nPix), sizeof(int32_t));
            out.write(reinterpret_cast<char*>(&nW), sizeof(int32_t));
            out.write(reinterpret_cast<char*>(&nT), sizeof(int32_t));

            for (auto& px : sparseViews[p])
                out.write(reinterpret_cast<char*>(&px.wire), sizeof(int16_t));

            for (auto& px : sparseViews[p])
                out.write(reinterpret_cast<char*>(&px.tick), sizeof(int16_t));

            for (auto& px : sparseViews[p])
                out.write(reinterpret_cast<char*>(&px.charge), sizeof(float));

            totalFileSize += 12 + nPix * (2 + 2 + 4);
            out.close();
        }
    }

    std::vector<TruthPixel> truthData;
    FillTruthData(fWirePlanes[2], hitX, hitY, hitZ, hitEnergy, hitTime, hitPDG, truthData);

    G4String truthFile = fileBase + "_truth.coo";
    G4String truthRelFile = relBase + "_truth.coo";

    {
        std::ofstream out(truthFile.c_str(), std::ios::binary);
        if (out.is_open()) {
            int32_t nEntries = truthData.size();
            int32_t nW = fWirePlanes[2].nWiresReal;
            int32_t nT = fNTickBins;
            out.write(reinterpret_cast<char*>(&nEntries), sizeof(int32_t));
            out.write(reinterpret_cast<char*>(&nW), sizeof(int32_t));
            out.write(reinterpret_cast<char*>(&nT), sizeof(int32_t));

            for (auto& tp : truthData)
                out.write(reinterpret_cast<char*>(&tp.wire), sizeof(int16_t));
            for (auto& tp : truthData)
                out.write(reinterpret_cast<char*>(&tp.tick), sizeof(int16_t));
            for (auto& tp : truthData)
                out.write(reinterpret_cast<char*>(&tp.pdg), sizeof(int32_t));
            for (auto& tp : truthData)
                out.write(reinterpret_cast<char*>(&tp.energy), sizeof(float));

            totalFileSize += 12 + nEntries * (2 + 2 + 4 + 4);
            out.close();
        }
    }

    {
        G4AutoLock lock(&sparseWriterMutex);
        if (fLabelsFileOpen && fLabelsFile.is_open()) {
            fLabelsFile << eventID << ","
                        << eventLabel << ","
                        << interactionType << ","
                        << viewFiles[0] << ","
                        << viewFiles[1] << ","
                        << viewFiles[2] << ","
                        << truthRelFile << ","
                        << hitX.size() << ","
                        << std::fixed << std::setprecision(2) << visEnergy << ","
                        << sparseViews[0].size() << ","
                        << sparseViews[1].size() << ","
                        << sparseViews[2].size() << ","
                        << truthData.size() << ","
                        << std::fixed << std::setprecision(1) << totalFileSize / 1024.0 << ","
                        << fWirePlanes[0].nWiresReal << ","
                        << fWirePlanes[1].nWiresReal << ","
                        << fWirePlanes[2].nWiresReal << ","
                        << fNTickBins
                        << "\n";
            fLabelsFile.flush();
        }
        fImagesWritten++;
    }
}

void SparseImageWriter::FillSparseView(
    const WirePlane& plane,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4double>& hitDeDx,
    std::vector<SparsePixel>& sparseData)
{
    sparseData.clear();

    G4double driftVelocity = fDetResponse->GetDriftVelocity();
    G4double samplingPeriod = fDetResponse->GetSamplingPeriod();
    G4double halfX = fTPC_X / 2.0;
    G4double halfZ = fTPC_Z / 2.0;
    G4double halfY = fTPC_Y / 2.0;
    G4double wireCoordMin, wireCoordMax;
    if (plane.planeID == 2) {
        wireCoordMin = -halfZ;
        wireCoordMax = halfZ;
    } else {
        G4double maxCoord = halfY * std::sin(60.0*M_PI/180.0) + halfZ * std::cos(60.0*M_PI/180.0);
        wireCoordMin = -maxCoord;
        wireCoordMax = maxCoord;
    }

    std::map<std::pair<int,int>, double> pixelMap;

    for (size_t i = 0; i < hitX.size(); i++) {
        G4double energy = hitEnergy[i];
        if (energy < 1e-6) continue;

        G4double dEdx = (i < hitDeDx.size() && hitDeDx[i] > 0) ? hitDeDx[i] : 2.1;
        G4double nElectrons;
        if (fApplyRecombination)
            nElectrons = fDetResponse->EnergyToElectrons(energy, dEdx);
        else
            nElectrons = energy / fDetResponse->GetWionization();

        if (nElectrons < 1.0) continue;

        G4double driftDist = halfX - hitX[i];
        if (driftDist < 0.0) driftDist = 0.0;
        if (driftDist > fTPC_X) driftDist = fTPC_X;
        G4double driftTime = driftDist / driftVelocity;

        G4double t0 = (i < hitTime.size() && hitTime[i] >= 0.0) ? hitTime[i] : 0.0;
        G4double arrivalTime = t0 + driftTime;
        G4int fullTick = static_cast<G4int>(arrivalTime / samplingPeriod);
        if (fullTick < kMicroBooNETrimFront || fullTick >= (kMicroBooNERawTicks - kMicroBooNETrimBack)) continue;
        G4int tickBin = fullTick - kMicroBooNETrimFront;
        if (tickBin < 0 || tickBin >= fNTickBins) continue;

        if (fApplyLifetime) {
            nElectrons *= fDetResponse->CalculateLifetimeAttenuation(driftTime);
        }

        G4double sigFrac = 1.0;
        if (!plane.isCollection) {
            sigFrac = (plane.planeID == 0) ? 0.75 : 0.80;
        }
        nElectrons *= sigFrac;

        if (nElectrons < 10.0) continue;

        G4double wireCoord = CalculateWireCoord(plane, hitY[i], hitZ[i]);
        G4double wireFrac = (wireCoord - wireCoordMin) / (wireCoordMax - wireCoordMin);
        int wireBin = static_cast<int>(wireFrac * plane.nWiresReal);
        if (wireBin < 0 || wireBin >= plane.nWiresReal) continue;

        if (fApplyDiffusion && driftTime > 0) {
            G4double sigW = fDetResponse->CalculateTransverseDiffusion(driftTime);
            G4double sigT = fDetResponse->CalculateLongitudinalDiffusion(driftTime);

            int rW = static_cast<int>(std::ceil(2.0 * sigW));
            int rT = static_cast<int>(std::ceil(2.0 * sigT));
            rW = std::max(0, std::min(rW, 5));
            rT = std::max(0, std::min(rT, 5));

            if (rW > 0 || rT > 0) {
                G4double totalWeight = 0;
                std::vector<std::tuple<int,int,double>> spread;

                for (int dw = -rW; dw <= rW; dw++) {
                    for (int dt = -rT; dt <= rT; dt++) {
                        int w = wireBin + dw;
                        int t = tickBin + dt;
                        if (w < 0 || w >= plane.nWiresReal || t < 0 || t >= fNTickBins) continue;

                        double weight = std::exp(-0.5 * (dw*dw/(sigW*sigW+0.01) + dt*dt/(sigT*sigT+0.01)));
                        spread.push_back({w, t, weight});
                        totalWeight += weight;
                    }
                }

                if (totalWeight > 0) {
                    for (auto& [w, t, weight] : spread) {
                        double charge = nElectrons * weight / totalWeight;
                        if (charge > 5.0) {
                            pixelMap[{w, t}] += charge;
                        }
                    }
                }
            } else {
                pixelMap[{wireBin, tickBin}] += nElectrons;
            }
        } else {
            pixelMap[{wireBin, tickBin}] += nElectrons;
        }
    }

    sparseData.reserve(pixelMap.size());
    for (auto& [coord, charge] : pixelMap) {
        if (charge > 10.0) {
            SparsePixel px;
            px.wire = static_cast<int16_t>(coord.first);
            px.tick = static_cast<int16_t>(coord.second);
            px.charge = static_cast<float>(charge);
            sparseData.push_back(px);
        }
    }
}

void SparseImageWriter::FillTruthData(
    const WirePlane& plane,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4int>& hitPDG,
    std::vector<TruthPixel>& truthData)
{
    truthData.clear();

    G4double samplingPeriod = fDetResponse->GetSamplingPeriod();
    G4double halfY = fTPC_Y / 2.0;

    G4double halfZ = fTPC_Z / 2.0;
    G4double wireCoordMin, wireCoordMax;
    if (plane.planeID == 2) {
        wireCoordMin = -halfZ; wireCoordMax = halfZ;
    } else {
        G4double mc = halfY*std::sin(60.0*M_PI/180.0) + halfZ*std::cos(60.0*M_PI/180.0);
        wireCoordMin = -mc; wireCoordMax = mc;
    }

    std::map<std::tuple<int,int,int>, double> truthMap;

    G4double driftVelocity = fDetResponse->GetDriftVelocity();
    G4double halfX = fTPC_X / 2.0;

    for (size_t i = 0; i < hitX.size(); i++) {
        if (hitEnergy[i] < 1e-6) continue;

        G4double wireCoord = CalculateWireCoord(plane, hitY[i], hitZ[i]);

        G4double driftDist = halfX - hitX[i];
        if (driftDist < 0.0) driftDist = 0.0;
        if (driftDist > fTPC_X) driftDist = fTPC_X;
        G4double driftTime = driftDist / driftVelocity;
        G4double t0 = (i < hitTime.size() && hitTime[i] >= 0.0) ? hitTime[i] : 0.0;
        G4double arrivalTime = t0 + driftTime;
        G4int fullTick = static_cast<G4int>(arrivalTime / samplingPeriod);
        if (fullTick < kMicroBooNETrimFront || fullTick >= (kMicroBooNERawTicks - kMicroBooNETrimBack)) continue;
        int tickBin = fullTick - kMicroBooNETrimFront;

        G4double wireFrac = (wireCoord - wireCoordMin) / (wireCoordMax - wireCoordMin);
        int wireBin = static_cast<int>(wireFrac * plane.nWiresReal);

        if (wireBin < 0 || wireBin >= plane.nWiresReal) continue;
        if (tickBin < 0 || tickBin >= fNTickBins) continue;

        G4int pdg = (i < hitPDG.size()) ? hitPDG[i] : 0;

        truthMap[{wireBin, tickBin, pdg}] += hitEnergy[i];
    }

    truthData.reserve(truthMap.size());
    for (auto& [key, energy] : truthMap) {
        auto& [w, t, pdg] = key;
        if (energy > 1e-4) {
            TruthPixel tp;
            tp.wire = static_cast<int16_t>(w);
            tp.tick = static_cast<int16_t>(t);
            tp.pdg = pdg;
            tp.energy = static_cast<float>(energy);
            truthData.push_back(tp);
        }
    }
}
