#include "ImageWriter.hh"
#include "G4AutoLock.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <map>

namespace {
    G4Mutex imageWriterMutex = G4MUTEX_INITIALIZER;
}

ImageWriter::ImageWriter()
    : fOutputDir("event_images"),
      fImageWidth(512),
      fImageHeight(512),
      fPixelContent(kEnergyDeposition),
      fUseLogScale(true),
      fXMin(-128.175),
      fXMax(128.175),
      fYMin(-116.5),
      fYMax(116.5),
      fZMin(-518.4),
      fZMax(518.4),
      fLabelsFileOpen(false),
      fImagesWritten(0)
{
}

ImageWriter::~ImageWriter()
{
    FinalizeLabelsFile();
}

void ImageWriter::SetImageResolution(G4int width, G4int height)
{
    fImageWidth = width;
    fImageHeight = height;
}

void ImageWriter::SetTPCBounds(G4double xMin, G4double xMax,
                                G4double yMin, G4double yMax,
                                G4double zMin, G4double zMax)
{
    fXMin = xMin; fXMax = xMax;
    fYMin = yMin; fYMax = yMax;
    fZMin = zMin; fZMax = zMax;
}

void ImageWriter::InitializeLabelsFile(const G4String& filename)
{
    G4AutoLock lock(&imageWriterMutex);

    mkdir(fOutputDir.c_str(), 0755);

    std::vector<G4String> subdirs = {
        "nue_cc", "numu_cc", "nutau_cc", "nc", "cosmic"
    };
    for (const auto& subdir : subdirs) {
        G4String fullPath = fOutputDir + "/" + subdir;
        mkdir(fullPath.c_str(), 0755);
    }

    G4String fullFilename = fOutputDir + "/" + filename;
    fLabelsFile.open(fullFilename.c_str());

    if (!fLabelsFile.is_open()) {
        G4cerr << "ERROR: Cannot open labels file: " << fullFilename << G4endl;
        return;
    }

    fLabelsFile << "# LAr TPC Event Image Dataset\n";
    fLabelsFile << "# Image format: raw binary float32 (energy deposition in MeV)\n";
    fLabelsFile << "# Mask format: raw binary int32 (PDG codes / track IDs)\n";
    fLabelsFile << "# Projections: XZ (top), YZ (side), XY (front)\n";
    fLabelsFile << "#\n";
    fLabelsFile << "event_id,"
                << "event_label,"
                << "interaction_type,"
                << "file_xz,"
                << "file_yz,"
                << "file_xy,"
                << "file_xz_pdg,"
                << "file_yz_pdg,"
                << "file_xy_pdg,"
                << "file_xz_track,"
                << "file_yz_track,"
                << "file_xy_track,"
                << "n_hits,"
                << "visible_energy_MeV,"
                << "image_width,"
                << "image_height"
                << "\n";

    fLabelsFile.flush();
    fLabelsFileOpen = true;
}

void ImageWriter::FinalizeLabelsFile()
{
    G4AutoLock lock(&imageWriterMutex);
    if (fLabelsFileOpen && fLabelsFile.is_open()) {
        fLabelsFile.close();
        fLabelsFileOpen = false;
    }
}

void ImageWriter::WriteEventImages(
    G4int eventID,
    const G4String& eventLabel,
    const G4String& interactionType,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    const std::vector<G4int>& hitPDG,
    const std::vector<G4int>& hitTrackID)
{
    if (hitX.empty()) return;

    if (hitX.size() < 5) return;

    G4int nPixels = fImageWidth * fImageHeight;

    std::ostringstream ssBase;
    ssBase << fOutputDir << "/" << eventLabel
           << "/evt_" << std::setfill('0') << std::setw(6) << eventID;
    G4String fileBase = ssBase.str();

    G4double visibleEnergy = 0.0;
    for (const auto& e : hitEnergy) visibleEnergy += e;

    std::ostringstream ssRelBase;
    ssRelBase << eventLabel << "/evt_" << std::setfill('0') << std::setw(6) << eventID;
    G4String relBase = ssRelBase.str();

    ProjectionType projections[] = {kXZ, kYZ, kXY};
    G4String projNames[] = {"xz", "yz", "xy"};

    G4String energyFiles[3], pdgFiles[3], trackFiles[3];

    for (G4int p = 0; p < 3; p++) {

        std::vector<float> pixelData(nPixels, 0.0f);
        FillProjection(projections[p], hitX, hitY, hitZ, hitEnergy, hitTime, pixelData);

        if (fUseLogScale) {
            ApplyLogScale(pixelData);
        }

        G4String energyFile = fileBase + "_" + projNames[p] + ".bin";
        WriteBinaryFloat(energyFile, pixelData, fImageWidth, fImageHeight);
        energyFiles[p] = relBase + "_" + projNames[p] + ".bin";

        std::vector<int32_t> pdgMask(nPixels, 0);
        FillPDGMask(projections[p], hitX, hitY, hitZ, hitEnergy, hitPDG, pdgMask);

        G4String pdgFile = fileBase + "_" + projNames[p] + "_pdg.bin";
        WriteBinaryInt(pdgFile, pdgMask, fImageWidth, fImageHeight);
        pdgFiles[p] = relBase + "_" + projNames[p] + "_pdg.bin";

        std::vector<int32_t> trackMask(nPixels, 0);
        FillTrackMask(projections[p], hitX, hitY, hitZ, hitEnergy, hitTrackID, trackMask);

        G4String trackFile = fileBase + "_" + projNames[p] + "_track.bin";
        WriteBinaryInt(trackFile, trackMask, fImageWidth, fImageHeight);
        trackFiles[p] = relBase + "_" + projNames[p] + "_track.bin";
    }

    {
        G4AutoLock lock(&imageWriterMutex);

        if (fLabelsFileOpen && fLabelsFile.is_open()) {
            fLabelsFile << eventID << ","
                        << eventLabel << ","
                        << interactionType << ","
                        << energyFiles[0] << ","
                        << energyFiles[1] << ","
                        << energyFiles[2] << ","
                        << pdgFiles[0] << ","
                        << pdgFiles[1] << ","
                        << pdgFiles[2] << ","
                        << trackFiles[0] << ","
                        << trackFiles[1] << ","
                        << trackFiles[2] << ","
                        << hitX.size() << ","
                        << std::fixed << std::setprecision(4) << visibleEnergy << ","
                        << fImageWidth << ","
                        << fImageHeight
                        << "\n";
            fLabelsFile.flush();
        }

        fImagesWritten++;
    }
}

void ImageWriter::FillProjection(
    ProjectionType proj,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4double>& hitTime,
    std::vector<float>& pixelData)
{
    for (size_t i = 0; i < hitX.size(); i++) {
        G4int row, col;
        if (!GetPixelCoords(proj, hitX[i], hitY[i], hitZ[i], row, col)) {
            continue;
        }

        G4int pixelIndex = row * fImageWidth + col;

        switch (fPixelContent) {
            case kEnergyDeposition:
                pixelData[pixelIndex] += static_cast<float>(hitEnergy[i]);
                break;
            case kHitCount:
                pixelData[pixelIndex] += 1.0f;
                break;
            case kTimeWeighted:
                if (i < hitTime.size() && hitTime[i] > 0) {
                    pixelData[pixelIndex] += static_cast<float>(hitEnergy[i] / hitTime[i]);
                }
                break;
        }
    }
}

void ImageWriter::FillPDGMask(
    ProjectionType proj,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4int>& hitPDG,
    std::vector<int32_t>& maskData)
{
    G4int nPixels = fImageWidth * fImageHeight;

    std::vector<std::map<G4int, G4double>> pixelPDGEnergy(nPixels);

    for (size_t i = 0; i < hitX.size(); i++) {
        G4int row, col;
        if (!GetPixelCoords(proj, hitX[i], hitY[i], hitZ[i], row, col)) {
            continue;
        }

        G4int pixelIndex = row * fImageWidth + col;
        G4int pdg = (i < hitPDG.size()) ? hitPDG[i] : 0;

        pixelPDGEnergy[pixelIndex][pdg] += hitEnergy[i];
    }

    for (G4int p = 0; p < nPixels; p++) {
        if (pixelPDGEnergy[p].empty()) {
            maskData[p] = 0;
            continue;
        }

        G4int dominantPDG = 0;
        G4double maxEnergy = 0.0;

        for (const auto& pair : pixelPDGEnergy[p]) {
            if (pair.second > maxEnergy) {
                maxEnergy = pair.second;
                dominantPDG = pair.first;
            }
        }

        maskData[p] = dominantPDG;
    }
}

void ImageWriter::FillTrackMask(
    ProjectionType proj,
    const std::vector<G4double>& hitX,
    const std::vector<G4double>& hitY,
    const std::vector<G4double>& hitZ,
    const std::vector<G4double>& hitEnergy,
    const std::vector<G4int>& hitTrackID,
    std::vector<int32_t>& maskData)
{
    G4int nPixels = fImageWidth * fImageHeight;
    std::vector<std::map<G4int, G4double>> pixelTrackEnergy(nPixels);

    for (size_t i = 0; i < hitX.size(); i++) {
        G4int row, col;
        if (!GetPixelCoords(proj, hitX[i], hitY[i], hitZ[i], row, col)) {
            continue;
        }

        G4int pixelIndex = row * fImageWidth + col;
        G4int trackID = (i < hitTrackID.size()) ? hitTrackID[i] : 0;

        pixelTrackEnergy[pixelIndex][trackID] += hitEnergy[i];
    }

    for (G4int p = 0; p < nPixels; p++) {
        if (pixelTrackEnergy[p].empty()) {
            maskData[p] = 0;
            continue;
        }

        G4int dominantTrack = 0;
        G4double maxEnergy = 0.0;

        for (const auto& pair : pixelTrackEnergy[p]) {
            if (pair.second > maxEnergy) {
                maxEnergy = pair.second;
                dominantTrack = pair.first;
            }
        }

        maskData[p] = dominantTrack;
    }
}

G4bool ImageWriter::GetPixelCoords(ProjectionType proj,
                                    G4double x, G4double y, G4double z,
                                    G4int& pixelRow, G4int& pixelCol) const
{
    G4double u = 0, v = 0;
    G4double uMin, uMax, vMin, vMax;

    switch (proj) {
        case kXZ:

            u = z;  v = x;
            uMin = fZMin; uMax = fZMax;
            vMin = fXMin; vMax = fXMax;
            break;

        case kYZ:

            u = z;  v = y;
            uMin = fZMin; uMax = fZMax;
            vMin = fYMin; vMax = fYMax;
            break;

        case kXY:

            u = x;  v = y;
            uMin = fXMin; uMax = fXMax;
            vMin = fYMin; vMax = fYMax;
            break;
    }

    if (u < uMin || u > uMax || v < vMin || v > vMax) {
        return false;
    }

    G4double uFrac = (u - uMin) / (uMax - uMin);
    G4double vFrac = (v - vMin) / (vMax - vMin);

    pixelCol = static_cast<G4int>(uFrac * (fImageWidth - 1));
    pixelRow = static_cast<G4int>((1.0 - vFrac) * (fImageHeight - 1));

    pixelCol = std::max(0, std::min(pixelCol, fImageWidth - 1));
    pixelRow = std::max(0, std::min(pixelRow, fImageHeight - 1));

    return true;
}

void ImageWriter::ApplyLogScale(std::vector<float>& data)
{

    float maxVal = 0.0f;
    for (const auto& v : data) {
        if (v > maxVal) maxVal = v;
    }

    if (maxVal <= 0.0f) return;

    float scaleFactor = maxVal * 0.01f;
    if (scaleFactor < 0.001f) scaleFactor = 0.001f;

    for (auto& v : data) {
        if (v > 0.0f) {
            v = std::log10(1.0f + v / scaleFactor);
        }
    }
}

void ImageWriter::WriteBinaryFloat(const G4String& filename,
                                    const std::vector<float>& data,
                                    G4int width, G4int height)
{
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file.is_open()) {
        G4cerr << "ERROR: Cannot open file for writing: " << filename << G4endl;
        return;
    }

    file.write(reinterpret_cast<const char*>(data.data()),
               data.size() * sizeof(float));
    file.close();
}

void ImageWriter::WriteBinaryInt(const G4String& filename,
                                  const std::vector<int32_t>& data,
                                  G4int width, G4int height)
{
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file.is_open()) {
        G4cerr << "ERROR: Cannot open file for writing: " << filename << G4endl;
        return;
    }

    file.write(reinterpret_cast<const char*>(data.data()),
               data.size() * sizeof(int32_t));
    file.close();
}

G4String ImageWriter::GetProjectionName(ProjectionType proj) const
{
    switch (proj) {
        case kXZ: return "xz";
        case kYZ: return "yz";
        case kXY: return "xy";
        default:  return "unknown";
    }
}
