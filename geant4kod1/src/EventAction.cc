#include "EventAction.hh"
#include "RunAction.hh"
#include "LArHit.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventInformation.hh"
#include "ImageWriter.hh"
#include "WireImageWriter.hh"
#include "SparseImageWriter.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "DetectorConstruction.hh"
#include <algorithm>
#include <cmath>
#include <random>
#include <map>
#include <set>

G4int EventAction::fGlobalEventCounter = 0;

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fImageWriter(nullptr),
      fWireImageWriter(nullptr),
      fSparseImageWriter(nullptr),
      fEnergyDeposit(0.),
      fTrackLength(0.),
      fEventID(0),
      fGlobalEventID(0),
      fEventType(0),
      fNeutrinoQ2(-1.),
      fNeutrinoW(-1.),
      fNeutrinoEnu(-1.),
      fInteractionType("Unknown")
{
    ResetEventData();

    fImageWriter = nullptr;

    fWireImageWriter = new WireImageWriter();
    fWireImageWriter->SetImageDimensions(800, 800);
    fWireImageWriter->SetOutputDirectory("wire_images");
    fWireImageWriter->SetApplyAllEffects(true);
    fWireImageWriter->InitializeLabelsFile("labels.csv");

    fSparseImageWriter = nullptr;
}

EventAction::~EventAction()
{
    if (fImageWriter) {
        fImageWriter->FinalizeLabelsFile();
        delete fImageWriter;
    }
    if (fWireImageWriter) {
        fWireImageWriter->FinalizeLabelsFile();
        G4cout << "Wire images written: " << fWireImageWriter->GetImagesWritten() << G4endl;
        delete fWireImageWriter;
    }
    if (fSparseImageWriter) {
        fSparseImageWriter->FinalizeLabelsFile();
        G4cout << "Sparse images written: " << fSparseImageWriter->GetImagesWritten() << G4endl;
        delete fSparseImageWriter;
    }
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    fEventID = event->GetEventID();
    fGlobalEventID = fGlobalEventCounter;
    fGlobalEventCounter++;
    ResetEventData();

    const PrimaryGeneratorAction* primaryGenerator =
        static_cast<const PrimaryGeneratorAction*>(
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    switch(primaryGenerator->GetGeneratorMode()) {
        case PrimaryGeneratorAction::kCosmicRayMode:
            fEventType = 0; break;
        case PrimaryGeneratorAction::kNeutrinoMode:
            fEventType = 1; break;
        case PrimaryGeneratorAction::kTestMode:
            fEventType = 2; break;
        default:
            fEventType = 0; break;
    }
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    EventInformation* eventInfo =
        static_cast<EventInformation*>(event->GetUserInformation());

    if (eventInfo) {
        fNeutrinoQ2  = eventInfo->GetQ2();
        fNeutrinoW   = eventInfo->GetW();
        fNeutrinoEnu = eventInfo->GetEnu();
        fInteractionType = eventInfo->GetInteractionType();

        G4ThreeVector vertex   = eventInfo->GetPrimaryVertex();
        G4ThreeVector momentum = eventInfo->GetPrimaryMomentum();

        fVertexX = vertex.x() / cm;
        fVertexY = vertex.y() / cm;
        fVertexZ = vertex.z() / cm;

        fPrimaryPx = momentum.x();
        fPrimaryPy = momentum.y();
        fPrimaryPz = momentum.z();
        fPrimaryP  = momentum.mag();

        fPrimaryE     = eventInfo->GetPrimaryEnergy();
        fPrimaryTheta = eventInfo->GetPrimaryTheta();
        fPrimaryPhi   = eventInfo->GetPrimaryPhi();
        fPrimaryPDG   = eventInfo->GetPrimaryPDG();
    } else {
        fNeutrinoQ2 = fNeutrinoW = fNeutrinoEnu = -1.0;
        fInteractionType = "Unknown";
        fVertexX = fVertexY = fVertexZ = 0.0;
        fPrimaryPx = fPrimaryPy = fPrimaryPz = 0.0;
        fPrimaryP = fPrimaryE = 0.0;
        fPrimaryTheta = fPrimaryPhi = 0.0;
        fPrimaryPDG = 0;
    }

    CollectHitData(event);

    CalculateGeometricFeatures();

    SaveEventData();

    SaveEventImages();

    SaveWireImages();

    SaveSparseImages();
}

void EventAction::CollectHitData(const G4Event* event)
{
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) return;

    G4int larHCID = G4SDManager::GetSDMpointer()->GetCollectionID("LArHitsCollection");
    if (larHCID < 0) return;

    LArHitsCollection* larHC = static_cast<LArHitsCollection*>(hce->GetHC(larHCID));
    if (!larHC) return;

    fHitX.clear(); fHitY.clear(); fHitZ.clear();
    fHitEnergy.clear(); fHitTime.clear();
    fWirePlane.clear(); fWireNumber.clear();
    fParticlePDG.clear();
    fTrackID.clear();
    fHitDeDx.clear();

    G4int nHits = larHC->entries();
    fNHits = nHits;

    fHitX.reserve(nHits);
    fHitY.reserve(nHits);
    fHitZ.reserve(nHits);
    fHitEnergy.reserve(nHits);
    fHitTime.reserve(nHits);
    fWirePlane.reserve(nHits);
    fWireNumber.reserve(nHits);
    fParticlePDG.reserve(nHits);
    fTrackID.reserve(nHits);
    fHitDeDx.reserve(nHits);

    for (G4int i = 0; i < nHits; i++) {
        LArHit* hit = (*larHC)[i];

        G4ThreeVector pos = hit->GetPosition();
        fHitX.push_back(pos.x()/cm);
        fHitY.push_back(pos.y()/cm);
        fHitZ.push_back(pos.z()/cm);
        fHitEnergy.push_back(hit->GetEnergyDeposit()/MeV);

        G4double hitTime = hit->GetTime()/microsecond;
        if (hitTime < 0.0 || hitTime >= 1e5) hitTime = 0.0;
        fHitTime.push_back(hitTime);

        fWirePlane.push_back(hit->GetWirePlaneID());
        fWireNumber.push_back(hit->GetWireNumber());
        fParticlePDG.push_back(hit->GetParticleType());
        fTrackID.push_back(hit->GetTrackID());

        G4double dEdx_real = hit->GetDeDx();
        if (dEdx_real <= 0.0) dEdx_real = 2.1;
        fHitDeDx.push_back(dEdx_real);

        fEnergyDeposit += hit->GetEnergyDeposit();
    }
}

void EventAction::SaveEventImages()
{
    if (!fImageWriter) return;
    if (fHitX.empty()) return;

    G4String eventLabel = DetermineEventLabel();

    fImageWriter->WriteEventImages(
        fGlobalEventID,
        eventLabel,
        fInteractionType,
        fHitX, fHitY, fHitZ,
        fHitEnergy,
        fHitTime,
        fParticlePDG,
        fTrackID
    );
}

void EventAction::SaveWireImages()
{
    if (!fWireImageWriter) return;
    if (fHitX.empty()) return;

    G4String eventLabel = DetermineEventLabel();

    fWireImageWriter->WriteEventImages(
        fGlobalEventID,
        eventLabel,
        fInteractionType,
        fHitX, fHitY, fHitZ,
        fHitEnergy,
        fHitTime,
        fHitDeDx,
        fParticlePDG,
        fTrackID
    );
}

void EventAction::SaveSparseImages()
{
    if (!fSparseImageWriter) return;
    if (fHitX.empty()) return;

    G4String eventLabel = DetermineEventLabel();

    fSparseImageWriter->WriteEventSparse(
        fGlobalEventID,
        eventLabel,
        fInteractionType,
        fHitX, fHitY, fHitZ,
        fHitEnergy,
        fHitTime,
        fHitDeDx,
        fParticlePDG,
        fTrackID
    );
}

G4String EventAction::DetermineEventLabel() const
{

    if (fInteractionType == "Cosmic") {
        return "cosmic";
    }
    else if (fInteractionType.find("NuE_") != std::string::npos) {
        return "nue_cc";
    }
    else if (fInteractionType.find("NuMu_") != std::string::npos) {
        return "numu_cc";
    }
    else if (fInteractionType.find("NuTau_") != std::string::npos) {
        return "nutau_cc";
    }
    else if (fInteractionType.find("NC_") != std::string::npos) {
        return "nc";
    }
    else if (fInteractionType == "CCQE" || fInteractionType == "CC1Pi" || fInteractionType == "DIS") {
        return "numu_cc";
    }
    else {
        return "unknown";
    }
}

void EventAction::CalculateGeometricFeatures()
{
    if (fHitEnergy.empty() || fNHits == 0) {
        fVisibleEnergy = 0; fAvgHitEnergy = 0;
        fHitLengthX = fHitLengthY = fHitLengthZ = 0;
        fHitAspectRatio = 1.0;
        fPCAEigenvalue1 = fPCAEigenvalue2 = fPCAEigenvalue3 = 0;
        fPCARatio12 = fPCARatio13 = 1.0;
        fOpeningAngle = 0; fTrackAngleWrtBeam = 0;
        fDeDxMean = fDeDxStd = 0;
        fEnergyFrontFraction = fEnergyBackFraction = 0;
        fHitDensity = 0; fHitTimeSpread = 0; fNIsolatedHits = 0;
        fNHitsU = fNHitsV = fNHitsY = 0; fHitPlaneRatioYUV = 0;
        return;
    }

    fVisibleEnergy = fEnergyDeposit / MeV;

    fAvgHitEnergy = 0;
    for (const auto& energy : fHitEnergy) fAvgHitEnergy += energy;
    fAvgHitEnergy /= fNHits;

    G4double minX = *std::min_element(fHitX.begin(), fHitX.end());
    G4double maxX = *std::max_element(fHitX.begin(), fHitX.end());
    G4double minY = *std::min_element(fHitY.begin(), fHitY.end());
    G4double maxY = *std::max_element(fHitY.begin(), fHitY.end());
    G4double minZ = *std::min_element(fHitZ.begin(), fHitZ.end());
    G4double maxZ = *std::max_element(fHitZ.begin(), fHitZ.end());

    fHitLengthX = maxX - minX;
    fHitLengthY = maxY - minY;
    fHitLengthZ = maxZ - minZ;

    G4double maxLength = std::max({fHitLengthX, fHitLengthY, fHitLengthZ});
    fHitAspectRatio = (maxLength > 0.1) ? maxLength / std::max(0.1, std::min({fHitLengthX, fHitLengthY, fHitLengthZ})) : 1.0;

    CalculatePCAFeatures();
    CalculateOpeningAngle();
    CalculateDeDxProfile();
    CalculateHitDensity();
    CalculateWirePlaneHits();
}

void EventAction::CalculatePCAFeatures()
{
    if (fNHits < 3) {
        fPCAEigenvalue1 = fPCAEigenvalue2 = fPCAEigenvalue3 = 0;
        fPCARatio12 = fPCARatio13 = 1.0;
        return;
    }

    G4double meanX = 0, meanY = 0, meanZ = 0;
    for (size_t i = 0; i < fHitX.size(); i++) {
        meanX += fHitX[i]; meanY += fHitY[i]; meanZ += fHitZ[i];
    }
    meanX /= fNHits; meanY /= fNHits; meanZ /= fNHits;

    G4double cov[3][3] = {{0}};
    for (size_t i = 0; i < fHitX.size(); i++) {
        G4double dx = fHitX[i] - meanX;
        G4double dy = fHitY[i] - meanY;
        G4double dz = fHitZ[i] - meanZ;
        cov[0][0] += dx*dx; cov[0][1] += dx*dy; cov[0][2] += dx*dz;
        cov[1][1] += dy*dy; cov[1][2] += dy*dz; cov[2][2] += dz*dz;
    }
    cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cov[i][j] /= fNHits;

    G4double trace = cov[0][0] + cov[1][1] + cov[2][2];

    G4double eigenvalues[3];
    eigenvalues[0] = trace / 3.0;
    eigenvalues[1] = trace / 6.0;
    eigenvalues[2] = trace / 12.0;

    for (int iter = 0; iter < 10; iter++) {
        for (int i = 0; i < 3; i++) {
            G4double sum = 0;
            for (int j = 0; j < 3; j++) sum += cov[i][j] * cov[j][i];
            eigenvalues[i] = std::sqrt(sum);
        }
        G4double norm = std::sqrt(eigenvalues[0]*eigenvalues[0] +
                                  eigenvalues[1]*eigenvalues[1] +
                                  eigenvalues[2]*eigenvalues[2]);
        if (norm > 0) {
            for (int i = 0; i < 3; i++) eigenvalues[i] /= norm;
        }
    }

    std::sort(eigenvalues, eigenvalues + 3, std::greater<G4double>());
    fPCAEigenvalue1 = eigenvalues[0] * trace;
    fPCAEigenvalue2 = eigenvalues[1] * trace;
    fPCAEigenvalue3 = eigenvalues[2] * trace;
    fPCARatio12 = (fPCAEigenvalue2 > 1e-6) ? fPCAEigenvalue1 / fPCAEigenvalue2 : 100.0;
    fPCARatio13 = (fPCAEigenvalue3 > 1e-6) ? fPCAEigenvalue1 / fPCAEigenvalue3 : 100.0;
}

void EventAction::CalculateOpeningAngle()
{
    if (fNHits < 2) {
        fOpeningAngle = 0; fTrackAngleWrtBeam = 0;
        return;
    }

    G4ThreeVector vertex(fVertexX, fVertexY, fVertexZ);
    const size_t maxSamples = std::min((size_t)100, fHitX.size());
    std::vector<size_t> indices(fHitX.size());
    for (size_t i = 0; i < fHitX.size(); i++) indices[i] = i;

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    G4double maxAngle = 0;
    for (size_t i = 0; i < maxSamples && i < indices.size(); i++) {
        size_t idx1 = indices[i];
        G4ThreeVector pos1(fHitX[idx1], fHitY[idx1], fHitZ[idx1]);
        G4ThreeVector dir1 = (pos1 - vertex).unit();

        for (size_t j = i + 1; j < maxSamples && j < indices.size(); j++) {
            size_t idx2 = indices[j];
            G4ThreeVector pos2(fHitX[idx2], fHitY[idx2], fHitZ[idx2]);
            G4ThreeVector dir2 = (pos2 - vertex).unit();

            G4double dotProduct = std::max(-1.0, std::min(1.0, dir1.dot(dir2)));
            G4double angle = std::acos(dotProduct) * 180.0 / M_PI;
            if (angle > maxAngle) maxAngle = angle;
        }
    }
    fOpeningAngle = maxAngle;

    G4ThreeVector beamDir(0, 0, 1);
    G4ThreeVector primaryMom(fPrimaryPx, fPrimaryPy, fPrimaryPz);
    if (primaryMom.mag() > 0) {
        G4double dotBeam = std::max(-1.0, std::min(1.0, primaryMom.unit().dot(beamDir)));
        fTrackAngleWrtBeam = std::acos(dotBeam) * 180.0 / M_PI;
    } else {
        fTrackAngleWrtBeam = 0;
    }
}

void EventAction::CalculateDeDxProfile()
{
    if (fNHits < 3) {
        fDeDxMean = fDeDxStd = 0;
        fEnergyFrontFraction = fEnergyBackFraction = 0;
        return;
    }

    std::vector<std::pair<G4double, size_t>> distanceIndex;
    G4ThreeVector vertex(fVertexX, fVertexY, fVertexZ);

    distanceIndex.reserve(fHitX.size());
    for (size_t i = 0; i < fHitX.size(); i++) {
        G4ThreeVector hitPos(fHitX[i], fHitY[i], fHitZ[i]);
        G4double dist = (hitPos - vertex).mag();
        distanceIndex.push_back(std::make_pair(dist, i));
    }
    std::sort(distanceIndex.begin(), distanceIndex.end());

    std::vector<G4double> dEdxValues;
    dEdxValues.reserve(distanceIndex.size());
    G4double totalPathLength = 0;

    const G4double energyThreshold = 0.1;
    for (size_t i = 0; i < distanceIndex.size() - 1; i++) {
        size_t idx1 = distanceIndex[i].second;
        size_t idx2 = distanceIndex[i+1].second;

        if (fHitEnergy[idx1] < energyThreshold || fHitEnergy[idx2] < energyThreshold) continue;

        G4ThreeVector pos1(fHitX[idx1], fHitY[idx1], fHitZ[idx1]);
        G4ThreeVector pos2(fHitX[idx2], fHitY[idx2], fHitZ[idx2]);
        G4double segmentLength = (pos2 - pos1).mag();

        if (segmentLength > 0.1) {
            G4double segmentEnergy = (fHitEnergy[idx1] + fHitEnergy[idx2]) / 2.0;
            dEdxValues.push_back(segmentEnergy / segmentLength);
            totalPathLength += segmentLength;
        }
    }

    if (!dEdxValues.empty()) {
        fDeDxMean = 0;
        for (const auto& dedx : dEdxValues) fDeDxMean += dedx;
        fDeDxMean /= dEdxValues.size();

        fDeDxStd = 0;
        for (const auto& dedx : dEdxValues)
            fDeDxStd += (dedx - fDeDxMean) * (dedx - fDeDxMean);
        fDeDxStd = sqrt(fDeDxStd / dEdxValues.size());
    } else {
        fDeDxMean = fDeDxStd = 0;
    }

    G4int frontIndex = distanceIndex.size() * 0.3;
    G4int backIndex = distanceIndex.size() * 0.7;
    G4double frontEnergy = 0, backEnergy = 0;

    for (size_t i = 0; i < distanceIndex.size(); i++) {
        size_t idx = distanceIndex[i].second;
        if ((G4int)i < frontIndex) frontEnergy += fHitEnergy[idx];
        else if ((G4int)i >= backIndex) backEnergy += fHitEnergy[idx];
    }

    fEnergyFrontFraction = (fVisibleEnergy > 0) ? frontEnergy / fVisibleEnergy : 0;
    fEnergyBackFraction = (fVisibleEnergy > 0) ? backEnergy / fVisibleEnergy : 0;
}

void EventAction::CalculateHitDensity()
{
    if (fNHits < 2) {
        fHitDensity = 0; fHitTimeSpread = 0; fNIsolatedHits = 0;
        return;
    }

    G4double volumeX = std::max(1.0, fHitLengthX);
    G4double volumeY = std::max(1.0, fHitLengthY);
    G4double volumeZ = std::max(1.0, fHitLengthZ);
    fHitDensity = fNHits / (volumeX * volumeY * volumeZ);

    if (!fHitTime.empty()) {
        auto minMaxTime = std::minmax_element(fHitTime.begin(), fHitTime.end());
        G4double timeSpread = *minMaxTime.second - *minMaxTime.first;
        fHitTimeSpread = (timeSpread < 1e5) ? timeSpread : 0;
    } else {
        fHitTimeSpread = 0;
    }

    G4double isolationThreshold = 10.0;

    if (fNHits < 100) {
        fNIsolatedHits = 0;
        for (size_t i = 0; i < fHitX.size(); i++) {
            G4ThreeVector pos1(fHitX[i], fHitY[i], fHitZ[i]);
            G4double minDist = 1e10;
            for (size_t j = 0; j < fHitX.size(); j++) {
                if (i == j) continue;
                G4ThreeVector pos2(fHitX[j], fHitY[j], fHitZ[j]);
                G4double dist = (pos2 - pos1).mag();
                if (dist < minDist) minDist = dist;
                if (minDist < isolationThreshold) break;
            }
            if (minDist >= isolationThreshold) fNIsolatedHits++;
        }
    } else {
        G4double cellSize = isolationThreshold;
        std::map<std::tuple<int,int,int>, std::vector<size_t>> grid;

        for (size_t i = 0; i < fHitX.size(); i++) {
            int cx = static_cast<int>(std::floor(fHitX[i] / cellSize));
            int cy = static_cast<int>(std::floor(fHitY[i] / cellSize));
            int cz = static_cast<int>(std::floor(fHitZ[i] / cellSize));
            grid[std::make_tuple(cx,cy,cz)].push_back(i);
        }

        fNIsolatedHits = 0;
        for (size_t i = 0; i < fHitX.size(); i++) {
            G4ThreeVector pos1(fHitX[i], fHitY[i], fHitZ[i]);
            int cx = static_cast<int>(std::floor(fHitX[i] / cellSize));
            int cy = static_cast<int>(std::floor(fHitY[i] / cellSize));
            int cz = static_cast<int>(std::floor(fHitZ[i] / cellSize));

            G4double minDist = 1e10;
            bool foundNearby = false;

            for (int dx = -1; dx <= 1 && !foundNearby; dx++)
                for (int dy = -1; dy <= 1 && !foundNearby; dy++)
                    for (int dz = -1; dz <= 1 && !foundNearby; dz++) {
                        auto key = std::make_tuple(cx+dx, cy+dy, cz+dz);
                        if (grid.find(key) != grid.end()) {
                            for (size_t j : grid[key]) {
                                if (i == j) continue;
                                G4ThreeVector pos2(fHitX[j], fHitY[j], fHitZ[j]);
                                G4double dist = (pos2 - pos1).mag();
                                if (dist < minDist) minDist = dist;
                                if (minDist < isolationThreshold) { foundNearby = true; break; }
                            }
                        }
                    }

            if (minDist >= isolationThreshold) fNIsolatedHits++;
        }
    }
}

void EventAction::CalculateWirePlaneHits()
{

    const G4double halfY = DetectorConstruction::fTPC_Y / 2.0;
    const G4double halfZ = DetectorConstruction::fTPC_Z / 2.0;
    const G4double wireSpacing = DetectorConstruction::fWireSpacing;
    const G4double sin60 = std::sin(60.0 * M_PI / 180.0);
    const G4double cos60 = std::cos(60.0 * M_PI / 180.0);
    const G4double maxUV = halfY * sin60 + halfZ * cos60;

    std::set<G4int> wiresU, wiresV, wiresY;
    for (size_t i = 0; i < fHitY.size(); i++) {
        const G4double y = fHitY[i];
        const G4double z = fHitZ[i];

        const G4double uCoord = y * sin60 + z * cos60;
        const G4int wU = static_cast<G4int>((uCoord + maxUV) / wireSpacing);
        if (wU >= 0 && wU < DetectorConstruction::fNWires_U) wiresU.insert(wU);

        const G4double vCoord = -y * sin60 + z * cos60;
        const G4int wV = static_cast<G4int>((vCoord + maxUV) / wireSpacing);
        if (wV >= 0 && wV < DetectorConstruction::fNWires_V) wiresV.insert(wV);

        const G4int wY = static_cast<G4int>((z + halfZ) / wireSpacing);
        if (wY >= 0 && wY < DetectorConstruction::fNWires_Y) wiresY.insert(wY);
    }

    fNHitsU = static_cast<G4int>(wiresU.size());
    fNHitsV = static_cast<G4int>(wiresV.size());
    fNHitsY = static_cast<G4int>(wiresY.size());
    const G4double uvTotal = fNHitsU + fNHitsV;
    fHitPlaneRatioYUV = (uvTotal > 0) ? static_cast<G4double>(fNHitsY) / uvTotal : 0.0;
}

void EventAction::SaveEventData()
{
    fRunAction->WriteEventData(
        fGlobalEventID, fEventType, fInteractionType,
        fNeutrinoEnu, fNeutrinoQ2, fNeutrinoW,
        fVertexX, fVertexY, fVertexZ,
        fPrimaryPx, fPrimaryPy, fPrimaryPz,
        fPrimaryP, fPrimaryE, fPrimaryTheta, fPrimaryPhi, fPrimaryPDG,
        fNHits, fVisibleEnergy, fAvgHitEnergy,
        fHitLengthX, fHitLengthY, fHitLengthZ, fHitAspectRatio,
        fPCAEigenvalue1, fPCAEigenvalue2, fPCAEigenvalue3,
        fPCARatio12, fPCARatio13,
        fOpeningAngle, fTrackAngleWrtBeam,
        fDeDxMean, fDeDxStd, fEnergyFrontFraction, fEnergyBackFraction,
        fHitDensity, fHitTimeSpread, fNIsolatedHits,
        fNHitsU, fNHitsV, fNHitsY, fHitPlaneRatioYUV
    );
}

void EventAction::ResetEventData()
{
    fEnergyDeposit = 0.; fTrackLength = 0.;
    fHitX.clear(); fHitY.clear(); fHitZ.clear();
    fHitEnergy.clear(); fHitTime.clear();
    fWirePlane.clear(); fWireNumber.clear();
    fParticlePDG.clear();
    fTrackID.clear();
    fHitDeDx.clear();

    fNeutrinoQ2 = -1.; fNeutrinoW = -1.; fNeutrinoEnu = -1.;
    fInteractionType = "Unknown";
    fVertexX = fVertexY = fVertexZ = 0.;
    fPrimaryPx = fPrimaryPy = fPrimaryPz = fPrimaryP = 0.;
    fPrimaryE = 0.; fPrimaryTheta = fPrimaryPhi = 0.; fPrimaryPDG = 0;
    fNHits = 0; fVisibleEnergy = 0; fAvgHitEnergy = 0;
    fHitLengthX = fHitLengthY = fHitLengthZ = 0; fHitAspectRatio = 1.0;
    fPCAEigenvalue1 = fPCAEigenvalue2 = fPCAEigenvalue3 = 0;
    fPCARatio12 = fPCARatio13 = 1.0;
    fOpeningAngle = 0; fTrackAngleWrtBeam = 0;
    fDeDxMean = fDeDxStd = 0;
    fEnergyFrontFraction = fEnergyBackFraction = 0;
    fHitDensity = 0; fHitTimeSpread = 0; fNIsolatedHits = 0;
    fNHitsU = fNHitsV = fNHitsY = 0; fHitPlaneRatioYUV = 0;
}

void EventAction::SetNeutrinoKinematics(G4double Q2, G4double W, G4double Enu, G4String type)
{
    fNeutrinoQ2 = Q2; fNeutrinoW = W; fNeutrinoEnu = Enu; fInteractionType = type;
}
