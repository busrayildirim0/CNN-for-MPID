#include "LArSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include <cmath>

LArSensitiveDetector::LArSensitiveDetector(const G4String& name, const G4String& hitsCollectionName)
    : G4VSensitiveDetector(name),
      fHitsCollection(nullptr),
      fHCID(-1)
{
    collectionName.insert(hitsCollectionName);
}

LArSensitiveDetector::~LArSensitiveDetector()
{}

void LArSensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new LArHitsCollection(SensitiveDetectorName, collectionName[0]);

    if (fHCID < 0) {
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }

    hce->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool LArSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return false;

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4ThreeVector position = preStepPoint->GetPosition();
    G4double time = preStepPoint->GetGlobalTime();

    G4Track* track = step->GetTrack();
    G4int pdgCode = track->GetDefinition()->GetPDGEncoding();
    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();

    G4String processName = "";
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (process) processName = process->GetProcessName();

    G4double stepLength = step->GetStepLength() / cm;
    G4double dEdx;
    if (stepLength > 0.001) {
        dEdx = (edep / MeV) / stepLength;
    } else {

        G4int absPDG = std::abs(pdgCode);
        if (absPDG == 2212)      dEdx = 10.0;
        else if (absPDG == 13)   dEdx = 2.1;
        else if (absPDG == 211)  dEdx = 2.1;
        else if (absPDG == 11)   dEdx = 2.0;
        else if (absPDG == 22)   dEdx = 2.0;
        else if (absPDG == 1000180400 || absPDG > 100000) dEdx = 20.0;
        else                     dEdx = 2.1;
    }

    dEdx = std::max(0.5, std::min(50.0, dEdx));

    LArHit* newHit = new LArHit();

    newHit->SetWirePlaneID(-1);
    newHit->SetWireNumber(-1);
    newHit->SetEnergyDeposit(edep);
    newHit->SetTime(time);
    newHit->SetPosition(position);
    newHit->SetParticleType(pdgCode);
    newHit->SetTrackID(trackID);
    newHit->SetParentTrackID(parentID);
    newHit->SetProcessName(processName);
    newHit->SetDeDx(dEdx);

    G4int photons = static_cast<G4int>(edep / MeV * 24000);
    newHit->SetScintillationPhotons(photons);

    fHitsCollection->insert(newHit);

    return true;
}

std::vector<LArSensitiveDetector::WireResponse>
LArSensitiveDetector::CalculateWireResponse(const G4ThreeVector& position,
                                             G4double energy, G4double time,
                                             G4double dEdx)
{
    std::vector<WireResponse> responses;

    for (G4int planeID = 0; planeID < 3; planeID++) {
        G4int wireNumber = CalculateWireNumber(position, planeID);
        G4double driftTime = CalculateDriftTime(position, planeID);
        G4double signalStrength = CalculateSignalStrength(position, energy, planeID, dEdx);

        G4int maxWires = (planeID == 0) ? fNWires_U :
                         (planeID == 1) ? fNWires_V :
                                          fNWires_Y;

        if (wireNumber >= 0 && wireNumber < maxWires) {
            WireResponse response;
            response.wireNumber = wireNumber;
            response.planeID = planeID;
            response.signal = signalStrength;
            response.time = time + driftTime;
            responses.push_back(response);
        }
    }

    return responses;
}

G4int LArSensitiveDetector::CalculateWireNumber(const G4ThreeVector& position, G4int planeID)
{
    G4double y = position.y() / cm;
    G4double z = position.z() / cm;

    G4double wireCoord = 0.0;
    G4double halfY = DetectorConstruction::fTPC_Y / 2.0;
    G4double halfZ = DetectorConstruction::fTPC_Z / 2.0;
    G4int wireNumber = -1;

    switch(planeID) {
        case 0:
        {

            wireCoord = y * std::sin(60.0 * M_PI / 180.0) + z * std::cos(60.0 * M_PI / 180.0);

            G4double minCoord = -halfY * std::sin(60.0 * M_PI / 180.0) - halfZ * std::cos(60.0 * M_PI / 180.0);
            wireNumber = static_cast<G4int>((wireCoord - minCoord) / fWireSpacing);
            break;
        }
        case 1:
        {

            wireCoord = -y * std::sin(60.0 * M_PI / 180.0) + z * std::cos(60.0 * M_PI / 180.0);
            G4double minCoord = -halfY * std::sin(60.0 * M_PI / 180.0) - halfZ * std::cos(60.0 * M_PI / 180.0);
            wireNumber = static_cast<G4int>((wireCoord - minCoord) / fWireSpacing);
            break;
        }
        case 2:
        {
            wireNumber = static_cast<G4int>((z + halfZ) / fWireSpacing);
            break;
        }
    }

    return wireNumber;
}

G4double LArSensitiveDetector::CalculateDriftTime(const G4ThreeVector& position, G4int planeID)
{
    G4double x = position.x() / cm;
    G4double halfX = DetectorConstruction::fTPC_X / 2.0;

    G4double planeX;
    switch(planeID) {
        case 0:  planeX = halfX - DetectorConstruction::fWirePlaneU_X; break;
        case 1:  planeX = halfX - DetectorConstruction::fWirePlaneV_X; break;
        case 2:  planeX = halfX - DetectorConstruction::fWirePlaneY_X; break;
        default: planeX = halfX; break;
    }

    G4double driftDistance = planeX - x;
    if (driftDistance < 0.0) driftDistance = 0.0;
    if (driftDistance > DetectorConstruction::fTPC_X) driftDistance = DetectorConstruction::fTPC_X;

    G4double driftTime = driftDistance / fDriftVelocity;
    return driftTime * microsecond;
}

G4double LArSensitiveDetector::CalculateSignalStrength(const G4ThreeVector& position,
                                                        G4double energy, G4int planeID,
                                                        G4double )
{
    G4double signalFactor = 1.0;

    switch(planeID) {
        case 0:  signalFactor = 0.75; break;
        case 1:  signalFactor = 0.80; break;
        case 2:  signalFactor = 1.00; break;
    }

    G4double x = position.x() / cm;
    G4double halfX = DetectorConstruction::fTPC_X / 2.0;
    G4double driftDistance = halfX - x;
    if (driftDistance < 0.0) driftDistance = 0.0;

    G4double driftTime = driftDistance / fDriftVelocity;
    G4double electronLifetime = 18000.0;
    G4double attenuation = std::exp(-driftTime / electronLifetime);

    return energy * signalFactor * attenuation;
}

void LArSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
    if (verboseLevel > 1) {
        G4int nofHits = fHitsCollection->entries();
        G4cout << G4endl
               << "-------->Hits Collection: " << nofHits
               << " hits in the LAr TPC:" << G4endl;

        G4int maxPrint = std::min(nofHits, 10);
        for (G4int i = 0; i < maxPrint; i++) {
            (*fHitsCollection)[i]->Print();
        }
        if (nofHits > maxPrint) {
            G4cout << "... and " << (nofHits - maxPrint) << " more hits" << G4endl;
        }
    }
}
