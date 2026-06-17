#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

class RunAction;
class ImageWriter;
class WireImageWriter;
class SparseImageWriter;
class G4Event;

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEnergyDeposit(G4double edep) { fEnergyDeposit += edep; }
    void AddTrackLength(G4double length) { fTrackLength += length; }

    void SetNeutrinoKinematics(G4double Q2, G4double W, G4double Enu, G4String type);
    G4double GetNeutrinoQ2() const { return fNeutrinoQ2; }
    G4double GetNeutrinoW() const { return fNeutrinoW; }
    G4double GetNeutrinoEnu() const { return fNeutrinoEnu; }
    G4String GetInteractionType() const { return fInteractionType; }

    ImageWriter* GetImageWriter() { return fImageWriter; }

private:
    void CollectHitData(const G4Event* event);
    void CalculateGeometricFeatures();
    void SaveEventData();
    void SaveEventImages();
    void SaveWireImages();
    void SaveSparseImages();
    void ResetEventData();

    void CalculatePCAFeatures();
    void CalculateDeDxProfile();
    void CalculateHitDensity();
    void CalculateWirePlaneHits();
    void CalculateOpeningAngle();

    G4String DetermineEventLabel() const;

    RunAction* fRunAction;
    ImageWriter* fImageWriter;
    WireImageWriter* fWireImageWriter;
    SparseImageWriter* fSparseImageWriter;

    G4double fEnergyDeposit;
    G4double fTrackLength;
    G4int fEventID;
    G4int fGlobalEventID;
    G4int fEventType;

    static G4int fGlobalEventCounter;

    G4double fNeutrinoQ2;
    G4double fNeutrinoW;
    G4double fNeutrinoEnu;
    G4String fInteractionType;

    std::vector<G4double> fHitX, fHitY, fHitZ;
    std::vector<G4double> fHitEnergy;
    std::vector<G4double> fHitTime;
    std::vector<G4int> fWirePlane;
    std::vector<G4int> fWireNumber;
    std::vector<G4int> fParticlePDG;
    std::vector<G4int> fTrackID;
    std::vector<G4double> fHitDeDx;

    G4double fVertexX, fVertexY, fVertexZ;
    G4double fPrimaryPx, fPrimaryPy, fPrimaryPz;
    G4double fPrimaryE;
    G4double fPrimaryP;
    G4double fPrimaryTheta;
    G4double fPrimaryPhi;
    G4int fPrimaryPDG;

    G4int fNHits;
    G4double fVisibleEnergy;
    G4double fAvgHitEnergy;

    G4double fHitLengthX, fHitLengthY, fHitLengthZ;
    G4double fHitAspectRatio;

    G4double fPCAEigenvalue1, fPCAEigenvalue2, fPCAEigenvalue3;
    G4double fPCARatio12, fPCARatio13;

    G4double fOpeningAngle;
    G4double fTrackAngleWrtBeam;

    G4double fDeDxMean, fDeDxStd;
    G4double fEnergyFrontFraction, fEnergyBackFraction;

    G4double fHitDensity;
    G4double fHitTimeSpread;
    G4int fNIsolatedHits;

    G4int fNHitsU, fNHitsV, fNHitsY;
    G4double fHitPlaneRatioYUV;
};

#endif
