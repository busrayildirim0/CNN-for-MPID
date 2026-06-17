#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class EventInformation;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

    enum GeneratorMode {
        kCosmicRayMode,
        kNeutrinoMode,
        kTestMode
    };

    enum NeutrinoFlavor {
        kNuMu,
        kNuE,
        kNuTau,
        kNC,
        kAllFlavors
    };

    enum SamplingProfile {
        kRealisticProfile,
        kMLBalancedProfile
    };

    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    void SetGeneratorMode(GeneratorMode mode) { fGeneratorMode = mode; }
    GeneratorMode GetGeneratorMode() const { return fGeneratorMode; }

    void SetNeutrinoFlavor(NeutrinoFlavor flavor) { fNeutrinoFlavor = flavor; }
    NeutrinoFlavor GetNeutrinoFlavor() const { return fNeutrinoFlavor; }

    void SetSamplingProfile(SamplingProfile profile) { fSamplingProfile = profile; }
    SamplingProfile GetSamplingProfile() const { return fSamplingProfile; }

    void SetCosmicEnergyRange(G4double minE, G4double maxE);
    void SetCosmicAngleRange(G4double minTheta, G4double maxTheta);
    void SetNeutrinoEnergyRange(G4double minE, G4double maxE);
    void SetEnableCosmicOverlay(G4bool enable) { fEnableCosmicOverlay = enable; }
    G4bool GetEnableCosmicOverlay() const { return fEnableCosmicOverlay; }

    void SetFlavorRatios(G4double nuMuFrac, G4double nuEFrac,
                         G4double nuTauFrac, G4double ncFrac);

private:

    void GenerateCosmicRay(G4Event* anEvent, EventInformation* eventInfo);
    void GenerateNeutrinoEvent(G4Event* anEvent, EventInformation* eventInfo);
    void GenerateTestParticle(G4Event* anEvent);
    void GenerateCosmicOverlay(G4Event* anEvent);

    void GenerateNuMuCCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNuMuCC1Pi(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNuMuDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);

    void GenerateNueCCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNueCC1Pi(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNueDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);

    void GenerateNuTauCCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNuTauDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);

    void GenerateNCQE(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNCRes(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);
    void GenerateNCDIS(G4Event* anEvent, EventInformation* eventInfo, G4double Enu);

    void Generate2p2hMEC(G4Event* anEvent, EventInformation* eventInfo,
                         G4double Enu, G4int leptonPDG);

    void SimulateTauDecay(G4Event* anEvent, const G4ThreeVector& position,
                          G4double tauEnergy, const G4ThreeVector& tauDirection);

    void ApplyFSI(G4Event* anEvent, const G4ThreeVector& vertex,
                  std::vector<G4int>& hadronPDGs,
                  std::vector<G4double>& hadronKEs,
                  std::vector<G4ThreeVector>& hadronDirs);

    G4bool PassesPauliBlocking(G4double nucleonMomentum) const;

    G4ThreeVector SampleForwardDirection(const G4ThreeVector& axis,
                                         G4double thetaScale);

    G4double SampleCosmicEnergy();
    G4double SampleCosmicAngle();
    G4ThreeVector SampleCosmicPosition();
    G4ThreeVector FindLArEntryPoint(G4ThreeVector position, G4ThreeVector direction);

    void InitializeBNBFlux();
    G4double SampleBNBEnergy();
    G4double SampleCCQE_Q2(G4double Enu);
    G4double SampleResonance_W();
    G4double SampleDIS_y(G4double Enu);

    G4ThreeVector SampleVertexInTPC();

    G4ParticleGun* fParticleGun;
    GeneratorMode fGeneratorMode;
    NeutrinoFlavor fNeutrinoFlavor;
    SamplingProfile fSamplingProfile;

    G4double fCosmicMinEnergy;
    G4double fCosmicMaxEnergy;
    G4double fCosmicMinTheta;
    G4double fCosmicMaxTheta;

    G4double fNeutrinoMinEnergy;
    G4double fNeutrinoMaxEnergy;

    G4int fCosmicEventCount;
    G4int fNeutrinoEventCount;
    G4int fNuMuCCCount;
    G4int fNueCCCount;
    G4int fNuTauCCCount;
    G4int fNCCount;
    G4int f2p2hCount;
    G4int fFSIAbsorptionCount;
    G4int fFSIChargeExchangeCount;

    static const G4int fNBNBBins = 60;
    G4double fBNBEnergyBins[61];
    G4double fBNBFluxWeights[60];
    G4double fBNBCumulativeWeights[60];

    G4double fNuMuFraction;
    G4double fNuEFraction;
    G4double fNuTauFraction;
    G4double fNCFraction;
    G4bool fEnableCosmicOverlay;
    G4double fCosmicOverlayMean;

    PrimaryGeneratorMessenger* fMessenger;
};

#endif
