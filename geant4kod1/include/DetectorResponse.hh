#ifndef DetectorResponse_h
#define DetectorResponse_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>
#include <cmath>

class DetectorResponse
{
public:
    DetectorResponse();
    ~DetectorResponse();

    G4double CalculateRecombination(G4double dEdx_MeV_per_cm) const;

    G4double EnergyToElectrons(G4double energyMeV, G4double dEdx_MeV_per_cm) const;

    G4double CalculateDriftTime(G4double distanceToPlaneCm) const;
    G4double CalculateLifetimeAttenuation(G4double driftTimeMicrosec) const;

    G4double CalculateLongitudinalDiffusion(G4double driftTimeMicrosec) const;
    G4double CalculateTransverseDiffusion(G4double driftTimeMicrosec) const;

    G4double GetElectronicsNoise() const;
    G4double GetADCConversion() const;

    std::vector<G4double> GetFieldResponse(G4int planeID) const;

    std::vector<G4double> GetElectronicsResponse() const;

    std::vector<G4double> GetSignalResponse(G4int planeID) const;

    void GenerateNoise(std::vector<float>& waveform, G4int planeID) const;

    std::vector<float> ApplyWienerDeconvolution(
        const std::vector<float>& rawWaveform,
        const std::vector<G4double>& signalResponse,
        G4double noiseRMS,
        G4int planeID = 2);

    struct ResponseDFTCache {
        std::vector<G4double> re;
        std::vector<G4double> im;
    };

    const ResponseDFTCache& GetCachedResponseDFT(
        const std::vector<G4double>& signalResponse,
        G4int planeID, G4int N);

    struct TrigTable {
        std::vector<G4double> cosTable;
        std::vector<G4double> sinTable;
    };
    const TrigTable& GetCachedTrigTable(G4int N);

    void SetElectricField(G4double field) { fElectricField = field; }
    G4double GetElectricField() const { return fElectricField; }

    G4double GetDriftVelocity() const { return fDriftVelocity; }

    void SetElectronLifetime(G4double lifetime) { fElectronLifetime = lifetime; }

    G4double GetSamplingPeriod() const { return fSamplingPeriod; }
    G4double GetSamplingRate() const { return 1.0 / fSamplingPeriod; }

    G4double GetWireSpacing() const { return fWireSpacing; }

    G4double GetLArDensity() const { return fLArDensity; }

    G4double GetWionization() const { return fWionization; }

private:

    G4double fElectricField;
    G4double fDriftVelocity;
    G4double fElectronLifetime;
    G4double fSamplingPeriod;
    G4double fWireSpacing;
    G4double fLArDensity;
    G4double fTemperature;
    G4double fWionization;

    G4double fRecombAlpha;
    G4double fRecombBeta;

    G4double fDiffusionL;
    G4double fDiffusionT;

    G4double fENC;
    G4double fADCPerElectron;
    G4double fShapingTime;

    G4double fFieldResponseTime;
    G4int fResponseNTicks;

    G4double fNoise1fKnee;
    G4double fCoherentNoiseRMS;

    G4double fWienerFilterCutoff;

    std::map<G4int, ResponseDFTCache> fResponseDFTCache;

    std::map<G4int, TrigTable> fTrigTableCache;
};

#endif
