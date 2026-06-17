#ifndef WireImageWriter_h
#define WireImageWriter_h 1

#include "globals.hh"
#include "DetectorResponse.hh"
#include <vector>
#include <string>
#include <fstream>
#include <map>

class WireImageWriter
{
public:
    WireImageWriter();
    ~WireImageWriter();

    void SetOutputDirectory(const G4String& dir) { fOutputDir = dir; }

    void SetImageDimensions(G4int nWires, G4int nTicks);

    void SetApplyRecombination(G4bool apply) { fApplyRecombination = apply; }
    void SetApplyDiffusion(G4bool apply) { fApplyDiffusion = apply; }
    void SetApplyLifetime(G4bool apply) { fApplyLifetime = apply; }
    void SetApplyNoise(G4bool apply) { fApplyNoise = apply; }
    void SetApplyConvolution(G4bool apply) { fApplyConvolution = apply; }
    void SetApplyDeconvolution(G4bool apply) { fApplyDeconvolution = apply; }
    void SetApplyAllEffects(G4bool apply);

    void WriteEventImages(
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
        const std::vector<G4int>& hitTrackID
    );

    void InitializeLabelsFile(const G4String& filename);
    void FinalizeLabelsFile();

    G4int GetImagesWritten() const { return fImagesWritten; }
    DetectorResponse* GetDetectorResponse() { return fDetResponse; }

private:

    struct WirePlane {
        G4int planeID;
        G4String name;
        G4double angle;
        G4int nWiresReal;
        G4double planeZ;
        G4bool isCollection;
    };

    std::vector<WirePlane> fWirePlanes;
    void InitializeWirePlanes();

    G4double CalculateWireCoord(const WirePlane& plane,
                                 G4double y, G4double z) const;

    G4int CalculateWireNumber(const WirePlane& plane,
                               G4double y, G4double z) const;

    G4double CalculateDriftDistance(const WirePlane& plane, G4double x) const;

    void FillWireTimeImage(
        const WirePlane& plane,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        const std::vector<G4double>& hitDeDx,
        std::vector<float>& imageData
    );

    void DepositChargeOnWires(
        const WirePlane& plane,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        const std::vector<G4double>& hitDeDx,
        std::map<G4int, std::vector<float>>& wireWaveforms,
        G4int nPhysicalTicks
    );

    void ConvolveWireWithResponse(
        std::vector<float>& waveform,
        const std::vector<G4double>& signalResponse
    );

    void DownsampleToImage(
        const std::map<G4int, std::vector<float>>& wireWaveforms,
        G4int nPhysicalTicks,
        G4int nWiresReal,
        std::vector<float>& imageData
    );

    void FillWireTimePDGMask(
        const WirePlane& plane,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4int>& hitPDG,
        std::vector<int32_t>& maskData
    );

    struct TruthEntry {
        int16_t wire;
        int16_t tick;
        int32_t pdg;
        float   energy;
    };

    void FillWireTimeTruthSparse(
        const WirePlane& plane,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        const std::vector<G4int>& hitPDG,
        std::vector<TruthEntry>& truthData
    );

    void WriteTruthCoo(const G4String& filename,
                       const std::vector<TruthEntry>& truthData);

    void AddNoise(std::vector<float>& imageData, G4int planeID);

    void ApplyThreshold(std::vector<float>& imageData, float threshold);

    void WriteSparseImage(const G4String& filename, const std::vector<float>& imageData);
    void WriteBinaryFloat(const G4String& filename, const std::vector<float>& data);
    void WriteBinaryInt(const G4String& filename, const std::vector<int32_t>& data);

    std::vector<G4double> PrecomputeNetKernel(G4int planeID);

    DetectorResponse* fDetResponse;
    G4String fOutputDir;

    G4int fNWireBins;
    G4int fNTickBins;

    G4double fTPC_X;
    G4double fTPC_Y;
    G4double fTPC_Z;

    std::vector<G4double> fWireCoordMin;
    std::vector<G4double> fWireCoordMax;

    G4bool fApplyRecombination;
    G4bool fApplyDiffusion;
    G4bool fApplyLifetime;
    G4bool fApplyNoise;
    G4bool fApplyConvolution;
    G4bool fApplyDeconvolution;

    std::vector<std::vector<G4double>> fSignalResponse;

    std::vector<std::vector<G4double>> fNetKernel;

    std::ofstream fLabelsFile;
    G4bool fLabelsFileOpen;

    G4int fImagesWritten;
};

#endif
