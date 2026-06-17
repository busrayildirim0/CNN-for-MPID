#ifndef SparseImageWriter_h
#define SparseImageWriter_h 1

#include "globals.hh"
#include "DetectorResponse.hh"
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <tuple>

class SparseImageWriter
{
public:
    SparseImageWriter();
    ~SparseImageWriter();

    void SetOutputDirectory(const G4String& dir) { fOutputDir = dir; }
    void SetImageDimensions(G4int nWires, G4int nTicks);

    void SetApplyRecombination(G4bool b) { fApplyRecombination = b; }
    void SetApplyDiffusion(G4bool b) { fApplyDiffusion = b; }
    void SetApplyLifetime(G4bool b) { fApplyLifetime = b; }
    void SetApplyNoise(G4bool b) { fApplyNoise = b; }
    void SetApplyAllEffects(G4bool b);

    void WriteEventSparse(
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
        G4bool isCollection;
    };

    std::vector<WirePlane> fWirePlanes;
    void InitializeWirePlanes();

    struct SparsePixel {
        int16_t wire;
        int16_t tick;
        float charge;
    };

    struct TruthPixel {
        int16_t wire;
        int16_t tick;
        G4int pdg;
        float energy;
    };

    G4double CalculateWireCoord(const WirePlane& plane, G4double y, G4double z) const;

    void FillSparseView(
        const WirePlane& plane,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        const std::vector<G4double>& hitDeDx,
        std::vector<SparsePixel>& sparseData
    );

    void FillTruthData(
        const WirePlane& plane,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        const std::vector<G4int>& hitPDG,
        std::vector<TruthPixel>& truthData
    );

    void WriteNPZ(
        const G4String& filename,
        const std::vector<SparsePixel> sparseViews[3],
        const std::vector<TruthPixel>& truthData,
        G4int eventID,
        const G4String& eventLabel,
        const G4String& interactionType,
        G4int nHits,
        G4double visibleEnergy
    );

    void WriteNpyArray_int16(std::ostream& os, const std::string& name,
                              const std::vector<int16_t>& data);
    void WriteNpyArray_int32(std::ostream& os, const std::string& name,
                              const std::vector<int32_t>& data);
    void WriteNpyArray_float(std::ostream& os, const std::string& name,
                              const std::vector<float>& data);

    DetectorResponse* fDetResponse;
    G4String fOutputDir;
    G4int fNWireBins;
    G4int fNTickBins;
    G4double fTPC_X, fTPC_Y, fTPC_Z;

    G4bool fApplyRecombination;
    G4bool fApplyDiffusion;
    G4bool fApplyLifetime;
    G4bool fApplyNoise;

    std::ofstream fLabelsFile;
    G4bool fLabelsFileOpen;
    G4int fImagesWritten;
};

#endif
