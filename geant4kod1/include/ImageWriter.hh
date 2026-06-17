#ifndef ImageWriter_h
#define ImageWriter_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <string>
#include <fstream>
#include <map>

class ImageWriter
{
public:

    enum ProjectionType {
        kXZ = 0,
        kYZ = 1,
        kXY = 2
    };

    enum PixelContent {
        kEnergyDeposition,
        kHitCount,
        kTimeWeighted
    };

    ImageWriter();
    ~ImageWriter();

    void SetOutputDirectory(const G4String& dir) { fOutputDir = dir; }
    void SetImageResolution(G4int width, G4int height);
    void SetPixelContent(PixelContent content) { fPixelContent = content; }
    void SetLogScale(G4bool useLog) { fUseLogScale = useLog; }

    void SetTPCBounds(G4double xMin, G4double xMax,
                      G4double yMin, G4double yMax,
                      G4double zMin, G4double zMax);

    void WriteEventImages(
        G4int eventID,
        const G4String& eventLabel,
        const G4String& interactionType,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        const std::vector<G4int>& hitPDG,
        const std::vector<G4int>& hitTrackID
    );

    void InitializeLabelsFile(const G4String& filename);
    void FinalizeLabelsFile();

    G4int GetImagesWritten() const { return fImagesWritten; }

private:

    void FillProjection(
        ProjectionType proj,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4double>& hitTime,
        std::vector<float>& pixelData
    );

    void FillPDGMask(
        ProjectionType proj,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4int>& hitPDG,
        std::vector<int32_t>& maskData
    );

    void FillTrackMask(
        ProjectionType proj,
        const std::vector<G4double>& hitX,
        const std::vector<G4double>& hitY,
        const std::vector<G4double>& hitZ,
        const std::vector<G4double>& hitEnergy,
        const std::vector<G4int>& hitTrackID,
        std::vector<int32_t>& maskData
    );

    G4bool GetPixelCoords(ProjectionType proj,
                          G4double x, G4double y, G4double z,
                          G4int& pixelRow, G4int& pixelCol) const;

    void WriteBinaryFloat(const G4String& filename,
                          const std::vector<float>& data,
                          G4int width, G4int height);
    void WriteBinaryInt(const G4String& filename,
                        const std::vector<int32_t>& data,
                        G4int width, G4int height);

    void ApplyLogScale(std::vector<float>& data);

    G4String GetProjectionName(ProjectionType proj) const;

    G4String fOutputDir;
    G4int fImageWidth;
    G4int fImageHeight;
    PixelContent fPixelContent;
    G4bool fUseLogScale;

    G4double fXMin, fXMax;
    G4double fYMin, fYMax;
    G4double fZMin, fZMax;

    std::ofstream fLabelsFile;
    G4bool fLabelsFileOpen;

    G4int fImagesWritten;
};

#endif
