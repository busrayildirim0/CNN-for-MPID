#include "DetectorConstruction.hh"
#include "LArSensitiveDetector.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SDManager.hh"

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fLAr(nullptr), fAir(nullptr), fSteel(nullptr),
      fWorldLogical(nullptr), fLArLogical(nullptr), fCryostatLogical(nullptr),
      fWorldPhysical(nullptr), fLArPhysical(nullptr)
{
    DefineMaterials();
}

DetectorConstruction::~DetectorConstruction()
{}

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* Ar = nist->FindOrBuildElement("Ar");

    fLAr = new G4Material("LiquidArgon", 1.396*g/cm3, 1, kStateLiquid, 87.*kelvin);
    fLAr->AddElement(Ar, 1);

    const G4int nEntries = 32;

    G4double photonEnergy[nEntries] = {
        2.034*eV,2.068*eV,2.103*eV,2.139*eV,2.177*eV,2.216*eV,2.256*eV,2.298*eV,
        2.341*eV,2.386*eV,2.433*eV,2.481*eV,2.532*eV,2.585*eV,2.640*eV,2.697*eV,
        2.757*eV,2.820*eV,2.885*eV,2.954*eV,3.026*eV,3.102*eV,3.181*eV,3.265*eV,
        3.353*eV,3.446*eV,3.545*eV,3.649*eV,3.760*eV,3.877*eV,4.002*eV,4.136*eV
    };

    G4double refractiveIndex[nEntries];
    G4double absorption[nEntries];
    G4double scintillation[nEntries];
    G4double rayleigh[nEntries];

    for(int i=0;i<nEntries;i++)
    {
        refractiveIndex[i]=1.38;
        absorption[i]=100.*m;
        scintillation[i]=1.0;
        rayleigh[i]=90.*cm;
    }

    G4MaterialPropertiesTable* larMPT = new G4MaterialPropertiesTable();

    larMPT->AddProperty("RINDEX",photonEnergy,refractiveIndex,nEntries);
    larMPT->AddProperty("ABSLENGTH",photonEnergy,absorption,nEntries);
    larMPT->AddProperty("RAYLEIGH",photonEnergy,rayleigh,nEntries);
    larMPT->AddProperty("SCINTILLATIONCOMPONENT1",photonEnergy,scintillation,nEntries);
    larMPT->AddProperty("SCINTILLATIONCOMPONENT2",photonEnergy,scintillation,nEntries);

    larMPT->AddConstProperty("SCINTILLATIONYIELD",24000./MeV,true);
    larMPT->AddConstProperty("RESOLUTIONSCALE",1.0,true);
    larMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",6.*ns,true);
    larMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2",1590.*ns,true);
    larMPT->AddConstProperty("SCINTILLATIONYIELD1",0.23,true);
    larMPT->AddConstProperty("SCINTILLATIONYIELD2",0.77,true);

    fLAr->SetMaterialPropertiesTable(larMPT);

    fAir   = nist->FindOrBuildMaterial("G4_AIR");
    fSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

    G4double cryostatThickness = 10.*cm;
    G4double cryostatHalfZ = 0.5 * (fTPC_Z * cm + 2 * cryostatThickness);
    G4double worldSize = 2.0 * (cryostatHalfZ + 1.0 * m);

    G4Box* solidWorld = new G4Box("World",
                                  0.5*worldSize,
                                  0.5*worldSize,
                                  0.5*worldSize);

    fWorldLogical =
        new G4LogicalVolume(solidWorld,
                            fAir,
                            "World");

    fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

    fWorldPhysical =
        new G4PVPlacement(0,
                          G4ThreeVector(),
                          fWorldLogical,
                          "World",
                          0,
                          false,
                          0,
                          true);

    G4Box* solidCryostat =
        new G4Box("Cryostat",
                  0.5*(fTPC_X*cm + 2*cryostatThickness),
                  0.5*(fTPC_Y*cm + 2*cryostatThickness),
                  0.5*(fTPC_Z*cm + 2*cryostatThickness));

    fCryostatLogical =
        new G4LogicalVolume(solidCryostat,
                            fSteel,
                            "Cryostat");

    new G4PVPlacement(0,
                      G4ThreeVector(),
                      fCryostatLogical,
                      "Cryostat",
                      fWorldLogical,
                      false,
                      0,
                      true);

    G4Box* solidLAr =
        new G4Box("LArTPC",
                  0.5*fTPC_X*cm,
                  0.5*fTPC_Y*cm,
                  0.5*fTPC_Z*cm);

    fLArLogical =
        new G4LogicalVolume(solidLAr,
                            fLAr,
                            "LArTPC");

    fLArPhysical =
        new G4PVPlacement(0,
                          G4ThreeVector(),
                          fLArLogical,
                          "LArTPC",
                          fCryostatLogical,
                          false,
                          0,
                          true);

    G4VisAttributes* cryoVis =
        new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.15));
    cryoVis->SetForceSolid(true);
    fCryostatLogical->SetVisAttributes(cryoVis);

    G4VisAttributes* larVis =
        new G4VisAttributes(G4Colour(0.0,0.6,1.0,0.45));
    larVis->SetForceSolid(true);
    fLArLogical->SetVisAttributes(larVis);

    return fWorldPhysical;
}

void DetectorConstruction::ConstructSDandField()
{
    G4String sdName = "LArTPC/LArSD";

    LArSensitiveDetector* larSD =
        new LArSensitiveDetector(sdName,
                                 "LArHitsCollection");

    G4SDManager::GetSDMpointer()->AddNewDetector(larSD);

    fLArLogical->SetSensitiveDetector(larSD);
}
