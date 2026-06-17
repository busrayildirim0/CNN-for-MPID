#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* genAction)
 : G4UImessenger(),
   fPrimaryGenerator(genAction),
   fGeneratorDir(nullptr),
   fSetModeCmd(nullptr),
   fSetFlavorCmd(nullptr),
   fSetProfileCmd(nullptr),
   fSetCosmicOverlayCmd(nullptr),
   fSetNeutrinoEnergyCmd(nullptr),
   fSetCosmicEnergyCmd(nullptr),
   fSetCosmicAngleCmd(nullptr)
{
    fGeneratorDir = new G4UIdirectory("/generator/");
    fGeneratorDir->SetGuidance("Primary generator control");

    fSetModeCmd = new G4UIcmdWithAString("/generator/setMode", this);
    fSetModeCmd->SetGuidance("Set generator mode: neutrino, cosmic, or test");
    fSetModeCmd->SetParameterName("mode", false);
    fSetModeCmd->SetCandidates("neutrino cosmic test");
    fSetModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetFlavorCmd = new G4UIcmdWithAString("/generator/setFlavor", this);
    fSetFlavorCmd->SetGuidance("Set neutrino flavor/category: numu, nue, nutau, nc, or all");
    fSetFlavorCmd->SetParameterName("flavor", false);
    fSetFlavorCmd->SetCandidates("numu nue nutau nc all");
    fSetFlavorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetProfileCmd = new G4UIcmdWithAString("/generator/setProfile", this);
    fSetProfileCmd->SetGuidance("Set sampling profile: realistic or ml");
    fSetProfileCmd->SetParameterName("profile", false);
    fSetProfileCmd->SetCandidates("realistic ml");
    fSetProfileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetCosmicOverlayCmd = new G4UIcmdWithAString("/generator/setCosmicOverlay", this);
    fSetCosmicOverlayCmd->SetGuidance("Enable or disable cosmic overlay in neutrino mode: on or off");
    fSetCosmicOverlayCmd->SetParameterName("overlay", false);
    fSetCosmicOverlayCmd->SetCandidates("on off");
    fSetCosmicOverlayCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetNeutrinoEnergyCmd = new G4UIcommand("/generator/setNeutrinoEnergyRange", this);
    fSetNeutrinoEnergyCmd->SetGuidance("Set neutrino energy range (min max unit)");

    G4UIparameter* minEnuParam = new G4UIparameter("minE", 'd', false);
    minEnuParam->SetDefaultValue(0.2);
    fSetNeutrinoEnergyCmd->SetParameter(minEnuParam);

    G4UIparameter* maxEnuParam = new G4UIparameter("maxE", 'd', false);
    maxEnuParam->SetDefaultValue(3.0);
    fSetNeutrinoEnergyCmd->SetParameter(maxEnuParam);

    G4UIparameter* energyUnitParam = new G4UIparameter("unit", 's', false);
    energyUnitParam->SetDefaultValue("GeV");
    fSetNeutrinoEnergyCmd->SetParameter(energyUnitParam);

    fSetNeutrinoEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetCosmicEnergyCmd = new G4UIcommand("/generator/setCosmicEnergyRange", this);
    fSetCosmicEnergyCmd->SetGuidance("Set cosmic ray energy range (min max unit)");

    G4UIparameter* minCosmicEParam = new G4UIparameter("minE", 'd', false);
    minCosmicEParam->SetDefaultValue(0.1);
    fSetCosmicEnergyCmd->SetParameter(minCosmicEParam);

    G4UIparameter* maxCosmicEParam = new G4UIparameter("maxE", 'd', false);
    maxCosmicEParam->SetDefaultValue(100.0);
    fSetCosmicEnergyCmd->SetParameter(maxCosmicEParam);

    G4UIparameter* cosmicUnitParam = new G4UIparameter("unit", 's', false);
    cosmicUnitParam->SetDefaultValue("GeV");
    fSetCosmicEnergyCmd->SetParameter(cosmicUnitParam);

    fSetCosmicEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSetCosmicAngleCmd = new G4UIcommand("/generator/setCosmicAngleRange", this);
    fSetCosmicAngleCmd->SetGuidance("Set cosmic ray angle range from zenith (min max unit)");

    G4UIparameter* minAngleParam = new G4UIparameter("minTheta", 'd', false);
    minAngleParam->SetDefaultValue(0.0);
    fSetCosmicAngleCmd->SetParameter(minAngleParam);

    G4UIparameter* maxAngleParam = new G4UIparameter("maxTheta", 'd', false);
    maxAngleParam->SetDefaultValue(85.0);
    fSetCosmicAngleCmd->SetParameter(maxAngleParam);

    G4UIparameter* angleUnitParam = new G4UIparameter("unit", 's', false);
    angleUnitParam->SetDefaultValue("deg");
    fSetCosmicAngleCmd->SetParameter(angleUnitParam);

    fSetCosmicAngleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete fSetModeCmd;
    delete fSetFlavorCmd;
    delete fSetProfileCmd;
    delete fSetCosmicOverlayCmd;
    delete fSetNeutrinoEnergyCmd;
    delete fSetCosmicEnergyCmd;
    delete fSetCosmicAngleCmd;
    delete fGeneratorDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fSetModeCmd) {
        if (newValue == "neutrino") {
            fPrimaryGenerator->SetGeneratorMode(PrimaryGeneratorAction::kNeutrinoMode);
        }
        else if (newValue == "cosmic") {
            fPrimaryGenerator->SetGeneratorMode(PrimaryGeneratorAction::kCosmicRayMode);
        }
        else if (newValue == "test") {
            fPrimaryGenerator->SetGeneratorMode(PrimaryGeneratorAction::kTestMode);
        }
    }
    else if (command == fSetFlavorCmd) {
        if (newValue == "numu") {
            fPrimaryGenerator->SetNeutrinoFlavor(PrimaryGeneratorAction::kNuMu);
        } else if (newValue == "nue") {
            fPrimaryGenerator->SetNeutrinoFlavor(PrimaryGeneratorAction::kNuE);
        } else if (newValue == "nutau") {
            fPrimaryGenerator->SetNeutrinoFlavor(PrimaryGeneratorAction::kNuTau);
        } else if (newValue == "nc") {
            fPrimaryGenerator->SetNeutrinoFlavor(PrimaryGeneratorAction::kNC);
        } else if (newValue == "all") {
            fPrimaryGenerator->SetNeutrinoFlavor(PrimaryGeneratorAction::kAllFlavors);
        }
    }
    else if (command == fSetProfileCmd) {
        if (newValue == "realistic") {
            fPrimaryGenerator->SetSamplingProfile(PrimaryGeneratorAction::kRealisticProfile);
        } else if (newValue == "ml") {
            fPrimaryGenerator->SetSamplingProfile(PrimaryGeneratorAction::kMLBalancedProfile);
        }
    }
    else if (command == fSetCosmicOverlayCmd) {
        fPrimaryGenerator->SetEnableCosmicOverlay(newValue == "on");
    }
    else if (command == fSetNeutrinoEnergyCmd) {
        G4double minE, maxE;
        G4String unit;
        std::istringstream is(newValue);
        is >> minE >> maxE >> unit;

        G4double unitValue = G4UIcommand::ValueOf(unit.c_str());
        fPrimaryGenerator->SetNeutrinoEnergyRange(minE * unitValue, maxE * unitValue);
    }
    else if (command == fSetCosmicEnergyCmd) {
        G4double minE, maxE;
        G4String unit;
        std::istringstream is(newValue);
        is >> minE >> maxE >> unit;

        G4double unitValue = G4UIcommand::ValueOf(unit.c_str());
        fPrimaryGenerator->SetCosmicEnergyRange(minE * unitValue, maxE * unitValue);
    }
    else if (command == fSetCosmicAngleCmd) {
        G4double minTheta, maxTheta;
        G4String unit;
        std::istringstream is(newValue);
        is >> minTheta >> maxTheta >> unit;

        G4double unitValue = G4UIcommand::ValueOf(unit.c_str());
        fPrimaryGenerator->SetCosmicAngleRange(minTheta * unitValue, maxTheta * unitValue);
    }
}
