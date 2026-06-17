#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
    : G4VModularPhysicsList(),
      fVerboseLevel(0)
{
    SetVerboseLevel(fVerboseLevel);
    defaultCutValue = 0.1*mm;

    RegisterPhysics(new G4EmStandardPhysics_option4(fVerboseLevel));
    RegisterPhysics(new G4EmExtraPhysics(fVerboseLevel));
    RegisterPhysics(new G4DecayPhysics(fVerboseLevel));
    RegisterPhysics(new G4RadioactiveDecayPhysics(fVerboseLevel));

    RegisterPhysics(new G4HadronElasticPhysics(fVerboseLevel));
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT(fVerboseLevel));
    RegisterPhysics(new G4IonPhysics(fVerboseLevel));
    RegisterPhysics(new G4StoppingPhysics(fVerboseLevel));
}

PhysicsList::~PhysicsList()
{}

void PhysicsList::SetCuts()
{
    G4VUserPhysicsList::SetCuts();

    SetCutValue(0.1*mm, "gamma");
    SetCutValue(0.1*mm, "e-");
    SetCutValue(0.1*mm, "e+");
    SetCutValue(0.1*mm, "proton");
}
