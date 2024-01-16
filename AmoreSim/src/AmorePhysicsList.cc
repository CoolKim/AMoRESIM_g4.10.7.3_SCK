#include <iomanip>

#include "AmoreSim/AmorePhysicsList.hh"
#include "AmoreSim/AmorePhysicsOp.hh"
#include "AmoreSim/PhysListEmStandardNR.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics.hh"

// Constructor /////////////////////////////////////////////////////////////
AmorePhysicsList::AmorePhysicsList() : CupPhysicsList() {}

// Destructor //////////////////////////////////////////////////////////////
AmorePhysicsList::~AmorePhysicsList() {}

void AmorePhysicsList::ConstructProcess() {

    AddTransportation();

    AddParameterisation();

    if (fEMName == "livermore") {
        auto a = new G4EmLivermorePhysics;
        a->ConstructProcess();
    } else if (fEMName == "emstandardNR") {
        auto a = new PhysListEmStandardNR();
        a->ConstructProcess();
    } else {
        ConstructEM();
        G4cout << "EM Physics is ConstructEM(): default!" << G4endl;
    }
    if (fOpName == "amorephysicsOp") {
        auto a = new AmorePhysicsOp();
        a->ConstructProcess();
    } else {
        ConstructOp();
	G4cout << "Op Physics is ConstructOp(): default!" << G4endl;
    }

    ConstructHad();

    ConstructGeneral();
}

// Add Physics List ////////////////////////////////////////////////////////
void AmorePhysicsList::AddPhysicsList(const G4String &name) {
    G4cout << "\n>>>   AmorePhysicsList::AddPhysicsList: <<<" << name << ">>> \n" << G4endl;

    // EM physics
    if (name == "livermore") {
        G4cout << "Physics : EmLivermore is selected \n";
        fEMName = name;
    } else if (name == "emstandardNR") {
        G4cout << "Physics : EmStandardNR is selected \n";
        fEMName = name;
    // Op physics
    } else if (name == "amorephysicsOp") {
        G4cout << "Physics : AmorePhysicsOp is selected \n";
        fOpName = name;
    } else {
        G4cout << "AmorePhysicsList::AddPhysicsList: <" << name << ">"
               << " is not defined" << G4endl;
    }
}
