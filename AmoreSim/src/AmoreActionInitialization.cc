#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1000

#include "AmoreSim/AmoreActionInitialization.hh"
#include "CupSim/CupPrimaryGeneratorAction.hh"
#include "CupSim/CupRunAction.hh"
#include "AmoreSim/AmoreSteppingAction.hh"
#include "AmoreSim/AmoreTrackingAction.hh"
#include "AmoreSim/AmoreEventAction.hh"

void AmoreActionInitialization::BuildForMaster() const {
    SetUserAction(new CupRunAction(fRecorders));
}

void AmoreActionInitialization::Build() const {
    auto p = new CupPrimaryGeneratorAction(fDetConstruction);
    SetUserAction(p);
    SetUserAction(new CupRunAction(fRecorders));
    SetUserAction(new AmoreEventAction(fRecorders));
    SetUserAction(new AmoreTrackingAction(fRecorders));
    SetUserAction(new AmoreSteppingAction(fRecorders, p));
}

#endif
