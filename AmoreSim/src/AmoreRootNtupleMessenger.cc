////////////////////////////////////////////////////////////////
// AmoreRootNtupleMessenger
////////////////////////////////////////////////////////////////

#include "AmoreSim/AmoreRootNtupleMessenger.hh"
#include "AmoreSim/AmoreRootNtuple.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <cstdlib> // for strtol
#include <fstream> // for file streams
#include <sstream>

AmoreRootNtupleMessenger::AmoreRootNtupleMessenger(AmoreRootNtuple *myntuple)
    : CupRootNtupleMessenger(myntuple), myNtuple(myntuple) {
    // the AmoreRootNtuple directory
    AmoreRootNtupleDir = new G4UIdirectory("/ntuple/");
    AmoreRootNtupleDir->SetGuidance("Control the detector geometry options.");

    CUTCmd = new G4UIcommand("/ntuple/recordWithCut", this);
    CUTCmd->SetGuidance("Select on/off of applying cut on recording of ntuples.");
    CUTCmd->AvailableForStates(G4State_PreInit);
    CUTCmd->SetParameter(new G4UIparameter("recordWithCut", 'b', true));

    PrimCmd = new G4UIcommand("/ntuple/recordPrimaries", this);
    PrimCmd->SetGuidance("Select on/off of recording primaries for neutron flux.");
    PrimCmd->AvailableForStates(G4State_PreInit);
    PrimCmd->SetParameter(new G4UIparameter("recordPrimaries", 'b', true));
}

AmoreRootNtupleMessenger::~AmoreRootNtupleMessenger() {
    delete CUTCmd;
    delete PrimCmd;

    delete AmoreRootNtupleDir;
}

void AmoreRootNtupleMessenger::SetNewValue(G4UIcommand *command, G4String newValues) {
    if (command == CUTCmd) {
        G4bool input = StoB(newValues);
        myNtuple->SetRecordCut(input);
    } else if (command == PrimCmd) {
        G4bool input = StoB(newValues);
        myNtuple->SetRecordPrim(input);
    } else {
        CupRootNtupleMessenger::SetNewValue(command, newValues);
    }
}

G4String AmoreRootNtupleMessenger::GetCurrentValue(G4UIcommand *command) {
    // CalDeviceCmd
    if (command == CUTCmd) {
        return BtoS(myNtuple->GetRecordCut());
    } else if (command == CUTCmd) {
        return BtoS(myNtuple->GetRecordPrim());
    } else { // invalid command
        return CupRootNtupleMessenger::GetCurrentValue(command);
    }
}
