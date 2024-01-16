//
//  Cup version by Glenn Horton-Smith December, 2004.
//  Based on earlier work by H. Ikeda and G. Horton-Smith
//
//  It was modified by E.J.Jeon

#include <fstream>

#include "AmoreSim/AmoreEventAction.hh"
#include "AmoreSim/AmoreTrajectory.hh"
#include "CupSim/CupVEventAction.hh"

#include "MCObjs/EvtStep.hh"
#include "MCObjs/EvtTrack.hh"
#include "MCObjs/TStep.hh"
#include "MCObjs/TTrack.hh"

#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4PrimaryVertex.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VHitsCollection.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "CupSim/CupRecorderBase.hh"

G4ThreeVector AmoreEventAction::fgPrimPosSkew = G4ThreeVector();

AmoreEventAction::AmoreEventAction(AmoreRootNtuple *r)
    : CupVEventAction(r), recorder(r), fPrimSkewEnable(false), fPrimSkewEnableCmd(nullptr),
      fPrimDirectory(nullptr) {
    fPrimDirectory     = new G4UIdirectory("/event/primary");
    fPrimSkewEnableCmd = new G4UIcommand("/event/primary/enablePrimarySkew", this);
    fPrimSkewEnableCmd->SetGuidance("Select enable/disable of primary skewing.");
    fPrimSkewEnableCmd->SetParameter(new G4UIparameter("enable", 'b', false));
    fPrimSkewEnableCmd->AvailableForStates(G4State_Idle);
}

AmoreEventAction::~AmoreEventAction() {
    delete fPrimSkewEnableCmd;
    delete fPrimDirectory;
}

G4String AmoreEventAction::GetCurrentValue(G4UIcommand *nowCommand) {
    if (nowCommand == fPrimSkewEnableCmd) {
        return fPrimSkewEnable ? "true" : "false";
    } else
        return CupVEventAction::GetCurrentValue(nowCommand);
}

// EJ
void AmoreEventAction::SetNewValue(G4UIcommand *nowCommand, G4String newValue) {
    if (nowCommand == fPrimSkewEnableCmd) {
        fPrimSkewEnable = StoB(newValue);
    } else
        CupVEventAction::SetNewValue(nowCommand, newValue);
}

void AmoreEventAction::BeginOfEventAction(const G4Event *evt) {
    recorder->ClearET();
    if (fPrimSkewEnable) {
        G4PrimaryVertex *temp = evt->GetPrimaryVertex();
        for (int i = 0; i < evt->GetNumberOfPrimaryVertex(); i++) {
            G4ThreeVector nowvec = temp->GetPosition();
            nowvec += fgPrimPosSkew;
            temp->SetPosition(nowvec.x(), nowvec.y(), nowvec.z());
            temp = temp->GetNext();
        }
    }
    recorder->RecordBeginOfEvent(evt);
    CupVEventAction::BeginOfEventAction(evt);
}

void AmoreEventAction::EndOfEventAction(const G4Event *evt) {
    CupVEventAction::EndOfEventAction(evt);
}
