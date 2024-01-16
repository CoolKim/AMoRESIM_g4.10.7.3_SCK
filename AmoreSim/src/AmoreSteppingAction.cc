//
//  Concrete implementation of G4UserSteppingAction
//
//  Current uses:
//    * Measure inter-step CPU time, broken down by process and particle type
//
//  Anticipated uses:
//    * Find PMT _fast_ when entering outer buffer
//
//  Author: Glenn Horton-Smith, April 7, 2000

#include "AmoreSim/AmoreSteppingAction.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CupSim/CupPrimaryGeneratorAction.hh"
#include "CupSim/CupRecorderBase.hh" // EJ
#include "CupSim/CupScintillation.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"
#include "G4VisExtent.hh"
#include "G4ios.hh"
#include "globals.hh"

AmoreSteppingAction::AmoreSteppingAction(AmoreRootNtuple *r) : CupSteppingAction(r){};
AmoreSteppingAction::AmoreSteppingAction(AmoreRootNtuple *r, CupPrimaryGeneratorAction *p)
    : CupSteppingAction(r, p){};

void AmoreSteppingAction::UserSteppingAction(const G4Step *aStep) {
    CupSteppingAction::UserSteppingAction(aStep);
}
