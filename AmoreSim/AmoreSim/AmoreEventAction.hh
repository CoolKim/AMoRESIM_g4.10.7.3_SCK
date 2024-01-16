// This file is part of the GenericLAND software library.
//
//
//  GenericLAND Simulation
//
//  Cup version by Glenn Horton-Smith December, 2004.
//  Based on earlier work by H. Ikeda, O. Tajima and G. Horton-Smith
//

#ifndef AmoreEventAction_h
#define AmoreEventAction_h 1

#include "G4UImessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "fstream"

#include "AmoreSim/AmoreRootNtuple.hh"
#include "CupSim/CupVEventAction.hh"

#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VHitsCollection.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

class G4Event;         // EJ
class CupRecorderBase; // EJ

class G4UIcmdWith3VectorAndUnit;

class AmoreEventAction : public CupVEventAction {
  public:
    // constructor, destructor
    //  AmoreEventAction();
    AmoreEventAction(AmoreRootNtuple *r = 0); // EJ
    ~AmoreEventAction();

    // overrides for G4UserEventAction methods
    virtual void BeginOfEventAction(const G4Event *);
    virtual void EndOfEventAction(const G4Event *);

    // overrides for G4UImessenger methods
    virtual void SetNewValue(G4UIcommand *command, G4String newValue);
    G4String GetCurrentValue(G4UIcommand *command);

    static const G4ThreeVector &GetPrimSkew() { return fgPrimPosSkew; }
    static void SetPrimSkew(const G4ThreeVector& a) { fgPrimPosSkew = a; }

  private: // EJ
    // Save the CupRecorderBase object to be called by the UserEventAction.
    AmoreRootNtuple *recorder; // EJ
    static G4ThreeVector fgPrimPosSkew;
    G4bool fPrimSkewEnable;
    G4UIcommand *fPrimSkewEnableCmd;
    G4UIdirectory *fPrimDirectory;

  protected:
};

#endif
