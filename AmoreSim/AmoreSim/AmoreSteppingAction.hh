// This file is part of the GenericLAND software library.
// $Id: AmoreSteppingAction.hh,v 1.1.1.1 2016/10/31 08:41:44 ejjeon Exp $
//
#ifndef __AmoreSteppingAction_H__
#define __AmoreSteppingAction_H__ 1

#include "AmoreSim/AmoreRootNtuple.hh"
#include "CupSim/CupPrimaryGeneratorAction.hh"
#include "CupSim/CupSteppingAction.hh"
#include "G4UserSteppingAction.hh"
#include "globals.hh"

class CupPrimaryGeneratorAction;
class CupRecorderBase; // EJ

class AmoreSteppingAction : public CupSteppingAction {
  public:
    AmoreSteppingAction(AmoreRootNtuple *r);
    AmoreSteppingAction(AmoreRootNtuple *r, CupPrimaryGeneratorAction *p);
    virtual ~AmoreSteppingAction(){};

    // void UserSteppingAction(const G4Step* aStep);
    virtual void UserSteppingAction(const G4Step *aStep);

  private:
};

#endif
