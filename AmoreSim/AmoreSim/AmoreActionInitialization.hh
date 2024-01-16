#ifndef __AmoreActionInitialization_h__
#define __AmoreActionInitialization_h__ 1

#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1000

#include "G4MTRunManager.hh"
#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class AmoreDetectorConstruction;
class AmoreRootNtuple;

class AmoreActionInitialization : public G4VUserActionInitialization {
  public:
    AmoreActionInitialization(AmoreRootNtuple *aRec, AmoreDetectorConstruction *aDet)
        : fDetConstruction(aDet), fRecorders(aRec){};
    virtual ~AmoreActionInitialization(){};

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    AmoreDetectorConstruction *fDetConstruction;

    AmoreRootNtuple *fRecorders;
};

#endif

#endif
