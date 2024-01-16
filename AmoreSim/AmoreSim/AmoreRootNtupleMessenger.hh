//
// AmoreRootNtupleMessenger.hh
//
#ifndef __AmoreRootNtupleMessenger_hh__
#define __AmoreRootNtupleMessenger_hh__ 1

#include "CupSim/CupRootNtupleMessenger.hh"
#include "G4UImessenger.hh"

class G4UIcommand;
class AmoreRootNtuple;

class AmoreRootNtupleMessenger : public CupRootNtupleMessenger {
  public:
    AmoreRootNtupleMessenger(AmoreRootNtuple *myNtuple);
    ~AmoreRootNtupleMessenger();

    void SetNewValue(G4UIcommand *command, G4String newValues);
    G4String GetCurrentValue(G4UIcommand *command);

  private:
    AmoreRootNtuple *myNtuple;

    G4UIdirectory *AmoreRootNtupleDir;
    G4UIcommand *CUTCmd;

    G4UIcommand *PrimCmd;
};

#endif
