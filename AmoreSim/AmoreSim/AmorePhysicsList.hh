#ifndef AmorePhysicsList_h
#define AmorePhysicsList_h 1

#include "CupSim/CupPhysicsList.hh"

#include "G4VUserPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AmorePhysicsList : public CupPhysicsList {
  public:
    AmorePhysicsList();
    ~AmorePhysicsList();
    virtual void AddPhysicsList(const G4String &name);
    virtual void ConstructProcess();

  private:
    G4String fOpName;
    G4String fEMName;
};

#endif
