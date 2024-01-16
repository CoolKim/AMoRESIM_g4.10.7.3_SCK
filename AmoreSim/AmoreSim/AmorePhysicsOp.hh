#ifndef AmorePhysicsOp_h
#define AmorePhysicsOp_h 1

#include "G4VPhysicsConstructor.hh"
//#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AmorePhysicsOp : public G4VPhysicsConstructor {
  public:
    AmorePhysicsOp(const G4String &name = "standardOp");
    ~AmorePhysicsOp();

  public:
    // This method is dummy for physics
    virtual void ConstructParticle(){};

    virtual void ConstructProcess();

  private:
    G4int OpVerbLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
