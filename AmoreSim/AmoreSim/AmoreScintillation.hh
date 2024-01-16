#ifndef AmoreScintillation_h
#define AmoreScintillation_h 1

#include "CupSim/CupScintillation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AmoreScintillation : public CupScintillation {
  public:
    AmoreScintillation(const G4String &processName = "Scintillation",
                       G4ProcessType type          = fElectromagnetic);
    ~AmoreScintillation();

  public:
    G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);

    static G4double GetTotEdepQuenched() { return TotalEnergyDepositQuenched; }

  private:
    static G4double TotalEnergyDepositQuenched;
};

#endif
