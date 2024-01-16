#ifndef AmoreScintSD_h
#define AmoreScintSD_h 1

#include "CupSim/CupScintSD.hh"

class G4TouchableHistory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AmoreScintSD : public CupScintSD {
  public:
    AmoreScintSD(G4String name, int max_tgs = 1000);
    ~AmoreScintSD();

  public:
    virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
};

#endif
