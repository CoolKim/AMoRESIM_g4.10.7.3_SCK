//
#ifndef AmoreModuleSD_h
#define AmoreModuleSD_h 1
#include <set>

#include "AmoreSim/AmoreDetectorStaticInfo.hh"
#include "AmoreSim/AmoreModuleHit.hh"
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class AmoreModuleSD : public G4VSensitiveDetector {
    // EJ:
  protected:
  public:
    using CrystalModuleInfo = AmoreDetectorStaticInfo::CrystalModuleInfo;
    AmoreModuleSD(G4String aName, std::set<AmoreModuleSDInfo> &aModuleSDInfoList);
    virtual ~AmoreModuleSD();

    virtual void Initialize(G4HCofThisEvent *HCE);
    virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *aTHist);
    virtual void EndOfEvent(G4HCofThisEvent *HCE);

  private:
    AmoreModuleHitsCollection *fHitsColl;
    G4int fHitCollID;
    G4int fModulesNumber;
    std::set<AmoreModuleSDInfo> &fModuleSDInfoList;
};

#endif
