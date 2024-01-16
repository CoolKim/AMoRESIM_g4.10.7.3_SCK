//
//  Created by E.J.Jeon 2015/06/16
//

#include "AmoreSim/AmoreModuleSD.hh"
#include "AmoreSim/AmoreModuleHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Version.hh"
#include "G4ios.hh"

#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

AmoreModuleSD::AmoreModuleSD(G4String name, std::set<AmoreModuleSDInfo> &aModuleSDInfoList)
    : G4VSensitiveDetector(name), fModuleSDInfoList(aModuleSDInfoList) {
    G4String HCname;
    collectionName.insert(HCname = "AmoreModuleSDColl");
    fModulesNumber = fModuleSDInfoList.size();
    fHitCollID     = -1;
}

AmoreModuleSD::~AmoreModuleSD() { ; }

void AmoreModuleSD::Initialize(G4HCofThisEvent *HCE) {
    fHitsColl = new AmoreModuleHitsCollection(SensitiveDetectorName, collectionName[0]);
    if (fHitCollID < 0) {
        fHitCollID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsColl);
    }
    HCE->AddHitsCollection(fHitCollID, fHitsColl);

    auto hitCollVector = fHitsColl->GetVector();
    hitCollVector->resize(fModulesNumber);
    for (auto &nowSDInfo : fModuleSDInfoList) {
        AmoreModuleHit *aHit = new AmoreModuleHit(&nowSDInfo);

        (*hitCollVector)[aHit->GetModuleID()] = aHit;
    }
}

G4bool AmoreModuleSD::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
    static CrystalModuleInfo aDummy;
    G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();

    G4int nowMotherLevel = 1;
    G4int envelopeCopyNo = 1;
		G4int nowCopyNo;
    G4double energyDeposit, qEnergyDeposit;
    G4ParticleDefinition *nowParticle;
    G4String particleName;

    G4VPhysicalVolume *theMotherPhysical, *nowPhysical;
    G4LogicalVolume *nowLogical;

    nowParticle   = aStep->GetTrack()->GetDefinition();
    energyDeposit = aStep->GetTotalEnergyDeposit();
    particleName  = nowParticle->GetParticleName();
    if (particleName == "opticalphoton" || energyDeposit == 0.) return true;

#if G4VERSION_NUMBER <= 1020
    qEnergyDeposit = emSaturation->VisibleEnergyDeposition(aStep);
#else
    qEnergyDeposit = emSaturation->VisibleEnergyDepositionAtAStep(aStep); // for geant4.10.4.2
#endif


    G4StepPoint *preStepPoint      = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();

    nowPhysical = theTouchable->GetVolume(0);
    nowLogical  = nowPhysical->GetLogicalVolume();
    nowCopyNo   = nowPhysical->GetCopyNo();

    theMotherPhysical = theTouchable->GetVolume(nowMotherLevel++);

    while (theMotherPhysical != nullptr) {
        if (theMotherPhysical->GetLogicalVolume()->IsRootRegion()) {
            if (theMotherPhysical->GetMotherLogical() == nullptr) { // Reached the end of the world
                G4Exception(__PRETTY_FUNCTION__, "MDSD_REGION_FAIL", FatalException,
                            "Finding a root region for SD has been failed.");
            }
            envelopeCopyNo = theMotherPhysical->GetCopyNo();
            break;
        }

        theMotherPhysical = theTouchable->GetVolume(nowMotherLevel++);
    }
    AmoreModuleHit *aHit = (*fHitsColl)[envelopeCopyNo];

    if (nowLogical == aHit->GetCrystalLogicalVolume()) {
        //aHit->SetCrystalEdep(energyDeposit);
        //aHit->SetCrystalQEdep(qEnergyDeposit);
        aHit->AddCrystalEdep(energyDeposit);   // JW modified
        aHit->AddCrystalQEdep(qEnergyDeposit); // JW modified
    } else if (nowLogical == aHit->GetGeWaferLogicalVolume()) {
        aHit->AddGeWaferEdep(energyDeposit);
        aHit->AddGeWaferQEdep(qEnergyDeposit);
    } else if (nowLogical == aHit->GetCrystalGoldFilmLogicalVolume()) {
        aHit->AddCrystalGoldFilmEdep(0, energyDeposit);
    } else if (nowLogical == aHit->GetGeWaferGoldFilmLogicalVolume()) {
        aHit->AddGeWaferGoldFilmEdep(nowCopyNo, energyDeposit);
    } else {
        G4Exception(__PRETTY_FUNCTION__, "MDSD_LV_NOTMATCH", G4ExceptionSeverity::JustWarning,
                    "SD has hitted but there is no LV which fits to this step.");
    }

    return true;
}

void AmoreModuleSD::EndOfEvent(G4HCofThisEvent * /*HCE*/) { ; }
