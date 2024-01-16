//
//  Created by E.J.Jeon 2015/06/16
//

#include "AmoreSim/AmoreModuleHit.hh"
#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Version.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#if G4VERSION_NUMBER <= 999
G4Allocator<AmoreModuleHit> *AmoreModuleHitAllocator = nullptr;
#else
G4ThreadLocal G4Allocator<AmoreModuleHit> *AmoreModuleHitAllocator = nullptr;
#endif

G4int AmoreModuleHit::fgCrystalGoldFilmNum = 1;
G4int AmoreModuleHit::fgGeWaferGoldFilmNum = 3;

AmoreModuleHit::AmoreModuleHit(const AmoreModuleSDInfo *aMSDInfo, G4int aCrystalGoldFilmNum,
                               G4int aGeWaferGoldFilmNum)
    : G4VHit(), fModuleSDInfo(aMSDInfo) {
    fEdepOnCrystal         = 0.;
    fQuenchedEdepOnCrystal = 0.;
    fEdepOnGeWafer         = 0.;
    fQuenchedEdepOnGeWafer = 0.;
    fEdepOnGeWaferGoldFilm.resize(aGeWaferGoldFilmNum, 0.);
    fEdepOnCrystalGoldFilm.resize(aCrystalGoldFilmNum, 0.);
}

AmoreModuleHit::~AmoreModuleHit() { ; }

AmoreModuleHit::AmoreModuleHit(const AmoreModuleHit &right)
    : G4VHit(), fModuleSDInfo(right.fModuleSDInfo) {
    fEdepOnCrystal         = right.fEdepOnCrystal;
    fQuenchedEdepOnCrystal = right.fQuenchedEdepOnCrystal;
    fEdepOnGeWafer         = right.fEdepOnGeWafer;
    fQuenchedEdepOnGeWafer = right.fQuenchedEdepOnGeWafer;
    fEdepOnGeWaferGoldFilm = right.fEdepOnGeWaferGoldFilm;
    fEdepOnCrystalGoldFilm = right.fEdepOnCrystalGoldFilm;
}

const AmoreModuleHit &AmoreModuleHit::operator=(const AmoreModuleHit &right) {
    fModuleSDInfo          = right.fModuleSDInfo;
    fEdepOnCrystal         = right.fEdepOnCrystal;
    fQuenchedEdepOnCrystal = right.fQuenchedEdepOnCrystal;
    fEdepOnGeWafer         = right.fEdepOnGeWafer;
    fQuenchedEdepOnGeWafer = right.fQuenchedEdepOnGeWafer;
    fEdepOnGeWaferGoldFilm = right.fEdepOnGeWaferGoldFilm;
    fEdepOnCrystalGoldFilm = right.fEdepOnCrystalGoldFilm;

    return *this;
}

int AmoreModuleHit::operator==(const AmoreModuleHit &right) const {
    return (fModuleSDInfo->fModuleID == right.fModuleSDInfo->fModuleID);
}

const std::map<G4String, G4AttDef> *AmoreModuleHit::GetAttDefs() const {
    G4bool isNew;
    std::map<G4String, G4AttDef> *store =
        G4AttDefStore::GetInstance("CupSim/AmoreModuleHit", isNew);
    if (isNew) {
        G4String nowTempStr;
        nowTempStr           = "ModuleID";
        (*store)[nowTempStr] = G4AttDef(nowTempStr, "Module ID", "Physics", "", "G4int");

        nowTempStr = G4String("PositionID");
        (*store)[nowTempStr] =
            G4AttDef(nowTempStr, "Position ID of module", "Physics", "", "G4int[2]");

        nowTempStr = G4String("EdepOnCrystal");
        (*store)[nowTempStr] =
            G4AttDef(nowTempStr, "Energy deposit on crystal", "Physics", "G4BestUnit", "G4double");

        nowTempStr = G4String("QEdepOnCrystal");
        (*store)[nowTempStr] =
            G4AttDef(nowTempStr, "Energy deposit on crystal including quenching effects", "Physics",
                     "G4BestUnit", "G4double");

        nowTempStr           = G4String("EdepOnGeWafer");
        (*store)[nowTempStr] = G4AttDef(nowTempStr, "Energy deposit on germanium wafer", "Physics",
                                        "G4BestUnit", "G4double");

        nowTempStr = G4String("QEdepOnGeWafer");
        (*store)[nowTempStr] =
            G4AttDef(nowTempStr, "Energy deposit on germanium wafer including quenching effects",
                     "Physics", "G4BestUnit", "G4double");

        nowTempStr = G4String("CrystalLV");
        (*store)[nowTempStr] =
            G4AttDef(nowTempStr, "Logical volume of crystal", "Physics", "", "G4String");

        nowTempStr = G4String("GeWaferLV");
        (*store)[nowTempStr] =
            G4AttDef(nowTempStr, "Logical volume of germanium wafer", "Physics", "", "G4String");
    }
    return store;
}

std::vector<G4AttValue> *AmoreModuleHit::CreateAttValues() const {
    std::vector<G4AttValue> *values = new std::vector<G4AttValue>;

    G4String nowTempStr;
    nowTempStr = "ModuleID";
    values->push_back(G4AttValue(nowTempStr, "AmoreModuleHit", ""));

    nowTempStr = G4String("PositionID");
    values->push_back(G4AttValue(nowTempStr, G4UIcommand::ConvertToString(fModuleSDInfo->fModuleID), ""));

    nowTempStr = G4String("EdepOnCrystal");
    values->push_back(G4AttValue(nowTempStr, G4BestUnit(fEdepOnCrystal, "Energy"), ""));

    nowTempStr = G4String("QEdepOnCrystal");
    values->push_back(G4AttValue(nowTempStr, G4BestUnit(fQuenchedEdepOnCrystal, "Energy"), ""));

    nowTempStr = G4String("EdepOnGeWafer");
    values->push_back(G4AttValue(nowTempStr, G4BestUnit(fEdepOnGeWafer, "Energy"), ""));

    nowTempStr = G4String("QEdepOnGeWafer");
    values->push_back(G4AttValue(nowTempStr, G4BestUnit(fQuenchedEdepOnGeWafer, "Energy"), ""));

    return values;
}

void AmoreModuleHit::Print() {
    G4cout << "  Cell[" << fModuleSDInfo->fModuleID << "]: Edep on crystal = " << fEdepOnCrystal / MeV << " (MeV)"
           << " | Edep on germanium wafer = " << fEdepOnGeWafer / MeV << " (MeV)" << G4endl;
}
