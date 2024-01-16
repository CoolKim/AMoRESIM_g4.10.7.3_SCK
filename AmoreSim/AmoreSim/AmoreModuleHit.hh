// AmoreModuleHit class
// A class for storing hit information of crystal detector module of AMoRE
// Author: BaseHardware (a.k.a. basehw) ; 19/03/2019
#ifndef __AMOREMODULEHIT__h_
#define __AMOREMODULEHIT__h_ 1

#include <vector>

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4VHit.hh"
#include "G4Version.hh"

#define MACRO_IS_IN_RANGE_OF(A, B, C) (((A) <= (B)) && ((B) < (C)))

class G4AttDef;
class G4AttValue;

using namespace CLHEP;

struct AmoreModuleSDInfo {
    G4String fModuleName;
    G4int fModuleID;
    G4int fCrystalPosIdx[2];
    G4VPhysicalVolume *fModulePV;
    G4LogicalVolume *fCrystalLV;
    G4LogicalVolume *fGeWaferLV;
    G4LogicalVolume *fCrystalGoldFilmLV;
    G4LogicalVolume *fGeWaferGoldFilmLV;
    // Comparer for fast searching
    inline G4bool operator<(const AmoreModuleSDInfo &rhs) const {
        return fModulePV < rhs.fModulePV;
    }
    inline G4bool operator==(const AmoreModuleSDInfo &rhs) const {
        return fModulePV == rhs.fModulePV;
    }
};

class AmoreModuleHit : public G4VHit {
  public:
    AmoreModuleHit(const AmoreModuleSDInfo *aMSDInfo,
                   G4int aCrystalGoldFilmNum = fgCrystalGoldFilmNum,
                   G4int aGeWaferGoldFilmNum = fgGeWaferGoldFilmNum);
    virtual ~AmoreModuleHit();
    AmoreModuleHit(const AmoreModuleHit &right);
    const AmoreModuleHit &operator=(const AmoreModuleHit &right);
    int operator==(const AmoreModuleHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    virtual const std::map<G4String, G4AttDef> *GetAttDefs() const;
    virtual std::vector<G4AttValue> *CreateAttValues() const;
    virtual void Print();

    inline void SetCrystalEdep(G4double aE) { fEdepOnCrystal = aE; }
    inline void SetCrystalQEdep(G4double aE) { fQuenchedEdepOnCrystal = aE; }
    inline void SetGeWaferEdep(G4double aE) { fEdepOnGeWafer = aE; }
    inline void SetGeWaferQEdep(G4double aE) { fQuenchedEdepOnGeWafer = aE; }

    inline G4bool SetGeWaferGoldFilmEdep(G4int aIdx, G4double aE);
    inline G4bool SetCrystalGoldFilmEdep(G4int aIdx, G4double aE);

    inline G4double GetCrystalEdep() const { return fEdepOnCrystal; }
    inline G4double GetCrystalQEdep() const { return fQuenchedEdepOnCrystal; }
    inline G4double GetGeWaferEdep() const { return fEdepOnGeWafer; }
    inline G4double GetGeWaferQEdep() const { return fQuenchedEdepOnGeWafer; }
    inline G4double GetGeWaferGoldFilmEdep(G4int) const;
    inline G4double GetCrystalGoldFilmEdep(G4int) const;

    inline void ResetCrystalEdep() { SetCrystalEdep(0); }
    inline void ResetCrystalQEdep() { SetCrystalQEdep(0); }
    inline void ResetGeWaferEdep() { SetGeWaferEdep(0); }
    inline void ResetGeWaferQEdep() { SetGeWaferQEdep(0); }
    inline void ResetGeWaferGoldFilmEdepAll();
    inline void ResetCrystalGoldFilmEdepAll();
    inline G4bool ResetGeWaferGoldFilmEdep(G4int aIdx) { return SetGeWaferGoldFilmEdep(aIdx, 0); };
    inline G4bool ResetCrystalGoldFilmEdep(G4int aIdx) { return SetCrystalGoldFilmEdep(aIdx, 0); };

    inline void AddCrystalEdep(G4double aE) { fEdepOnCrystal += aE; }
    inline void AddCrystalQEdep(G4double aE) { fQuenchedEdepOnCrystal += aE; }
    inline void AddGeWaferEdep(G4double aE) { fEdepOnGeWafer += aE; }
    inline void AddGeWaferQEdep(G4double aE) { fQuenchedEdepOnGeWafer += aE; }
    inline G4bool AddGeWaferGoldFilmEdep(G4int aIdx, G4double aE);
    inline G4bool AddCrystalGoldFilmEdep(G4int aIdx, G4double aE);

    inline void SetModuleSDInfo(const AmoreModuleSDInfo *aMSDInfo) { fModuleSDInfo = aMSDInfo; }
    inline const AmoreModuleSDInfo *GetModuleSDInfo() const { return fModuleSDInfo; }

    inline const G4LogicalVolume *GetCrystalLogicalVolume() const {
        return fModuleSDInfo->fCrystalLV;
    }
    inline const G4LogicalVolume *GetGeWaferLogicalVolume() const {
        return fModuleSDInfo->fGeWaferLV;
    }
    inline const G4LogicalVolume *GetCrystalGoldFilmLogicalVolume() const {
        return fModuleSDInfo->fCrystalGoldFilmLV;
    }
    inline const G4LogicalVolume *GetGeWaferGoldFilmLogicalVolume() const {
        return fModuleSDInfo->fGeWaferGoldFilmLV;
    }

    inline G4String GetModuleName() const { return fModuleSDInfo->fModuleName; }

    inline G4int GetModuleID() const { return fModuleSDInfo->fModuleID; }
    inline G4int GetXPositionIdx() const { return fModuleSDInfo->fCrystalPosIdx[0]; }
    inline G4int GetYPositionIdx() const { return fModuleSDInfo->fCrystalPosIdx[1]; }
    inline void GetPositionIdx(G4double &aXIdx, G4double &aYIdx) const {
        aXIdx = fModuleSDInfo->fCrystalPosIdx[0];
        aYIdx = fModuleSDInfo->fCrystalPosIdx[1];
    }

  private:
    G4double fEdepOnCrystal;
    G4double fQuenchedEdepOnCrystal;
    G4double fEdepOnGeWafer;
    G4double fQuenchedEdepOnGeWafer;
    std::vector<G4double> fEdepOnGeWaferGoldFilm;
    std::vector<G4double> fEdepOnCrystalGoldFilm;

    const AmoreModuleSDInfo *fModuleSDInfo;

  protected:
    static G4int fgGeWaferGoldFilmNum;
    static G4int fgCrystalGoldFilmNum;

  public:
    static inline void SetGeWaferGoldFilmNum(G4int aNum);
    static inline void SetCrystalGoldFilmNum(G4int aNum);
    static inline G4int GetGeWaferGoldFilmNum() { return fgGeWaferGoldFilmNum; }
    static inline G4int GetCrystalGoldFilmNum() { return fgCrystalGoldFilmNum; }
};

typedef G4THitsCollection<AmoreModuleHit> AmoreModuleHitsCollection;

#if G4VERSION_NUMBER <= 999
extern G4Allocator<AmoreModuleHit> *AmoreModuleHitAllocator;
#else
extern G4ThreadLocal G4Allocator<AmoreModuleHit> *AmoreModuleHitAllocator;
#endif

inline void *AmoreModuleHit::operator new(size_t) {
    void *aHit;
    if (!AmoreModuleHitAllocator) AmoreModuleHitAllocator = new G4Allocator<AmoreModuleHit>;
    aHit = (void *)AmoreModuleHitAllocator->MallocSingle();
    return aHit;
}

inline void AmoreModuleHit::operator delete(void *aHit) {
    AmoreModuleHitAllocator->FreeSingle((AmoreModuleHit *)aHit);
}

inline void AmoreModuleHit::SetCrystalGoldFilmNum(G4int aNum) {
    if (aNum >= 1)
        fgCrystalGoldFilmNum = aNum;
    else
        fgCrystalGoldFilmNum = 1;
}

inline void AmoreModuleHit::ResetCrystalGoldFilmEdepAll() {
    for (auto &nowVal : fEdepOnCrystalGoldFilm)
        nowVal = 0;
}

inline G4bool AmoreModuleHit::SetCrystalGoldFilmEdep(G4int aIdx, G4double aE) {
    if (fEdepOnCrystalGoldFilm.size() > 0 &&
        MACRO_IS_IN_RANGE_OF(0, aIdx, (G4int)fEdepOnCrystalGoldFilm.size())) {
        fEdepOnCrystalGoldFilm[aIdx] = aE;
        return true;
    } else
        return false;
}

inline G4bool AmoreModuleHit::AddCrystalGoldFilmEdep(G4int aIdx, G4double aE) {
    if (fEdepOnCrystalGoldFilm.size() > 0 &&
        MACRO_IS_IN_RANGE_OF(0, aIdx, (G4int)fEdepOnCrystalGoldFilm.size())) {
        fEdepOnCrystalGoldFilm[aIdx] += aE;
        return true;
    } else
        return false;
}

inline void AmoreModuleHit::SetGeWaferGoldFilmNum(G4int aNum) {
    if (aNum >= 1)
        fgGeWaferGoldFilmNum = aNum;
    else
        fgGeWaferGoldFilmNum = 1;
}

inline void AmoreModuleHit::ResetGeWaferGoldFilmEdepAll() {
    for (auto &nowVal : fEdepOnGeWaferGoldFilm)
        nowVal = 0;
}

inline G4bool AmoreModuleHit::SetGeWaferGoldFilmEdep(G4int aIdx, G4double aE) {
    if (fEdepOnGeWaferGoldFilm.size() > 0 &&
        MACRO_IS_IN_RANGE_OF(0, aIdx, (G4int)fEdepOnGeWaferGoldFilm.size())) {
        fEdepOnGeWaferGoldFilm[aIdx] = aE;
        return true;
    } else
        return false;
}

inline G4bool AmoreModuleHit::AddGeWaferGoldFilmEdep(G4int aIdx, G4double aE) {
    if (fEdepOnGeWaferGoldFilm.size() > 0 &&
        MACRO_IS_IN_RANGE_OF(0, aIdx, (G4int)fEdepOnGeWaferGoldFilm.size())) {
        fEdepOnGeWaferGoldFilm[aIdx] += aE;
        return true;
    } else
        return false;
}

inline G4double AmoreModuleHit::GetCrystalGoldFilmEdep(G4int aIdx) const {
    if (fEdepOnGeWaferGoldFilm.size() > 0 &&
        MACRO_IS_IN_RANGE_OF(0, aIdx, (G4int)fEdepOnCrystalGoldFilm.size())) {
        return fEdepOnCrystalGoldFilm[aIdx];
    } else
        return -1;
}

inline G4double AmoreModuleHit::GetGeWaferGoldFilmEdep(G4int aIdx) const {
    if (fEdepOnGeWaferGoldFilm.size() > 0 &&
        MACRO_IS_IN_RANGE_OF(0, aIdx, (G4int)fEdepOnGeWaferGoldFilm.size())) {
        return fEdepOnGeWaferGoldFilm[aIdx];
    } else
        return -1;
}

#endif
