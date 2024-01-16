// AmoreCMOParameterisation header file
// 2017/9/11 by BaseHardware

#ifndef AmoreCMOParameterisation_h

#define AmoreCMOParameterisation_h 1

#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4VPVParameterisation.hh"
#include "globals.hh"

#include <vector>

class G4VPhysicalVolume;
class G4Box;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

///  A parameterisation that describes a series of boxes along Z.
///
///  The boxes have equal width, & their lengths are a linear equation.
///  They are spaced an equal distance apart, starting from given location.

class AmoreCMOParameterisation : public G4VPVParameterisation {
  public:
    typedef enum { kChuhMethod, kBasehwMethod } eTranslationMethod;

    AmoreCMOParameterisation(G4double aCellsize_R, G4double aCellsize_H, G4double aRSpacing,
                             G4double aZSpacing, G4int aFloorNum, G4int aLayerNum,
                             eTranslationMethod aMethod = kBasehwMethod);
    AmoreCMOParameterisation(G4double aCellsize_R, G4double aCellsize_H, G4double aRSpacing,
                             G4double aZSpacing, eTranslationMethod aMethod = kBasehwMethod);

    virtual ~AmoreCMOParameterisation();

    void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const;
    G4bool Optimize(G4int, G4double);
    G4bool Initialize();

    G4double GetEnvelopeLength() { return fInitialized ? fEnvelL : -1; };
    G4double GetTotalHeight() { return fInitialized ? fTotH : -1; };
    G4int GetNumberOfTotalCellNum() { return fInitialized ? fTotCellNum : -1; };
    G4int GetNumberOfFloorCellNum() { return fInitialized ? fFloorCellNum : -1; };
    G4VSolid *GetEnvelopeSolid() { return fInitialized ? fEnvelopeSolid : nullptr; };
    G4VSolid *GetCellSolid() { return fInitialized ? fCellSolid : nullptr; };

    void SetCellsize_R(G4double a) { fCSR = a; };
    void SetCellsize_H(G4double a) { fCSH = a; };
    void SetSpacing_R(G4double a) { fRS = a; };
    void SetSpacing_H(G4double a) { fHS = a; };
    void SetFloorNumber(G4int a) { fFN = a; };
    void SetLayerNumber(G4int a) { fLN = a; };
    void SetTranslationMethod(eTranslationMethod a) { fM = a; };

  private: // Dummy declarations to get rid of warnings ...
    void ComputeDimensions(G4Tubs &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Box &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Trd &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Trap &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Cons &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Sphere &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Orb &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Ellipsoid &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Torus &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Para &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Hype &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Polycone &, const G4int, const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Polyhedra &, const G4int, const G4VPhysicalVolume *) const {}

  private:
    G4bool CheckSettings();
    G4bool GenerateTranslationTable();
    G4bool BuildEnvelopeSolid();

    G4bool fInitialized;
    G4double fCSR, fCSH, fRS, fHS;
    G4int fFN, fLN;
    eTranslationMethod fM;

    G4double fEnvelL;
    G4double fTotH;
    G4int fFloorCellNum;
    G4int fTotCellNum;

    G4VSolid *fEnvelopeSolid;
    G4VSolid *fCellSolid;

    G4ThreeVector *fTranslationTable;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
