#ifndef __AmoreDetectorStaticType_HH__
#define __AmoreDetectorStaticType_HH__ 1

#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class AmoreModuleHit;

#define IS_APPROX_SAME_MACRO(X, Y)                                                                 \
    (fabs((X) - (Y)) < (G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()))

namespace AmoreDetectorStaticInfo {
    /** @struct CrystalModuleInfo
     *  @brief This struct contains dimension and meta information of a single crystal detector
     * module
     * */
    struct CrystalModuleInfo {
        G4String fName;                  //< Name of detector module
        G4String fFrameMaterialName;     //< Name of material for frame of module
        G4String fCrystalMaterialName;   //< Name of material for crystal in the module
        G4String fReflectorMaterialName; //< Name of material for reflector around the crystal
				G4String fPhononCollectorMaterialName; //< Name of material for phonon collector 
        G4double fCrystalHeight;         //< Crystal height
        G4double fCrystalMajorDiameter;  //< Crystal major radius (elliptical)
        G4double fCrystalMinorDiameter;  //< Crystal minor radius (elliptical)
				G4double fCrystalMass;           //< Crystal mass

        /// = operator for std::map storing
        G4bool operator=(const CrystalModuleInfo &rhs) const {
            return fFrameMaterialName == rhs.fFrameMaterialName &&
                   fCrystalMaterialName == rhs.fCrystalMaterialName &&
                   fReflectorMaterialName == rhs.fReflectorMaterialName &&
                   IS_APPROX_SAME_MACRO(fCrystalHeight, rhs.fCrystalHeight) &&
                   IS_APPROX_SAME_MACRO(fCrystalMajorDiameter, rhs.fCrystalMajorDiameter) &&
                   IS_APPROX_SAME_MACRO(fCrystalMinorDiameter, rhs.fCrystalMinorDiameter) &&
									 IS_APPROX_SAME_MACRO(fCrystalMass, rhs.fCrystalMass) ;
        };
    };

		struct AMoRE200CrystalModuleInfo {
			G4String fName;
			G4double fCrystalRadius;
			G4double fCrystalHeight;
			G4double fTowerX;
			G4double fTowerY;

      G4bool operator=(const AMoRE200CrystalModuleInfo &rhs) const {
				return 
					IS_APPROX_SAME_MACRO(fCrystalRadius, rhs.fCrystalRadius) &&
					IS_APPROX_SAME_MACRO(fCrystalHeight, rhs.fCrystalHeight) &&
					IS_APPROX_SAME_MACRO(fTowerX, rhs.fTowerX) &&
					IS_APPROX_SAME_MACRO(fTowerY, rhs.fTowerY);
				
			};
		};
} // namespace AmoreDetectorStaticInfo

#endif
