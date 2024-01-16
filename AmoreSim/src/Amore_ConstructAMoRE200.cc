#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "CupSim/CupPMTSD.hh"
#include "CupSim/CupParam.hh"
#include "CupSim/CupScintSD.hh"            // for making sensitive photocathodes
#include "CupSim/CupTorusStack.hh"         // for making the cavern and LSC envelope
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh" // for making spheres
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"

//#include "CLHEP/Matrix/Matrix.h"
#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include "G4ios.hh"

#include "G4Types.hh"

#include <sstream>

//#include "Randomize.hh"                     // for G4UniformRand()

//#include "fstream"
//#include "strstream"
//#include <stdio.h>

// == Construct Geometry for AMoRE-200 ======================================
// uses parameters from database or file
void AmoreDetectorConstruction::ConstructAMoRE200() {
    // --- put in the Inner Detector tanks, PMTs, and details
    G4LogicalVolume* odLV = ConstructAMoRE200_OD();

    // --- put in the Inner Detector tanks, PMTs, and details
    ConstructAMoRE200_ID(odLV);
}

void AmoreDetectorConstruction::ConstructAMoRE200_SDandField() {
    //////////////////////////////
    // --- TG sensitive detector
    //////////////////////////////
    // get pointer to logical volume

    CupScintSD *TGSD;
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    G4String SDname, VetoActiveMatPVName;

    TGSD = new CupScintSD(SDname = "/CupDet/TGSD", f200_TotCrystalNum);
    SDman->AddNewDetector(TGSD);
		//f200_logiCrystalCell->SetSensitiveDetector(TGSD);
		for(int i = 0; i < f200_TotCrystalNum; i++){
			f200_logiCrystalCell[i]->SetSensitiveDetector(TGSD);
		}

		/*
    G4LogicalVolume *logiGeWafer = f200_physGeWafer->GetLogicalVolume();
    G4LogicalVolume *logiVacDisk = f200_physVacDisk->GetLogicalVolume();
    CupPMTSD *pmtSD;
    pmtSD = new CupPMTSD("/cupdet/pmt/MLCS", f200_CMOTotCNum + f200_CMOFloorCNum);
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD);
    logiGeWafer->SetSensitiveDetector(pmtSD);
    logiVacDisk->SetSensitiveDetector(pmtSD);
		*/
}
