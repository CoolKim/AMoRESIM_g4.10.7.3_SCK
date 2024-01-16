//
//  Original by G. Horton-Smith 2004/12/02
//
//  Modified by E.J.Jeon 2007/06/14

#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "CupSim/CupBoxSD.hh"                    // for making sensitive boxes
#include "CupSim/CupPMTOpticalModel.hh"          // for same PMT optical model as main sim
#include "CupSim/CupPMTSD.hh"                    // for making sensitive photocathodes
#include "CupSim/CupParam.hh"
#include "CupSim/CupScintSD.hh"            // for making sensitive photocathodes
#include "CupSim/CupTorusStack.hh"         // for making the balloon
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Element.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh" // for making spheres
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"

////#include "CLHEP/Matrix/Matrix.h"
#include "G4LogicalVolume.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include "Randomize.hh" // for G4UniformRand()

#include "G4PVDivision.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"

#include "G4VSensitiveDetector.hh"

#include <cstdio>
#include <fstream>
#include <iostream>

#include "G4Region.hh"
#include "G4RegionStore.hh"

//using namespace CLHEP;
using namespace std;

void AmoreDetectorConstruction::ConstructMyDetector() {
    /////////////////////////////////////////
    // Construct your detector world here.
    // You must at least do the following:
    // * Define a world physical volume.
    // * Set the class variable "world_phys" to point to the world phys volume.

    // -- Make a 20m x 20m x 20m "Experimental Hall" for the world volume
    G4double bounding_size    = 20. * meter / 2.0;
    G4Box *boxHall            = new G4Box("hallbox", bounding_size, bounding_size, bounding_size);
    G4LogicalVolume *logiHall = new G4LogicalVolume(boxHall, _air, "logiHall");
    logiHall->SetVisAttributes(G4VisAttributes::Invisible);

    G4VPhysicalVolume *physHall = new G4PVPlacement(0,                      // rotation
                                                    G4ThreeVector(0, 0, 0), // translation
                                                    logiHall,   // associated logical vol
                                                    "physHall", // name
                                                    NULL,       // parent
                                                    false,      // no "Many"
                                                    0);         // copy number

		G4double LMOCell_radius      = 50. * mm / 2.;
		G4double LMOCell_half_height = 50. * mm / 2.;

		G4Tubs *LMOCell = new G4Tubs("LMOCell",0, LMOCell_radius, LMOCell_half_height, 0, 360.);
		G4LogicalVolume *LMOCell_LV = new G4LogicalVolume(LMOCell, _Li2MoO4, "LMOCell_lv");

		G4Colour yellow(1.0, 1.0, 0.0);
		G4VisAttributes *LMOCell_Vis = new G4VisAttributes(yellow);
		LMOCell_Vis->SetForceSolid(true);

		new G4PVPlacement(nullptr, {0,0,0}, LMOCell_LV, "LMOCell_PV", logiHall, false, 0, false);

    // the class variable "world_phys" must be set to the world phys volume.
    world_phys = physHall;

    CupScintSD *TGSD;
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    TGSD = new CupScintSD(SDname = "/CupDet/TGSD", 1);
    SDman->AddNewDetector(TGSD);
    LMOCell_LV->SetSensitiveDetector(TGSD);
}
