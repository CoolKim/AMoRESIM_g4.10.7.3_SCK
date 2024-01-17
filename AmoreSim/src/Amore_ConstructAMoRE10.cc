//
//  Original by G. Horton-Smith 2004/12/02
//
//  Modified by E.J.Jeon 2007/06/14
//  Modified by Y.S.Yoon 2015/06/15

#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "CupSim/CupBoxSD.hh"                    // for making sensitive boxes
#include "CupSim/CupPMTOpticalModel.hh"          // for same PMT optical model as main sim
#include "CupSim/CupPMTSD.hh"                    // for making sensitive photocathodes
#include "CupSim/CupParam.hh"
#include "AmoreSim/AmoreScintSD.hh"            // for making sensitive photocathodes
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Element.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"
//#include "CLHEP/Matrix/Matrix.h"

#include <fstream>
//#include <sstream>
#include <cstdio>

#include "G4Region.hh"
#include "G4Types.hh"

#include "G4NistManager.hh"

// -- database
// == Construct Geometry for GenericLAND ======================================
// uses parameters from database or file
void AmoreDetectorConstruction::ConstructAMoRE10() {
    /*
      // -- database
      CupParam &db ( CupParam::GetDB() );
      // add spherical-geometry-specific parameters to parameter list
      if ( getenv("AmoreDATA") != NULL )
        db.ReadFile( (G4String(getenv("AmoreDATA"))+"/settings.dat").c_str() );
      else
        //db.ReadFile("../Cupsim/data/settings_cylindrical.dat");
        //db.ReadFile("./data/settings_cylindrical.dat");	// EJ
        db.ReadFile((G4String(getenv("CUPSOFT_DIR"))+"/CupSim/data/settings.dat").c_str());


      // -- retrieve some parameters from database we'll need right away
      G4double tunnel_arch_radius= db["tunnel_arch_radius"];
      G4double pit_depth= db["pit_depth"];
      G4double pit_radius= db["pit_radius"];
      G4double rock_shell_thickness= db["rock_shell_thickness"];

      G4double LS_tank_x= db["LS_tank_x"];
      G4double LS_tank_y= db["LS_tank_y"];
      G4double LS_tank_z= db["LS_tank_z"];
      G4double LS_tank_thickness= db["LS_tank_thickness"];

      G4double pb_tank_x= db["pb_tank_x"];
      G4double pb_tank_y= db["pb_tank_y"];
      G4double pb_tank_z= db["pb_tank_z"];
      G4double pb_tank_thickness= db["pb_tank_thickness"];

      G4double buffer_tank_x= db["buffer_tank_x"];
      G4double buffer_tank_y= db["buffer_tank_y"];
      G4double buffer_tank_z= db["buffer_tank_z"];
      G4double buffer_tank_thickness= db["buffer_tank_thickness"];

      G4double gamma_catcher_vessel_height= db["gamma_catcher_vessel_height"];
      G4double gamma_catcher_vessel_radius= db["gamma_catcher_vessel_radius"];
      G4double gamma_catcher_vessel_thickness=db["gamma_catcher_vessel_thickness"];

      G4double target_vessel_height= db["target_vessel_height"];
      G4double target_vessel_radius= db["target_vessel_radius"];
      G4double target_vessel_thickness= db["target_vessel_thickness"];
    */
    G4int flagInvisible = 0;
    G4int flagOneCell   = 0;

    ///////////////////////////////////////////////////////
    // -- make nested cylindrical tanks and fill volumes
    // make colours
    G4Colour white(1.0, 1.0, 1.0);
    G4Colour grey(0.5, 0.5, 0.5);
    G4Colour greyl(0.5, 0.5, 0.5, 0.3);
    G4Colour lgrey(.85, .85, .85);
    G4Colour red(1.0, 0.0, 0.0);
    G4Colour blue(0.0, 0.0, 1.0);
    G4Colour cyan(0.0, 1.0, 1.0);
    G4Colour cyanl(0.0, 1.0, 1.0, 0.5);
    G4Colour magenta(1.0, 0.0, 1.0);
    G4Colour yellow(1.0, 1.0, 0.0);
    G4Colour orange(.75, .55, 0.0);
    G4Colour orangel(.75, .55, 0.0, 0.5);
    G4Colour lblue(0.0, 0.0, .75);
    G4Colour bluel(0.0, 0.0, 1.0, 0.5);
    G4Colour lgreen(0.0, .75, 0.0);
    G4Colour green(0.0, 1.0, 0.0);
    G4Colour greenl(0.0, 1.0, 0.0, 0.3);
    G4Colour brown(0.7, 0.4, 0.1);
    G4Colour brownl(0.7, 0.4, 0.1, 0.5);

    G4double deg0 = 0. * deg;
    //  G4double  twopi	= 360.* deg;

    ///////////////////////////////////////////////
    // * Define a world physical volume.
    // * Set the class variable "world_phys" to point to the world phys volume.
    // -- Make a 20m x 20m x 20m "Experimental Hall" for the world volume
    G4double bounding_size = 20. * meter / 2.0;
    G4Box *boxHall;
    G4LogicalVolume *logiHall;
    G4VisAttributes *logiHallVis;
    G4VPhysicalVolume *physHall;

    /////////////////////////////////////////////////////////////////
    // Rock Shell
    // Sizes: 11m x 11m x 5.5m
    // Air Room inside Rock Shell
    // Sizes: 10m x 10m x 4.5m
    G4double RockShellSize   = 11.0 * meter / 2.0;
    G4double RockShellHeight = 5.5 * meter / 2.0;
    G4Box *RockSolid;
    G4LogicalVolume *logiRockShell;
    G4VisAttributes *logiRockShellVis;
    G4VPhysicalVolume *physRockShell;

    /////////////////////////////////////////
    // Work Area   10 m x 10 m x 4.5 m
    G4double workbox  = 10. * m / 2; // 3000. ;
    G4double workboxH = 4.5 * m / 2; // 3000.;
    G4Box *WorkArea;
    G4LogicalVolume *logiWorkArea;
    G4VisAttributes *logiWorkAreaVis;
    G4VPhysicalVolume *physWorkArea;

    /////////////////////////////////////////
    //  Pb Top Box
    G4double pbtopbox_size  = 150. * cm / 2.0;
    G4double pbtopbox_zsize = 15. * cm / 2.0;
    G4double vgaptoTop      = 49.4 * cm;   // (46.4 cm  + 3 cm sstop thickness)
    G4double PbTBoxVSize    = 175.45 * cm; //  114.05 cm (center to top of Stainless Steel cover)
                                           //   + 46.4 cm (gap)
                                           //   + 15.0 cm (Lead top thickmess)
    G4double CenterPbTZ    = PbTBoxVSize - pbtopbox_zsize;
    G4ThreeVector PosTopPb = G4ThreeVector(0., 0., CenterPbTZ);

    G4Box *TopPbBox;
    G4LogicalVolume *logiTopPbBox;
    G4VisAttributes *logiTopPbBoxVis;
    G4VPhysicalVolume *physTopPbBox;

    /////////////////////////////////////////////////////////////////
    //  Pb  Box
    G4double pbbox_thick = 15. * cm;
    G4double pbbox_size  = (75. * cm + 2 * pbbox_thick) / 2.0;
    G4double pbbox_zsize = (162.1 * cm + pbbox_thick) / 2.0;
    G4Box *PbBoxOut;
    G4Box *PbBoxIn;
    G4SubtractionSolid *PbBox;
    G4LogicalVolume *logiPbBox;
    G4VisAttributes *logiPbBoxVis;
    G4VPhysicalVolume *physPbBox;

    G4double PbboxTopZ     = PbTBoxVSize - vgaptoTop - pbtopbox_zsize * 2;
    G4double CenterPbZ     = PbboxTopZ - pbbox_zsize;
    G4ThreeVector PosPbBox = G4ThreeVector(0., 0., CenterPbZ);

    /////////////////////////////////////////////////////////////////
    // Layer5_OVC Stainless Steel (SS)
    G4double ss_radius      = 64. * cm / 2;
    G4double ss_thick       = 0.5 * cm;
    G4double ss_height      = 153.8 * cm / 2; // 153.8
    G4double sstopbox_xsize = 74.8 * cm / 2.0;
    G4double sstopbox_ysize = 74.8 * cm / 2.0;
    G4double sst_zsize      = 3. * cm / 2.;
    G4double ssb_zsize      = 2. * cm / 2;
    G4double plateGap       = 12. * cm;

    G4double CenterTSSZ         = ss_height + sst_zsize;
    G4double CenterSSZ          = CenterPbZ + pbbox_zsize - ss_height;
    G4ThreeVector PosSS         = G4ThreeVector(0., 0., CenterSSZ);
    G4ThreeVector SSTopShift    = G4ThreeVector(0., 0., ss_height + sst_zsize);
    G4ThreeVector SSBottomShift = G4ThreeVector(0., 0., -ss_height - ssb_zsize);

    G4Tubs *SSCylinder;
    G4Box *SSTop;
    G4Tubs *SSBottom;
    G4VSolid *SSOVC0;
    G4VSolid *SSOVC1;
    G4LogicalVolume *logiSSOVC;
    G4VisAttributes *logiSSOVCVis;
    G4VPhysicalVolume *physSSOVC;

    /////////////////////////////////////////////////////////////////
    // Layer4_ 50K-SHIELD Copper
    G4double cu4_radius = 57.6 * cm / 2;
    G4double cu4_thick  = 0.3 * cm;
    G4double cu4_height = 134.2 * cm / 2; // 136.9-1.0-1.7
    G4double cu4t_zsize = 1.7 * cm / 2;
    G4double cu4b_zsize = 1.0 * cm / 2;

    G4double CenterCu4Z          = CenterSSZ + ss_height - cu4_height - plateGap - cu4t_zsize * 2;
    G4ThreeVector PosCu4         = G4ThreeVector(0., 0., CenterCu4Z);
    G4ThreeVector Cu4TopShift    = G4ThreeVector(0., 0., cu4_height + cu4t_zsize);
    G4ThreeVector Cu4BottomShift = G4ThreeVector(0., 0., -cu4_height - cu4b_zsize);

    G4Tubs *Cu4Cylinder;
    G4Tubs *Cu4Top;
    G4Tubs *Cu4Bottom;
    G4VSolid *Cu4Shield0;
    G4VSolid *Cu4Shield1;
    G4LogicalVolume *logiCu4;
    G4VisAttributes *logiCu4Vis;
    G4VPhysicalVolume *physCu4;

    ///////////////////////////////////////////////////////
    // Layer3_Cu IVC

    G4double cu3_radius = 51.4 * cm / 2;
    G4double cu3_thick  = 0.5 * cm;
    G4double cu3_height = 116.18 * cm / 2; // 119.88 -1.8 - 1.9
    G4double cu3t_zsize = 1.8 * cm / 2;
    G4double cu3b_zsize = 1.9 * cm / 2;

    G4double CenterCu3Z          = CenterCu4Z + cu4_height - cu3_height - plateGap - cu3t_zsize * 2;
    G4ThreeVector PosCu3         = G4ThreeVector(0., 0., CenterCu3Z);
    G4ThreeVector Cu3TopShift    = G4ThreeVector(0., 0., cu3_height + cu3t_zsize);
    G4ThreeVector Cu3BottomShift = G4ThreeVector(0., 0., -cu3_height - cu3b_zsize);

    G4Tubs *Cu3Cylinder;
    G4Tubs *Cu3Top;
    G4Tubs *Cu3Bottom;
    G4VSolid *Cu3Shield0;
    G4VSolid *Cu3Shield1;
    G4LogicalVolume *logiCu3;
    G4VisAttributes *logiCu3Vis;
    G4VPhysicalVolume *physCu3;

    ///////////////////////////////////////////////////////
    // Layer2_Cu SHIELD-STILL

    G4double cu2_radius = 45.5 * cm / 2;
    G4double cu2_thick  = 0.1 * cm;
    G4double cu2_height = 100.2 * cm / 2;
    G4double cu2t_zsize = 2.0 * cm / 2;
    G4double cu2b_zsize = 0.1 * cm / 2;

    G4double CenterCu2Z          = CenterCu3Z + cu3_height - cu2_height - plateGap - cu2t_zsize * 2;
    G4ThreeVector PosCu2         = G4ThreeVector(0., 0., CenterCu2Z);
    G4ThreeVector Cu2TopShift    = G4ThreeVector(0., 0., cu2_height + cu2t_zsize);
    G4ThreeVector Cu2BottomShift = G4ThreeVector(0., 0., -cu2_height - cu2b_zsize);

    G4Tubs *Cu2Cylinder;
    G4Tubs *Cu2Top;
    G4Tubs *Cu2Bottom;
    G4VSolid *Cu2Shield0;
    G4VSolid *Cu2Shield1;
    G4LogicalVolume *logiCu2;
    G4VisAttributes *logiCu2Vis;
    G4VPhysicalVolume *physCu2;

    ///////////////////////////////////////////////////////
    // Layer1_Cu 50mK-SHIELD

    G4double cu1_radius = 43.5 * cm / 2;
    G4double cu1_thick  = 0.1 * cm;
    G4double cu1_height = 85.1 * cm / 2;
    G4double cu1t_zsize = 2.0 * cm / 2;
    G4double cu1b_zsize = 0.1 * cm / 2;

    G4double CenterCu1Z          = CenterCu2Z + cu2_height - cu1_height - plateGap - cu1t_zsize * 2;
    G4ThreeVector PosCu1         = G4ThreeVector(0., 0., CenterCu1Z);
    G4ThreeVector Cu1TopShift    = G4ThreeVector(0., 0., cu1_height + cu1t_zsize);
    G4ThreeVector Cu1BottomShift = G4ThreeVector(0., 0., -cu1_height - cu1b_zsize);

    G4Tubs *Cu1Cylinder;
    G4Tubs *Cu1Top;
    G4Tubs *Cu1Bottom;
    G4VSolid *Cu1Shield0;
    G4VSolid *Cu1Shield1;
    G4LogicalVolume *logiCu1;
    G4VisAttributes *logiCu1Vis;
    G4VPhysicalVolume *physCu1;

    ///////////////////////////////////////////////////////
    // Cu Plate: Mixing Chamber SHIELD
    G4double cumcp_radius = 40.8 * cm / 2;
    G4double cumcp_height = 3. * cm / 2;

    G4double CenterCuMCPZ  = CenterCu1Z + cu1_height - cumcp_height - plateGap;
    G4ThreeVector PosCuMCP = G4ThreeVector(0., 0., CenterCuMCPZ);

    G4Tubs *CuMCPlate;
    G4LogicalVolume *logiCuMCP;
    G4VisAttributes *logiCuMCPVis;
    G4VPhysicalVolume *physCuMCP;

    ///////////////////////////////////////////////////////
    // Plate1_Cu
    G4double cup1_zsize   = 0.5 * cm / 2;
    G4double plateGap2    = 1.0 * cm;
    G4double CenterCuP1Z  = CenterCuMCPZ - cumcp_height - cup1_zsize - plateGap2;
    G4ThreeVector PosCuP1 = G4ThreeVector(0., 0., CenterCuP1Z);

    G4Tubs *CuPlate1;
    G4LogicalVolume *logiCuP1;
    G4VisAttributes *logiCuP1Vis;
    G4VPhysicalVolume *physCuP1;

    ///////////////////////////////////////////////////////
    // Plate2_Pb Shield
    G4double pbp2_zsize   = 5. * cm / 2;
    G4double CenterPbP2Z  = CenterCuP1Z - cup1_zsize - pbp2_zsize;
    G4ThreeVector PosPbP2 = G4ThreeVector(0., 0., CenterPbP2Z);

    G4Tubs *PbPlate2;
    G4LogicalVolume *logiPbP2;
    G4VisAttributes *logiPbP2Vis;
    G4VPhysicalVolume *physPbP2;

    ///////////////////////////////////////////////////////
    // Plate3_Cu
    G4double cup3_zsize   = 1. * cm / 2;
    G4double CenterCuP3Z  = CenterPbP2Z - pbp2_zsize - cup3_zsize;
    G4ThreeVector PosCuP3 = G4ThreeVector(0., 0., CenterCuP3Z);

    G4Tubs *CuPlate3;
    G4LogicalVolume *logiCuP3;
    G4VisAttributes *logiCuP3Vis;
    G4VPhysicalVolume *physCuP3;

    ///////////////////////////////////////////////////////
    // Plate4_Cu
    G4double CenterCuP4Z  = CenterCuP3Z - cup3_zsize - cup1_zsize - plateGap2;
    G4ThreeVector PosCuP4 = G4ThreeVector(0., 0., CenterCuP4Z);
    G4VPhysicalVolume *physCuP4;

    ///////////////////////////////////////////////////////
    // Plate5_Pb Shield
    G4double CenterPbP5Z  = CenterCuP4Z - cup1_zsize - pbp2_zsize;
    G4ThreeVector PosPbP5 = G4ThreeVector(0., 0., CenterPbP5Z);
    G4VPhysicalVolume *physPbP5;

    ///////////////////////////////////////////////////////
    // Plate6_Cu
    G4double CenterCuP6Z  = CenterPbP5Z - pbp2_zsize - cup3_zsize;
    G4ThreeVector PosCuP6 = G4ThreeVector(0., 0., CenterCuP6Z);
    G4VPhysicalVolume *physCuP6;

    ///////////////////////////////////////////////////////
    //  Target Room
    G4double TR_radius = 30. * cm / 2;
    G4double TR_height = 45.1 * cm / 2; // Nb shield size needs to be confirmed
    G4double GapCuNb   = 7 * cm;
    G4double CenterTR = cu1_height - plateGap - cumcp_height * 2 - pbp2_zsize * 2 - cup3_zsize * 2 -
                        GapCuNb - TR_height;
    //  G4ThreeVector		PosTR = G4ThreeVector( 0.,0., CenterTR);
    G4ThreeVector PosTR = G4ThreeVector(0., 0., 0.);
    G4Tubs *TargetRoom;
    G4LogicalVolume *logiTargetRoom;
    G4VisAttributes *logiTargetRoomVis;
    G4VPhysicalVolume *physTargetRoom;

    ///////////////////////////////////////////////////////
    // CMO cell
    G4double cmocell_r     = 4.5 * cm / 2.;
    G4double cmocell_zsize = 4.5 * cm / 2.;
    G4double CMOSupportGap = 1 * cm;
    G4double CMOSupportGapV= 1*cm; //EJ added
    G4double CMOSupport_W  = 0.5 * cm;
    G4double CMOSupport_H  = 0.5 * cm / 2.;
    G4double copperRods_r  = 0.5 * cm / 2.;
    G4double copperRods_H  = 6. * cm / 2.;

    G4Tubs *CMOCell;
    G4LogicalVolume *logiCMOCell;
    G4VisAttributes *logiCMOVis;
    G4VPhysicalVolume *physCMOCell;

    G4double xCell[7][6];
    G4double yCell[7][6];
    G4double zCell[7][6];
    G4double zCup[7][6];
    G4double zCdown[7][6];
    G4int copyNo;

    G4ThreeVector PosCMOCell;

    ///////////////////////////////////////////////////////
    // Photon Detector: Ge Wafer + Vaccum disk

    G4double PhotonDet_r     = 2.54 * 2 * cm / 2. + 1. * cm;
    G4double PhotonDet_zsize = 1.0 * cm / 2.;

    G4double PhotonDetZ = cmocell_zsize + PhotonDet_zsize * 2;
    G4ThreeVector PosPhotonDet;

    G4Tubs *PhotonDet;
    G4LogicalVolume *logiPhotonDet;
    G4VisAttributes *logiPhotonDetVis;
    G4VPhysicalVolume *physPhotonDet;

    ///////////////////////////////////////////////////////
    // Ge wafer on the bottom of the Crystal

    G4double GeWafer_r     = 2.54 * 2 * cm / 2.;
    G4double GeWafer_zsize = 0.00004 * cm / 2.; // 300 nm + 100 nm Vac
    G4double GeVacGap      = 0.2 * cm;

    G4double cusupportx      = 7.5 * cm / 2.;
    G4double cusupportz      = 0.2 * cm / 2.;
    G4double cusupportinnerx = 6.5 * cm / 2.;

    G4ThreeVector PosGeWafer = G4ThreeVector(0., 0., 0); // Center of PhotonDet

    G4Tubs *GeWafer;
    G4LogicalVolume *logiGeWafer;
    G4VisAttributes *logiGeWaferVis;
    G4VPhysicalVolume *physGeWafer;

    ///////////////////////////////////////////////////////
    // Vacuum Disks on Ge Wafer

    G4double rxDisk;
    G4double ryDisk;
    G4int NDisk            = 0;
    G4double VacDisk_r     = GeWafer_r;
    G4double VacDisk_zsize = 0.00001 * cm / 2.; // 100 nm
    G4double VacDiskZ      = -GeWafer_zsize + VacDisk_zsize;

    G4Tubs *VacDisk;
    G4LogicalVolume *logiVacDisk;
    G4VisAttributes *logiVacDiskVis;
    G4VPhysicalVolume *physVacDisk;
    G4ThreeVector PosVacDisk;

    ///////////////////////////////////////////////////////
    // Copper frame
    //
    G4Tubs *CopperframeSolid;
    G4ThreeVector PosCopperframeUp;
    G4ThreeVector PosCopperframeDown;

    G4LogicalVolume *logiCopperframe;
    G4VisAttributes *logiCopperframeVis;
    G4VPhysicalVolume *physCopperUpperframe;
    G4VPhysicalVolume *physCopperLowerframe;

    ///////////////////////////////////////////////////////
    // Copper Rod

    G4int NRod = 0;
    G4double rxRod;
    G4double ryRod;

    G4Tubs *CopperframeRod;
    G4LogicalVolume *logiCopperframeRod;
    G4VisAttributes *logiCopperframeRodVis;
    G4VPhysicalVolume *physCopperframeRod;
    G4ThreeVector PosCuRod;

    G4double zero = 0. * cm;

  // EJ added
  ///////////////////////////////////////////////////////  
  // Envelop for Reflector
  //
  G4double              EnvelopReflector_r = cmocell_r + 0.2*cm; // 64 micro-meter
  G4double              EnvelopReflector_thick = 0.1*cm; // 64 micro-meter
  G4double              EnvelopReflector_height = cmocell_zsize + 0.3*cm;
  G4Tubs*               EnvelopReflector;
  G4ThreeVector         PosEnvRef;

  G4LogicalVolume*      logiEnvelopReflector;
  G4VisAttributes*      logiEnvelopReflectorVis;
  G4VPhysicalVolume*    physEnvelopReflector;

    ///////////////////////////////////////////////////////
    // Reflector around the CMO (Vm2000 shield)
    //
    G4double reflector_r     = cmocell_r + 0.2 * cm; // 64 micro-meter
    //G4double reflector_thick = 0.0064 * cm;          // 64 micro-meter
    G4double reflector_thick = 0.0164 * cm; // EJ corrected
    G4double reflector_height = cmocell_zsize ; // EJ added
    G4Tubs *Reflector;
    G4ThreeVector PosRef;

    G4LogicalVolume *logiReflector;
    G4VisAttributes *logiReflectorVis;
    G4VPhysicalVolume *physReflector;

  // EJ added
  ///////////////////////////////////////////////////////  
  // Vacuum inside Reflector 
  //
  G4double              reflectorVacuum_r = cmocell_r + 0.2*cm; // 64 micro-meter
  G4double              reflectorVacuum_thick = 0.01*cm; // 64 micro-meter
  G4double              reflectorVacuum_height = cmocell_zsize ; // 64 micro-meter
  G4Tubs*               ReflectorVacuum;
  G4ThreeVector         PosRefVacuum = G4ThreeVector(0.,0.,0);

  G4LogicalVolume*      logiReflectorVacuum;
  G4VisAttributes*      logiReflectorVacuumVis;
  G4VPhysicalVolume*    physReflectorVacuum;

    ////////////////////////////////////////////////////////////////
    /*
      PrintAMoRE10Coordinate("Top Pb Shield",
                            pbtopbox_size, PbTBoxVSize,
                            zero, CenterPbTZ, pbtopbox_zsize );

      PrintAMoRE10Coordinate("Pb Shield Box",
                            pbbox_size, pbbox_zsize-sst_zsize,
                            zero, CenterPbZ, pbbox_zsize);

      PrintAMoRE10Coordinate("SS TopCover", sstopbox_xsize, sst_zsize,
                            zero, CenterTSSZ, sst_zsize );

      PrintAMoRE10Coordinate("SS OVC    ",ss_radius, ss_height,
                            zero, CenterSSZ, ss_height );

      PrintAMoRE10Coordinate("Cu 50K-Shield", cu4_radius, cu4_height,
                            zero, CenterCu4Z, cu4_height );

      PrintAMoRE10Coordinate("Cu IVC    ", cu3_radius, cu3_height,
                            zero, CenterCu3Z, cu3_height );

      PrintAMoRE10Coordinate("Cu STILLShield", cu2_radius, cu2_height,
                            zero, CenterCu2Z, cu2_height );

      PrintAMoRE10Coordinate("Cu 50mK-STILL", cu1_radius, cu1_height,
                            zero, CenterCu1Z, cu1_height );

      PrintAMoRE10Coordinate("Cu MixingPlate", cumcp_radius, cumcp_height,
                            zero, CenterCuMCPZ, cumcp_height );
    */
    ////////////////////////////////////////////////////////////////
    // Physics Hall
    ////////////////////////////////////////////////////////////////

    boxHall  = new G4Box("hallbox", bounding_size, bounding_size, bounding_size);
    logiHall = new G4LogicalVolume(boxHall, _air, "logiHall");
    physHall = new G4PVPlacement(0,                      // rotation
                                 G4ThreeVector(0, 0, 0), // translation
                                 logiHall,               // associated logical vol
                                 "physHall",             // name
                                 NULL,                   // parent
                                 false,                  // no "Many"
                                 0);                     // copy number

    // the class variable "world_phys" must be set to the world phys volume.
    world_phys = physHall;

    /////////////////////////////////////////////////////////////////
    // Rock Shell
    ////////////////////////////////////////////////////////////////

    RockSolid     = new G4Box("RockShell", RockShellSize, RockShellSize, RockShellHeight);
    logiRockShell = new G4LogicalVolume(RockSolid, _rock, "logiRockShell");
    physRockShell = new G4PVPlacement(0,                      // rotation
                                      G4ThreeVector(0, 0, 0), // translation
                                      "physRockShell",        // name
                                      logiRockShell,          // associated logical vol
                                      physHall,               // parent
                                      false,                  // no "Many"
                                      0);                     // copy number
    (void)physRockShell;

    ///////////////////////////////////////////////////////
    // Work Area   3 m x 3 m x 3 m
    //
    WorkArea     = new G4Box("WorkArea", workbox, workbox, workboxH);
    logiWorkArea = new G4LogicalVolume(WorkArea, _air, "logiWorkArea");
    physWorkArea = new G4PVPlacement(0, // no rotation
                                     G4ThreeVector(0., 0., 0), logiWorkArea, "physWorkArea",
                                     logiRockShell, false, 0);
    (void)physWorkArea;

    ///////////////////////////////////////////////////////
    // Layer6_Pb Top Box 1500x1500x1674.5 mm
    //           Air box under the top plate volume

    TopPbBox     = new G4Box("TopPbBox", pbtopbox_size, pbtopbox_size, pbtopbox_zsize);
    logiTopPbBox = new G4LogicalVolume(TopPbBox, _lead, "logiTopPbBox");
    physTopPbBox =
        new G4PVPlacement(0, // no rotation
                          PosTopPb, logiTopPbBox, "physTopPbBox", logiWorkArea, false, 0);
    (void)physTopPbBox;

    ///////////////////////////////////////////////////////
    // Layer6_Pb body of Lead Box  1050
    PbBoxOut  = new G4Box("PbBoxOut", pbbox_size, pbbox_size, pbbox_zsize);
    PbBoxIn   = new G4Box("PbBoxIn", pbbox_size - pbbox_thick, pbbox_size - pbbox_thick,
                        pbbox_zsize - pbbox_thick / 2.);
    PbBox     = new G4SubtractionSolid("PbBox", PbBoxOut, PbBoxIn, 0,
                                   G4ThreeVector(0, 0, pbbox_thick / 2.));
    logiPbBox = new G4LogicalVolume(PbBox, _lead, "logiPbBox");
    physPbBox = new G4PVPlacement(0, // no rotation
                                  PosPbBox, logiPbBox, "physPbBox", logiWorkArea, false, 0);
    (void)physPbBox;

    ///////////////////////////////////////////////////////
    // Layer5_OVC Stainless Steel (SS)
    //
    SSCylinder = new G4Tubs("SSCylinder", ss_radius - ss_thick, ss_radius, ss_height, deg0, twopi);
    SSTop      = new G4Box("SSTop", sstopbox_xsize, sstopbox_ysize, sst_zsize);
    SSBottom   = new G4Tubs("SSBottom", 0, ss_radius, ssb_zsize, deg0, twopi);
    SSOVC0     = new G4UnionSolid("SSOVC0", SSCylinder, SSTop, 0, SSTopShift);
    SSOVC1     = new G4UnionSolid("SSOVC1", SSOVC0, SSBottom, 0, SSBottomShift);
    logiSSOVC  = new G4LogicalVolume(SSOVC1, _stainless, "logiSSOVC");
    physSSOVC  = new G4PVPlacement(0, // no rotation
                                  PosSS, logiSSOVC, "physSSOVC", logiWorkArea, false, 0);
    (void)physSSOVC;

    ///////////////////////////////////////////////////////
    // Layer4_ 50K-SHIELD Copper

    Cu4Cylinder =
        new G4Tubs("Cu4Cylinder", cu4_radius - cu4_thick, cu4_radius, cu4_height, deg0, twopi);
    Cu4Top     = new G4Tubs("Cu4Top", 0, cu4_radius, cu4t_zsize, deg0, twopi);
    Cu4Bottom  = new G4Tubs("Cu4Bottom", 0, cu4_radius, cu4b_zsize, deg0, twopi);
    Cu4Shield0 = new G4UnionSolid("Cu4Shield0", Cu4Cylinder, Cu4Top, 0, Cu4TopShift);
    Cu4Shield1 = new G4UnionSolid("Cu4Shield1", Cu4Shield0, Cu4Bottom, 0, Cu4BottomShift);
    logiCu4    = new G4LogicalVolume(Cu4Shield1, _copper, "logiCu4", 0, 0, 0);
    physCu4    = new G4PVPlacement(0, // no rotation
                                PosCu4, logiCu4, "physCu4", logiWorkArea, false, 0);
    (void)physCu4;

    ///////////////////////////////////////////////////////
    // Layer3_Cu IVC

    Cu3Cylinder =
        new G4Tubs("Cu3Cylinder", cu3_radius - cu3_thick, cu3_radius, cu3_height, deg0, twopi);
    Cu3Top     = new G4Tubs("Cu3Top", 0, cu3_radius, cu3t_zsize, deg0, twopi);
    Cu3Bottom  = new G4Tubs("Cu3Bottom", 0, cu3_radius, cu3b_zsize, deg0, twopi);
    Cu3Shield0 = new G4UnionSolid("Cu3Shield0", Cu3Cylinder, Cu3Top, 0, Cu3TopShift);
    Cu3Shield1 = new G4UnionSolid("Cu3Shield1", Cu3Shield0, Cu3Bottom, 0, Cu3BottomShift);
    logiCu3    = new G4LogicalVolume(Cu3Shield1, _copper, "logiCu3", 0, 0, 0);
    physCu3    = new G4PVPlacement(0, // no rotation
                                PosCu3, logiCu3, "physCu3", logiWorkArea, false, 0);
    (void)physCu3;

    ///////////////////////////////////////////////////////
    // Layer2_Cu3 SHIELD-STILL

    Cu2Cylinder =
        new G4Tubs("Cu2Cylinder", cu2_radius - cu2_thick, cu2_radius, cu2_height, deg0, twopi);
    Cu2Top     = new G4Tubs("Cu2Top", 0, cu2_radius, cu2t_zsize, deg0, twopi);
    Cu2Bottom  = new G4Tubs("Cu2Bottom", 0, cu2_radius, cu2b_zsize, deg0, twopi);
    Cu2Shield0 = new G4UnionSolid("Cu2Shield0", Cu2Cylinder, Cu2Top, 0, Cu2TopShift);
    Cu2Shield1 = new G4UnionSolid("Cu2Shield1", Cu2Shield0, Cu2Bottom, 0, Cu2BottomShift);
    logiCu2    = new G4LogicalVolume(Cu2Shield1, _copper, "logiCu2", 0, 0, 0);
    physCu2    = new G4PVPlacement(0, // no rotation
                                PosCu2, logiCu2, "physCu2", logiWorkArea, false, 0);
    (void)physCu2;

    ///////////////////////////////////////////////////////
    // Layer1_Cu4 50mk-SHIELD

    Cu1Cylinder =
        new G4Tubs("Cu1Cylinder", cu1_radius - cu1_thick, cu1_radius, cu1_height, deg0, twopi);
    Cu1Top     = new G4Tubs("Cu1Top", 0, cu1_radius, cu1t_zsize, deg0, twopi);
    Cu1Bottom  = new G4Tubs("Cu1Bottom", 0, cu1_radius, cu1b_zsize, deg0, twopi);
    Cu1Shield0 = new G4UnionSolid("Cu1Shield0", Cu1Cylinder, Cu1Top, 0, Cu1TopShift);
    Cu1Shield1 = new G4UnionSolid("Cu1Shield1", Cu1Shield0, Cu1Bottom, 0, Cu1BottomShift);
    logiCu1    = new G4LogicalVolume(Cu1Shield1, _copper, "logiCu1", 0, 0, 0);
    physCu1    = new G4PVPlacement(0, // no rotation
                                PosCu1, logiCu1, "physCu1", logiWorkArea, false, 0);
    (void)physCu1;

    ///////////////////////////////////////////////////////
    // Mixing Chamber Cu Plate

    CuMCPlate = new G4Tubs("CuMCPlate", 0, cumcp_radius, cumcp_height, deg0, twopi);
    logiCuMCP = new G4LogicalVolume(CuMCPlate, _copper, "logiCuMCP", 0, 0, 0);
    physCuMCP = new G4PVPlacement(0, // no rotation
                                  PosCuMCP, logiCuMCP, "physCuMCP", logiWorkArea, false, 0);
    (void)physCuMCP;

    ///////////////////////////////////////////////////////
    // Plate1_Cu
    CuPlate1 = new G4Tubs("CuPlate1", 0, cumcp_radius, cup1_zsize, deg0, twopi);
    logiCuP1 = new G4LogicalVolume(CuPlate1, _copper, "logiCuP1", 0, 0, 0);
    physCuP1 = new G4PVPlacement(0, // no rotation
                                 PosCuP1, logiCuP1, "physCuP1", logiWorkArea, false, 0);
    (void)physCuP1;

    ///////////////////////////////////////////////////////
    // Plate2_Pb Shield
    PbPlate2 = new G4Tubs("PbPlate2", 0, cumcp_radius, pbp2_zsize, deg0, twopi);
    logiPbP2 = new G4LogicalVolume(PbPlate2, _lead, "logiPbP2", 0, 0, 0);
    physPbP2 = new G4PVPlacement(0, // no rotation
                                 PosPbP2, logiPbP2, "physPbP2", logiWorkArea, false, 0);
    (void)physPbP2;

    ///////////////////////////////////////////////////////
    CuPlate3 = new G4Tubs("CuPlate3", 0, cumcp_radius, cup3_zsize, deg0, twopi);
    logiCuP3 = new G4LogicalVolume(CuPlate3, _copper, "logiCuP3", 0, 0, 0);
    physCuP3 = new G4PVPlacement(0, // no rotation
                                 PosCuP3, logiCuP3, "physCuP3", logiWorkArea, false, 0);

    physCuP4 = new G4PVPlacement(0, // no rotation
                                 PosCuP4, logiCuP1, "physCuP4", logiWorkArea, false, 0);
    physPbP5 = new G4PVPlacement(0, // no rotation
                                 PosPbP5, logiPbP2, "physPbP5", logiWorkArea, false, 0);
    physCuP6 = new G4PVPlacement(0, // no rotation
                                 PosCuP6, logiCuP3, "physCuP6", logiWorkArea, false, 0);
    (void)physCuP3;
    (void)physCuP4;
    (void)physPbP5;
    (void)physCuP6;

    ///////////////////////////////////////////////////////
    // CMO  Nb SHIELD

    //  G4double Nbthickness= 0.1*cm ;
    //  G4double Nbthickness= 0*cm ;
    /*
      G4Tubs* NbShieldSolid=
        new G4Tubs("NbShieldSolid",0,radius0,halfNbHeight,0.*deg, 360.*deg);
      G4LogicalVolume* NbShieldLogi=
        new G4LogicalVolume(NbShieldSolid, Nb,"NbShieldLogical",0,0,0 );
      G4VisAttributes* NbShieldVisAtt = new  G4VisAttributes(orange);
      if ( flagOneCell)
        NbShieldLogi->SetVisAttributes (G4VisAttributes::GetInvisible);
      else
        NbShieldLogi->SetVisAttributes (NbShieldVisAtt);
      G4VPhysicalVolume* NbShieldPhys= new G4PVPlacement(0, // no rotation
                                           G4ThreeVector(0,0,0),
                                            "NbShieldPhys",
                                            NbShieldLogi,
                                            Lyr1CuPhys,
                                            false,
                                            0);
    */
    //  TargetRoomLogi->SetVisAttributes (G4VisAttributes::GetInvisible);

    ///////////////////////////////////////////////////////
    // CMO Target area

    TargetRoom     = new G4Tubs("TargetRoom", 0, TR_radius, TR_height, deg0, twopi);
    logiTargetRoom = new G4LogicalVolume(TargetRoom, _vacuum, "logiTargetRoom", 0, 0, 0);
    physTargetRoom =
        new G4PVPlacement(0, // no rotation
                          PosTR, logiTargetRoom, "physTargetRoom", logiWorkArea, false, 0);

    (void)physTargetRoom;

// EJ added
  for(int NLayer=0;NLayer<6;NLayer++) { //  5 layers + additional Photon Det.
    for(int NinL=0;NinL<7;NinL++)  {// 7 CMOs on a layer
        if ( NinL == 0 ) {
           xCell[NinL][NLayer] = 0;
           yCell[NinL][NLayer] = 0;
        }
        else {
           xCell[NinL][NLayer] = ((cmocell_r+CMOSupport_W+CMOSupportGap*2) * 2.
                                * cos((180.-NinL*60.)*deg));
           yCell[NinL][NLayer] = ((cmocell_r+CMOSupport_W+CMOSupportGap*2) * 2.
                                * sin((180.-NinL*60.)*deg));
        }
        zCell[NinL][NLayer] = ((2-NLayer)
                           * (cmocell_zsize+CMOSupport_H*2+CMOSupportGapV)*2);
        if ( NLayer == 5 ) continue;
        zCup[NinL][NLayer] = ((2-NLayer)
                           * (cmocell_zsize+CMOSupport_H*2+CMOSupportGapV)*2+(cmocell_zsize));
        zCdown[NinL][NLayer] = ((2-NLayer)
                           * (cmocell_zsize+CMOSupport_H*2+CMOSupportGapV)*2-(cmocell_zsize));
     }
  }

  ///////////////////////////////////////////////////////  
  // Envelop for Reflector
  // Reflector around the CMO (Vm2000 shield)
  // Reflector inside vacuum shell for Optical Surface
  //
  EnvelopReflector      = new G4Tubs("EnvelopReflector",
                        0, EnvelopReflector_r + EnvelopReflector_thick,
                        EnvelopReflector_height, deg0, twopi);
  logiEnvelopReflector  = new G4LogicalVolume(EnvelopReflector,
                        _vacuum, "logiEnvelopReflector");

  Reflector             = new G4Tubs("Reflector",
                        0, reflector_r + reflector_thick,
                        reflector_height, deg0, twopi);
  logiReflector         = new G4LogicalVolume(Reflector,
                        _vm2000, "logiReflector");

  physReflector         = new G4PVPlacement(0, // no rotation              
                                     PosRef,
                                     logiReflector,
                                     "physReflector",
                                     logiEnvelopReflector,
                                     false,
                                     0);

  ReflectorVacuum       = new G4Tubs("ReflectorVacuum",
                        0, reflectorVacuum_r + reflectorVacuum_thick,
                        reflectorVacuum_height, deg0, twopi);
  logiReflectorVacuum   = new G4LogicalVolume(ReflectorVacuum,
                        _vacuum, "logiReflectorVacuum");

   physReflectorVacuum  = new G4PVPlacement(0, // no rotation              
                                     PosRefVacuum,
                                     logiReflectorVacuum,
                                     "physReflectorVacuum",
                                     logiReflector,
                                     false,
                                     0);
  CMOCell       = new G4Tubs("CMOCell",  // name
                        0, cmocell_r, cmocell_zsize, deg0, twopi);
  logiCMOCell   = new G4LogicalVolume(CMOCell,
                        CaMoO4, "logiCMOCell", 0,0,0 );

  physCMOCell   = new G4PVPlacement(0, // no rotation              
                                     PosCMOCell,
                                     logiCMOCell,
                                     "physCMOCell",
                                     logiReflectorVacuum,
                                     false,
                                     0);

        PosEnvRef       = G4ThreeVector(0, 0, 0);
        physEnvelopReflector    = new G4PVPlacement(0, // no rotation              
                                          PosEnvRef,
                                          logiEnvelopReflector,
                                          "physEnvelopReflector",
                                          logiTargetRoom,
                                          false,
                                          0);

/* // EJ start
    ///////////////////////////////////////////////////////
    // CMO cell
    CMOCell     = new G4Tubs("CMOCell", // name
                         0, cmocell_r, cmocell_zsize, deg0, twopi);
    logiCMOCell = new G4LogicalVolume(CMOCell, CaMoO4, "logiCMOCell", 0, 0, 0);

    for (int NLayer = 0; NLayer < 5; NLayer++) { //  5 layers
        for (int NinL = 0; NinL < 7; NinL++) {   // 7 CMOs on a layer
            copyNo = NLayer * 7 + NinL;
            if (NinL == 0) {
                xCell[NinL][NLayer] = 0;
                yCell[NinL][NLayer] = 0;
            } else {
                xCell[NinL][NLayer] = ((cmocell_r + CMOSupport_W + CMOSupportGap / 2.) * 2. *
                                       cos((180. - NinL * 60.) * deg));
                yCell[NinL][NLayer] = ((cmocell_r + CMOSupport_W + CMOSupportGap / 2.) * 2. *
                                       sin((180. - NinL * 60.) * deg));
            }
            zCell[NinL][NLayer] =
                ((2 - NLayer) * (cmocell_zsize + CMOSupport_H * 2 + CMOSupportGap) * 2);
            zCup[NinL][NLayer] =
                ((2 - NLayer) * (cmocell_zsize + CMOSupport_H * 2 + CMOSupportGap) * 2 +
                 (cmocell_zsize));
            zCdown[NinL][NLayer] =
                ((2 - NLayer) * (cmocell_zsize + CMOSupport_H * 2 + CMOSupportGap) * 2 -
                 (cmocell_zsize));

            PosCMOCell =
                G4ThreeVector(xCell[NinL][NLayer], yCell[NinL][NLayer], zCell[NinL][NLayer]);

            physCMOCell = new G4PVPlacement(0, // no rotation
                                            PosCMOCell, logiCMOCell, "physCMOCell", logiTargetRoom,
                                            false, copyNo);

            if (flagOneCell) {
                NinL   = 7;
                NLayer = 5;
            }
        }
    }
    (void)physCMOCell;
*/ // EJ end

    ///////////////////////////////////////////////////////
    // Ge wafer on the bottom of the Crystal
    //

    PhotonDet     = new G4Tubs("PhotonDet", 0, PhotonDet_r, PhotonDet_zsize, deg0, twopi);
    logiPhotonDet = new G4LogicalVolume(PhotonDet, _vacuum, "logiPhotonDet");

    GeWafer     = new G4Tubs("GeWafer", 0, GeWafer_r, GeWafer_zsize, deg0, twopi);
    logiGeWafer = new G4LogicalVolume(GeWafer, _gewafer, "logiGeWafer", 0, 0, 0);

    physGeWafer =
        new G4PVPlacement(0, // no rotation
                          PosGeWafer, logiGeWafer, "physGeWafer", logiPhotonDet, false, 0);

    VacDisk     = new G4Tubs("VacDisk", 0, VacDisk_r, VacDisk_zsize, deg0, twopi);
    logiVacDisk = new G4LogicalVolume(VacDisk, _vacuum, "logiVacDisk", 0, 0, 0);

    PosVacDisk  = G4ThreeVector(0, 0, VacDiskZ);
    physVacDisk = new G4PVPlacement(0, // no rotation
                                    PosVacDisk, logiVacDisk, "physVacDisk", logiGeWafer, false, 0);
    (void)physVacDisk;

    ///////////////////////////////////////////////////////
    // Duplicate Photon Detectors
    //

    for (int NLayer = 0; NLayer < 6; NLayer++) { //  5 layers
        for (int NinL = 0; NinL < 7; NinL++) {   // 7 CMOs on a layer
            copyNo = NLayer * 7 + NinL;
            if (flagOneCell) {
                NLayer = 1;
            }
            if (NLayer == 5) {
                if (NinL == 0) {
                    xCell[NinL][NLayer] = 0;
                    yCell[NinL][NLayer] = 0;
                } else {
                    xCell[NinL][NLayer] = ((cmocell_r + CMOSupport_W + CMOSupportGap / 2.) * 2. *
                                           cos((180. - NinL * 60.) * deg));
                    yCell[NinL][NLayer] = ((cmocell_r + CMOSupport_W + CMOSupportGap / 2.) * 2. *
                                           sin((180. - NinL * 60.) * deg));
                }
                zCell[NinL][NLayer] =
                    ((2 - NLayer) * (cmocell_zsize + CMOSupport_H * 2 + CMOSupportGap) * 2);
            }
            PosPhotonDet  = G4ThreeVector(xCell[NinL][NLayer], yCell[NinL][NLayer],
                                         zCell[NinL][NLayer] + PhotonDetZ);
            physPhotonDet = new G4PVPlacement(0, // no rotation
                                              PosPhotonDet, logiPhotonDet, "physPhotonDet",
                                              logiTargetRoom, false, copyNo);
            if (flagOneCell) {
                NinL   = 7;
                NLayer = 6;
            }
        }
    }
    (void)physPhotonDet;

    ///////////////////////////////////////////////////////
    // Copper frame
    CopperframeSolid = new G4Tubs("CopperframeSolid", cmocell_r, cmocell_r + CMOSupport_W,
                                  CMOSupport_H, deg0, twopi);
    logiCopperframe  = new G4LogicalVolume(CopperframeSolid, _copper, "logiCopperframe");
    for (int NLayer = 0; NLayer < 5; NLayer++) { //  5 layers
        for (int NinL = 0; NinL < 7; NinL++) {   // 7 CMOs on a layer
            PosCopperframeUp   = G4ThreeVector(xCell[NinL][NLayer], yCell[NinL][NLayer],
                                             zCup[NinL][NLayer] + CMOSupportGap / 2.);
            PosCopperframeDown = G4ThreeVector(xCell[NinL][NLayer], yCell[NinL][NLayer],
                                               zCdown[NinL][NLayer] - CMOSupportGap / 2.);
            physCopperUpperframe =
                new G4PVPlacement(0, // no rotation
                                  PosCopperframeUp, logiCopperframe, "physCopperUpperframe",
                                  logiTargetRoom, false, 0);
            physCopperLowerframe =
                new G4PVPlacement(0, // no rotation
                                  PosCopperframeDown, logiCopperframe, "physCopperLowerframe",
                                  logiTargetRoom, false, 0);
            if (flagOneCell) {
                NinL   = 7;
                NLayer = 5;
            }
        }
    }
    (void)physCopperUpperframe;
    (void)physCopperLowerframe;

    ///////////////////////////////////////////////////////
    // Copper Rod
    CopperframeRod     = new G4Tubs("CopperframeRod", 0, copperRods_r, copperRods_H, deg0, twopi);
    logiCopperframeRod = new G4LogicalVolume(CopperframeRod, _copper, "logiCuSopportRod");

    for (int NLayer = 0; NLayer < 5; NLayer++) { //  5 layers
        for (int NinL = 0; NinL < 7; NinL++) {   // 7 CMOs on a layer
            for (NRod = 0; NRod < 3; NRod++) {   //
                rxRod = xCell[NinL][NLayer] +
                        ((cmocell_r + CMOSupport_W / 2.) * cos((180 - NRod * 120) * deg));
                ryRod = yCell[NinL][NLayer] +
                        ((cmocell_r + CMOSupport_W / 2.) * sin((180 - NRod * 120) * deg));
                PosCuRod           = G4ThreeVector(rxRod, ryRod, zCell[NinL][NLayer]);
                physCopperframeRod = new G4PVPlacement(
                    0, // no rotation
                    PosCuRod, logiCopperframeRod, "physCopperframeRod", logiTargetRoom, false, 0);
            }
            if (flagOneCell) {
                NinL   = 7;
                NLayer = 5;
            }
        }
    }
    (void)physCopperframeRod;

/* // EJ start
    ///////////////////////////////////////////////////////
    // Reflector around the CMO (Vm2000 shield)
    //
    Reflector = new G4Tubs("Reflector", reflector_r, reflector_r + reflector_thick, cmocell_zsize,
                           deg0, twopi);
    logiReflector = new G4LogicalVolume(Reflector, _vm2000, "logiReflector");

    for (int NLayer = 0; NLayer < 5; NLayer++) { //  5 layers
        for (int NinL = 0; NinL < 7; NinL++) {   // 7 CMOs on a layer
            PosRef = G4ThreeVector(xCell[NinL][NLayer], yCell[NinL][NLayer], zCell[NinL][NLayer]);
            physReflector =
                new G4PVPlacement(0, // no rotation
                                  PosRef, logiReflector, "physReflector", logiTargetRoom, false, 0);
            if (flagOneCell) {
                NinL   = 7;
                NLayer = 5;
            }
        }
    }
    (void)physReflector;
*/ // EJ end

    //////////////////////////////
    // Set Attributes
    //////////////////////////////

    logiHallVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.1));
    logiHall->SetVisAttributes(logiHallVis);
    logiHall->SetVisAttributes(G4VisAttributes::GetInvisible);

    logiTargetRoomVis = new G4VisAttributes(orange);
    logiTargetRoom->SetVisAttributes(logiTargetRoomVis);
    logiTargetRoom->SetVisAttributes(G4VisAttributes::GetInvisible);

    logiRockShellVis = new G4VisAttributes(red);
    logiWorkAreaVis  = new G4VisAttributes(grey);
    logiTopPbBoxVis  = new G4VisAttributes(grey);
    logiPbBoxVis     = new G4VisAttributes(greyl);
    logiSSOVCVis     = new G4VisAttributes(cyanl);
    logiCu4Vis       = new G4VisAttributes(bluel);
    logiCu3Vis       = new G4VisAttributes(green);
    logiCu2Vis       = new G4VisAttributes(cyan);
    logiCu1Vis       = new G4VisAttributes(brown);
    logiCuMCPVis     = new G4VisAttributes(lblue);
    logiCuP1Vis      = new G4VisAttributes(lgreen);
    logiPbP2Vis      = new G4VisAttributes(grey);
    logiCuP3Vis      = new G4VisAttributes(lgreen);

    logiCMOVis            = new G4VisAttributes(yellow);
    logiCopperframeVis    = new G4VisAttributes(brownl);
    logiCopperframeRodVis = new G4VisAttributes(brownl);
    logiReflectorVis      = new G4VisAttributes(orangel);

    logiPhotonDetVis = new G4VisAttributes(white);
    logiGeWaferVis   = new G4VisAttributes(blue);
    logiVacDiskVis   = new G4VisAttributes(orange);

    if (flagInvisible == 999) {
        logiPhotonDetVis->SetVisibility(true);
        logiGeWaferVis->SetVisibility(true);
        logiGeWaferVis->SetForceSolid(true);
        logiVacDiskVis->SetVisibility(true);
        logiVacDiskVis->SetForceSolid(true);
        logiPhotonDet->SetVisAttributes(logiPhotonDetVis);
        logiGeWafer->SetVisAttributes(logiGeWaferVis);
        logiVacDisk->SetVisAttributes(logiVacDiskVis);
    } else {
        logiPhotonDet->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiVacDisk->SetVisAttributes(G4VisAttributes::GetInvisible);
    }
    if (flagInvisible || flagOneCell) {
        logiRockShell->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiWorkArea->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiTopPbBox->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiPbBox->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiSSOVC->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu4->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu3->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu2->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu1->SetVisAttributes(G4VisAttributes::GetInvisible);

        if (flagInvisible == 999) {
            logiCMOCell->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiCopperframe->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiCopperframeRod->SetVisAttributes(G4VisAttributes::GetInvisible);

            logiCuMCP->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiCuP1->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiPbP2->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiCuP3->SetVisAttributes(G4VisAttributes::GetInvisible);
        }
    } else {
        logiCuP1Vis->SetForceSolid(true);
        logiPbP2Vis->SetForceSolid(true);
        logiCuP3Vis->SetForceSolid(true);
        logiCMOVis->SetVisibility(true);
        logiCMOVis->SetForceSolid(true);
        logiCopperframeVis->SetVisibility(true);
        logiCopperframeVis->SetForceSolid(true);
        logiCopperframeRodVis->SetVisibility(true);
        logiCopperframeRodVis->SetForceSolid(true);
        logiReflectorVis->SetVisibility(true);
        logiGeWaferVis->SetVisibility(true);
        logiGeWaferVis->SetForceSolid(true);

        logiRockShell->SetVisAttributes(logiRockShellVis);
        logiWorkArea->SetVisAttributes(logiWorkAreaVis);
        logiTopPbBox->SetVisAttributes(logiTopPbBoxVis);
        logiPbBox->SetVisAttributes(logiPbBoxVis);
        logiSSOVC->SetVisAttributes(logiSSOVCVis);
        logiCu4->SetVisAttributes(logiCu4Vis);
        logiCu3->SetVisAttributes(logiCu3Vis);
        logiCu2->SetVisAttributes(logiCu2Vis);
        logiCu1->SetVisAttributes(logiCu1Vis);
        logiCuMCP->SetVisAttributes(logiCuMCPVis);
        logiCuP1->SetVisAttributes(logiCuP1Vis);
        logiPbP2->SetVisAttributes(logiPbP2Vis);
        logiCuP3->SetVisAttributes(logiCuP3Vis);
        logiCMOCell->SetVisAttributes(logiCMOVis);
        logiCopperframe->SetVisAttributes(logiCopperframeVis);
        logiCopperframeRod->SetVisAttributes(logiCopperframeRodVis);
        logiGeWafer->SetVisAttributes(logiGeWaferVis);
        logiReflector->SetVisAttributes(logiReflectorVis);
    }

    //////////////////////////////
    // --- TG sensitive detector
    //////////////////////////////

    AmoreScintSD *TGSD;
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    TGSD = new AmoreScintSD(SDname = "/CupDet/TGSD", 35);
    SDman->AddNewDetector(TGSD);
    logiCMOCell->SetSensitiveDetector(TGSD);

    CupPMTSD *pmtSD;
    pmtSD = new CupPMTSD("/cupdet/pmt/MLCS", 35);
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD);
    logiGeWafer->SetSensitiveDetector(pmtSD);
    logiVacDisk->SetSensitiveDetector(pmtSD);

    new G4LogicalBorderSurface("MLCS_GeWafer_logsurf",
                               physGeWafer, // exiting glass into vac.
                               physVacDisk, GeWafer_opsurf);

    G4Region *PmtRegion = new G4Region("MLCS");
    PmtRegion->AddRootLogicalVolume(logiGeWafer);
    CupPMTOpticalModel *pmtOpticalModel = new CupPMTOpticalModel("MLCS_optical_model", physGeWafer);

  // EJ added
  //////////////////////////////
  // --- Logical Border Surfaces
  //////////////////////////////
  G4cout <<"Logical Border Surfaces between Reflector and Vacuum\n";

  new G4LogicalBorderSurface("MLCS_VacuumReflector_logsurf",
                       physReflectorVacuum,     // daughter
                       physReflector,           // mother
                       Vikuiti_opsurf);

    /////////////////////////////////////////////////////////////////
    //  Set Region
    //

    G4Region *crystalsRegion = new G4Region("crystals");
    crystalsRegion->AddRootLogicalVolume(logiCMOCell);
}
/*
void AmoreDetectorConstruction::PrintAMoRE10Coordinate
//void PrintAMoRE10Coordinate
        (char * Name, G4double xsize, G4double zsize,
        G4double xpos, G4double zpos, G4double height)
{
   printf (" %s\t\t%.1f\t%1.f\t%1.f\t%1.f\t%1.f\t%1.f\n",
           Name, (float) xsize, (float) zsize, (float) xpos, (float) zpos,
           (float) zpos + height, (float) zpos - height);
}
*/
