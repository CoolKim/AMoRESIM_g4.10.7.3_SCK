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
#include "CupSim/CupScintSD.hh"            // for making sensitive photocathodes
#include "CupSim/CupVetoSD.hh"             // for making sensitive photocathodes
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Element.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Version.hh"
#include "G4VisAttributes.hh"

#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"
//#include "CLHEP/Matrix/Matrix.h"

#include <fstream>
#include <string>
//#include <strstream>
#include <cmath>
#include <cstdio>
#include <stdlib.h>

#include "G4Region.hh"
#include "G4Types.hh"

#include "G4NistManager.hh"

// -- database
// == Construct Geometry for GenericLAND ======================================
// uses parameters from database or file
void AmoreDetectorConstruction::ConstructAMoREPilot() {
    // =-=-= CUT HERE =-=-=
    G4int flagInvisible = 0;

    ///////////////////////////////////////////////////////
    // -- make nested cylindrical tanks and fill volumes
    // make colours
    G4Colour white(1.0, 1.0, 1.0);
    G4Colour grey(0.5, 0.5, 0.5);
    G4Colour greyl(0.5, 0.5, 0.5, 0.3);
    G4Colour lgrey(.85, .85, .85);
    G4Colour red(1.0, 0.0, 0.0);
    G4Colour redl(1.0, .0, 0.3);
    G4Colour blue(0.0, 0.0, 1.0);
    G4Colour cyan(0.0, 1.0, 1.0);
    G4Colour cyanl(0.0, 1.0, 1.0, 0.5);
    G4Colour magenta(1.0, 0.0, 1.0);
    G4Colour yellow(1.0, 1.0, 0.0);
    G4Colour yellowl(1.0, 1.0, 0.0, 0.5);
    G4Colour orange(.75, .55, 0.0);
    G4Colour orangel(.75, .55, 0.0, 0.5);
    G4Colour lblue(0.0, 0.0, .75);
    G4Colour bluel(0.0, 0.0, 1.0, 0.5);
    G4Colour lgreen(0.0, .75, 0.0);
    G4Colour green(0.0, 1.0, 0.0);
    G4Colour green2(0.0, 1.0, 0.0, 0.5);
    G4Colour greenl(0.0, 1.0, 0.0, 0.3);
    G4Colour brown(0.7, 0.4, 0.1);
    G4Colour brown2(0.9, 0.6, 0.3);
    G4Colour brownl(0.7, 0.4, 0.1, 0.5);

    ///////////////////////////////////////////////
    // Primitive values
    ///////////////////////////////////////////////

    // WorldPV
    G4double bounding_size = 30. * meter / 2.0;

    // Rock Shell
    // Sizes: Orb, Radius: 20*meter
    G4double RockRadius = 10. * m;

    // Cavity in the Rock Shell
    //  : Air Room inside Rock Shell
    // Sizes: Hemisphere, Radius: 10*meter
    G4double RockCavityRadius  = 5.0 * m;
    G4double NEUT_CavityRadius = 3.5 * m + 5 * mm;

    // Rock Disk solid   5m x 5m x 5m
    G4double RockDiskThick_half = 1.0 * m;

    // Work Area   5m x 5m x 5m
    G4double workboxX = 3.7 * m / 2; // 3000. ;
    G4double workboxY = 2.8 * m / 2; // 3000. ;
    G4double workboxH = 4.5 * m / 2; // 4000. ;

    // Mounting structure of H-beam
    G4double Horizontal_beamX_length      = 3000. / 2;
    G4double Horizontal_beamXLarge_length = 3000. / 2;
    G4double Horizontal_beamY_length      = 1300. / 2;
    G4double Horizontal_beamYLarge_length = 1300. / 2;

    // H Beam points
    std::vector<G4TwoVector> large_points; // Large beam points.
    large_points.push_back(G4TwoVector(100.0, 100.0));
    large_points.push_back(G4TwoVector(100.0, -100.0));
    large_points.push_back(G4TwoVector(88.0, -100.0));
    large_points.push_back(G4TwoVector(88.0, -4.0));
    large_points.push_back(G4TwoVector(-88.0, -4.0));
    large_points.push_back(G4TwoVector(-88.0, -100.0));
    large_points.push_back(G4TwoVector(-100.0, -100.0));
    large_points.push_back(G4TwoVector(-100.0, 100.0));
    large_points.push_back(G4TwoVector(-88.0, 100.0));
    large_points.push_back(G4TwoVector(-88.0, 4.0));
    large_points.push_back(G4TwoVector(88.0, 4.0));
    large_points.push_back(G4TwoVector(88.0, 100.0));

    std::vector<G4TwoVector> small_points; // Small beam points.
    small_points.push_back(G4TwoVector(75.0, 75.0));
    small_points.push_back(G4TwoVector(75.0, -75.0));
    small_points.push_back(G4TwoVector(65.0, -75.0));
    small_points.push_back(G4TwoVector(65.0, -3.5));
    small_points.push_back(G4TwoVector(-65.0, -3.5));
    small_points.push_back(G4TwoVector(-65.0, -75.0));
    small_points.push_back(G4TwoVector(-75.0, -75.0));
    small_points.push_back(G4TwoVector(-75.0, 75.0));
    small_points.push_back(G4TwoVector(-65.0, 75.0));
    small_points.push_back(G4TwoVector(-65.0, 3.5));
    small_points.push_back(G4TwoVector(65.0, 3.5));
    small_points.push_back(G4TwoVector(65.0, 75.0));

    G4double Bplate_size       = 300. / 2;
    G4double Bplate_thick_half = 40. / 2.;
    G4double Hbeam1_length     = 3460.1 / 2;
    G4double Hbeam2_length     = 800. / 2;

    double x_beam[4] = {1.0, -1.0, -1.0, 1.0};
    double y_beam[4] = {1.0, 1.0, -1.0, -1.0};

    double HoLayer[3] = {1000., 2925., 3480.1};
    double HoLine[2]  = {1, -1};

    double HoLayerY[3] = {1000., 2925., 3480.1};
    double HoLineY[2]  = {1, -1};

    // PMT
    G4double AMoRE_unit_inch = 2.54 * cm;

    G4double PMTBacksideThickness = 1 * mm;
    G4double PMTTubeThickness     = 2 * mm;
    G4double PMTWindowThickness   = 4 * mm;
    G4double PMTSizeR             = 51.5 * mm / 2;
    G4double PMTSizeZ             = 109 * mm;
    G4double PMTTopSizeRmin       = 0;
    G4double PMTVacuTopSizeRmin   = 0;
    G4double PMTGreaseShieldThick = 5 * um;
    G4double PMTGreaseThickness   = 2 * mm;

    // Scintillator
    G4double scint_Thick_half        = 50.000 * mm / 2.;
    G4double scint_FlatL_ratio       = 0.7;
    G4double scint_reflector_thick   = 1 * mm;
    G4double scint_boolean_tolerence = 5 * um;

    G4double topScint_trapZ_half     = 368.855 * mm / 2.;
    G4double topScint_trapX1_half    = 37.000 * mm / 2.;
    G4double topScint_trapX2_half    = 381.000 * mm / 2.;
    G4double topScint_boxHeight_half = 1725.000 * mm / 2.;
    G4double topScint_boxWidth_half  = 762.063 * mm / 2.;

    G4double sideFBScint_trapZ_half     = 250.000 * mm / 2.;
    G4double sideFBScint_trapX1_half    = 37.000 * mm / 2.;
    G4double sideFBScint_trapX2_half    = 275.000 * mm / 2.;
    G4double sideFBScint_boxHeight_half = 1700.000 * mm / 2.;
    G4double sideFBScint_boxWidth_half  = 550.000 * mm / 2.;

    G4double sideLRScint_trapZ_half     = 500.350 * mm / 2.;
    G4double sideLRScint_trapX1_half    = 46.699 * mm / 2.;
    G4double sideLRScint_trapX2_half    = 513.333 * mm / 2.;
    G4double sideLRScint_boxHeight_half = 1700.000 * mm / 2.;
    G4double sideLRScint_boxWidth_half  = 600.000 * mm / 2.;

    //  Pb Top Box
    G4double pbtopbox_size              = 150. * cm / 2.0;
    G4double pbtopbox_zsize             = 15. * cm / 2.0;
    G4double pbtopbox_housing_thickness = 3 * mm;
    G4double lifting_PbTopBox           = 22.0 * cm;

    //  Pb  Box
    G4double pbbox_thick             = 15. * cm;
    G4double pbbox_housing_thickness = 3 * mm;
    G4double pbbox_XYinnerlength     = 80. * cm;
    G4double pbbox_Zinnerlength      = 162.1 * cm;
    G4double vgaptoMuVeto            = 885 * mm;

    // Layer5_OVC Stainless Steel (SS)
    G4double ss_radius      = 64. * cm / 2;
    G4double ss_thick       = 0.5 * cm;
    G4double ss_height      = 153.8 * cm / 2; // 153.8
    G4double sstopbox_xsize = 74.8 * cm / 2.0;
    G4double sstopbox_ysize = 74.8 * cm / 2.0;
    G4double sst_zsize      = 3. * cm / 2.;
    G4double ssb_zsize      = 2. * cm / 2;
    G4double plateGap       = 12. * cm;

    //  Mu Metal shield
    G4double mu_radius = 65.5 * cm / 2;
    G4double mu_thick  = 0.1 * cm;
    G4double mu_height = 156. * cm / 2; // 153.8
    G4double mub_zsize = 0.1 * cm / 2;

    // Layer4_ 50K-SHIELD Copper
    G4double cu4_radius       = 57.6 * cm / 2;
    G4double cu4_thick        = 0.3 * cm;
    G4double cu4_inner_height = 134.2 * cm / 2; // 136.9-1.0-1.7
    G4double cu4t_zsize       = 1.7 * cm / 2;
    G4double cu4b_zsize       = 1.0 * cm / 2;
    G4double cu4_gap_fromTop  = 120 * mm;

    // Layer3_Cu IVC
    G4double cu3_radius       = 51.4 * cm / 2;
    G4double cu3_thick        = 0.5 * cm;
    G4double cu3_inner_height = 116.18 * cm / 2; // 119.88 -1.8 - 1.9
    G4double cu3t_zsize       = 1.8 * cm / 2;
    G4double cu3b_zsize       = 1.9 * cm / 2;
    G4double cu3_gap_fromTop  = 120 * mm;

    // Layer2_Cu SHIELD-STILL
    G4double cu2_radius       = 45.5 * cm / 2;
    G4double cu2_thick        = 0.1 * cm;
    G4double cu2_inner_height = 100.2 * cm / 2;
    G4double cu2t_zsize       = 2.0 * cm / 2;
    G4double cu2b_zsize       = 0.1 * cm / 2;
    G4double cu2_gap_fromTop  = 120 * mm;

    // Layer1_Cu 50mK-SHIELD
    G4double cu1_radius       = 43.5 * cm / 2;
    G4double cu1_thick        = 0.1 * cm;
    G4double cu1_inner_height = 85.1 * cm / 2;
    G4double cu1t_zsize       = 2.0 * cm / 2;
    G4double cu1b_zsize       = 0.1 * cm / 2;
    G4double cu1_gap_fromTop  = 120 * mm;

    // Cu Plate: Mixing Chamber SHIELD
    G4double cumcp_radius = 40.8 * cm / 2;
    G4double cumcp_height = 3. * cm / 2;

    // Plate1_Cu
    G4double cup1_zsize = 0.5 * cm / 2;
    G4double plateGap2  = 1.0 * cm;

    // Plate2_Pb Shield
    G4double pbp2_zsize = 5. * cm / 2;

    // Plate3_Cu
    G4double cup3_zsize = 1. * cm / 2;

    // superconducting magnetic shielding tube
    G4double SMS_radius         = 15.9 * cm / 2;
    G4double SMS_height         = 40. * cm / 2; // Nb shield size needs to be confirmed
    G4double SMS_thickness      = 2. * mm;
    G4double SMS_gap_fromBottom = 41. * mm;

    ///////////////////////////////////////////////////////
    //  Target Room
    G4double TR_radius         = (15.4) * cm / 2;
    G4double TR_height         = (40. - 0.2) * cm / 2; // Nb shield size needs to be confirmed
    G4double TR_gap_fromBottom = 0. * mm;

    // CMO cell
    G4double cmocell_sr[7] = {4.30 * cm / 2., 4.09 * cm / 2,  4.79 * cm / 2, 4.27 * cm / 2,
                              4.37 * cm / 2,  4.043 * cm / 2, 4.5 * cm / 2};
    G4double cmocell_lr[7] = {5.25 * cm / 2., 4.60 * cm / 2,  5.38 * cm / 2, 4.68 * cm / 2,
                              5.14 * cm / 2,  4.835 * cm / 2, 4.5 * cm / 2};
    G4double cmocell_h[7]  = {
        2.6 * cm / 2, 4.05 * cm / 2, 4.15 * cm / 2, 5.2 * cm / 2,
        5.1 * cm / 2, 5.1 * cm / 2,  3. * cm / 2}; /// Heights of SB28, NSB29, S35, SS68,
    /// and 4.5
    /// The fifth height need to be updated to SE01's height
    /// The sixth height is dummy for the sixth Photon Detector
    G4double CMOSupportGap = 0.15 * cm;
    G4double CMOSupport_H  = 0.5 * cm / 2.;

    // Photon Detector: Ge Wafer + Vaccum disk
    G4double PhotonDet_r     = 2.54 * 2 * cm / 2. + 0.05 * cm;
    G4double PhotonDet_zsize = 0.6 * cm / 2.;
    G4double PhotonFrame_x   = 5.7 * cm / 2;
    G4double PhotonFrame_y   = 7.38 * cm / 2;
    G4double PhotonFrame_z   = 6 / 2;

    // Ge wafer on the bottom of the Crystal
    G4double GeWafer_r       = 2.54 * 2 * cm / 2.;
    G4double GeWafer_zsize   = 0.004 * cm / 2.;           // 300 nm + 100 nm Vac
    G4ThreeVector PosGeWafer = G4ThreeVector(0., 0., 0.); // Center of PhotonDet

    // Vacuum Disks on Ge Wafer
    G4double VacDisk_zsize = 0.00001 * cm / 2.; // 100 nm

    // Top Assembly
    G4double Assem_px = 0, Assem_py = 0;
    G4double top_pz      = 47.;
    G4double pillarH     = 40.;
    G4int kk             = 0;
    double x_smallbox[4] = {1.0, -1.0, -1.0, 1.0};
    double y_smallbox[4] = {1.0, 1.0, -1.0, -1.0};

    // CopperFrame
    G4double CopperFx = 7.7 * cm;
    G4double CopperFy = 7.4 * cm;
    G4double CopperFz = 4.64 * mm;

    // Copper Rod
    G4double copperRods_r = 0.7 * cm / 2.;

    // Reflector around the CMO (Vm2000 shield)
    G4double reflector_thick  = 0.0064 * cm; // 64 micro-meter
    G4double reflector_thick2 = 0.1 * cm;

    // Polyethylene(+Borated) geometry
    G4double toyPE_thickness         = 30. * cm;
    G4double realPE_thickness1       = 10. * cm;
    G4double realPE_thickness2       = 30. * cm;
    G4double realPE_LRCuttingLength1 = 60. * cm;
    G4double realPE_LRCuttingLength2 = 13. * cm;
    G4double boratedPE_thickness     = 2.5 * cm;
    G4double boratedPE_BottomGap     = 19. * cm;

    // Boric acid panel geometry
    G4double thickness_BoricAcid        = 3. * mm;
    G4double thickness_BoricAcidHousing = 1. * mm;
    G4double thickness_BoricAcidSpacing = 10. * mm;

    // G10 material holders between SS and 50K-SHIELD Copper
    G4double Cuholder_R     = 2.540 * cm / 2; // outer radius
    G4double G10holderXY[5] = {162.5 * mm, 190 * mm, 190 * mm, 230 * mm, 260 * mm};

    // BOLTS in Copper Frame
    // SSCB BOLTS: M4-12
    G4double dBoltM4Top = 8.22 * mm / 2;
    G4double zBoltM4Top = 5.0 * mm / 2;
    G4double dBoltM4    = 9.2 * mm / 2;
    G4double zBoltM4    = 4.0 * mm / 2;

    // SSCB BOLTS: M3-6
    G4double dBoltM3 = 6.5 * mm / 2;
    G4double zBoltM3 = 3.0 * mm / 2;

    // SSCB BOLTS: M3-4
    G4double dBoltM3F = 5.85 * mm / 2;
    G4double zBoltM3F = 3.0 * mm / 2;

    G4double dBoltM3S = 6.96 * mm / 2;
    G4double zBoltM3S = 3.0 * mm / 2;

    // additional Bolts : AB
    // beside of pilar
    G4double dBoltAB = 7.525 * mm / 2;
    G4double zBoltAB = 3.5 * mm / 2;

    // GOLD FILM
    // a gold film below crystal
    G4double dGOLDb = 2 * cm / 2;
    G4double zGOLDb = 0.0003 * mm / 2;

    // 3 gold film above crystal
    G4double dGOLDa = 5 * mm / 2;
    G4double zGOLDa = 0.0003 * mm / 2;

    // B4C rubber shield
    G4double rubber_thickness      = 4. * mm;
    G4double rubber_UnderPb_XYsize = 100. * cm;
    G4double rubber_GapFromSSOVC   = 10. * cm;

    ///////////////////////////////////////////////
    // Composite values
    ///////////////////////////////////////////////

    //  Mounting structure of H-beam
    G4double posZ_Bplate1    = Hbeam1_length + Bplate_thick_half;
    G4double posZ_Bplate2    = Hbeam2_length + Bplate_thick_half;
    G4double posMountOrigin1 = -workboxH + Hbeam1_length + 2. * Bplate_thick_half;
    G4double posMountOrigin2 = -workboxH + Hbeam2_length + 2. * Bplate_thick_half - 1 * mm;
    G4ThreeVector posBplate1 = G4ThreeVector(0., 0., posZ_Bplate1);
    G4ThreeVector posBplate2 = G4ThreeVector(0., 0., posZ_Bplate2 - 1 * mm);

    /// Vertical H beam Rotatrion Matrix
    G4RotationMatrix *VeBeamRotation = new G4RotationMatrix();
    VeBeamRotation->rotateX(0. * deg);
    VeBeamRotation->rotateY(0. * deg);

    /// Horizontal H beam Rotation Matrix X
    G4RotationMatrix *HoBXRotation = new G4RotationMatrix();
    HoBXRotation->rotateX(90. * deg);
    HoBXRotation->rotateZ(90. * deg);
    HoBXRotation->rotateX(90. * deg);

    /// Horizontal H beam Rotation Matrix Y
    G4RotationMatrix *HoBYRotation = new G4RotationMatrix();
    HoBYRotation->rotateY(90. * deg);
    HoBYRotation->rotateX(90. * deg);

    // PMT
    G4double PMTWindowSizeR     = AMoRE_unit_inch * 3. / 2.;
    G4double PMTTopSizeRmax     = std::max(PMTWindowSizeR, PMTSizeR);
    G4double PMTVacuSizeR       = PMTSizeR - PMTTubeThickness;
    G4double PMTVacuSizeZ       = PMTSizeZ - PMTBacksideThickness;
    G4double PMTVacuTopSizeRmax = PMTWindowSizeR - PMTTubeThickness;
    G4double PMTGreaseSizeR     = PMTWindowSizeR - PMTGreaseShieldThick;
    G4Material *material_Grease = _grease;

    // Scintillator
    G4double topScint_trapFlatL_half  = topScint_trapZ_half * scint_FlatL_ratio;
    G4double topScint_trapThick1_half = scint_Thick_half;
    G4double topScint_trapThick2_half = topScint_trapX1_half;
    G4double topScint_boxThick_half   = scint_Thick_half;
    G4double topScint_TrapSlope_double =
        (topScint_trapX2_half - topScint_trapX1_half) / topScint_trapZ_half;

    G4double sideFBScint_trapFlatL_half  = sideFBScint_trapZ_half * scint_FlatL_ratio;
    G4double sideFBScint_trapThick1_half = scint_Thick_half;
    G4double sideFBScint_trapThick2_half = sideFBScint_trapX1_half;
    G4double sideFBScint_boxThick_half   = scint_Thick_half;
    G4double sideFBScint_TrapSlope_double =
        (sideFBScint_trapX2_half - sideFBScint_trapX1_half) / sideFBScint_trapZ_half;

    G4double sideLRScint_trapFlatL_half  = sideLRScint_trapZ_half * scint_FlatL_ratio;
    G4double sideLRScint_trapThick1_half = scint_Thick_half;
    G4double sideLRScint_trapThick2_half = sideLRScint_trapX1_half;
    G4double sideLRScint_boxThick_half   = scint_Thick_half;
    G4double sideLRScint_TrapSlope_double =
        (sideLRScint_trapX2_half - sideLRScint_trapX1_half) / sideLRScint_trapZ_half;

    G4RotationMatrix *scint_alignMtx_Uside = new G4RotationMatrix;
    G4RotationMatrix *scint_alignMtx_Dside = new G4RotationMatrix;
    G4RotationMatrix *scint_alignMtx_Rside = new G4RotationMatrix;

    G4RotationMatrix *scintMother_alignMtx_FBSide  = new G4RotationMatrix();
    G4RotationMatrix *scintMother_alignMtx_LRSide1 = new G4RotationMatrix();
    G4RotationMatrix *scintMother_alignMtx_LRSide2 = new G4RotationMatrix();

    scint_alignMtx_Uside->rotateX(-90 * deg);
    scint_alignMtx_Dside->rotateX(90 * deg);
    scint_alignMtx_Rside->rotateX(90 * deg);
    scint_alignMtx_Rside->rotateY(-90 * deg);

    scintMother_alignMtx_FBSide->rotateY(-90 * deg);
    scintMother_alignMtx_FBSide->rotateZ(-90 * deg);

    scintMother_alignMtx_LRSide1->rotateX(-90 * deg);

    scintMother_alignMtx_LRSide2->rotateX(-90 * deg);
    scintMother_alignMtx_LRSide2->rotateY(-180 * deg);

    //  Pb Top Box
    G4double CenterPbTZ = Hbeam1_length + posMountOrigin1 + pbtopbox_zsize + 1 + lifting_PbTopBox;

    //  Pb  Box
    G4double pbbox_size       = (pbbox_XYinnerlength + 2 * pbbox_thick) / 2.0;
    G4double pbbox_zsize      = (pbbox_Zinnerlength + pbbox_thick) / 2.0;
    G4double vgaptoTop        = -pbtopbox_zsize + vgaptoMuVeto;
    G4double PbboxTopZ        = CenterPbTZ;
    G4double CenterPbZ        = PbboxTopZ - pbbox_zsize - vgaptoTop;
    G4ThreeVector PosPbBox    = G4ThreeVector(0., 0., CenterPbZ);
    G4ThreeVector PosInnerDet = PosPbBox;

    // SSOVC
    G4ThreeVector SSTopShift = G4ThreeVector(0., 0., ss_height + ssb_zsize + sst_zsize);

    //  Mu Metal shield
    G4ThreeVector Posmu         = G4ThreeVector(0., 0., pbbox_zsize - mu_height);
    G4ThreeVector muBottomShift = G4ThreeVector(0., 0., -mu_height - mub_zsize);

    // Layer4_ 50K-SHIELD Copper
    G4double cu4_total_height_half = cu4_inner_height + cu4b_zsize + cu4t_zsize;
    G4ThreeVector PosCu4Outer =
        G4ThreeVector(0., 0., ss_height - cu4_total_height_half - cu4_gap_fromTop);
    G4ThreeVector PosCu4Inner = G4ThreeVector(0., 0., cu4b_zsize - cu4t_zsize);

    // Layer3_Cu IVC
    G4double cu3_total_height_half = cu3_inner_height + cu3b_zsize + cu3t_zsize;
    G4ThreeVector PosCu3Outer(0., 0., cu4_inner_height - cu3_total_height_half - cu3_gap_fromTop);
    G4ThreeVector PosCu3Inner(0, 0, cu3b_zsize - cu3t_zsize);

    // Layer2_Cu SHIELD-STILL
    G4double cu2_total_height_half = cu2_inner_height + cu2b_zsize + cu2t_zsize;
    G4ThreeVector PosCu2Outer =
        G4ThreeVector(0., 0., cu3_inner_height - cu2_total_height_half - cu2_gap_fromTop);
    G4ThreeVector PosCu2Inner = G4ThreeVector(0., 0., cu2b_zsize - cu2t_zsize);

    // Layer1_Cu 50mK-SHIELD
    G4double cu1_total_height_half = cu1_inner_height + cu1b_zsize + cu1t_zsize;
    G4ThreeVector PosCu1Outer(0, 0, cu2_inner_height - cu1_total_height_half - cu1_gap_fromTop);
    G4ThreeVector PosCu1Inner(0, 0, cu1b_zsize - cu1t_zsize);

    // Cu Plate: Mixing Chamber SHIELD
    G4double CenterCuMCPZ  = cu1_inner_height - cumcp_height - plateGap;
    G4ThreeVector PosCuMCP = G4ThreeVector(0., 0., CenterCuMCPZ);

    // Plate1_Cu
    G4double CenterCuP1Z  = CenterCuMCPZ - cumcp_height - cup1_zsize - plateGap2;
    G4ThreeVector PosCuP1 = G4ThreeVector(0., 0., CenterCuP1Z);

    // Plate2_Pb Shield
    G4double CenterPbP2Z  = CenterCuP1Z - cup1_zsize - pbp2_zsize;
    G4ThreeVector PosPbP2 = G4ThreeVector(0., 0., CenterPbP2Z);

    // Plate3_Cu
    G4double CenterCuP3Z  = CenterPbP2Z - pbp2_zsize - cup3_zsize;
    G4ThreeVector PosCuP3 = G4ThreeVector(0., 0., CenterCuP3Z);

    // Plate4_Cu
    G4double CenterCuP4Z  = CenterCuP3Z - cup3_zsize - cup1_zsize - plateGap2;
    G4ThreeVector PosCuP4 = G4ThreeVector(0., 0., CenterCuP4Z);

    // Plate5_Pb Shield
    G4double CenterPbP5Z  = CenterCuP4Z - cup1_zsize - pbp2_zsize;
    G4ThreeVector PosPbP5 = G4ThreeVector(0., 0., CenterPbP5Z);

    // Plate6_Cu
    G4double CenterCuP6Z  = CenterPbP5Z - pbp2_zsize - cup3_zsize;
    G4ThreeVector PosCuP6 = G4ThreeVector(0., 0., CenterCuP6Z);

    // superconducting magnetic shielding tube
    G4ThreeVector PosSMS =
        G4ThreeVector(0., 0., -cu1_inner_height + SMS_height + SMS_gap_fromBottom);

    // Target Room
    G4ThreeVector PosTR = G4ThreeVector(0., 0.,
                                        -cu1_inner_height + TR_height + SMS_gap_fromBottom +
                                            SMS_thickness + TR_gap_fromBottom);

    // Reflector around the CMO (Vm2000 shield)
    G4double reflector_x    = CopperFx / 2. - 1.15 * cm;     // 64 micro-meter
    G4double reflector_y    = CopperFy / 2 - 1.1 * cm;       // 64 micro-meter
    G4double reflector_inx  = reflector_x - reflector_thick; // 64 micro-meter
    G4double reflector_iny  = reflector_y - reflector_thick; // 64 micro-meter
    G4double reflector_inx2 = reflector_x - reflector_thick2;
    G4double reflector_iny2 = reflector_y - reflector_thick2;

    // Boric acid panel
    G4double totalThickness_BoricAcidPanel = thickness_BoricAcid + thickness_BoricAcidHousing * 2.;

    G4double height_LR = pbbox_zsize - pbbox_thick / 2.;
    G4double width_LR  = pbbox_size - pbbox_thick;

    G4double height_FB = pbbox_zsize - pbbox_thick / 2.;
    G4double width_FB  = pbbox_size - pbbox_thick - totalThickness_BoricAcidPanel;

    G4double length_Bottom = width_LR - totalThickness_BoricAcidPanel;
    G4double width_Bottom  = width_FB - totalThickness_BoricAcidPanel;

    G4double length_Top = width_LR;
    G4double width_Top  = width_FB;

    G4ThreeVector sideLRBoricAcid_Pos[2];
    sideLRBoricAcid_Pos[0] =
        G4ThreeVector(0, pbbox_size - pbbox_thick - totalThickness_BoricAcidPanel / 2.,
                      -pbbox_zsize + pbbox_thick + height_LR);
    sideLRBoricAcid_Pos[1] =
        G4ThreeVector(0, -pbbox_size + pbbox_thick + totalThickness_BoricAcidPanel / 2.,
                      -pbbox_zsize + pbbox_thick + height_LR);

    G4ThreeVector sideFBBoricAcid_Pos[2];
    sideFBBoricAcid_Pos[0] =
        G4ThreeVector(pbbox_size - pbbox_thick - totalThickness_BoricAcidPanel / 2., 0,
                      -pbbox_zsize + pbbox_thick + height_FB);
    sideFBBoricAcid_Pos[1] =
        G4ThreeVector(-pbbox_size + pbbox_thick + totalThickness_BoricAcidPanel / 2., 0,
                      -pbbox_zsize + pbbox_thick + height_FB);

    G4ThreeVector bottomBoricAcid_Pos =
        G4ThreeVector(0, 0, -pbbox_zsize + pbbox_thick + totalThickness_BoricAcidPanel / 2.);
    G4ThreeVector topBoricAcid_Pos =
        G4ThreeVector(0, 0, +pbbox_zsize + sst_zsize * 2. + totalThickness_BoricAcidPanel / 2.);

    // G10 material holders between SS and 50K-SHIELD Copper
    G4double Cuholder_r = Cuholder_R - 1.9038 * mm;
    G4double CuholderH  = plateGap; // height

    G4double xG10[6][6];
    G4double yG10[6][6];
    G4double CenterPlateGap4VolZ[5] = {cu1_inner_height - CuholderH / 2,
                                       cu2_inner_height - CuholderH / 2,
                                       cu3_inner_height - CuholderH / 2,
                                       cu4_inner_height - CuholderH / 2, ss_height - CuholderH / 2};

    G4double PlateGapR[5] = {cu1_radius - cu1_thick, cu2_radius - cu2_thick, cu3_radius - cu3_thick,
                             cu4_radius - cu4_thick, ss_radius - ss_thick};
    ///////////////////////////////////////////////
    // Variables
    ///////////////////////////////////////////////

    // * Define a world physical volume.
    // * Set the class variable "world_phys" to point to the world phys volume.
    // -- Make a 20m x 20m x 20m "Experimental Hall" for the world volume
    G4Box *boxHall;
    G4LogicalVolume *logiHall;
    G4VisAttributes *logiHallVis;
    G4VPhysicalVolume *physHall;

    /////////////////////////////////////////////////////////////////
    // Rock Shell
    G4Sphere *RockSolid;
    G4LogicalVolume *logiRock;
    G4VisAttributes *logiRockVis;

    /////////////////////////////////////////////////////////////////
    // Cavity in the Rock Shell
    //  : Air Room inside Rock Shell
    /////////////////////////////////////////////////////////////////
    G4Tubs *RockCavitySolid;
    G4LogicalVolume *logiRockCavity;
    G4VisAttributes *logiRockCavityVis;

    /////////////////////////////////////////
    // Rock Disk solid   5m x 5m x 5m
    /////////////////////////////////////////
    G4Tubs *RockDiskSolid;
    G4LogicalVolume *logiRockDisk;
    G4VisAttributes *logiRockDiskVis;

    /////////////////////////////////////////
    // Work Area   5m x 5m x 5m
    /////////////////////////////////////////
    G4Box *WorkArea;
    G4LogicalVolume *logiWorkArea;
    G4VisAttributes *logiWorkAreaVis;

    ////////////////////////////////////////////////////////////
    //  Mounting structure of H-beam properties & variables
    ////////////////////////////////////////////////////////////
    G4Box *Bplate;
    G4UnionSolid *UnionHbeam1;
    G4UnionSolid *UnionHbeam2;
    G4ExtrudedSolid *Hbeam1;
    G4ExtrudedSolid *Hbeam2;
    G4ExtrudedSolid *Ho_beamX;
    G4ExtrudedSolid *Ho_beamXLar;
    G4ExtrudedSolid *Ho_beamY;
    G4ExtrudedSolid *Ho_beamYLar;

    G4LogicalVolume *logicHbeam1;
    G4LogicalVolume *logicHbeam2;
    G4LogicalVolume *logicHo_beamX;
    G4LogicalVolume *logicHo_beamXLar;
    G4LogicalVolume *logicHo_beamY;
    G4LogicalVolume *logicHo_beamYLar;

    G4VisAttributes *Hb1VisAtt;
    G4VisAttributes *Hb2VisAtt;
    G4VisAttributes *HoBXVisAtt;
    G4VisAttributes *HoBXLVisAtt;
    G4VisAttributes *HoBYLVisAtt;
    G4VisAttributes *HoBYVisAtt;

    /////////////////////////////////////////
    // G4double PbTBoxVSize = 175.45 * cm; //  114.05 cm (center to top of Stainless Steel cover)
    //   + 46.4 cm (gap)
    //   + 15.0 cm (Lead top thickmess)

    //  Pb Top Box
    G4Box *TopPbBox;
    G4Box *TopPbBlock;
    G4LogicalVolume *logiTopPbBox;
    G4VPhysicalVolume *physTopPbBox1 = nullptr;
    G4LogicalVolume *logiTopPbBlock;

    //  Pb  Box
    G4Box *PbBoxOut;
    G4Box *PbBoxIn;
    G4Box *PbBlockOut;
    G4Box *PbBlockIn;
    G4SubtractionSolid *PbBox;
    G4SubtractionSolid *PbBlock;

    G4LogicalVolume *logiPbBox;
    G4LogicalVolume *logiPbBlock;
    G4VisAttributes *logiPbBoxVis;

    ////////////////
    // Inner detector
    ////////////////

    G4Box *InnerDetectorBox;
    G4Box *SSTopBoxForInnerDet;
    G4VSolid *InnerDetector;
    G4LogicalVolume *logiInnerDetector;
    G4VisAttributes *logiInnerDetectorVis;

    G4ThreeVector PosSS;

    G4Tubs *SSCylinder;
    G4Box *SSTop;
    G4VSolid *SSOVC0;
    G4Tubs *SSOVCInnerSolid;
    G4LogicalVolume *logiSSOVCOuter;
    G4LogicalVolume *logiSSOVCInner;
    G4VisAttributes *logiSSOVCVis;

    ///////////////////////////////////////////////////////////////////
    //  Mu Metal shield

    G4Tubs *muPipe;
    G4LogicalVolume *logimu;
    G4VisAttributes *logimuVis;

    /////////////////////////////////////////////////////////////////
    // Layer4_ 50K-SHIELD Copper

    G4Tubs *Cu4OuterCylinder;
    G4Tubs *Cu4InnerCylinder;
    G4LogicalVolume *logiCu4Outer;
    G4LogicalVolume *logiCu4Inner;
    G4VisAttributes *logiCu4Vis;

    ///////////////////////////////////////////////////////
    // Layer3_Cu IVC

    G4Tubs *Cu3OuterCylinder;
    G4Tubs *Cu3InnerCylinder;
    G4LogicalVolume *logiCu3Outer;
    G4LogicalVolume *logiCu3Inner;
    G4VisAttributes *logiCu3Vis;

    ///////////////////////////////////////////////////////
    // Layer2_Cu SHIELD-STILL

    G4Tubs *Cu2OuterCylinder;
    G4Tubs *Cu2InnerCylinder;
    G4LogicalVolume *logiCu2Outer;
    G4LogicalVolume *logiCu2Inner;
    G4VisAttributes *logiCu2Vis;

    ///////////////////////////////////////////////////////
    // Layer1_Cu 50mK-SHIELD

    G4Tubs *Cu1OuterCylinder;
    G4Tubs *Cu1InnerCylinder;
    G4LogicalVolume *logiCu1Outer;
    G4LogicalVolume *logiCu1Inner;
    G4VisAttributes *logiCu1Vis;

    ///////////////////////////////////////////////////////
    // Cu Plate: Mixing Chamber SHIELD

    G4Tubs *CuMCPlate;
    G4LogicalVolume *logiCuMCP;
    G4VisAttributes *logiCuMCPVis;

    ///////////////////////////////////////////////////////
    // Plate1_Cu

    G4Tubs *CuPlate1;
    G4LogicalVolume *logiCuP1;
    G4VisAttributes *logiCuP1Vis;

    ///////////////////////////////////////////////////////
    // Plate2_Pb Shield

    G4Tubs *PbPlate2;
    G4LogicalVolume *logiPbP2;
    G4VisAttributes *logiPbP2Vis;

    ///////////////////////////////////////////////////////
    // Plate3_Cu

    G4Tubs *CuPlate3;
    G4LogicalVolume *logiCuP3;
    G4VisAttributes *logiCuP3Vis;

    ///////////////////////////////////////////////////////
    // Plate4_Cu

    ///////////////////////////////////////////////////////
    // Plate5_Pb Shield

    ///////////////////////////////////////////////////////
    // Plate6_Cu

    ///////////////////////////////////////////////////////
    // G10
    G4ThreeVector PosPlateGap4Vol[5];
    G4LogicalVolume *PlateGap4VolLogi[5];

    G4ThreeVector PosG10holder;

    //////////////////////////////////////////////////////
    // superconducting magnetic shielding tube

    G4Tubs *SuperMS0;
    G4Tubs *SuperMS1;
    G4VSolid *SuperMS;
    G4LogicalVolume *logiSuperMS;
    G4VisAttributes *logiSuperMSVis;

    //  Target Room
    G4Tubs *TargetRoom;
    G4LogicalVolume *logiTargetRoom;
    G4VisAttributes *logiTargetRoomVis;

    ///////////////////////////////////////////////////////
    // CMO cell
    ///////////////////////////////////////////////////////
    char namecmo[20];
    char namecuframe[20];

    G4EllipticalTube *CMOCell[7];
    G4LogicalVolume *logiCMOCell[7];
    G4VisAttributes *logiCMOVis[7];

    G4double xCell[8] = {0};
    G4double yCell[8] = {0};
    G4double zCell[8] = {0};
    G4double zCup[8];
    G4double zCdown[8];
    G4int copyNo;

    G4ThreeVector PosCMOCell;

    ///////////////////////////////////////////////////////
    // Photon Detector: Ge Wafer + Vaccum disk

    G4double PhotonDetZ;
    G4ThreeVector PosPhotonDet;

    G4Tubs *PhotonDet;
    G4LogicalVolume *logiPhotonDet;
    G4VisAttributes *logiPhotonDetVis;

    G4Box *PhotonFramePlate;
    G4VSolid *PhotonFrame;
    G4VSolid *PhotonFrameX;
    G4VSolid *PhotonFrameY;
    G4VSolid *PhotonFrameW;
    G4Box *sol_CopPlatehole2;
    G4VSolid *PhotonFrame2;
    G4LogicalVolume *lv_PhotonFrame;
    G4VisAttributes *lv_PhotonFrameVis;

    //////////////////////////////////////////////////////
    // staycaste & PinConnector
    char namePin[20];

    ///////////////////////////////////////////////////////
    // Ge wafer on the bottom of the Crystal

    G4Tubs *GeWafer;
    G4LogicalVolume *logiGeWafer;
    G4VisAttributes *logiGeWaferVis;

    ///////////////////////////////////////////////////////
    // Vacuum Disks on Ge Wafer
    G4double VacDisk_r = GeWafer_r;
    G4double VacDiskZ  = -GeWafer_zsize + VacDisk_zsize;

    G4Tubs *VacDisk;
    G4LogicalVolume *logiVacDisk;
    G4VisAttributes *logiVacDiskVis;
    G4ThreeVector PosVacDisk;

    ///////////////////////////////////////////////////////
    // Top Assembly

    G4VSolid *TopCopper[4];
    G4VSolid *Top_co;
    G4VSolid *smallbox;
    G4VSolid *smallboxB;
    G4Tubs *hole;
    G4Tubs *pillar;
    G4LogicalVolume *lv_smallbox;
    G4LogicalVolume *lv_smallboxB;
    G4LogicalVolume *lv_copillar;

    G4ThreeVector PosCopperframeUp;
    G4ThreeVector PosCopperframeDown;

    G4LogicalVolume *logiSquareDiskUp[7];
    G4VisAttributes *logiSquareDiskUpVis;

    /////CopperFrame by Hoon /////////////////////////

    G4VSolid *sol_CopperFrame[4];
    G4VSolid *sol_CopperFrameb[4];
    G4VSolid *sol_CopperFrame2[4];
    G4VSolid *sol_CopperFrameb2[4];
    G4VSolid *sol_CopperFrame3[4];
    G4VSolid *sol_CopperFrameb3[4];
    G4VSolid *sol_CopperFrame4[4];
    G4VSolid *sol_CopperFrame5[4];
    G4Tubs *sol_Cophole;
    G4Box *sol_CopPlatehole;
    G4Box *sol_CopPlate;
    G4VSolid *sol_Copframe;
    G4Box *sol_CopBoxA;
    G4Box *sol_CopBoxB;
    G4Box *sol_CopDownDiskA;
    G4VSolid *sol_CopDownDiskA_Re;
    G4Box *sol_CopDownDiskB;
    G4VSolid *sol_CopDownDiskB_Re;

    ///////////////////////////////////////////////////////
    // Copper Rod
    G4double rxRod;
    G4double ryRod;
    G4double copperRods_H;

    G4Tubs *CopperframeRodX[7];
    G4VSolid *CopperframeRod[7];
    G4ThreeVector PosCuRod;

    ///////////////////////////////////////////////////////
    // Reflector around the CMO (Vm2000 shield)
    //

    // Scintillator
    G4VSolid *topScint_MotherSolid;
    G4VSolid *sideFBScint_MotherSolid;
    G4VSolid *sideLRScint_MotherSolid;

    G4Box *topScint_MotherBox;
    G4Box *topScint_Box;
    G4VSolid *topScint_MotherTrap;
    G4Trd *topScint_MotherPMTTrap;
    G4Trd *topScint_FlatTrap;
    G4Trd *topScint_PMTTrap;

    G4Box *sideFBScint_MotherBox;
    G4VSolid *sideFBScint_MotherTrap;
    G4Trd *sideFBScint_MotherPMTTrap;
    G4Box *sideFBScint_Box;
    G4Trd *sideFBScint_FlatTrap;
    G4Trd *sideFBScint_PMTTrap;

    G4Box *sideLRScint_MotherBox;
    G4VSolid *sideLRScint_MotherTrap;
    G4Trd *sideLRScint_MotherPMTTrap;
    G4Box *sideLRScint_Box;
    G4Trd *sideLRScint_FlatTrap;
    G4Trd *sideLRScint_PMTTrap;

    G4LogicalVolume *topScint_MotherLogical;
    G4LogicalVolume *sideFBScint_MotherLogical;
    G4LogicalVolume *sideLRScint_MotherLogical;

    ///////////////////////////////////////////////////////
    // B4C 24% rubber geometry
    ///////////////////////////////////////////////////////
    G4Box *toyRubberB4C_Box;
    G4Box *realPbTopRubberB4C_Box;
    G4Box *realPbUnderRubberB4C_Box;

    G4LogicalVolume *toyRubberB4C_BoxLogical;
    G4LogicalVolume *realPbTopRubberB4C_BoxLogical;
    G4LogicalVolume *realPbUnderRubberB4C_BoxLogical;

    G4VisAttributes *rubberB4C_BoxVisAttr;

    ///////////////////////////////////////////////////////
    // Polyethylene geometry
    ///////////////////////////////////////////////////////
    G4Box *toyPE_Box;
    G4Box *realPE_Top_Box;
    G4Box *realPE_Bottom_Box;
    G4Box *realPE_SideFB_Box;
    G4Box *realPE_SideLR1_Box;
    G4Box *realPE_SideLR2_Box;

    G4LogicalVolume *toyPE_BoxLogical;
    G4LogicalVolume *realPE_Top_BoxLogical;
    G4LogicalVolume *realPE_Bottom_BoxLogical;
    G4LogicalVolume *realPE_SideFB_BoxLogical;
    G4LogicalVolume *realPE_SideLR1_BoxLogical;
    G4LogicalVolume *realPE_SideLR2_BoxLogical;

    G4VisAttributes *pe_BoxVisAttr;

    ///////////////////////////////////////////////////////
    // Borated PE geometry
    ///////////////////////////////////////////////////////
    G4Box *boratedPE_Bottom_Box;
    G4Box *boratedPE_SideFB_Box;
    G4Box *boratedPE_SideLR_Box;

    G4LogicalVolume *boratedPE_Bottom_BoxLogical;
    G4LogicalVolume *boratedPE_SideFB_BoxLogical;
    G4LogicalVolume *boratedPE_SideLR_BoxLogical;

    G4VisAttributes *boratedPE_BoxVisAttr;

    ///////////////////////////////////////////////////////
    // Boric acid panel geometry
    ///////////////////////////////////////////////////////

    G4Box *sideLRBoricAcid_HousingBox; // X-Axis: FB | Y-Axis: LR
    G4Box *sideLRBoricAcid_Block;
    G4Box *sideFBBoricAcid_HousingBox;
    G4Box *sideFBBoricAcid_Block;
    G4Box *bottomBoricAcid_HousingBox;
    G4Box *bottomBoricAcid_Block;
    G4Box *topBoricAcid_HousingBox;
    G4Box *topBoricAcid_Block;

    G4LogicalVolume *sideLRBoricAcid_MotherLogical;
    G4LogicalVolume *sideLRBoricAcid_BlockLogical;
    G4LogicalVolume *sideFBBoricAcid_MotherLogical;
    G4LogicalVolume *sideFBBoricAcid_BlockLogical;
    G4LogicalVolume *bottomBoricAcid_MotherLogical;
    G4LogicalVolume *bottomBoricAcid_BlockLogical;
    G4LogicalVolume *topBoricAcid_MotherLogical;
    G4LogicalVolume *topBoricAcid_BlockLogical;
    G4VisAttributes *boricAcid_VisAttr;
    G4VisAttributes *boricAcidHousing_VisAttr;

    ////////////////////////////////////////////////////////////////
    // Physics Hall
    ////////////////////////////////////////////////////////////////

    boxHall  = new G4Box("hallbox", bounding_size, bounding_size, bounding_size);
    logiHall = new G4LogicalVolume(boxHall, _air, "logiHall");
    physHall = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), // translation
                                 logiHall,                        // associated logical vol
                                 "physHall",                      // name
                                 NULL,                            // parent
                                 false,                           // no "Many"
                                 0);                              // copy number

    // the class variable "world_phys" must be set to the world phys volume.
    world_phys = physHall;

    /////////////////////////////////////////////////////////////////
    // Rock·
    ////////////////////////////////////////////////////////////////

    RockSolid =
        new G4Sphere("Rock Hemisphere", 0, RockRadius, 0 * deg, 360 * deg, 0 * deg, 90 * deg);
    logiRock = new G4LogicalVolume(RockSolid, _rock, "logiRock");

    /////////////////////////////////////////////////////////////////
    // Cavity in the Rock
    ////////////////////////////////////////////////////////////////

    G4RotationMatrix *tube_alignMtx = new G4RotationMatrix();
    tube_alignMtx->rotateX(-90 * deg);
    RockCavitySolid =
        new G4Tubs("Inner tunnel", 0, RockCavityRadius, RockCavityRadius, 0, 180 * deg);
    logiRockCavity    = new G4LogicalVolume(RockCavitySolid, _air, "Rock Cavity Logical");
    logiRockCavityVis = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));
    logiRockCavityVis->SetForceSolid(false);

    RockDiskSolid = new G4Tubs("Rock Disk", 0, RockRadius, RockDiskThick_half, 0 * deg, 360 * deg);
    logiRockDisk  = new G4LogicalVolume(RockDiskSolid, _rock, "Rock Disk Logical");

    ///////////////////////////////////////////////////////
    // Work Area   3 m x 3 m x 3 m
    //
    // G4cout << "WorkArea Positon : " << workboxH << G4endl;
    WorkArea = new G4Box("WorkArea", workboxX, workboxY, workboxH);

    G4RotationMatrix *WA_alignMtx = new G4RotationMatrix();
    WA_alignMtx->rotateX(90 * deg);
    logiWorkArea = new G4LogicalVolume(WorkArea, _air, "logiWorkArea");

    ///////////////////////////////////////////////////////
    // Layer6_Pb Top Box 1500x1500x1674.5 mm
    //           Air box under the top plate volume
    TopPbBox   = new G4Box("TopPbBox", pbtopbox_size, pbtopbox_size / 2, pbtopbox_zsize);
    TopPbBlock = new G4Box("TopPbBlock", pbtopbox_size - pbtopbox_housing_thickness,
                           (pbtopbox_size / 2.) - pbtopbox_housing_thickness,
                           pbtopbox_zsize - pbtopbox_housing_thickness);

    logiTopPbBox = new G4LogicalVolume(TopPbBox, _stainless, "logiTopPbBox");

    logiTopPbBlock = new G4LogicalVolume(TopPbBlock, _lead, "logiTopPbBlock");

    new G4PVPlacement(nullptr, G4ThreeVector(), logiTopPbBlock, "physTopPbBlock", logiTopPbBox,
                      false, 0, false);

    ///////////////////////////////////////////////////////
    // Layer6_Pb body of Lead Box  1050
    PbBoxOut  = new G4Box("PbBoxOut", pbbox_size / 2., pbbox_size, pbbox_zsize);
    PbBoxIn   = new G4Box("PbBoxIn", (pbbox_size - pbbox_thick) / 2., pbbox_size - pbbox_thick,
                        pbbox_zsize - pbbox_thick / 2.);
    PbBox     = new G4SubtractionSolid("PbBox", PbBoxOut, PbBoxIn, 0,
                                   G4ThreeVector(pbbox_thick / 2., 0, pbbox_thick / 2.));
    logiPbBox = new G4LogicalVolume(PbBox, _stainless, "logiPbBox");

    PbBlockOut =
        new G4Box("PbBlockOut", pbbox_size / 2. - pbbox_housing_thickness,
                  pbbox_size - pbbox_housing_thickness, pbbox_zsize - pbbox_housing_thickness);

    PbBlockIn = new G4Box("PbBlockIn", (pbbox_size - pbbox_thick + pbbox_housing_thickness) / 2.,
                          pbbox_size - pbbox_thick + pbbox_housing_thickness,
                          pbbox_zsize - pbbox_thick / 2. + pbbox_housing_thickness / 2.);

    PbBlock = new G4SubtractionSolid("PbBlock", PbBlockOut, PbBlockIn, nullptr,
                                     G4ThreeVector((pbbox_thick - pbbox_housing_thickness) / 2., 0,
                                                   (pbbox_thick - pbbox_housing_thickness) / 2.));

    logiPbBlock = new G4LogicalVolume(PbBlock, _lead, "logiPbBlock");
    new G4PVPlacement(nullptr, G4ThreeVector(), logiPbBlock, "physPbBlock", logiPbBox, false, 0,
                      false);

    ///////////////////////////////////////////////////////
    // InnerDetector
    ///////////////////////////////////////////////////////
    InnerDetectorBox = new G4Box("InnerDetectorBox", pbbox_size, pbbox_size, pbbox_zsize);
    SSTopBoxForInnerDet =
        new G4Box("SSTopBoxForInnerDet", sstopbox_xsize, sstopbox_ysize, sst_zsize);
    InnerDetector = new G4UnionSolid("InnerDetectorSolid", InnerDetectorBox, SSTopBoxForInnerDet,
                                     nullptr, G4ThreeVector(0, 0, pbbox_zsize + sst_zsize));

    logiInnerDetector = new G4LogicalVolume(InnerDetector, _air, "logiInnerDetector");

    ///////////////////////////////////////////////////////
    // Layer5_OVC Stainless Steel (SS)
    //
    SSCylinder     = new G4Tubs("SSCylinder", 0., ss_radius, ss_height + ssb_zsize, 0, 360. * deg);
    SSTop          = new G4Box("SSTop", sstopbox_xsize, sstopbox_ysize, sst_zsize);
    SSOVC0         = new G4UnionSolid("SSOVC0", SSCylinder, SSTop, 0, SSTopShift);
    logiSSOVCOuter = new G4LogicalVolume(SSOVC0, _stainless, "logiSSOVCOuterLV");

    SSOVCInnerSolid =
        new G4Tubs("SSOVCInnerCylinder", 0, ss_radius - ss_thick, ss_height, 0, 360 * deg);
    logiSSOVCInner = new G4LogicalVolume(SSOVCInnerSolid, _vacuum, "SSOVCInnerLV");

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, ssb_zsize), logiSSOVCInner, "physSSOVCInner",
                      logiSSOVCOuter, false, 0, false);

    ///////////////////////////////////////////////
    // Mu-metal
    muPipe = new G4Tubs("muPipe", mu_radius - mu_thick, mu_radius, mu_height, 0, 360. * deg);
    logimu = new G4LogicalVolume(muPipe, _mumetal, "logimu");

    //////////////////////////////////////////////////////////
    // OVC welding parts

    G4double OVCwe_thickness = 3 * mm / 2.;
    G4ThreeVector OVCwePos   = G4ThreeVector(0, 0, pbbox_zsize - ss_height * 2.);
    G4Tubs *OVCwRing =
        new G4Tubs("SSRing", ss_radius, ss_radius + ss_thick, OVCwe_thickness, 0, 360 * deg);

    G4Box *OVCwBox   = new G4Box("SSweldingbox", OVCwe_thickness, ss_thick / 2., ss_height);
    G4VSolid *OVCwel = new G4UnionSolid("SS_Welding", OVCwRing, OVCwBox, 0,
                                        G4ThreeVector(0, ss_radius + ss_thick / 2., ss_height));

    G4LogicalVolume *lv_OVCwe = new G4LogicalVolume(OVCwel, _aluminum, "OVCwelding");

    new G4PVPlacement(nullptr, OVCwePos, lv_OVCwe, "physOVCwelding", logiInnerDetector, false, 0);

    ///////////////////////////////////////////////////////
    // Layer4_ 50K-SHIELD Copper

    Cu4OuterCylinder = new G4Tubs("Cu4OuterCylinder", 0., cu4_radius,
                                  cu4_inner_height + cu4b_zsize + cu4t_zsize, 0, 360. * deg);
    Cu4InnerCylinder =
        new G4Tubs("Cu4InnerCylinder", 0, cu4_radius - cu4_thick, cu4_inner_height, 0., 360. * deg);
    logiCu4Outer = new G4LogicalVolume(Cu4OuterCylinder, _copper, "logiCu4", 0, 0, 0);
    logiCu4Inner = new G4LogicalVolume(Cu4InnerCylinder, _vacuum, "logiCu4Inner");

    new G4PVPlacement(nullptr, PosCu4Inner, logiCu4Inner, "physCu4Inner", logiCu4Outer, false, 0,
                      false);

    ///////////////////////////////////////////////////////
    // Layer3_Cu IVC

    Cu3OuterCylinder = new G4Tubs("Cu3OuterCylinder", 0, cu3_radius,
                                  cu3t_zsize + cu3_inner_height + cu3b_zsize, 0, 360 * deg);
    Cu3InnerCylinder =
        new G4Tubs("Cu3InnerCylinder", 0, cu3_radius - cu3_thick, cu3_inner_height, 0, 360 * deg);
    logiCu3Outer = new G4LogicalVolume(Cu3OuterCylinder, _stainless, "logiCu3", 0, 0, 0);
    logiCu3Inner = new G4LogicalVolume(Cu3InnerCylinder, _vacuum, "logiCu3Inner");

    new G4PVPlacement(nullptr, PosCu3Outer, logiCu3Outer, "physIVCOuter", logiCu4Inner, false, 0);
    new G4PVPlacement(nullptr, PosCu3Inner, logiCu3Inner, "physIVCInner", logiCu3Outer, false, 0,
                      false);

    ////////////////////////////////////////////////////////
    // IVC welding part

    G4double IVCwe_thickness = 3 * mm / 2.;
    G4ThreeVector IVCwePos   = G4ThreeVector(
        0, 0, cu4_inner_height - cu3_gap_fromTop - (cu3_inner_height + cu3t_zsize) * 2.);
    G4Tubs *IVCwRing =
        new G4Tubs("Cu3Ring", cu3_radius, cu3_radius + cu3_thick, IVCwe_thickness, 0, 360 * deg);

    G4Box *IVCwBox = new G4Box("Cu3weldingbox", IVCwe_thickness, cu3_thick / 2., cu3_inner_height);
    G4VSolid *IVCwel =
        new G4UnionSolid("IVC_Welding", IVCwRing, IVCwBox, 0,
                         G4ThreeVector(0, cu3_radius + cu3_thick / 2., cu3_inner_height));

    G4LogicalVolume *lv_IVCwe = new G4LogicalVolume(IVCwel, _copper, "IVCwelding");

    new G4PVPlacement(nullptr, IVCwePos, lv_IVCwe, "physIVCwelding", logiCu4Inner, false, 0);

    ///////////////////////////////////////////////////////
    // Layer2_Cu3 SHIELD-STILL

    Cu2OuterCylinder = new G4Tubs("Cu2OuterCylinder", 0, cu2_radius,
                                  cu2_inner_height + cu2b_zsize + cu2t_zsize, 0, 360. * deg);
    Cu2InnerCylinder =
        new G4Tubs("Cu2InnerCylinder", 0, cu2_radius - cu2_thick, cu2_inner_height, 0, 360. * deg);
    logiCu2Outer = new G4LogicalVolume(Cu2OuterCylinder, _copper, "logiCu2Outer");
    logiCu2Inner = new G4LogicalVolume(Cu2InnerCylinder, _vacuum, "logiCu2Inner");

    new G4PVPlacement(nullptr, PosCu2Outer, logiCu2Outer, "physCu2Outer", logiCu3Inner, false, 0);
    new G4PVPlacement(nullptr, PosCu2Inner, logiCu2Inner, "physCu2Inner", logiCu2Outer, false, 0);

    ///////////////////////////////////////////////////////
    // Layer1_Cu1 50mk-SHIELD

    Cu1OuterCylinder = new G4Tubs("Cu1OuterCylinder", 0., cu1_radius,
                                  cu1_inner_height + cu1b_zsize + cu1t_zsize, 0., 360. * deg);
    Cu1InnerCylinder = new G4Tubs("Cu1InnerCylinder", 0., cu1_radius - cu1_thick, cu1_inner_height,
                                  0., 360. * deg);

    logiCu1Outer = new G4LogicalVolume(Cu1OuterCylinder, _copper, "logiCu1Outer");
    logiCu1Inner = new G4LogicalVolume(Cu1InnerCylinder, _vacuum, "logiCu1Inner");

    new G4PVPlacement(nullptr, PosCu1Inner, logiCu1Inner, "physCu1Inner", logiCu1Outer, false, 0);
    new G4PVPlacement(nullptr, PosCu1Outer, logiCu1Outer, "physCu1Outer", logiCu2Inner, false, 0);

    //////////////////////////////////////////////////////
    // Super Conducting magnetic shielding

    SuperMS0 = new G4Tubs("SuperMS0", 0, SMS_radius, SMS_height, 0, 360. * deg);
    SuperMS1 = new G4Tubs("SuperMS1", 0, SMS_radius - SMS_thickness, SMS_height, 0, 360. * deg);
    SuperMS  = new G4SubtractionSolid("Top-co-hole[4]", SuperMS0, SuperMS1, 0,
                                     G4ThreeVector(0, 0, SMS_thickness));

    logiSuperMS = new G4LogicalVolume(SuperMS, _lead, "logiSuperMS", 0, 0, 0);

    ///////////////////////////////////////////////////////
    // CMO Target area

    TargetRoom     = new G4Tubs("TargetRoom", 0, TR_radius, TR_height, 0, 360. * deg);
    logiTargetRoom = new G4LogicalVolume(TargetRoom, _vacuum, "logiTargetRoom", 0, 0, 0);

    ///////////////////////////////////////////////////////
    // Mixing Chamber Cu Plate

    CuMCPlate = new G4Tubs("CuMCPlate", 0, cumcp_radius, cumcp_height, 0, 360. * deg);
    logiCuMCP = new G4LogicalVolume(CuMCPlate, _copper, "logiCuMCP", 0, 0, 0);

    ///////////////////////////////////////////////////////
    // Plate1_Cu
    CuPlate1 = new G4Tubs("CuPlate1", 0, cumcp_radius, cup1_zsize, 0, 360. * deg);
    logiCuP1 = new G4LogicalVolume(CuPlate1, _copper, "logiCuP1", 0, 0, 0);

    ///////////////////////////////////////////////////////
    // Plate2_Pb Shield
    PbPlate2 = new G4Tubs("PbPlate2", 0, cumcp_radius, pbp2_zsize, 0, 360. * deg);
    logiPbP2 = new G4LogicalVolume(PbPlate2, _lead, "logiPbP2", 0, 0, 0);

    ///////////////////////////////////////////////////////
    CuPlate3 = new G4Tubs("CuPlate3", 0, cumcp_radius, cup3_zsize, 0, 360. * deg);
    logiCuP3 = new G4LogicalVolume(CuPlate3, _copper, "logiCuP3", 0, 0, 0);

    /////////////////////////////////////////////////////////
    // G10 material holders between SS and 50K-SHIELD Copper

    for (int g10Lyr = 0; g10Lyr < 5; g10Lyr++) {
        G4Tubs *PlateGap4Vol = new G4Tubs("PlateGap4Vol", 0, PlateGapR[g10Lyr] - 3. * mm,
                                          CuholderH / 2., 0. * deg, 360. * deg);
        PlateGap4VolLogi[g10Lyr] =
            new G4LogicalVolume(PlateGap4Vol, _vacuum, "PlateGap4VolLogical", 0, 0, 0);

        G4Tubs *G10holderSolid =
            new G4Tubs("G10holderSolid", // name
                       Cuholder_r, Cuholder_R, CuholderH / 2, 0. * deg, 360. * deg);
        G4LogicalVolume *G10holderLogi =
            new G4LogicalVolume(G10holderSolid, g10material, "G10holderLogical", 0, 0, 0);

        if (g10Lyr < 4) {
            for (int NinLyr = 0; NinLyr < 3; NinLyr++) {
                copyNo               = NinLyr;
                xG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * cos((NinLyr * 120.) * deg));
                yG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * sin((NinLyr * 120.) * deg));
                PosG10holder         = G4ThreeVector(xG10[NinLyr][g10Lyr], yG10[NinLyr][g10Lyr], 0);
                new G4PVPlacement(nullptr, PosG10holder, G10holderLogi, "G10holderPhys",
                                  PlateGap4VolLogi[g10Lyr], false, copyNo);
            }
        } else {
            for (int NinLyr = 0; NinLyr < 6; NinLyr++) {
                copyNo               = NinLyr;
                xG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * cos((NinLyr * 60.) * deg));
                yG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * sin((NinLyr * 60.) * deg));
                PosG10holder         = G4ThreeVector(xG10[NinLyr][g10Lyr], yG10[NinLyr][g10Lyr], 0);
                new G4PVPlacement(nullptr, PosG10holder, G10holderLogi, "G10holderPhys",
                                  PlateGap4VolLogi[g10Lyr], false, copyNo);
            }
        }
        G4VisAttributes *G10holderVisAtt = new G4VisAttributes(lgreen);
        G10holderVisAtt->SetVisibility(true);
        G10holderVisAtt->SetForceSolid(true);
        G10holderLogi->SetVisAttributes(G10holderVisAtt);
        G4cout << "###   " << G10holderLogi->GetName() << "  " << G10holderLogi->GetMass() / kg
               << G4endl;
        G4cout << "###   " << PlateGap4VolLogi[g10Lyr]->GetName() << "  "
               << PlateGap4VolLogi[g10Lyr]->GetMass() / kg << G4endl;
    }

    ///////////////////////////////////////////////////////
    // CMO cell
    //
    constexpr const int CMOCellNum = 6;
    for (int nn = 1; nn < CMOCellNum + 1; nn++) { //  5 layers
        zCell[nn] = zCell[nn - 1] + cmocell_h[nn - 1] + CMOSupport_H * 2 + cmocell_h[nn] +
                    CMOSupport_H * 2 + CMOSupportGap * 2 - 0.65 * mm;
    }
    G4double Zcen = zCell[2] + 26;
    for (int nn = 0; nn < CMOCellNum + 1; nn++) { //  5 layers
        zCell[nn]  = Zcen - zCell[nn] - .55 * cm; //-5.9*cm;
        zCup[nn]   = zCell[nn] + cmocell_h[nn];
        zCdown[nn] = zCell[nn] - cmocell_h[nn];
    }
    fPilot_logiCMOCell = new G4LogicalVolume *[CMOCellNum];
    for (int nn = 0; nn < CMOCellNum; nn++) { //  5 layers
        sprintf(namecmo, "CMOCell%d", nn);
        CMOCell[nn] = new G4EllipticalTube(namecmo, // name
                                           cmocell_lr[nn], cmocell_sr[nn], cmocell_h[nn]);
        sprintf(namecmo, "logiCMOCell%d", nn);
        logiCMOCell[nn]        = new G4LogicalVolume(CMOCell[nn], CaMoO4, namecmo, 0, 0, 0);
        fPilot_logiCMOCell[nn] = logiCMOCell[nn];

        copyNo = nn;

        PosCMOCell = G4ThreeVector(0, 0, zCell[nn]);

        sprintf(namecmo, "physCMOCell%d", nn);
        new G4PVPlacement(nullptr, PosCMOCell, logiCMOCell[nn], namecmo, logiTargetRoom, false,
                          copyNo);
        G4cout << "###   " << logiCMOCell[nn]->GetName() << "  " << logiCMOCell[nn]->GetMass() / kg
               << "  " << zCell[nn] << "  " << PosTR << "  " << PosInnerDet << G4endl;
    }

    ///////////////////////////////////////////////////////
    // New Photon Detector : Run-6
    //
    //
    G4EllipticalTube *circle1;
    G4EllipticalTube *circle2;

    G4RotationMatrix *aRotation = new G4RotationMatrix();

    PhotonFramePlate = new G4Box("PhotonFrame", PhotonFrame_x, PhotonFrame_y, PhotonFrame_z);
    hole = new G4Tubs("Photonhole", 0, 49 / 2, PhotonFrame_z + 10, -20 * deg, 225 * deg);
    PhotonFrameX =
        new G4SubtractionSolid("Top-co-hole[4]", PhotonFramePlate, hole, 0, G4ThreeVector(0, 0, 0));
    circle1 = new G4EllipticalTube("Photoncircle1", 8 * mm, 15.5 * mm, PhotonFrame_z + 10);
    aRotation->rotateZ(135 * deg);
    PhotonFrameW = new G4SubtractionSolid("Top-co-hole[4]", PhotonFrameX, circle1, aRotation,
                                          G4ThreeVector(10, 45 / 4 * cos(120 * deg) + 6 * mm, 0));

    circle2 = new G4EllipticalTube("Photoncircle1", 8 * mm, 15.5 * mm, PhotonFrame_z + 10);
    aRotation->rotateZ(270 * deg);
    PhotonFrameY = new G4SubtractionSolid("Top-co-hole[4]", PhotonFrameW, circle2, aRotation,
                                          G4ThreeVector(-10, 45 / 4 * cos(120 * deg) + 6 * mm, 0));

    // hole for peek
    G4VSolid *PEEKa;
    G4Tubs *aPEEK = new G4Tubs("aPEEK", 0, 2.9 * mm, 2 * mm / 2, 0, 360 * deg);
    PEEKa = new G4SubtractionSolid("aPEEK", PhotonFrameY, aPEEK, 0, G4ThreeVector(15, 21, 1));

    G4VSolid *PEEKb;
    G4Tubs *bPEEK = new G4Tubs("bPEEK", 0, 2.9 * mm, 2 * mm / 2, 0, 360 * deg);
    PEEKb         = new G4SubtractionSolid("bPEEK", PEEKa, bPEEK, 0, G4ThreeVector(-15, 21, 1));

    G4VSolid *PEEKc;
    G4Tubs *cPEEK = new G4Tubs("cPEEK", 0, 2.9 * mm, 2 * mm / 2, 0, 360 * deg);
    PEEKc         = new G4SubtractionSolid("cPEEK", PEEKb, cPEEK, 0, G4ThreeVector(15, -21, 1));

    G4VSolid *PEEKd;
    G4Tubs *dPEEK = new G4Tubs("dPEEK", 0, 2.9 * mm, 2 * mm / 2, 0, 360 * deg);
    PEEKd         = new G4SubtractionSolid("dPEEK", PEEKc, dPEEK, 0, G4ThreeVector(-15, -21, 1));

    hole        = new G4Tubs("Photonhole", 0, GeWafer_r + 1 * mm, 2 * mm, 0 * deg, 360 * deg);
    PhotonFrame = new G4SubtractionSolid("Top-co-hole[4]", PEEKd, hole, 0, G4ThreeVector(0, 0, 1.));

    G4Box *sol_CopPlatehole3 =
        new G4Box("CopperFrame", PhotonFrame_x - 7.5 * mm, PhotonFrame_y + 10 * mm, 5. / 2 * mm);
    G4VSolid *PhotonFrame3 =
        new G4SubtractionSolid("Top-co-hole[4]", PhotonFrame, sol_CopPlatehole3, 0,
                               G4ThreeVector(0, 0, -PhotonFrame_z + 1));

    sol_CopPlatehole2 =
        new G4Box("CopperFrame", PhotonFrame_x + 5 * mm, PhotonFrame_y - 6.5 * mm, 5. / 2 * mm);
    PhotonFrame2 = new G4SubtractionSolid("CopperFrameDetail", PhotonFrame3, sol_CopPlatehole2, 0,
                                          G4ThreeVector(0, 0, -2.1));

    lv_PhotonFrame = new G4LogicalVolume(PhotonFrame2, _copper2, "logiPhotonFrame");

    G4cout << "###   " << lv_PhotonFrame->GetName() << "  " << lv_PhotonFrame->GetMass() / kg
           << G4endl;

    PhotonDet     = new G4Tubs("PhotonDet", 0, PhotonDet_r, 0.2 * mm, 0, 360 * deg);
    logiPhotonDet = new G4LogicalVolume(PhotonDet, _vacuum, "logiPhotonDet");

    GeWafer     = new G4Tubs("GeWafer", 0, GeWafer_r, GeWafer_zsize, 0, 360 * deg);
    logiGeWafer = new G4LogicalVolume(GeWafer, _gewafer, "logiGeWafer", 0, 0, 0);

    new G4PVPlacement(nullptr, PosGeWafer, logiGeWafer, "physGeWafer", logiPhotonDet, false, 0);

    VacDisk     = new G4Tubs("VacDisk", 0, VacDisk_r, VacDisk_zsize, 0, 360 * deg);
    logiVacDisk = new G4LogicalVolume(VacDisk, _vacuum, "logiVacDisk", 0, 0, 0);

    PosVacDisk = G4ThreeVector(0, 0, VacDiskZ);
    new G4PVPlacement(nullptr, PosVacDisk, logiVacDisk, "physVacDisk", logiGeWafer, false, 0);

    // Springs of light detector : PEEK4

    G4Tubs *PEEK4           = new G4Tubs("PEEK4", 0, 2.775 * mm, 1.25 * mm / 2, 0, 360 * deg);
    G4LogicalVolume *lPEEK4 = new G4LogicalVolume(PEEK4, _peek1, "logiPEEK4", 0, 0, 0);
    char namePEEK4[20];

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        sprintf(namePEEK4, "PEEK4_%d", nn);
        new G4PVPlacement(0, G4ThreeVector(15, 21, PhotonDetZ + 1.05 * mm), lPEEK4, namePEEK4,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(-15, 21, PhotonDetZ + 1.05 * mm), lPEEK4, namePEEK4,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(15, -21, PhotonDetZ + 1.05 * mm), lPEEK4, namePEEK4,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(-15, -21, PhotonDetZ + 1.05 * mm), lPEEK4, namePEEK4,
                          logiTargetRoom, false, copyNo);
    }

    G4VisAttributes *logiPEEK4_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiPEEK4_VisAtt->SetVisibility(true);
    logiPEEK4_VisAtt->SetForceSolid(true);
    lPEEK4->SetVisAttributes(logiPEEK4_VisAtt);

    G4cout << "###   " << lPEEK4->GetName() << "  " << lPEEK4->GetMass() / kg << G4endl;

    // Bolts at the Clamp Teflon
    G4Tubs *BCT           = new G4Tubs("BCT", 0, 3. * mm / 2, 2. * mm / 2, 0, 360 * deg);
    G4LogicalVolume *lBCT = new G4LogicalVolume(BCT, _stainless, "logiBCT", 0, 0, 0);
    char nameBCT[20];

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        sprintf(nameBCT, "BCT_%d", nn);
        new G4PVPlacement(0, G4ThreeVector(19.5, 26.5, PhotonDetZ - 1. * mm), lBCT, nameBCT,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(19.5, -26.5, PhotonDetZ - 1. * mm), lBCT, nameBCT,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(-19.5, 26.5, PhotonDetZ - 1. * mm), lBCT, nameBCT,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(-19.5, -26.5, PhotonDetZ - 1. * mm), lBCT, nameBCT,
                          logiTargetRoom, false, copyNo);
    }

    G4VisAttributes *logiBCT_VisAtt = new G4VisAttributes(G4Colour::Brown());
    logiBCT_VisAtt->SetVisibility(true);
    logiBCT_VisAtt->SetForceSolid(true);
    lBCT->SetVisAttributes(logiBCT_VisAtt);

    G4cout << "###   " << lBCT->GetName() << "  " << lBCT->GetMass() / kg << G4endl;

    ////////////////////////////////////////////////////////////////////////
    // PEEK below CMO

    G4Tubs *PEEKCr           = new G4Tubs("PEEKCr", 0, 4.787 * mm, 0.42 * mm / 2, 0, 360 * deg);
    G4LogicalVolume *lPEEKCr = new G4LogicalVolume(PEEKCr, _peek1, "logiPEEKCr", 0, 0, 0);
    char namePEEKCr[20];

    for (int nn = 1; nn < 7; nn++) {
        copyNo = nn;
        sprintf(namePEEKCr, "PEEKCr_%d", nn);

        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);

        new G4PVPlacement(0, G4ThreeVector(2, -15, PhotonDetZ + 9.127 * mm), lPEEKCr, namePEEKCr,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(-15, 0, PhotonDetZ + 9.127 * mm), lPEEKCr, namePEEKCr,
                          logiTargetRoom, false, copyNo);
        new G4PVPlacement(0, G4ThreeVector(15, 2, PhotonDetZ + 9.127 * mm), lPEEKCr, namePEEKCr,
                          logiTargetRoom, false, copyNo);
    }

    G4VisAttributes *logiPEEKCr_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiPEEKCr_VisAtt->SetVisibility(true);
    logiPEEKCr_VisAtt->SetForceSolid(true);
    lPEEKCr->SetVisAttributes(logiPEEKCr_VisAtt);

    G4cout << "###   " << lPEEKCr->GetName() << "  " << lPEEKCr->GetMass() / kg << G4endl;

    /////////////////////////////////////////////////////////////////////////
    // PCB + Lead Tin Solder
    //

    G4Box *PCB       = new G4Box("PCB", 16.8 * mm / 2, 29.26 * mm / 2, 0.115 * mm / 2);
    G4Box *MothEpoxy = new G4Box("MothEpoxy", 15 * mm, 45 * mm / 2, 1 * mm / 2);
    G4Box *Epoxy     = new G4Box("Epoxy", 0.5624 * mm, 0.5624 * mm, 0.46 * mm / 2);
    G4Box *Solder    = new G4Box("Solder", 0.5 * mm, 2.53 * mm, 0.115 * mm / 2);
    G4Box *Solder2   = new G4Box("Solder2", 2 * mm, 0.633 * mm, 0.115 * mm / 2);
    G4Box *PCBa      = new G4Box("PCBa", 11.8 * mm / 2, 11.8 * mm / 2, 0.115 * mm / 2);
    G4Box *PCBb      = new G4Box("PCBb", 8.5 * mm / 2, 5.5 * mm / 2, 0.115 * mm / 2 + 2 * mm);

    char nameEpoxy[30];
    char namePCB[20];
    char nameSolder[30];
    char nameSolder2[30];
    char nameMothEpoxy[20];

    G4VSolid *UPCB;
    G4VSolid *SubPCB;
    G4VPhysicalVolume *physMothEpoxy;

    G4RotationMatrix *cRotation = new G4RotationMatrix();
    cRotation->rotateZ(135 * deg);

    UPCB   = new G4UnionSolid("PCB+PCBa", PCB, PCBa, cRotation, G4ThreeVector(0, 14.7 * mm, 0));
    SubPCB = new G4SubtractionSolid("SubPCB", UPCB, PCBb, cRotation, G4ThreeVector(-2, 21., 0));

    G4LogicalVolume *lPCB       = new G4LogicalVolume(SubPCB, _copper, "logiPCB", 0, 0, 0);
    G4LogicalVolume *lMothEpoxy = new G4LogicalVolume(MothEpoxy, _vacuum, "logiMothEpoxy", 0, 0, 0);
    G4LogicalVolume *lEpoxy     = new G4LogicalVolume(Epoxy, _solder, "logiEpoxy", 0, 0, 0);
    G4LogicalVolume *lSolder    = new G4LogicalVolume(Solder, _solder, "logiSolder", 0, 0, 0);
    G4LogicalVolume *lSolder2   = new G4LogicalVolume(Solder2, _solder, "logiSolder2", 0, 0, 0);

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);

        sprintf(nameMothEpoxy, "MothEpoxy_%d", nn);
        physMothEpoxy = new G4PVPlacement(0, G4ThreeVector(0, -20, PhotonDetZ + 3.775), lMothEpoxy,
                                          nameMothEpoxy, logiTargetRoom, false, copyNo);
        (void)physMothEpoxy;
    }

    sprintf(namePCB, "PCB_");
    new G4PVPlacement(0, G4ThreeVector(0, -5, 0.0225), lPCB, namePCB, lMothEpoxy, false, copyNo);

    sprintf(nameEpoxy, "Epoxy_");
    new G4PVPlacement(cRotation, G4ThreeVector(-5, 10, -0.2655), lEpoxy, nameEpoxy, lMothEpoxy,
                      false, copyNo);

    new G4PVPlacement(cRotation, G4ThreeVector(2, 13, -0.2665), lEpoxy, nameEpoxy, lMothEpoxy,
                      false, copyNo);

    sprintf(nameSolder, "Solder_");
    new G4PVPlacement(0, G4ThreeVector(3, 0, 0.14), lSolder, nameSolder, lMothEpoxy, false, copyNo);

    sprintf(nameSolder2, "Solder2_");
    new G4PVPlacement(0, G4ThreeVector(-3, -3, 0.14), lSolder2, nameSolder2, lMothEpoxy, false,
                      copyNo);

    G4VisAttributes *logiPCB_VisAtt = new G4VisAttributes(G4Colour::White());
    logiPCB_VisAtt->SetVisibility(true);
    logiPCB_VisAtt->SetForceSolid(true);
    lPCB->SetVisAttributes(logiPCB_VisAtt);

    G4VisAttributes *logiMothEpoxy_VisAtt = new G4VisAttributes(G4Colour::Green());
    logiMothEpoxy_VisAtt->SetVisibility(true);
    logiMothEpoxy_VisAtt->SetForceWireframe(true);
    lMothEpoxy->SetVisAttributes(logiMothEpoxy_VisAtt);

    G4VisAttributes *logiEpoxy_VisAtt = new G4VisAttributes(G4Colour::Gray());
    logiEpoxy_VisAtt->SetVisibility(true);
    logiEpoxy_VisAtt->SetForceSolid(true);
    lEpoxy->SetVisAttributes(logiEpoxy_VisAtt);

    G4VisAttributes *logiSolder_VisAtt = new G4VisAttributes(G4Colour::Gray());
    logiSolder_VisAtt->SetVisibility(true);
    logiSolder_VisAtt->SetForceSolid(true);
    lSolder->SetVisAttributes(logiSolder_VisAtt);

    G4VisAttributes *logiSolder2_VisAtt = new G4VisAttributes(G4Colour::Gray());
    logiSolder2_VisAtt->SetVisibility(true);
    logiSolder2_VisAtt->SetForceSolid(true);
    lSolder2->SetVisAttributes(logiSolder2_VisAtt);

    G4cout << "###   " << lSolder->GetName() << "  " << lSolder->GetMass() / kg << G4endl;

    G4cout << "###   " << lSolder2->GetName() << "  " << lSolder2->GetMass() / kg << G4endl;

    G4cout << "###   " << lEpoxy->GetName() << "  " << lEpoxy->GetMass() / kg << G4endl;

    G4cout << "###   " << lPCB->GetName() << "  " << lPCB->GetMass() / kg << G4endl;
    ///////////////////////////////////////////////////////
    // Ge wafer on the bottom of the Crystal
    //
    ///////////////////////////////////////////////////////
    // Duplicate Photon Detectors
    //

    for (int nn = 0; nn < 6; nn++) { //  5 layers
        copyNo = nn;

        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, PhotonDetZ + 0.22 * mm), logiPhotonDet,
                          "physPhotonDet", logiTargetRoom, false, copyNo);
        new G4PVPlacement(nullptr, PosPhotonDet, lv_PhotonFrame, "physPhotonFrame", logiTargetRoom,
                          false, copyNo);
    }

    ///////////////////////////////////////////////////////
    // TopAssembly frame

    Top_co    = new G4Box("topAss", 78. / 2, 73. / 2, 4. / 2);
    smallbox  = new G4Box("smallBox", 10 / 2, 14 / 2, 6 / 2);
    smallboxB = new G4Box("smallBoxB", 14 / 2, 14 / 2, 6 / 2);
    hole      = new G4Tubs("hole", 0, 5, 5, 0, 360 * deg);

    pillar = new G4Tubs("pillar", 0, 5, pillarH / 2., 0, 360 * deg);
    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            TopCopper[kk] = new G4SubtractionSolid(
                "Top-co-hole[4]", Top_co, hole, 0,
                G4ThreeVector(78 / 2 * cos(90 * kk * deg), 73 / 2 * sin(90 * kk * deg), 0));
        } else {
            TopCopper[kk] = new G4SubtractionSolid(
                "Top-co-hole[4]", TopCopper[kk - 1], hole, 0,
                G4ThreeVector(78 / 2 * cos(90 * kk * deg), 73 / 2 * sin(90 * kk * deg), 0));
        }
    }

    lv_smallbox  = new G4LogicalVolume(smallbox, _copper, "smallbox");
    lv_smallboxB = new G4LogicalVolume(smallboxB, _copper, "smallbox");
    lv_copillar  = new G4LogicalVolume(pillar, _copper, "copillar");

    G4cout << "###   " << lv_smallbox->GetName() << "  " << lv_smallbox->GetMass() / kg * 4
           << G4endl;
    G4cout << "###   " << lv_smallboxB->GetName() << "  " << lv_smallboxB->GetMass() / kg * 4
           << G4endl;
    G4cout << "###   " << lv_copillar->GetName() << "  " << lv_copillar->GetMass() / kg * 4
           << G4endl;

    // assembly locate information///////////

    for (int i = 0; i < 4; i++) {
        new G4PVPlacement(nullptr,
                          G4ThreeVector(Assem_px + x_smallbox[i] * (78 / 2 - 10 / 2),
                                        Assem_py + y_smallbox[i] * (73 / 2 - 8 - 14 / 2),
                                        top_pz + 5 + TR_height - 55.),
                          lv_smallbox, "smallbox", logiTargetRoom, false, 0);
        new G4PVPlacement(nullptr,
                          G4ThreeVector(Assem_px + x_smallbox[i] * (78 / 2 - 14 / 2),
                                        Assem_py + y_smallbox[i] * (73 / 2 - 8 - 14 / 2),
                                        -(top_pz + 5) + TR_height - 45. + pillarH),
                          lv_smallboxB, "smallboxB", logiTargetRoom, false, 0);
        new G4PVPlacement(nullptr,
                          G4ThreeVector(Assem_px + x_smallbox[i] * (78 / 2 - 14 - 5 / 2),
                                        Assem_py + y_smallbox[i] * (73 / 2 - 5),
                                        0. + TR_height - 50. + pillarH / 2.),
                          lv_copillar, "copillar", logiTargetRoom, false, 0);
    }

    G4LogicalVolume *lv_TopPlate  = new G4LogicalVolume(TopCopper[3], _copper, "topAss");
    G4LogicalVolume *lv_TopPlateB = new G4LogicalVolume(TopCopper[3], _copper, "topAss");
    new G4PVPlacement(nullptr, G4ThreeVector(Assem_px, Assem_py, top_pz + TR_height - 55),
                      lv_TopPlate, "TopAss", logiTargetRoom, false, 0);
    new G4PVPlacement(0, G4ThreeVector(Assem_px, Assem_py, -top_pz + TR_height - 45 + pillarH),
                      lv_TopPlate, "TopAss", logiTargetRoom, false, 0);

    G4cout << "###   " << lv_TopPlate->GetName() << "  " << lv_TopPlate->GetMass() / kg * 2
           << G4endl;

    //////////////////////////////////////////////////////
    // MASS-Spring
    //
    G4double MASSposZ     = CenterCuP6Z - cup3_zsize - 25. * cm;
    G4ThreeVector MASSPos = G4ThreeVector(0, 0., MASSposZ);

    G4Tubs *mass_tube =
        new G4Tubs("mass_tube", 186. / 2 * mm, 398 / 2 * mm, 18 / 2. * mm, 0. * deg, 360. * deg);
    G4Tubs *mass_tubeT =
        new G4Tubs("mass_tubeT", 0., 226. / 2 * mm, 11 / 2. * mm, 0. * deg, 360. * deg);

    G4Box *t_box          = new G4Box("t_box", 10 / 2 * mm, 10 / 2 * mm, 132 / 2. * mm);
    G4Box *t_box2         = new G4Box("t_box2", 37.13 / 2 * mm, 10 / 2 * mm, 5 / 2. * mm);
    G4UnionSolid *mass_ta = new G4UnionSolid("mass_ta", t_box, t_box2, 0,
                                             G4ThreeVector(0., 0., 132 / 2. * mm + 2.5 * mm));
    G4UnionSolid *mass_tb = new G4UnionSolid("mass_tb", mass_ta, t_box2, 0,
                                             G4ThreeVector(0., 0., -(132 / 2. * mm + 2.5 * mm)));

    G4UnionSolid *mass_t          = nullptr;
    G4RotationMatrix *myRotation2 = new G4RotationMatrix();
    for (int i = 0; i < 4; i++) {
        myRotation2->rotateZ(90. * deg);

        if (i == 0) {
            mass_t = new G4UnionSolid("mass_t", mass_tube, mass_tb, myRotation2,
                                      G4ThreeVector(196 / 2. * mm * cos(90 * i * deg),
                                                    196 / 2. * mm * sin(90 * i * deg),
                                                    (142. / 2 * mm + 17.9 / 2. * mm)));
        }

        else {
            mass_t = new G4UnionSolid("mass_t", mass_t, mass_tb, myRotation2,
                                      G4ThreeVector(196 / 2. * mm * cos(90 * i * deg),
                                                    196 / 2. * mm * sin(90 * i * deg),
                                                    (142. / 2 * mm + 17.9 / 2. * mm)));
        }
    }
    G4UnionSolid *mass_sp = new G4UnionSolid("mass_sp", mass_t, mass_tubeT, 0,
                                             G4ThreeVector(0., 0., (142. * mm + 18.95 / 2. * mm)));

    G4LogicalVolume *logicMASS = new G4LogicalVolume(mass_sp,        // its solid
                                                     _copper2,       // its material
                                                     "MASS-Spring"); // its name

    G4cout << "###   " << logicMASS->GetName() << "  " << logicMASS->GetMass() / kg << G4endl;
    G4Tubs *spring =
        new G4Tubs("spring", 3. / 2. * cm, 3.3 / 2. * cm, 24. / 2. * cm, 0. * deg, 365. * deg);
    G4LogicalVolume *lv_spring = new G4LogicalVolume(spring, _vacuum, "lv_spring");

    G4cout << "###   " << lv_spring->GetName() << "  " << lv_spring->GetMass() / kg << G4endl;

    ///////////////////////////////////////////////////////
    // Copper frame
    //

    sol_Cophole = new G4Tubs("hole", 0, 5.5 / 2, 10, 0, 360 * deg);
    sol_CopPlatehole =
        new G4Box("CopperFrame", CopperFx / 2 - 1 * cm, 7.39 * cm / 2 - 1 * cm, 5.65 / 2 * mm);
    sol_CopPlate = new G4Box("CopperFrame", CopperFx / 2, CopperFy / 2, CopperFz / 2);
    sol_Copframe =
        new G4SubtractionSolid("CopperFrame_plate-vaccum", sol_CopPlate, sol_CopPlatehole);
    sol_CopBoxA      = new G4Box("CopperBoxA", 1. / 2 * cm, 1. / 2. * cm, 6. / 2 * mm);
    sol_CopBoxB      = new G4Box("CopperBoxB", 1. / 2 * cm, 1. / 2. * cm, 4. / 2 * mm);
    sol_CopDownDiskA = new G4Box("CopperPlateA", 1.2 * cm / 2, CopperFy / 2., CopperFz / 4.);
    sol_CopDownDiskA_Re =
        new G4SubtractionSolid("CopperPlateA_Re", sol_CopDownDiskA, sol_CopPlatehole, 0,
                               G4ThreeVector(0, 0, 5.65 / 2 * mm + 4.64 / 4 * mm - 1 * mm));
    sol_CopDownDiskB = new G4Box("CopperPlateB", 1.95 * cm / 2, CopperFy / 2., CopperFz / 4.);
    sol_CopDownDiskB_Re =
        new G4SubtractionSolid("CopperPlateB_Re", sol_CopDownDiskB, sol_CopPlatehole, 0,
                               G4ThreeVector(0, 0, 5.65 / 2 * mm + 4.64 / 4 * mm - 1 * mm));

    new G4LogicalVolume(sol_Copframe, _copper, "logiSquareDisk", 0, 0, 0);

    ////////Make Upper Disk//////////////////////////////////////////////
    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            sol_CopperFrame[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxA", sol_Copframe, sol_CopBoxA, 0,
                G4ThreeVector((CopperFx / 2 - 5 * mm) * x_smallbox[kk],
                              (CopperFy / 2 - 15 * mm) * y_smallbox[kk], CopperFz / 2 + 3 * mm));
        } else {
            sol_CopperFrame[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxA", sol_CopperFrame[kk - 1], sol_CopBoxA, 0,
                G4ThreeVector((CopperFx / 2 - 5 * mm) * x_smallbox[kk],
                              (CopperFy / 2 - 15 * mm) * y_smallbox[kk], CopperFz / 2 + 3 * mm));
        }
    }
    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            sol_CopperFrame2[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrame[3], sol_CopBoxB, 0,
                G4ThreeVector((CopperFx / 2 - 5 * mm) * x_smallbox[kk],
                              (CopperFy / 2 - 5 * mm) * y_smallbox[kk], -CopperFz / 2 - 2 * mm));
        } else {
            sol_CopperFrame2[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrame2[kk - 1], sol_CopBoxB, 0,
                G4ThreeVector((CopperFx / 2 - 5 * mm) * x_smallbox[kk],
                              (CopperFy / 2 - 5 * mm) * y_smallbox[kk], -CopperFz / 2 - 2 * mm));
        }
    }
    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            sol_CopperFrame3[kk] =
                new G4SubtractionSolid("CopperFrame - CopBoxB", sol_CopperFrame2[3], sol_Cophole, 0,
                                       G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                                                     (CopperFy / 2 - 5) * y_smallbox[kk], -3));
        } else {
            sol_CopperFrame3[kk] = new G4SubtractionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrame3[kk - 1], sol_Cophole, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 5) * y_smallbox[kk], -3));
        }
    }

    ////////Make down Disk //////////////////////////////////////////////

    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            sol_CopperFrameb[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxA", sol_Copframe, sol_CopBoxA, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 15) * y_smallbox[kk], -CopperFz / 2 - 3));
        } else {
            sol_CopperFrameb[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxA", sol_CopperFrameb[kk - 1], sol_CopBoxA, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 15) * y_smallbox[kk], -CopperFz / 2 - 3));
        }
    }
    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            sol_CopperFrameb2[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrameb[3], sol_CopBoxB, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 5) * y_smallbox[kk], +CopperFz / 2 + 2));
        } else {
            sol_CopperFrameb2[kk] = new G4UnionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrameb2[kk - 1], sol_CopBoxB, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 5) * y_smallbox[kk], +CopperFz / 2 + 2));
        }
    }
    for (kk = 0; kk < 4; kk++) {
        if (kk == 0) {
            sol_CopperFrameb3[kk] = new G4SubtractionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrameb2[3], sol_Cophole, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 5) * y_smallbox[kk], +3));
        } else {
            sol_CopperFrameb3[kk] = new G4SubtractionSolid(
                "CopperFrame - CopBoxB", sol_CopperFrameb3[kk - 1], sol_Cophole, 0,
                G4ThreeVector((CopperFx / 2 - 5) * x_smallbox[kk],
                              (CopperFy / 2 - 5) * y_smallbox[kk], +3));
        }
    }

    // bottom part plate make (bottom plate for phonon senser)

    sol_CopperFrame4[0] = new G4UnionSolid(
        "CopperFrame + CoppperPlateA", sol_CopperFrameb3[3], sol_CopDownDiskA_Re, 0,
        G4ThreeVector(1.5 * cm / 2 + 4 * mm, 0, -CopperFz / 4. - CopperFz / 2.));
    sol_CopperFrame4[1] =
        new G4UnionSolid("CopperFrame + CoppperPlateB", sol_CopperFrame4[0], sol_CopDownDiskB_Re, 0,
                         G4ThreeVector(-1.9 * cm / 2 - 2 * mm, 0, -CopperFz / 4. - CopperFz / 2.));

    for (int nn = 0; nn < 6; nn++) { //  5 layers

        sol_CopperFrame4[2] =
            new G4UnionSolid("CopperTop + CoppperBottom", sol_CopperFrame3[3], sol_CopperFrame4[1],
                             0, G4ThreeVector(0, 0, -2 * cmocell_h[nn] + CopperFz));

        copperRods_H = cmocell_h[nn] - (CopperFz);
        CopperframeRodX[nn] =
            new G4Tubs("CopperframeRod", 0, copperRods_r, copperRods_H, 0, 360. * deg);
        sol_CopPlatehole =
            new G4Box("CopperFrame", copperRods_r, copperRods_r, copperRods_H - 4 * mm);

        CopperframeRod[nn] =
            new G4SubtractionSolid("CopperFrame - xB", CopperframeRodX[nn], sol_CopPlatehole, 0,
                                   G4ThreeVector(0, copperRods_r, 0));

        for (int NRod = 0; NRod < 4; NRod++) { //
            rxRod = (CopperFx / 2 - 5) * x_smallbox[NRod];
            ryRod = (CopperFy / 2 - 5) * y_smallbox[NRod];
            if (NRod == 0)
                sol_CopperFrame5[NRod] = new G4UnionSolid(
                    "CopperFrame + CoppperRod", sol_CopperFrame4[2], CopperframeRod[nn], 0,
                    G4ThreeVector(rxRod, ryRod, -cmocell_h[nn] + CopperFz / 2.));
            else
                sol_CopperFrame5[NRod] = new G4UnionSolid(
                    "CopperFrame + CoppperRod", sol_CopperFrame5[NRod - 1], CopperframeRod[nn], 0,
                    G4ThreeVector(rxRod, ryRod, -cmocell_h[nn] + CopperFz / 2.));
        }
        logiSquareDiskUp[nn] = new G4LogicalVolume(sol_CopperFrame5[3], _copper, "CopFrame");

        copyNo           = nn;
        PosCopperframeUp = G4ThreeVector(0, 0, zCup[nn] - CopperFz / 2);
        sprintf(namecuframe, "physCopperFrame%d", nn);
        new G4PVPlacement(nullptr, PosCopperframeUp, logiSquareDiskUp[nn], namecuframe,
                          logiTargetRoom, false, copyNo);
        G4cout << "###   " << logiSquareDiskUp[nn]->GetName() << "  "
               << logiSquareDiskUp[nn]->GetMass() / kg << G4endl;
    }

    //////////////////////////////////////////////////////////
    // Run6 Updated
    //
    G4RotationMatrix *zRotation = new G4RotationMatrix();
    G4Box *PCB2                 = new G4Box("PCB2", 14.9 * mm / 2, 23.2 * mm / 2, 0.115 * mm);
    G4Box *MotEpoxy             = new G4Box("MotEpoxy", 27 * mm / 2, 50 * mm / 2, 1 * mm / 2);
    G4Box *Epoxy2               = new G4Box("Epoxy2", 0.643 * mm, 0.643 * mm, 0.46 * mm / 2);
    G4Box *Epoxy3               = new G4Box("Epoxy3", 0.87 * mm, 0.475 * mm, 0.46 * mm / 2);
    G4Box *Solder3              = new G4Box("Solder3", 0.7 * mm, 2.37 * mm, 0.115 * mm / 2);
    G4Box *Solder4              = new G4Box("Solder4", 1.84 * mm, 0.9 * mm, 0.115 * mm / 2);
    G4VSolid *SPCB2;
    G4VSolid *UPCB2;
    G4VSolid *UPCB2a;
    G4LogicalVolume *lEpoxy2   = new G4LogicalVolume(Epoxy2, _solder, "Epoxy2");
    G4LogicalVolume *lEpoxy3   = new G4LogicalVolume(Epoxy3, _solder, "Epoxy3");
    G4LogicalVolume *lMotEpoxy = new G4LogicalVolume(MotEpoxy, _vacuum, "MotEpoxy");
    G4LogicalVolume *lSolder3  = new G4LogicalVolume(Solder3, _solder, "Solder3");
    G4LogicalVolume *lSolder4  = new G4LogicalVolume(Solder4, _solder, "Solder4");

    G4VPhysicalVolume *physMotEpoxy;

    char namePCB2[20];
    char nameEpoxy2[20];
    char nameEpoxy3[20];
    char nameSolder3[20];
    char nameSolder4[20];
    char nameMotEpoxy[20];

    G4Box *S1 = new G4Box("S1", 9 * mm / 2, 4.5 * mm / 2, 0.115 * mm + 2 * mm);
    SPCB2     = new G4SubtractionSolid("SPCB2", PCB2, S1, 0, G4ThreeVector(3.3 * mm, -9.5 * mm, 0));

    G4Box *PCB2_b = new G4Box("PCB2_b", 16.3 * mm / 2, 12.4 * mm / 2, 0.115 * mm);
    zRotation->rotateZ(125 * deg);
    UPCB2 = new G4UnionSolid("PCB2+PCB2b", SPCB2, PCB2_b, zRotation,
                             G4ThreeVector(2.34 * mm, 14.7 * mm, 0));

    G4Box *PCB2_a = new G4Box("PCB2_a", 14.9 * mm / 2, 11.2 * mm / 2, 0.115 * mm);
    UPCB2a =
        new G4UnionSolid("UPCB2+PCB2a", UPCB2, PCB2_a, 0, G4ThreeVector(4.67 * mm, 23.4 * mm, 0));
    G4LogicalVolume *lPCB2 = new G4LogicalVolume(UPCB2a, _copper, "logiPCB", 0, 0, 0);

    for (int nn = 1; nn < 7; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);

        sprintf(nameMotEpoxy, "MotEpoxy_%d", nn);
        physMotEpoxy = new G4PVPlacement(
            0, G4ThreeVector(-1.45 * cm / 2 - 2 * mm, 15.3, PhotonDetZ + 7.1 * mm - 0.058 * cm),
            lMotEpoxy, nameMotEpoxy, logiTargetRoom, false, copyNo);
        (void)physMotEpoxy;
    }

    sprintf(namePCB2, "PCB2");
    new G4PVPlacement(0, G4ThreeVector(0, -5, -0.079), lPCB2, namePCB2, lMotEpoxy, false, copyNo);

    sprintf(nameEpoxy3, "Epoxy3");
    new G4PVPlacement(0, G4ThreeVector(1.5, -11, 0.268 * mm), lEpoxy3, nameEpoxy3, lMotEpoxy, false,
                      copyNo);

    sprintf(nameEpoxy2, "Epoxy2");
    new G4PVPlacement(0, G4ThreeVector(-4.1, -15, 0.268 * mm), lEpoxy2, nameEpoxy2, lMotEpoxy,
                      false, copyNo);

    sprintf(nameSolder3, "Solder3");
    new G4PVPlacement(0, G4ThreeVector(3.5, 0, -0.256 * mm), lSolder3, nameSolder3, lMotEpoxy,
                      false, copyNo);

    sprintf(nameSolder4, "Solder4");
    new G4PVPlacement(0, G4ThreeVector(-3.5, 3, -0.256 * mm), lSolder4, nameSolder4, lMotEpoxy,
                      false, copyNo);

    G4VisAttributes *logiPCB2_VisAtt = new G4VisAttributes(G4Colour::Gray());
    logiPCB2_VisAtt->SetVisibility(true);
    logiPCB2_VisAtt->SetForceSolid(true);
    lPCB2->SetVisAttributes(logiPCB2_VisAtt);

    G4VisAttributes *logiMotEpoxy_VisAtt = new G4VisAttributes(G4Colour::Green());
    logiMotEpoxy_VisAtt->SetVisibility(true);
    logiMotEpoxy_VisAtt->SetForceWireframe(true);
    lMotEpoxy->SetVisAttributes(logiMotEpoxy_VisAtt);

    G4VisAttributes *logiEpoxy2_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiEpoxy2_VisAtt->SetVisibility(true);
    logiEpoxy2_VisAtt->SetForceSolid(true);
    lEpoxy2->SetVisAttributes(logiEpoxy2_VisAtt);

    G4VisAttributes *logiEpoxy3_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiEpoxy3_VisAtt->SetVisibility(true);
    logiEpoxy3_VisAtt->SetForceSolid(true);
    lEpoxy3->SetVisAttributes(logiEpoxy3_VisAtt);

    G4VisAttributes *logiSolder3_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiSolder3_VisAtt->SetVisibility(true);
    logiSolder3_VisAtt->SetForceSolid(true);
    lSolder3->SetVisAttributes(logiSolder3_VisAtt);

    G4VisAttributes *logiSolder4_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiSolder4_VisAtt->SetVisibility(true);
    logiSolder4_VisAtt->SetForceSolid(true);
    lSolder4->SetVisAttributes(logiSolder4_VisAtt);

    G4cout << "###   " << lEpoxy2->GetName() << "  " << lEpoxy2->GetMass() / kg << G4endl;

    G4cout << "###   " << lEpoxy3->GetName() << "  " << lEpoxy3->GetMass() / kg << G4endl;

    G4cout << "###   " << lSolder3->GetName() << "  " << lSolder3->GetMass() / kg << G4endl;

    G4cout << "###   " << lSolder4->GetName() << "  " << lSolder4->GetMass() / kg << G4endl;

    G4cout << "###   " << lPCB2->GetName() << "  " << lPCB2->GetMass() / kg << G4endl;

    ///////////////////////////////////////////////////////
    // Stabilization heater
    //

    G4Box *mSH = new G4Box("mSH", 5 * mm / 2, 5 * mm / 2, 0.49 * mm / 2);
    G4Box *SH  = new G4Box("SH", 3.7 * mm / 2, 3.8 * mm / 2, 0.38 * mm / 2);
    G4Box *eSH = new G4Box("eSH", 0.04 * mm / 2, 2.695 * mm / 2, 100 * nm / 2);

    G4VPhysicalVolume *physmSH;

    G4LogicalVolume *lmSH = new G4LogicalVolume(mSH, _vacuum, "logimSH");
    G4LogicalVolume *leSH = new G4LogicalVolume(eSH, _AuPd, "logieSH");
    G4LogicalVolume *lSH  = new G4LogicalVolume(SH, _SiWafer, "logiSH");

    char nameMSH[10];
    char nameSH[10];
    char nameESH[200];

    for (int nn = 1; nn < 7; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);

        sprintf(nameMSH, "mSH_%d", nn);
        physmSH = new G4PVPlacement(0, G4ThreeVector(0, 15.5, PhotonDetZ + 9.1 * mm), lmSH, nameMSH,
                                    logiTargetRoom, false, copyNo);
        (void)physmSH;
    }

    sprintf(nameSH, "SH");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0.054), lSH, nameSH, lmSH, false, copyNo);

    sprintf(nameESH, "eSH");
    new G4PVPlacement(0, G4ThreeVector(0.26 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(-0.26 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(0.78 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(-0.78 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(1.3 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(-1.3 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(1.82 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(-1.82 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(2.34 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(-2.34 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(2.96 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);
    new G4PVPlacement(0, G4ThreeVector(-2.86 * mm / 2, 0, -0.137), leSH, nameESH, lmSH, false,
                      copyNo);

    G4VisAttributes *logiSH_VisAtt = new G4VisAttributes(G4Colour::Gray());
    logiSH_VisAtt->SetVisibility(true);
    logiSH_VisAtt->SetForceSolid(true);
    lSH->SetVisAttributes(logiSH_VisAtt);

    G4VisAttributes *logimSH_VisAtt = new G4VisAttributes(G4Colour::Green());
    logimSH_VisAtt->SetVisibility(true);
    logimSH_VisAtt->SetForceWireframe(true);
    lmSH->SetVisAttributes(logimSH_VisAtt);

    G4VisAttributes *logieSH_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logieSH_VisAtt->SetVisibility(true);
    logieSH_VisAtt->SetForceSolid(true);
    leSH->SetVisAttributes(logieSH_VisAtt);

    G4cout << "###   " << lSH->GetName() << "  " << lSH->GetMass() / kg << G4endl;

    G4cout << "###   " << leSH->GetName() << "  " << leSH->GetMass() / kg << G4endl;

    ///////////////////////////////////////////////////////
    //  //Springs of light detector : Clamp Teflon
    //
    G4RotationMatrix *ClampRotation = new G4RotationMatrix();
    ClampRotation->rotateZ(145 * deg);
    G4Box *Clamp0    = new G4Box("Clamp0", 5.5 * mm / 2, 5.5 * mm / 2, 2 * mm / 2);
    G4Tubs *ClampH   = new G4Tubs("ClampHole", 0, 3.2 * mm / 2, 2 * mm / 2 + 2, 0, 360 * deg);
    G4VSolid *Clamp  = new G4SubtractionSolid("Clamp", Clamp0, ClampH, 0, G4ThreeVector(0, 0, 0));
    G4Box *Clamp1    = new G4Box("Clamp0", 2 * mm / 2, 7 * mm / 2, 2 * mm / 2);
    G4VSolid *Clamp2 = new G4UnionSolid("Clamp0+Clamp1", Clamp, Clamp1, 0, G4ThreeVector(0, -5, 0));
    G4LogicalVolume *lClamp = new G4LogicalVolume(Clamp2, _teflon, "logiClamp", 0, 0, 0);

    G4Box *Clamp0_a  = new G4Box("Clamp0_a", 5.5 * mm / 2, 5.5 * mm / 2, 2 * mm / 2);
    G4Tubs *ClampH_a = new G4Tubs("ClampHole_a", 0, 3.2 * mm / 2, 2 * mm / 2 + 2, 0, 360 * deg);
    G4VSolid *Clamp_a =
        new G4SubtractionSolid("Clamp_a", Clamp0_a, ClampH_a, 0, G4ThreeVector(0, 0, 0));
    G4Box *Clamp1_a = new G4Box("Clamp0_a", 2 * mm / 2, 7 * mm / 2, 2 * mm / 2);
    G4VSolid *Clamp2_a =
        new G4UnionSolid("Clamp0_a+Clamp1_a", Clamp_a, Clamp1_a, 0, G4ThreeVector(0, 5, 0));
    G4LogicalVolume *lClamp_a = new G4LogicalVolume(Clamp2_a, _teflon, "logiClamp", 0, 0, 0);

    char nameClamp[20];
    char nameClamp_a[20];

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        sprintf(nameClamp, "Clamp_%d", nn);
        new G4PVPlacement(ClampRotation, G4ThreeVector(19.5, -26.5, PhotonDetZ - 1.), lClamp,
                          nameClamp, logiTargetRoom, false, copyNo);
    }

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        sprintf(nameClamp_a, "ClampA_%d", nn);
        new G4PVPlacement(ClampRotation, G4ThreeVector(-19.5, 26.5, PhotonDetZ - 1.), lClamp_a,
                          nameClamp_a, logiTargetRoom, false, copyNo);
    }

    G4VisAttributes *logiClamp_VisAtt = new G4VisAttributes(G4Colour::White());
    logiClamp_VisAtt->SetVisibility(true);
    logiClamp_VisAtt->SetForceSolid(true);
    lClamp->SetVisAttributes(logiClamp_VisAtt);

    G4VisAttributes *logiClampA_VisAtt = new G4VisAttributes(G4Colour::White());
    logiClampA_VisAtt->SetVisibility(true);
    logiClampA_VisAtt->SetForceSolid(true);
    lClamp_a->SetVisAttributes(logiClampA_VisAtt);

    G4cout << "###   " << lClamp->GetName() << "  " << lClamp->GetMass() / kg << G4endl;
    G4cout << "###   " << lClamp_a->GetName() << "  " << lClamp_a->GetMass() / kg << G4endl;

    G4RotationMatrix *TClampRotation = new G4RotationMatrix();
    TClampRotation->rotateZ(130 * deg);
    G4Box *TClamp0  = new G4Box("TClamp0", 5.5 * mm / 2, 5.5 * mm / 2, 2 * mm / 2);
    G4Tubs *TClampH = new G4Tubs("TClampHole", 0, 3.2 * mm / 2, 2 * mm / 2 + 2, 0, 360 * deg);
    G4VSolid *TClamp =
        new G4SubtractionSolid("TClamp", TClamp0, TClampH, 0, G4ThreeVector(0, 0, 0));
    G4Box *TClamp1 = new G4Box("TClamp0", 7 * mm / 2, 2 * mm / 2, 2 * mm / 2);
    G4VSolid *TClamp2 =
        new G4UnionSolid("TClamp0+TClamp1", TClamp, TClamp1, 0, G4ThreeVector(5, 0, 0));
    G4LogicalVolume *lTClamp = new G4LogicalVolume(TClamp2, _teflon, "logiTClamp", 0, 0, 0);

    G4Box *TClamp0_a  = new G4Box("TClamp0_a", 5.5 * mm / 2, 5.5 * mm / 2, 2 * mm / 2);
    G4Tubs *TClampH_a = new G4Tubs("TClampHole_a", 0, 3.2 * mm / 2, 2 * mm / 2 + 2, 0, 360 * deg);
    G4VSolid *TClamp_a =
        new G4SubtractionSolid("TClamp_a", TClamp0_a, TClampH_a, 0, G4ThreeVector(0, 0, 0));
    G4Box *TClamp1_a = new G4Box("TClamp0_a", 7 * mm / 2, 2 * mm / 2, 2 * mm / 2);
    G4VSolid *TClamp2_a =
        new G4UnionSolid("TClamp0_a+TClamp1_a", TClamp_a, TClamp1_a, 0, G4ThreeVector(-5, 0, 0));
    G4LogicalVolume *lTClamp_a = new G4LogicalVolume(TClamp2_a, _teflon, "logiTClamp", 0, 0, 0);

    char nameTClamp[20];
    char nameTClamp_a[20];

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        sprintf(nameTClamp, "TClamp_%d", nn);
        new G4PVPlacement(TClampRotation, G4ThreeVector(19.5, 26.5, PhotonDetZ - 1.), lTClamp,
                          nameTClamp, logiTargetRoom, false, copyNo);
    }

    for (int nn = 0; nn < 6; nn++) {
        copyNo       = nn;
        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);
        sprintf(nameTClamp_a, "TClampA_%d", nn);
        new G4PVPlacement(TClampRotation, G4ThreeVector(-19.5, -26.5, PhotonDetZ - 1.), lTClamp_a,
                          nameTClamp_a, logiTargetRoom, false, copyNo);
    }

    G4VisAttributes *logiTClamp_VisAtt = new G4VisAttributes(G4Colour::White());
    logiTClamp_VisAtt->SetVisibility(true);
    logiTClamp_VisAtt->SetForceSolid(true);
    lTClamp->SetVisAttributes(logiTClamp_VisAtt);

    G4VisAttributes *logiTClampA_VisAtt = new G4VisAttributes(G4Colour::White());
    logiTClampA_VisAtt->SetVisibility(true);
    logiTClampA_VisAtt->SetForceSolid(true);
    lTClamp_a->SetVisAttributes(logiTClampA_VisAtt);

    G4cout << "###   " << lTClamp->GetName() << "  " << lTClamp->GetMass() / kg << G4endl;

    G4cout << "###   " << lTClamp_a->GetName() << "  " << lTClamp_a->GetMass() / kg << G4endl;

    ///////////////////////////////////////////////////////
    // GOLD FILM
    //

    /* CopyNo Scheme (TGSD)
     * 0 ~ 5: for CMO
     * 5 ~ 23: for GOLDa
     * 24 ~ 29: for GOLDb
     * */
    // one gold film below crystal
    char nameGOLDsb[20];

    G4Tubs *GOLDb;
    G4LogicalVolume *logiGOLDb;

    GOLDb            = new G4Tubs("GOLDsb", 0, dGOLDb, zGOLDb, 0, 360 * deg);
    logiGOLDb        = new G4LogicalVolume(GOLDb, _gold, "GOLDb");
    fPilot_logiGOLDb = logiGOLDb;

    for (int nn = 1; nn < 7; nn++) {
        copyNo = nn + 23;
        sprintf(nameGOLDsb, "GOLDsb_%d", nn);

        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);

        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, PhotonDetZ + 9.24 * mm), logiGOLDb,
                          nameGOLDsb, logiTargetRoom, false, copyNo);
    }

    G4VisAttributes *logiGOLDb_VisAtt = new G4VisAttributes(G4Colour::Yellow());
    logiGOLDb_VisAtt->SetVisibility(true);
    logiGOLDb_VisAtt->SetForceSolid(true);
    logiGOLDb->SetVisAttributes(logiGOLDb_VisAtt);

    // 3 gold film above crystal
    char nameGOLDsa[20];

    G4Tubs *GOLDa;
    G4LogicalVolume *logiGOLDa;

    GOLDa            = new G4Tubs("GOLsDa", 0, dGOLDa, zGOLDa, 0, 360 * deg);
    logiGOLDa        = new G4LogicalVolume(GOLDa, _gold, "GOLDa");
    fPilot_logiGOLDa = logiGOLDa;

    for (int nn = 0; nn < 6; nn++) {
        copyNo = 6 + nn;
        sprintf(nameGOLDsa, "GOLDsa_%d", nn);

        PhotonDetZ   = cmocell_h[nn] + PhotonDet_zsize + zCell[nn];
        PosPhotonDet = G4ThreeVector(0, 0, PhotonDetZ);

        new G4PVPlacement(nullptr, G4ThreeVector(0, 7 * mm, 0.11 * mm), logiGOLDa, nameGOLDsa,
                          logiPhotonDet, false, copyNo);

        new G4PVPlacement(nullptr, G4ThreeVector(-4.85 * mm, 0 * mm, 0.11 * mm), logiGOLDa,
                          nameGOLDsa, logiPhotonDet, false, copyNo + 1);

        new G4PVPlacement(nullptr, G4ThreeVector(4.85 * mm, 0 * mm, 0.11 * mm), logiGOLDa,
                          nameGOLDsa, logiPhotonDet, false, copyNo + 2);
    }

    G4VisAttributes *logiGOLDa_VisAtt = new G4VisAttributes(orange);
    logiGOLDa_VisAtt->SetVisibility(true);
    logiGOLDa_VisAtt->SetForceSolid(true);
    logiGOLDa->SetVisAttributes(logiGOLDa_VisAtt);

    ///////////////////////////////////////////////////////
    // Reflector around the CMO (Vm2000 shield)
    //
    G4VSolid *Reflector[7];
    G4Box *Reflector_in[7];
    G4Box *Reflector_out[7];
    G4Box *Reflector_c;
    G4VSolid *Reflector_b[7];
    G4ThreeVector PosRef;
    char namerefl[20];

    G4LogicalVolume *logiReflector[7];
    G4VisAttributes *logiReflectorVis[7];
    G4VPhysicalVolume *physReflector[7];

    G4VSolid *Reflector_d;
    G4Tubs *Reflector_h;

    for (int nn = 0; nn < 6; nn++) {
        sprintf(namerefl, "Reflector%d", nn);

        Reflector_in[nn] = new G4Box(namerefl, reflector_inx, reflector_iny,
                                     cmocell_h[nn] + 0.7 * mm); // not Top cover Ver.

        Reflector_out[nn] = new G4Box(namerefl, reflector_x, reflector_y, cmocell_h[nn] + 0.5 * mm);
        Reflector_b[nn] = new G4SubtractionSolid("Reflector_out - Reflector_in", Reflector_out[nn],
                                                 Reflector_in[nn]);

        Reflector_h = new G4Tubs(namerefl, 0, 1 * cm, cmocell_h[nn], 0, 360. * deg);
        Reflector_c = new G4Box(namerefl, reflector_x, reflector_y, reflector_thick / 2);
        Reflector_d =
            new G4SubtractionSolid("Reflector_out - Reflector_in", Reflector_c, Reflector_h);

        Reflector[nn] = new G4UnionSolid(
            "Reflector_out - Reflector_in", Reflector_b[nn], Reflector_d, 0,
            G4ThreeVector(0, 0, -cmocell_h[nn] - 0.5 * mm - reflector_thick / 2 + 0.0001 * cm));

        sprintf(namerefl, "logiReflector%d", nn);
        logiReflector[nn]    = new G4LogicalVolume(Reflector[nn], _vm2000, namerefl);
        logiReflectorVis[nn] = new G4VisAttributes(G4Colour(0.75, 0.0, 0.75, 0.3));

        if (nn == 1) {
            sprintf(namerefl, "Reflector%d", nn);

            Reflector_in[nn] = new G4Box(namerefl, reflector_inx2, reflector_iny2,
                                         cmocell_h[nn] + 1.5 * mm); // not Top cover Ver.

            Reflector_out[nn] =
                new G4Box(namerefl, reflector_x, reflector_y, cmocell_h[nn] + 0.5 * mm);
            Reflector_b[nn] = new G4SubtractionSolid("Reflector_out - Reflector_in",
                                                     Reflector_out[nn], Reflector_in[nn]);

            Reflector_h = new G4Tubs(namerefl, 0, 1 * cm, cmocell_h[nn], 0, 360. * deg);
            Reflector_c = new G4Box(namerefl, reflector_x, reflector_y, reflector_thick2 / 4);
            Reflector_d =
                new G4SubtractionSolid("Reflector_out - Reflector_in", Reflector_c, Reflector_h);

            Reflector[nn] = new G4UnionSolid(
                "Reflector_out - Reflector_in", Reflector_b[nn], Reflector_d, 0,
                G4ThreeVector(0, 0,
                              -cmocell_h[nn] - 0.5 * mm - reflector_thick2 / 4 + 0.0001 * cm));

            sprintf(namerefl, "logiReflector%d", nn);
            logiReflector[nn]    = new G4LogicalVolume(Reflector[nn], _teflon2, namerefl);
            logiReflectorVis[nn] = new G4VisAttributes(G4Colour::White());
        }

        PosRef = G4ThreeVector(xCell[nn], yCell[nn], zCell[nn]);
        sprintf(namerefl, "physReflector%d", nn);
        new G4PVPlacement(nullptr, PosRef, logiReflector[nn], namerefl, logiTargetRoom, false, 0);
        G4cout << "###   " << logiReflector[nn]->GetName() << "  "
               << logiReflector[nn]->GetMass() / kg << G4endl;

        logiReflectorVis[nn]->SetVisibility(true);
        logiReflectorVis[nn]->SetForceSolid(true);
        logiReflector[nn]->SetVisAttributes(logiReflectorVis[nn]);
    }

    // -----------------------------------
    // H beam construction
    // -----------------------------------
    //
    /// Vertical H beam
    //

    Bplate = new G4Box("Bplate", Bplate_size, Bplate_size, Bplate_thick_half);
    Hbeam1 = new G4ExtrudedSolid("Hbeam1", large_points, Hbeam1_length, 0.0, 1.0, 0.0, 1.0);
    Hbeam2 = new G4ExtrudedSolid("Hbeam2", large_points, Hbeam2_length, 1.0, 1.0, 1.0, 1.0);

    UnionHbeam1 = new G4UnionSolid("UnionHbeam1", Hbeam1, Bplate, 0, -posBplate1);
    UnionHbeam2 = new G4UnionSolid("UnionHbeam2", Hbeam2, Bplate, 0, -posBplate2);

    logicHbeam1 = new G4LogicalVolume(UnionHbeam1, _steel, "logicHbeam1", 0, 0, 0);
    logicHbeam2 = new G4LogicalVolume(UnionHbeam2, _steel, "logicHbeam2", 0, 0, 0);

    Hb1VisAtt = new G4VisAttributes(yellow);
    Hb2VisAtt = new G4VisAttributes(yellow);

    logicHbeam1->SetVisAttributes(Hb1VisAtt);
    logicHbeam2->SetVisAttributes(Hb2VisAtt);

    /// Horizontal H -beam X
    Ho_beamX =
        new G4ExtrudedSolid("Ho_beamX", small_points, Horizontal_beamX_length, 1.0, 1.0, 1.0, 1.0);

    logicHo_beamX = new G4LogicalVolume(Ho_beamX, _steel, "logicHo_beamX");

    HoBXVisAtt = new G4VisAttributes(green);
    logicHo_beamX->SetVisAttributes(HoBXVisAtt);

    Ho_beamXLar = new G4ExtrudedSolid("Ho_beamXLar", large_points, Horizontal_beamXLarge_length,
                                      1.0, 1.0, 1.0, 1.0);
    logicHo_beamXLar = new G4LogicalVolume(Ho_beamXLar, _steel, "logicHo_beamXLar");

    HoBXLVisAtt = new G4VisAttributes(yellow);
    logicHo_beamXLar->SetVisAttributes(HoBXLVisAtt);
    // HoBXLVisAtt->SetForceSolid(true);

    G4ThreeVector Gantry_Top_Trans(0, 0, -Hbeam1_length - 20. + HoLayer[2] - 75. + posMountOrigin1);

    /// Horizontal H -beam Y

    Ho_beamY =
        new G4ExtrudedSolid("Ho_beamY", small_points, Horizontal_beamY_length, 1.0, 1.0, 1.0, 1.0);

    logicHo_beamY = new G4LogicalVolume(Ho_beamY, _steel, "logicHo_beamY");

    HoBYVisAtt = new G4VisAttributes(green);
    logicHo_beamY->SetVisAttributes(HoBYVisAtt);
    // HoBYVisAtt->SetForceSolid(true);

    Ho_beamYLar = new G4ExtrudedSolid("Ho_beamYLar", large_points, Horizontal_beamYLarge_length,
                                      1.0, 1.0, 1.0, 1.0);

    logicHo_beamYLar = new G4LogicalVolume(Ho_beamYLar, _steel, "logicHo_beamYLar");

    HoBYLVisAtt = new G4VisAttributes(yellow);
    logicHo_beamYLar->SetVisAttributes(HoBYLVisAtt);

    //////////////////////////////
    // Build Scintillator
    //////////////////////////////

    //////////////////////////////
    // Build Solid, LV and PV
    //////////////////////////////

    topScint_MotherBox =
        new G4Box("Muon_Top_Scint_MotherBox", topScint_boxWidth_half + scint_reflector_thick,
                  topScint_boxHeight_half + scint_reflector_thick,
                  topScint_boxThick_half + scint_reflector_thick);

    topScint_MotherTrap =
        new G4Trd("Muon_Top_Scintillator_MS_PreTrap1",
                  topScint_trapX2_half - topScint_TrapSlope_double * (topScint_trapFlatL_half +
                                                                      scint_boolean_tolerence / 2.),
                  topScint_trapX2_half + topScint_TrapSlope_double * scint_boolean_tolerence / 2.,
                  topScint_trapThick1_half + scint_reflector_thick,
                  topScint_trapThick1_half + scint_reflector_thick,
                  topScint_trapFlatL_half + scint_boolean_tolerence);
    topScint_MotherPMTTrap =
        new G4Trd("Muon_Top_Scintillator_MS_PreTrap2",
                  topScint_trapX1_half + topScint_TrapSlope_double * scint_reflector_thick / 2.,
                  topScint_trapX2_half - topScint_TrapSlope_double * topScint_trapFlatL_half,
                  topScint_trapThick2_half + scint_reflector_thick,
                  topScint_trapThick1_half + scint_reflector_thick,
                  topScint_trapZ_half - topScint_trapFlatL_half - scint_reflector_thick / 2.);
    topScint_MotherTrap = new G4UnionSolid(
        "Muon_Top_Scintillator_MS_Trap", topScint_MotherTrap, topScint_MotherPMTTrap, 0,
        G4ThreeVector(0, 0, -topScint_trapZ_half + scint_reflector_thick / 2.));

    topScint_MotherSolid = new G4UnionSolid(
        "Muon_Top_Mother_Solid_Temp1", topScint_MotherBox, topScint_MotherTrap,
        scint_alignMtx_Uside,
        G4ThreeVector(-topScint_trapX2_half,
                      topScint_boxHeight_half + scint_reflector_thick + topScint_trapFlatL_half,
                      0));
    topScint_MotherSolid = new G4UnionSolid(
        "Muon_Top_Mother_Solid_Temp2", topScint_MotherSolid, topScint_MotherTrap,
        scint_alignMtx_Uside,
        G4ThreeVector(+topScint_trapX2_half,
                      topScint_boxHeight_half + scint_reflector_thick + topScint_trapFlatL_half,
                      0));
    topScint_MotherSolid = new G4UnionSolid(
        "Muon_Top_Mother_Mother_Solid_Temp3", topScint_MotherSolid, topScint_MotherTrap,
        scint_alignMtx_Dside,
        G4ThreeVector(-topScint_trapX2_half,
                      -topScint_boxHeight_half - scint_reflector_thick - topScint_trapFlatL_half,
                      0));
    topScint_MotherSolid = new G4UnionSolid(
        "Muon_Top_Mother_Solid", topScint_MotherSolid, topScint_MotherTrap, scint_alignMtx_Dside,
        G4ThreeVector(+topScint_trapX2_half,
                      -topScint_boxHeight_half - scint_reflector_thick - topScint_trapFlatL_half,
                      0));

    topScint_Box = new G4Box("Muon_Top_Scintillator Box", topScint_boxWidth_half,
                             topScint_boxHeight_half, topScint_boxThick_half);
    topScint_FlatTrap =
        new G4Trd("Muon_Top_Scintillator Flat Trd",
                  topScint_trapX2_half - topScint_TrapSlope_double * topScint_trapFlatL_half,
                  topScint_trapX2_half, topScint_trapThick1_half, topScint_trapThick1_half,
                  topScint_trapFlatL_half);
    topScint_PMTTrap =
        new G4Trd("Muon_Top_Scintillator PMT Trd", topScint_trapX1_half,
                  topScint_trapX2_half - (topScint_trapX2_half - topScint_trapX1_half) /
                                             topScint_trapZ_half * topScint_trapFlatL_half,
                  topScint_trapThick2_half, topScint_trapThick1_half,
                  topScint_trapZ_half - topScint_trapFlatL_half);

    topScint_MotherLogical =
        new G4LogicalVolume(topScint_MotherSolid, _air, "Muon_Top_Scintillator_Logical");
    new G4LogicalSkinSurface("Muon_Top_Scintillator_reflector_opsurf", topScint_MotherLogical,
                             Tyvek_opsurf);

    fPilot_TopScint_BoxLogical =
        new G4LogicalVolume(topScint_Box, _vinylt, "Muon_Top_Scintillator Box Logical");
    fPilot_TopScint_PMTTrapLogical =
        new G4LogicalVolume(topScint_PMTTrap, _vinylt, "Muon_Top_Scintillator PMT Trd Logical");
    fPilot_TopScint_FlatTrapLogical =
        new G4LogicalVolume(topScint_FlatTrap, _vinylt, "Muon_Top_Scintillator Flat Trd Logical");

    new G4PVPlacement(nullptr, G4ThreeVector(), fPilot_TopScint_BoxLogical,
                      "Muon_Top_Scintillator Box PV", topScint_MotherLogical, false, 0, false);

    new G4PVPlacement(
        scint_alignMtx_Uside,
        G4ThreeVector(-topScint_trapX2_half, topScint_boxHeight_half + topScint_trapFlatL_half, 0),
        fPilot_TopScint_FlatTrapLogical, "Muon_Top_Scintillator FlatTrd PV", topScint_MotherLogical,
        false, 0, false);
    new G4PVPlacement(
        scint_alignMtx_Uside,
        G4ThreeVector(-topScint_trapX2_half,
                      topScint_boxHeight_half + topScint_trapZ_half + topScint_trapFlatL_half, 0),
        fPilot_TopScint_PMTTrapLogical, "Muon_Top_Scintillator PMT Trd PV", topScint_MotherLogical,
        false, 0, false);

    new G4PVPlacement(
        scint_alignMtx_Uside,
        G4ThreeVector(topScint_trapX2_half, topScint_boxHeight_half + topScint_trapFlatL_half, 0),
        fPilot_TopScint_FlatTrapLogical, "Muon_Top_Scintillator FlatTrd PV", topScint_MotherLogical,
        false, 1, false);
    new G4PVPlacement(
        scint_alignMtx_Uside,
        G4ThreeVector(topScint_trapX2_half,
                      topScint_boxHeight_half + topScint_trapZ_half + topScint_trapFlatL_half, 0),
        fPilot_TopScint_PMTTrapLogical, "Muon_Top_Scintillator PMT Trd PV", topScint_MotherLogical,
        false, 1, false);

    new G4PVPlacement(
        scint_alignMtx_Dside,
        G4ThreeVector(-topScint_trapX2_half, -topScint_boxHeight_half - topScint_trapFlatL_half, 0),
        fPilot_TopScint_FlatTrapLogical, "Muon_Top_Scintillator FlatTrd PV", topScint_MotherLogical,
        false, 2, false);
    new G4PVPlacement(
        scint_alignMtx_Dside,
        G4ThreeVector(-topScint_trapX2_half,
                      -topScint_boxHeight_half - topScint_trapZ_half - topScint_trapFlatL_half, 0),
        fPilot_TopScint_PMTTrapLogical, "Muon_Top_Scintillator PMT Trd PV", topScint_MotherLogical,
        false, 2, false);

    new G4PVPlacement(
        scint_alignMtx_Dside,
        G4ThreeVector(topScint_trapX2_half, -topScint_boxHeight_half - topScint_trapFlatL_half, 0),
        fPilot_TopScint_FlatTrapLogical, "Muon_Top_Scintillator FlatTrd PV", topScint_MotherLogical,
        false, 3, false);
    new G4PVPlacement(
        scint_alignMtx_Dside,
        G4ThreeVector(topScint_trapX2_half,
                      -topScint_boxHeight_half - topScint_trapZ_half - topScint_trapFlatL_half, 0),
        fPilot_TopScint_PMTTrapLogical, "Muon_Top_Scintillator PMT Trd PV", topScint_MotherLogical,
        false, 3, false);

    sideFBScint_MotherBox  = new G4Box("Muon_FBSide_Scintillator_MotherBox",
                                      sideFBScint_boxWidth_half + scint_reflector_thick,
                                      sideFBScint_boxHeight_half + scint_reflector_thick,
                                      sideFBScint_boxThick_half + scint_reflector_thick);
    sideFBScint_MotherTrap = new G4Trd(
        "Muon_FBSide_Scintillator_MS_PreTrap1",
        sideFBScint_trapX2_half - sideFBScint_TrapSlope_double *
                                      (sideFBScint_trapFlatL_half + scint_boolean_tolerence / 2.),
        sideFBScint_trapX2_half + sideFBScint_TrapSlope_double * scint_boolean_tolerence / 2.,
        sideFBScint_trapThick1_half + scint_reflector_thick,
        sideFBScint_trapThick1_half + scint_reflector_thick,
        sideFBScint_trapFlatL_half + scint_boolean_tolerence);
    sideFBScint_MotherPMTTrap = new G4Trd(
        "Muon_FBSide_Scintillator_MS_PreTrap2",
        sideFBScint_trapX1_half + sideFBScint_TrapSlope_double * scint_reflector_thick / 2.,
        sideFBScint_trapX2_half - sideFBScint_TrapSlope_double * sideFBScint_trapFlatL_half,
        sideFBScint_trapThick2_half + scint_reflector_thick,
        sideFBScint_trapThick1_half + scint_reflector_thick,
        sideFBScint_trapZ_half - sideFBScint_trapFlatL_half - scint_reflector_thick / 2.);
    sideFBScint_MotherTrap = new G4UnionSolid(
        "Muon_FBSide_Scintillator_MS_Trap", sideFBScint_MotherTrap, sideFBScint_MotherPMTTrap, 0,
        G4ThreeVector(0, 0, -sideFBScint_trapZ_half + scint_reflector_thick / 2.));

    sideFBScint_MotherSolid = new G4UnionSolid(
        "Muon_FBSide_Mother_Solid_Temp", sideFBScint_MotherBox, sideFBScint_MotherTrap,
        scint_alignMtx_Uside,
        G4ThreeVector(
            sideFBScint_trapX2_half,
            sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half + scint_reflector_thick, 0));
    sideFBScint_MotherSolid = new G4UnionSolid(
        "Muon_FBSide_Mother_Solid", sideFBScint_MotherSolid, sideFBScint_MotherTrap,
        scint_alignMtx_Uside,
        G4ThreeVector(
            -sideFBScint_trapX2_half,
            sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half + scint_reflector_thick, 0));

    sideFBScint_Box      = new G4Box("Muon Side Scintillator FB Box", sideFBScint_boxWidth_half,
                                sideFBScint_boxHeight_half, sideFBScint_boxThick_half);
    sideFBScint_FlatTrap = new G4Trd("Muon Side Scintillator FB Flat Trd",
                                     sideFBScint_trapX2_half -
                                         sideFBScint_TrapSlope_double * sideFBScint_trapFlatL_half,
                                     sideFBScint_trapX2_half, sideFBScint_trapThick1_half,
                                     sideFBScint_trapThick1_half, sideFBScint_trapFlatL_half);
    sideFBScint_PMTTrap  = new G4Trd("Muon Side Scintillator FB PMT Trd", sideFBScint_trapX1_half,
                                    sideFBScint_trapX2_half -
                                        sideFBScint_TrapSlope_double * sideFBScint_trapFlatL_half,
                                    sideFBScint_trapThick2_half, sideFBScint_trapThick1_half,
                                    sideFBScint_trapZ_half - sideFBScint_trapFlatL_half);

    sideFBScint_MotherLogical =
        new G4LogicalVolume(sideFBScint_MotherSolid, _tyvek, "Muon Side Scintillator FB Logical");
    new G4LogicalSkinSurface("Muon_FBSide_Scintillator_reflector_opsurf", sideFBScint_MotherLogical,
                             Tyvek_opsurf);

    fPilot_SideFBScint_BoxLogical =
        new G4LogicalVolume(sideFBScint_Box, _vinylt, "Muon Side Scintillator FB Box Logical");
    fPilot_SideFBScint_PMTTrapLogical = new G4LogicalVolume(
        sideFBScint_PMTTrap, _vinylt, "Muon Side Scintillator FB PMT Trd Logical");
    fPilot_SideFBScint_FlatTrapLogical = new G4LogicalVolume(
        sideFBScint_FlatTrap, _vinylt, "Muon Side Scintillator FB Flat Trd Logical");

    new G4PVPlacement(nullptr, G4ThreeVector(), fPilot_SideFBScint_BoxLogical,
                      "Muon Scintillator FB Box PV", sideFBScint_MotherLogical, false, 0, false);

    new G4PVPlacement(scint_alignMtx_Uside,
                      G4ThreeVector(-sideFBScint_trapX2_half,
                                    sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half, 0),
                      fPilot_SideFBScint_FlatTrapLogical, "Muon Scintillator FB FlatTrd PV",
                      sideFBScint_MotherLogical, false, 0, false);
    new G4PVPlacement(scint_alignMtx_Uside,
                      G4ThreeVector(-sideFBScint_trapX2_half,
                                    sideFBScint_boxHeight_half + sideFBScint_trapZ_half +
                                        sideFBScint_trapFlatL_half,
                                    0),
                      fPilot_SideFBScint_PMTTrapLogical, "Muon Scintillator FB PMT Trd PV",
                      sideFBScint_MotherLogical, false, 0, false);

    new G4PVPlacement(scint_alignMtx_Uside,
                      G4ThreeVector(sideFBScint_trapX2_half,
                                    sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half, 0),
                      fPilot_SideFBScint_FlatTrapLogical, "Muon Scintillator FB FlatTrd PV",
                      sideFBScint_MotherLogical, false, 1, false);
    new G4PVPlacement(scint_alignMtx_Uside,
                      G4ThreeVector(sideFBScint_trapX2_half,
                                    sideFBScint_boxHeight_half + sideFBScint_trapZ_half +
                                        sideFBScint_trapFlatL_half,
                                    0),
                      fPilot_SideFBScint_PMTTrapLogical, "Muon Scintillator FB PMT Trd PV",
                      sideFBScint_MotherLogical, false, 1, false);

    sideLRScint_MotherBox  = new G4Box("Muon_LRSide_Scintillator_MotherBox",
                                      sideLRScint_boxWidth_half + scint_reflector_thick,
                                      sideLRScint_boxHeight_half + scint_reflector_thick,
                                      sideLRScint_boxThick_half + scint_reflector_thick);
    sideLRScint_MotherTrap = new G4Trd(
        "Muon_LRSide_Scintillator_MS_PreTrap1",
        sideLRScint_trapX2_half - sideLRScint_TrapSlope_double *
                                      (sideLRScint_trapFlatL_half + scint_boolean_tolerence / 2.),
        sideLRScint_trapX2_half + sideLRScint_TrapSlope_double * scint_boolean_tolerence / 2.,
        sideLRScint_trapThick1_half + scint_reflector_thick,
        sideLRScint_trapThick1_half + scint_reflector_thick,
        sideLRScint_trapFlatL_half + scint_boolean_tolerence);
    sideLRScint_MotherPMTTrap = new G4Trd(
        "Muon_LRSide_Scintillator_MS_PreTrap2",
        sideLRScint_trapX1_half + sideLRScint_TrapSlope_double * scint_reflector_thick / 2.,
        sideLRScint_trapX2_half - sideLRScint_TrapSlope_double * sideLRScint_trapFlatL_half,
        sideLRScint_trapThick2_half + scint_reflector_thick,
        sideLRScint_trapThick1_half + scint_reflector_thick,
        sideLRScint_trapZ_half - sideLRScint_trapFlatL_half - scint_reflector_thick / 2.);
    sideLRScint_MotherTrap = new G4UnionSolid(
        "Muon_LRSide_Scintillator_MS_Trap", sideLRScint_MotherTrap, sideLRScint_MotherPMTTrap, 0,
        G4ThreeVector(0, 0, -sideLRScint_trapZ_half + scint_reflector_thick / 2.));

    sideLRScint_MotherSolid = new G4UnionSolid(
        "Muon_LRSide_Scintillator_Solid_Temp1", sideLRScint_MotherBox, sideLRScint_MotherTrap,
        scint_alignMtx_Rside,
        G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                          scint_reflector_thick,
                      -sideLRScint_boxHeight_half + 1. * sideLRScint_trapX2_half, 0));
    sideLRScint_MotherSolid = new G4UnionSolid(
        "Muon_LRSide_Scintillator_Solid_Temp2", sideLRScint_MotherSolid, sideLRScint_MotherTrap,
        scint_alignMtx_Rside,
        G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                          scint_reflector_thick,
                      -sideLRScint_boxHeight_half + 3. * sideLRScint_trapX2_half, 0));
    sideLRScint_MotherSolid = new G4UnionSolid(
        "Muon_LRSide_Scintillator_Solid", sideLRScint_MotherSolid, sideLRScint_MotherTrap,
        scint_alignMtx_Rside,
        G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                          scint_reflector_thick,
                      -sideLRScint_boxHeight_half + 5. * sideLRScint_trapX2_half, 0));

    sideLRScint_Box      = new G4Box("Muon Side Scintillator LR Box", sideLRScint_boxWidth_half,
                                sideLRScint_boxHeight_half, sideLRScint_boxThick_half);
    sideLRScint_FlatTrap = new G4Trd("Muon Side Scintillator LR Flat Trd",
                                     sideLRScint_trapX2_half -
                                         sideLRScint_TrapSlope_double * sideLRScint_trapFlatL_half,
                                     sideLRScint_trapX2_half, sideLRScint_trapThick1_half,
                                     sideLRScint_trapThick1_half, sideLRScint_trapFlatL_half);
    sideLRScint_PMTTrap  = new G4Trd("Muon Side Scintillator LR PMT Trd", sideLRScint_trapX1_half,
                                    sideLRScint_trapX2_half -
                                        sideLRScint_TrapSlope_double * sideLRScint_trapFlatL_half,
                                    sideLRScint_trapThick2_half, sideLRScint_trapThick1_half,
                                    sideLRScint_trapZ_half - sideLRScint_trapFlatL_half);

    sideLRScint_MotherLogical =
        new G4LogicalVolume(sideLRScint_MotherSolid, _tyvek, "Muon Side Scintillator LR Logical");
    new G4LogicalSkinSurface("Muon_LRSide_Scintillator_reflector_opsurf", sideLRScint_MotherLogical,
                             Tyvek_opsurf);

    fPilot_SideLRScint_BoxLogical =
        new G4LogicalVolume(sideLRScint_Box, _vinylt, "Muon Side Scintillator LR Box Logical");
    fPilot_SideLRScint_PMTTrapLogical = new G4LogicalVolume(
        sideLRScint_PMTTrap, _vinylt, "Muon Side Scintillator LR PMT Trd Logical");
    fPilot_SideLRScint_FlatTrapLogical = new G4LogicalVolume(
        sideLRScint_FlatTrap, _vinylt, "Muon Side Scintillator LR Flat Trd Logical");

    new G4PVPlacement(nullptr, G4ThreeVector(), fPilot_SideLRScint_BoxLogical,
                      "Muon Scintillator LR Box PV", sideLRScint_MotherLogical, false, 0, false);

    new G4PVPlacement(scint_alignMtx_Rside,
                      G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half,
                                    -sideLRScint_boxHeight_half + sideLRScint_trapX2_half, 0),
                      fPilot_SideLRScint_FlatTrapLogical, "Muon Scintillator LR FlatTrd PV",
                      sideLRScint_MotherLogical, false, 0, false);
    new G4PVPlacement(scint_alignMtx_Rside,
                      G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                                        sideLRScint_trapZ_half,
                                    -sideLRScint_boxHeight_half + sideLRScint_trapX2_half, 0),
                      fPilot_SideLRScint_PMTTrapLogical, "Muon Scintillator LR PMT Trd PV",
                      sideLRScint_MotherLogical, false, 0, false);

    new G4PVPlacement(scint_alignMtx_Rside,
                      G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half,
                                    -sideLRScint_boxHeight_half + 3. * sideLRScint_trapX2_half, 0),
                      fPilot_SideLRScint_FlatTrapLogical, "Muon Scintillator LR FlatTrd PV",
                      sideLRScint_MotherLogical, false, 1, false);
    new G4PVPlacement(scint_alignMtx_Rside,
                      G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                                        sideLRScint_trapZ_half,
                                    -sideLRScint_boxHeight_half + 3. * sideLRScint_trapX2_half, 0),
                      fPilot_SideLRScint_PMTTrapLogical, "Muon Scintillator LR PMT Trd PV",
                      sideLRScint_MotherLogical, false, 1, false);

    new G4PVPlacement(scint_alignMtx_Rside,
                      G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half,
                                    -sideLRScint_boxHeight_half + 5. * sideLRScint_trapX2_half, 0),
                      fPilot_SideLRScint_FlatTrapLogical, "Muon Scintillator LR FlatTrd PV",
                      sideLRScint_MotherLogical, false, 2, false);
    new G4PVPlacement(scint_alignMtx_Rside,
                      G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                                        sideLRScint_trapZ_half,
                                    -sideLRScint_boxHeight_half + 5. * sideLRScint_trapX2_half, 0),
                      fPilot_SideLRScint_PMTTrapLogical, "Muon Scintillator LR PMT Trd PV",
                      sideLRScint_MotherLogical, false, 2, false);

    /////////////////////////////////////////////////////////
    // PMT
    /////////////////////////////////////////////////////////

    G4cout << "now making & placing PMTs...\n";
    G4VSolid *PMTTube = new G4Tubs("PMTbody0", 0, PMTSizeR, PMTSizeZ / 2., 0, 360. * deg);
    G4VSolid *PMTTop  = new G4Sphere("PMTTop", PMTTopSizeRmin, PMTTopSizeRmax, 0, 360 * degree,
                                    90 * degree, 180 * degree);
    G4VSolid *PMTWindow =
        new G4Tubs("PMTWindow", 0, PMTTopSizeRmax, PMTWindowThickness / 2., 0, 360. * deg);
    G4VSolid *PMTUnion0 =
        new G4UnionSolid("PMTUnion0", PMTTube, PMTTop, 0, G4ThreeVector(0, 0, PMTSizeZ / 2));
    G4VSolid *PMTUnion1 =
        new G4UnionSolid("PMTUnion1", PMTUnion0, PMTWindow, 0,
                         G4ThreeVector(0, 0, PMTSizeZ / 2 + PMTWindowThickness / 2));

    G4VSolid *PMTVacuTube =
        new G4Tubs("PMTVacubody0", 0, PMTVacuSizeR, PMTVacuSizeZ / 2, 0, 360. * deg);
    G4VSolid *PMTVacuTop   = new G4Sphere("PMTVacuTop", PMTVacuTopSizeRmin, PMTVacuTopSizeRmax, 0,
                                        360 * degree, 90 * degree, 180 * degree);
    G4VSolid *solidPMTVacu = new G4UnionSolid("PMTVacu", PMTVacuTube, PMTVacuTop, 0,
                                              G4ThreeVector(0, 0, PMTVacuSizeZ / 2));
    G4VSolid *solidPMT     = PMTUnion1;
    G4VSolid *solidPMTGrease =
        new G4Tubs("PMTGrease", 0, PMTGreaseSizeR, PMTGreaseThickness / 2, 0, 360. * deg);
    G4VSolid *solidPMTGreaseShield =
        new G4Tubs("PMTGS", 0, PMTWindowSizeR, PMTGreaseThickness / 2., 0, 360. * deg);
    G4LogicalVolume *logicPMTGrease =
        new G4LogicalVolume(solidPMTGrease, material_Grease, "PMTGrease");
    G4LogicalVolume *logicPMTGS = new G4LogicalVolume(solidPMTGreaseShield, _aluminum, "PMTGrease");
    // for visual
    G4VisAttributes *PMTGreaseVisAtt = new G4VisAttributes(G4Colour::Gray());
    PMTGreaseVisAtt->SetVisibility(true);
    PMTGreaseVisAtt->SetForceWireframe(true);
    logicPMTGrease->SetVisAttributes(PMTGreaseVisAtt);

    G4LogicalVolume *logicPMT = new G4LogicalVolume(solidPMT, _glass, "Pmt");
    fPilot_logicPMT           = logicPMT;
    // for visual
    G4VisAttributes *PMTVisAtt = new G4VisAttributes(G4Colour::Magenta());
    PMTVisAtt->SetVisibility(true);
    PMTVisAtt->SetForceWireframe(true);
    logicPMT->SetVisAttributes(PMTVisAtt);

    G4LogicalVolume *logicPMTVacu = new G4LogicalVolume(solidPMTVacu, _vacuum, "PMTVacu");
    fPilot_logicPMTVacu           = logicPMTVacu;
    // for visual
    G4VisAttributes *PMTVacuVisAtt = new G4VisAttributes(G4Colour::Brown());
    PMTVacuVisAtt->SetVisibility(true);
    PMTVacuVisAtt->SetForceWireframe(true);
    logicPMTVacu->SetVisAttributes(PMTVacuVisAtt);

    G4double PMTEnvelopeSizeZ = PMTSizeZ + PMTWindowThickness + PMTGreaseThickness;
    G4VSolid *solidPMTEnvelope =
        new G4Tubs("PMTEnvelope", 0, PMTTopSizeRmax, PMTEnvelopeSizeZ / 2., 0, 2 * M_PI);
    G4LogicalVolume *logicPMTEnvelope =
        new G4LogicalVolume(solidPMTEnvelope, _N2_Gas, "PMTEnvelope");
    // for visual
    G4VisAttributes *PMTEnvelopeVisAtt = new G4VisAttributes(G4Colour::Gray());
    PMTEnvelopeVisAtt->SetVisibility(true);
    logicPMTEnvelope->SetVisAttributes(PMTEnvelopeVisAtt);

    G4double zz = -PMTEnvelopeSizeZ / 2. + PMTSizeZ / 2.;
    G4VPhysicalVolume *physiPMT =
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zz), // translation position
                          logicPMT,                         // its logical volume
                          "phys_pmt",                       // its name
                          logicPMTEnvelope,                 // its mother logical volume
                          false,                            // no boolean operation
                          0                                 // its copy number
        );
    zz = -PMTBacksideThickness / 2;
    G4VPhysicalVolume *physiPMTVacu =
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zz), // translation position
                          logicPMTVacu,                     // its logical volume
                          "phys_pmtvacu",                   // its name
                          logicPMT,                         // its mother logical volume
                          false,                            // no boolean operation
                          0                                 // its copy number
        );
    zz = PMTEnvelopeSizeZ / 2 - PMTGreaseThickness / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zz), logicPMTGS, "phys_pmtgs", logicPMTEnvelope,
                      false, 0);

    new G4PVPlacement(nullptr, G4ThreeVector(), logicPMTGrease, "phys_pmtgrease", logicPMTGS, false,
                      0, false);

    /////////////////////////////////////////////////////////////////
    // B4C Rubber (24%)
    ////////////////////////////////////////////////////////////////
    toyRubberB4C_Box = new G4Box("ToyRubberB4C_Box", workboxX + rubber_thickness,
                                 workboxY + rubber_thickness, workboxH + rubber_thickness);
    toyRubberB4C_BoxLogical =
        new G4LogicalVolume(toyRubberB4C_Box, _B4CRubber24perCent, "ToyRubberB4C_LV");

    realPbTopRubberB4C_Box =
        new G4Box("RealPbTopRubberB4C_Box", pbtopbox_size, pbtopbox_size, rubber_thickness / 2.);
    realPbTopRubberB4C_BoxLogical =
        new G4LogicalVolume(realPbTopRubberB4C_Box, _B4CRubber24perCent, "RealPbTopRubberB4C_LV");

    realPbUnderRubberB4C_Box = new G4Box("RealPbUnderRubberB4C_Box", rubber_UnderPb_XYsize / 2.,
                                         rubber_UnderPb_XYsize / 2., rubber_thickness / 2.);
    realPbUnderRubberB4C_BoxLogical = new G4LogicalVolume(
        realPbUnderRubberB4C_Box, _B4CRubber24perCent, "RealPbUnderRubberB4C_LV");

    rubberB4C_BoxVisAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.8));
    rubberB4C_BoxVisAttr->SetForceSolid(true);
    toyRubberB4C_BoxLogical->SetVisAttributes(rubberB4C_BoxVisAttr);
    realPbTopRubberB4C_BoxLogical->SetVisAttributes(rubberB4C_BoxVisAttr);
    realPbUnderRubberB4C_BoxLogical->SetVisAttributes(rubberB4C_BoxVisAttr);

    /////////////////////////////////////////////////////////////////
    // Polyethelene
    ////////////////////////////////////////////////////////////////
    toyPE_Box = new G4Box("ToyPolyEthylene_Box", workboxX + rubber_thickness + toyPE_thickness,
                          workboxY + rubber_thickness + toyPE_thickness,
                          workboxH + rubber_thickness + toyPE_thickness);
    toyPE_BoxLogical = new G4LogicalVolume(toyPE_Box, _polyethylene, "PolyEthylene_LV");
    pe_BoxVisAttr    = new G4VisAttributes(G4Colour(129 / 255., 193 / 255., 71 / 255., 0.9));
    pe_BoxVisAttr->SetForceSolid(true);

    realPE_SideFB_Box = new G4Box("RealPolyEthylene_FB_Box", realPE_thickness1 / 2.,
                                  2. * sideFBScint_boxWidth_half, sideFBScint_boxHeight_half);
    realPE_SideFB_BoxLogical =
        new G4LogicalVolume(realPE_SideFB_Box, _polyethylene, "RealPolyEthylene_FB_LV");

    realPE_SideLR1_Box = new G4Box("RealPolyEthylene_LR1_Box", 2. * sideLRScint_boxWidth_half,
                                   realPE_thickness1 / 2.,
                                   sideLRScint_boxHeight_half - realPE_LRCuttingLength1 / 2.);
    realPE_SideLR1_BoxLogical =
        new G4LogicalVolume(realPE_SideLR1_Box, _polyethylene, "RealPolyEthylene_LR1_LV");

    realPE_SideLR2_Box = new G4Box("RealPolyEthylene_LR2_Box", 2. * sideLRScint_boxWidth_half,
                                   realPE_thickness2 / 2.,
                                   sideLRScint_boxHeight_half - realPE_LRCuttingLength2 / 2.);
    realPE_SideLR2_BoxLogical =
        new G4LogicalVolume(realPE_SideLR2_Box, _polyethylene, "RealPolyEthylene_LR2_LV");

    realPE_Bottom_Box =
        new G4Box("BoratedPE_Bottom_Box", (pbbox_thick * 2. + pbbox_XYinnerlength) / 2.,
                  (pbbox_thick * 2. + pbbox_XYinnerlength) / 2., realPE_thickness1 / 2.);
    realPE_Bottom_BoxLogical =
        new G4LogicalVolume(realPE_Bottom_Box, _polyethylene, "RealPolyEthylene_Bottom_LV");

    realPE_Top_Box =
        new G4Box("RealPolyEthylene_Top_Box", pbtopbox_size, pbtopbox_size, realPE_thickness1 / 2.);
    realPE_Top_BoxLogical =
        new G4LogicalVolume(realPE_Top_Box, _polyethylene, "RealPolyEthylene_Top_LV");

    toyPE_BoxLogical->SetVisAttributes(pe_BoxVisAttr);
    realPE_SideFB_BoxLogical->SetVisAttributes(pe_BoxVisAttr);
    realPE_SideLR1_BoxLogical->SetVisAttributes(pe_BoxVisAttr);
    realPE_SideLR2_BoxLogical->SetVisAttributes(pe_BoxVisAttr);
    realPE_Bottom_BoxLogical->SetVisAttributes(pe_BoxVisAttr);
    realPE_Top_BoxLogical->SetVisAttributes(pe_BoxVisAttr);

    //////////////////////////////
    // Borated PE
    //////////////////////////////
    boratedPE_SideFB_Box = new G4Box("BoratedPE_SideFB_Box", boratedPE_thickness / 2.,
                                     2. * sideFBScint_boxWidth_half, sideFBScint_boxHeight_half);
    boratedPE_SideFB_BoxLogical =
        new G4LogicalVolume(boratedPE_SideFB_Box, _BoratedPE_5perCent, "BoratedPE_SideFB_LV");

    boratedPE_SideLR_Box = new G4Box("BoratedPE_SideLR_Box", 2. * sideLRScint_boxWidth_half,
                                     boratedPE_thickness / 2., sideLRScint_boxHeight_half);
    boratedPE_SideLR_BoxLogical =
        new G4LogicalVolume(boratedPE_SideLR_Box, _BoratedPE_5perCent, "BoratedPE_SideLR_LV");

    boratedPE_Bottom_Box =
        new G4Box("BoratedPE_Bottom_Box", (pbbox_thick * 2. + pbbox_XYinnerlength) / 2.,
                  (pbbox_thick * 2. + pbbox_XYinnerlength) / 2., boratedPE_thickness / 2.);
    boratedPE_Bottom_BoxLogical =
        new G4LogicalVolume(boratedPE_Bottom_Box, _BoratedPE_5perCent, "BoratedPE_Bottom_LV");

    boratedPE_BoxVisAttr = new G4VisAttributes(redl);
    boratedPE_BoxVisAttr->SetForceSolid(true);
    boratedPE_SideFB_BoxLogical->SetVisAttributes(boratedPE_BoxVisAttr);
    boratedPE_SideLR_BoxLogical->SetVisAttributes(boratedPE_BoxVisAttr);
    boratedPE_Bottom_BoxLogical->SetVisAttributes(boratedPE_BoxVisAttr);

    //////////////////////////////
    // Boric acid panel
    //////////////////////////////
    sideLRBoricAcid_HousingBox = new G4Box("BoricAcid_LR_MotherBox", width_LR,
                                           totalThickness_BoricAcidPanel / 2., height_LR);
    sideFBBoricAcid_HousingBox = new G4Box("BoricAcid_FB_MotherBox",
                                           totalThickness_BoricAcidPanel / 2., width_FB, height_FB);
    bottomBoricAcid_HousingBox = new G4Box("BoricAcid_Bottom_MotherBox", length_Bottom,
                                           width_Bottom, totalThickness_BoricAcidPanel / 2.);
    topBoricAcid_HousingBox    = new G4Box("BoricAcid_Top_MotherBox", length_Top, width_Top,
                                        totalThickness_BoricAcidPanel / 2.);
    sideLRBoricAcid_Block =
        new G4Box("BoricAcid_LR_Block", width_LR - thickness_BoricAcidSpacing,
                  thickness_BoricAcid / 2., height_LR - thickness_BoricAcidSpacing);
    sideFBBoricAcid_Block =
        new G4Box("BoricAcid_FB_Block", thickness_BoricAcid / 2.,
                  width_FB - thickness_BoricAcidSpacing, height_FB - thickness_BoricAcidSpacing);
    bottomBoricAcid_Block =
        new G4Box("BoricAcid_Bottom_Block", length_Bottom - thickness_BoricAcidSpacing,
                  width_Bottom - thickness_BoricAcidSpacing, thickness_BoricAcid / 2.);
    topBoricAcid_Block =
        new G4Box("BoricAcid_Top_Block", length_Top - thickness_BoricAcidSpacing,
                  width_Top - thickness_BoricAcidSpacing, thickness_BoricAcid / 2.);

    sideLRBoricAcid_MotherLogical =
        new G4LogicalVolume(sideLRBoricAcid_HousingBox, _polyethylene, "BoricAcid_LR_MotherLV");
    sideFBBoricAcid_MotherLogical =
        new G4LogicalVolume(sideFBBoricAcid_HousingBox, _polyethylene, "BoricAcid_FB_MotherLV");
    bottomBoricAcid_MotherLogical =
        new G4LogicalVolume(bottomBoricAcid_HousingBox, _polyethylene, "BoricAcid_Bottom_MotherLV");
    topBoricAcid_MotherLogical =
        new G4LogicalVolume(topBoricAcid_HousingBox, _polyethylene, "BoricAcid_Top_MotherLV");

    sideLRBoricAcid_BlockLogical =
        new G4LogicalVolume(sideLRBoricAcid_Block, _BoricAcidPowder, "BoricAcid_LR_BlockLV");
    sideFBBoricAcid_BlockLogical =
        new G4LogicalVolume(sideFBBoricAcid_Block, _BoricAcidPowder, "BoricAcid_FB_BlockLV");
    bottomBoricAcid_BlockLogical =
        new G4LogicalVolume(bottomBoricAcid_Block, _BoricAcidPowder, "BoricAcid_Bottom_BlockLV");
    topBoricAcid_BlockLogical =
        new G4LogicalVolume(topBoricAcid_Block, _BoricAcidPowder, "BoricAcid_Top_BlockLV");

    new G4PVPlacement(nullptr, G4ThreeVector(), sideLRBoricAcid_BlockLogical,
                      "BoricAcid_LR_BlockPV", sideLRBoricAcid_MotherLogical, false, 0);
    new G4PVPlacement(nullptr, G4ThreeVector(), sideFBBoricAcid_BlockLogical,
                      "BoricAcid_FB_BlockPV", sideFBBoricAcid_MotherLogical, false, 0);
    new G4PVPlacement(nullptr, G4ThreeVector(), bottomBoricAcid_BlockLogical,
                      "BoricAcid_Bottom_BlockPV", bottomBoricAcid_MotherLogical, false, 0);
    new G4PVPlacement(nullptr, G4ThreeVector(), topBoricAcid_BlockLogical, "BoricAcid_Top_BlockPV",
                      topBoricAcid_MotherLogical, false, 0);
    boricAcid_VisAttr        = new G4VisAttributes(G4Colour(1, 1, 1, 0.8));
    boricAcidHousing_VisAttr = new G4VisAttributes(G4Colour(1, 1, 0.9, 0.4));
    boricAcid_VisAttr->SetForceSolid(true);
    boricAcidHousing_VisAttr->SetForceSolid(true);

    bottomBoricAcid_MotherLogical->SetVisAttributes(boricAcidHousing_VisAttr);
    bottomBoricAcid_BlockLogical->SetVisAttributes(boricAcid_VisAttr);

    sideLRBoricAcid_MotherLogical->SetVisAttributes(boricAcidHousing_VisAttr);
    sideLRBoricAcid_BlockLogical->SetVisAttributes(boricAcid_VisAttr);

    sideFBBoricAcid_MotherLogical->SetVisAttributes(boricAcidHousing_VisAttr);
    sideFBBoricAcid_BlockLogical->SetVisAttributes(boricAcid_VisAttr);

    topBoricAcid_MotherLogical->SetVisAttributes(boricAcidHousing_VisAttr);
    topBoricAcid_BlockLogical->SetVisAttributes(boricAcid_VisAttr);

    //////////////////////////////
    // Place all the geometries at Work Area
    //////////////////////////////

    // Calculation of position parameters
    G4ThreeVector FB_Trans =
        G4ThreeVector(0, 0, -PbBoxOut->GetZHalfLength() + sideFBScint_MotherBox->GetYHalfLength());
    G4ThreeVector LR_Trans =
        G4ThreeVector(0, 0, -PbBoxOut->GetZHalfLength() + sideLRScint_MotherBox->GetYHalfLength());
    G4double sideLRScint_YPos =
        std::fmax(2 * sideFBScint_MotherBox->GetXHalfLength(), PbBoxOut->GetXHalfLength() * 2.);

    if (fEnable_OriginalGeom) {
        physTopPbBox1 =
            new G4PVPlacement(nullptr, G4ThreeVector(0., -pbtopbox_size / 2., CenterPbTZ),
                              logiTopPbBox, "physTopPbBox1", logiWorkArea, false, 0);
        new G4PVPlacement(nullptr, G4ThreeVector(0., pbtopbox_size / 2., CenterPbTZ), logiTopPbBox,
                          "physTopPbBox2", logiWorkArea, false, 0);
        new G4PVPlacement(nullptr, PosInnerDet, logiInnerDetector, "physInnerDetector",
                          logiWorkArea, false, 0);

        if (fEnable_InnerDetector) {
            G4ThreeVector PbBoxXHalfTrans1(-pbbox_size / 2., 0, 0);
            G4ThreeVector PbBoxXHalfTrans2(pbbox_size / 2., 0, 0);
            G4ThreeVector ShiftPbBoxInInnerDet(0, 0, 0.);
            G4RotationMatrix *PbBox2_alignMtx = new G4RotationMatrix();
            PbBox2_alignMtx->rotateZ(180 * deg);
            new G4PVPlacement(nullptr, PbBoxXHalfTrans1 + ShiftPbBoxInInnerDet, logiPbBox,
                              "physPbBox1", logiInnerDetector, false, 0);
            new G4PVPlacement(PbBox2_alignMtx, PbBoxXHalfTrans2 + ShiftPbBoxInInnerDet, logiPbBox,
                              "physPbBox2", logiInnerDetector, false, 0);
            PosSS = G4ThreeVector(0, 0, pbbox_zsize - ss_height - ssb_zsize);
            new G4PVPlacement(nullptr, PosSS, logiSSOVCOuter, "physSSOVC", logiInnerDetector, false,
                              0);
        }

        if (fPilot_Enable_Mumetal) {
            new G4PVPlacement(nullptr, Posmu, logimu, "physmu", logiInnerDetector, false, 0);
        }

        if (fPilot_Enable_TargetRoom) {
            new G4PVPlacement(nullptr, PosTR, logiTargetRoom, "physTargetRoom", logiCu1Inner, false,
                              0);
            G4cout << " Target Room position : " << PosTR << G4endl;
        }

        if (fEnable_Innermost) {
            new G4PVPlacement(nullptr, PosCu4Outer, logiCu4Outer, "physCu4Outer", logiSSOVCInner,
                              false, 0);
            new G4PVPlacement(nullptr, PosSMS, logiSuperMS, "physSuperMS", logiCu1Inner, false, 0);
            G4cout << "###   " << logiSuperMS->GetName() << "  " << logiSuperMS->GetMass() / kg
                   << G4endl;

            new G4PVPlacement(nullptr, PosCuMCP, logiCuMCP, "physCuMCP", logiCu1Inner, false, 0);
            new G4PVPlacement(nullptr, PosCuP1, logiCuP1, "physCuP1", logiCu1Inner, false, 0);
            new G4PVPlacement(nullptr, PosPbP2, logiPbP2, "physPbP2", logiCu1Inner, false, 0);
            new G4PVPlacement(nullptr, PosCuP3, logiCuP3, "physCuP3", logiCu1Inner, false, 0);
            new G4PVPlacement(nullptr, PosCuP4, logiCuP1, "physCuP4", logiCu1Inner, false, 0);
            new G4PVPlacement(nullptr, PosPbP5, logiPbP2, "physPbP5", logiCu1Inner, false, 0);
            new G4PVPlacement(nullptr, PosCuP6, logiCuP3, "physCuP6", logiCu1Inner, false, 0);

            G4LogicalVolume *G10MotherLV[5] = {logiCu1Inner, logiCu2Inner, logiCu3Inner,
                                               logiCu4Inner, logiSSOVCInner};

            const char G10PVName[5][30] = {"PlateGap4VolPhys_Cu1", "PlateGap4VolPhys_Cu2",
                                           "PlateGap4VolPhys_Cu3", "PlateGap4VolPhys_Cu4",
                                           "PlateGap4VolPhys_SSOVC"};
            for (int tmp_i = 0; tmp_i < 5; tmp_i++) {
                new G4PVPlacement(nullptr, G4ThreeVector(0, 0, CenterPlateGap4VolZ[tmp_i]),
                                  PlateGap4VolLogi[tmp_i], G10PVName[tmp_i], G10MotherLV[tmp_i],
                                  false, 0);
            }

            new G4PVPlacement(nullptr,
                              MASSPos,       // at position
                              logicMASS,     // its logical volume
                              "MASS-Spring", // its name
                              logiCu1Inner,  // its mother  volume
                              false,         // no boolean operation
                              0              // copy number
            );
            for (int i = 0; i < 4; i++) {
                new G4PVPlacement(nullptr,
                                  G4ThreeVector(300. / 2 * mm * cos(45 * deg + 90 * (i + 1) * deg),
                                                300. / 2 * mm * sin(45 * deg + 90 * (i + 1) * deg),
                                                MASSposZ + 24. / 2 * cm + 9.05 * mm),
                                  lv_spring, "Spring", logiCu1Inner, false, i);
            }
        }

        if (fEnable_Gantry) {
            static constexpr char beamX_PVname[][50] = {"HoriXBeamLyr1PV", "HoriXBeamLyr2PV",
                                                        "HoriXBeamLyr3PV"};
            static constexpr char beamY_PVname[][50] = {"HoriYBeamLyr1PV", "HoriYBeamLyr2PV",
                                                        "HoriYBeamLyr3PV"};
            for (int i = 0; i < 4; i++) {
                new G4PVPlacement(
                    VeBeamRotation,
                    G4ThreeVector(x_beam[i] * 750., y_beam[i] * 1600., posMountOrigin1)
                        .rotateZ(90. * deg),
                    logicHbeam1, "VerticalBeamPV", logiWorkArea, false, i);
            }

            for (int j = 0; j < 2; j++) {
                new G4PVPlacement(
                    nullptr,
                    G4ThreeVector(HoLineY[j] * 750., 0, posMountOrigin2).rotateZ(90. * deg),
                    logicHbeam2, "VerticalBottomBeamPV", logiWorkArea, false, j);
                for (int i = 0; i < 3; i++) {
                    new G4PVPlacement(HoBXRotation,
                                      G4ThreeVector(HoLine[j] * 750., 0.,
                                                    -Hbeam1_length - 20. + HoLayer[i] -
                                                        ((i == 0) ? 80. : 75.) + posMountOrigin1)
                                          .rotateZ(90. * deg),
                                      (i == 0) ? logicHo_beamXLar : logicHo_beamX, beamX_PVname[i],
                                      logiWorkArea, false, j);
                }
            }
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 3; i++) {
                    new G4PVPlacement(HoBYRotation,
                                      G4ThreeVector(0., HoLineY[j] * 1600.,
                                                    -Hbeam1_length - 20. + HoLayerY[i] -
                                                        ((i == 0) ? 80. : 75.) + posMountOrigin1)
                                          .rotateZ(90. * deg),
                                      (i == 0) ? logicHo_beamYLar : logicHo_beamY, beamY_PVname[i],
                                      logiWorkArea, false, j);
                }
            }
        }
        if (fEnable_Scintillator) {
            G4ThreeVector Top_Trans =
                G4ThreeVector(0, 0,
                              physTopPbBox1->GetTranslation().z() + TopPbBox->GetZHalfLength() +
                                  topScint_MotherBox->GetZHalfLength());
            G4ThreeVector Top_PMTTrans = G4ThreeVector(
                0, topScint_boxHeight_half + topScint_trapZ_half * 2. + PMTEnvelopeSizeZ / 2., 0);

            G4RotationMatrix *PMT_alignMtx_Top1 = new G4RotationMatrix();
            PMT_alignMtx_Top1->rotateX(-90 * deg);
            G4RotationMatrix *PMT_alignMtx_Top2 = new G4RotationMatrix(*PMT_alignMtx_Top1);
            PMT_alignMtx_Top2->rotateX(180 * deg);

            G4VPhysicalVolume *T1_Phys = new G4PVPlacement(
                nullptr, G4ThreeVector(-topScint_MotherBox->GetXHalfLength(), 0, 0) + Top_Trans,
                topScint_MotherLogical, "Envelope_MuonTopScintillator", logiWorkArea, false, 0,
                false);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T1_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 0);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T1_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 1);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T1_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 2);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T1_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 3);

            G4VPhysicalVolume *T2_Phys = new G4PVPlacement(
                0, G4ThreeVector(topScint_MotherBox->GetXHalfLength(), 0, 0) + Top_Trans,
                topScint_MotherLogical, "Envelope_MuonTopScintillator", logiWorkArea, false, 1,
                false);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T2_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 4);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T2_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 5);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T2_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 6);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T2_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 7);

            G4ThreeVector LR_PMTTrans = G4ThreeVector(
                sideLRScint_trapZ_half * 2 + sideLRScint_boxWidth_half + PMTEnvelopeSizeZ / 2., 0,
                0);
            G4ThreeVector LR_Mod_PMTTrans =
                G4ThreeVector(0, 0, sideLRScint_boxHeight_half - sideLRScint_trapX2_half * 3);
            G4RotationMatrix *PMT_alignMtx_LRSide1 = new G4RotationMatrix();
            PMT_alignMtx_LRSide1->rotateY(90 * deg);
            G4RotationMatrix *PMT_alignMtx_LRSide2 = new G4RotationMatrix();
            PMT_alignMtx_LRSide2->rotateY(-90 * deg);

            G4VPhysicalVolume *L1_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide2,
                PosPbBox +
                    G4ThreeVector(sideLRScint_MotherBox->GetXHalfLength(),
                                  -sideLRScint_YPos - sideLRScint_MotherBox->GetZHalfLength(), 0) +
                    LR_Trans,
                sideLRScint_MotherLogical, "Envelope_MuonSideLRScintillator", logiWorkArea, false,
                2 /*9*/, false);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              L1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 8);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              L1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 9);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              L1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 10);

            G4VPhysicalVolume *L2_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide1,
                PosPbBox +
                    G4ThreeVector(-sideLRScint_MotherBox->GetXHalfLength(),
                                  -sideLRScint_YPos - sideLRScint_MotherBox->GetZHalfLength(), 0) +
                    LR_Trans,
                sideLRScint_MotherLogical, "Envelope_MuonSideLRScintillator", logiWorkArea, false,
                3 /*8*/, false);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              L2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 11);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              L2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 12);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              L2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 13);

            G4VPhysicalVolume *R1_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide2,
                PosPbBox +
                    G4ThreeVector(+sideLRScint_MotherBox->GetXHalfLength(),
                                  sideLRScint_YPos + sideLRScint_MotherBox->GetZHalfLength(), 0) +
                    LR_Trans,
                sideLRScint_MotherLogical, "Envelope_MuonSideLRScintillator", logiWorkArea, false,
                4 /*7*/, false);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              R1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 14);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              R1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 15);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              R1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 16);

            G4VPhysicalVolume *R2_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide1,
                PosPbBox +
                    G4ThreeVector(-sideLRScint_MotherBox->GetXHalfLength(),
                                  sideLRScint_YPos + sideLRScint_MotherBox->GetZHalfLength(), 0) +
                    LR_Trans,
                sideLRScint_MotherLogical, "Envelope_MuonSideLRScintillator", logiWorkArea, false,
                5 /*6*/, false);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              R2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 17);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              R2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 18);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              R2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 19);

            G4ThreeVector FB_PMTTrans = G4ThreeVector(
                0, 0,
                sideFBScint_boxHeight_half + sideFBScint_trapZ_half * 2. + PMTEnvelopeSizeZ / 2.);
            G4RotationMatrix *PMT_alignMtx_FBSide = new G4RotationMatrix();
            PMT_alignMtx_FBSide->rotateX(180 * deg);

            G4VPhysicalVolume *B1_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(-pbbox_size - sideFBScint_MotherBox->GetZHalfLength(),
                                  -sideFBScint_MotherBox->GetXHalfLength(), 0) +
                    FB_Trans,
                sideFBScint_MotherLogical, "Envelope_MuonSideFBScintillator", logiWorkArea, false,
                6 /*4*/, false);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 20);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 21);

            G4VPhysicalVolume *B2_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(-pbbox_size - sideFBScint_MotherBox->GetZHalfLength(),
                                  sideFBScint_MotherBox->GetXHalfLength(), 0) +
                    FB_Trans,
                sideFBScint_MotherLogical, "Envelope_MuonSideFBScintillator", logiWorkArea, false,
                7 /*5*/, false);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 22);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 23);

            G4VPhysicalVolume *F1_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength(),
                                  -sideFBScint_MotherBox->GetXHalfLength(), 0) +
                    FB_Trans,
                sideFBScint_MotherLogical, "Envelope_MuonSideFBScintillator", logiWorkArea, false,
                8 /*2*/, false);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 24);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 25);

            G4VPhysicalVolume *F2_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength(),
                                  sideFBScint_MotherBox->GetXHalfLength(), 0) +
                    FB_Trans,
                sideFBScint_MotherLogical, "Envelope_MuonSideFBScintillator", logiWorkArea, false,
                9 /*3*/, false);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 26);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT", logiWorkArea, false, 27);
        }
        if (fEnable_NeutronShield) {
            new G4PVPlacement(
                nullptr,
                PosPbBox +
                    G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                      boratedPE_thickness / 2.,
                                  0, 0) +
                    FB_Trans,
                boratedPE_SideFB_BoxLogical, "physBoratedPE_SideFB", logiWorkArea, false, 0);
            new G4PVPlacement(
                nullptr,
                PosPbBox -
                    G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                      boratedPE_thickness / 2.,
                                  0, 0) +
                    FB_Trans,
                boratedPE_SideFB_BoxLogical, "physBoratedPE_SideFB", logiWorkArea, false, 1);

            new G4PVPlacement(
                nullptr,
                PosPbBox +
                    G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                      boratedPE_thickness + realPE_thickness1 / 2.,
                                  0, 0) +
                    FB_Trans,
                realPE_SideFB_BoxLogical, "physPolyEthylene_SideFB", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                PosPbBox -
                    G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                      boratedPE_thickness + realPE_thickness1 / 2.,
                                  0, 0) +
                    FB_Trans,
                realPE_SideFB_BoxLogical, "physPolyEthylene_SideFB", logiWorkArea, false, 1);

            new G4PVPlacement(
                nullptr,
                PosPbBox +
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos - sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                      boratedPE_thickness / 2.,
                                  0) +
                    LR_Trans,
                boratedPE_SideLR_BoxLogical, "physBoratedPE_SideLR", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                PosPbBox -
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos - sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                      boratedPE_thickness / 2.,
                                  0) +
                    LR_Trans,
                boratedPE_SideLR_BoxLogical, "physBoratedPE_SideLR", logiWorkArea, false, 1);

            new G4PVPlacement(
                nullptr,
                PosPbBox -
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos - sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                      boratedPE_thickness - realPE_thickness2 / 2.,
                                  0) +
                    G4ThreeVector(0, 0, -realPE_LRCuttingLength2 / 2.) + LR_Trans,
                realPE_SideLR2_BoxLogical, "physPolyEthylene_SideLR2", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                PosPbBox +
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos - sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                      boratedPE_thickness - realPE_thickness1 / 2.,
                                  0) +
                    G4ThreeVector(0, 0, -realPE_LRCuttingLength1 / 2.) + LR_Trans,
                realPE_SideLR1_BoxLogical, "physPolyEthylene_SideLR1", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                PosPbBox - G4ThreeVector(0, 0,
                                         (pbbox_Zinnerlength + pbbox_thick) / 2. +
                                             boratedPE_thickness / 2. + boratedPE_BottomGap),
                boratedPE_Bottom_BoxLogical, "physBoratedPE_Bottom", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                PosPbBox -
                    G4ThreeVector(0, 0,
                                  (pbbox_Zinnerlength + pbbox_thick) / 2. + boratedPE_thickness +
                                      boratedPE_BottomGap + realPE_thickness1 / 2.),
                realPE_Bottom_BoxLogical, "physPolyEthylene_Bottom", logiWorkArea, false, 0);

            new G4PVPlacement(nullptr,
                              PosPbBox - G4ThreeVector(0, 0,
                                                       (pbbox_Zinnerlength + pbbox_thick) / 2. +
                                                           boratedPE_thickness * 1.5 +
                                                           boratedPE_BottomGap + realPE_thickness1),
                              boratedPE_Bottom_BoxLogical, "physBoratedPE_Bottom", logiWorkArea,
                              false, 1);

            new G4PVPlacement(
                nullptr,
                PosPbBox + G4ThreeVector(0, 0, pbbox_zsize + sst_zsize * 2. + rubber_GapFromSSOVC),
                realPbUnderRubberB4C_BoxLogical, "B4CRubber_PbUnder", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                G4ThreeVector(0, 0,
                              physTopPbBox1->GetTranslation().z() + TopPbBox->GetZHalfLength() +
                                  topScint_MotherBox->GetZHalfLength() * 2. +
                                  rubber_thickness / 2.),
                realPbTopRubberB4C_BoxLogical, "B4CRubber_PbTop", logiWorkArea, false, 0);

            new G4PVPlacement(
                nullptr,
                G4ThreeVector(0, 0,
                              physTopPbBox1->GetTranslation().z() + TopPbBox->GetZHalfLength() +
                                  topScint_MotherBox->GetZHalfLength() * 2. + rubber_thickness +
                                  realPE_thickness1 / 2.),
                realPE_Top_BoxLogical, "physPolyEthylene_Top", logiWorkArea, false, 0);

            new G4PVPlacement(nullptr, sideLRBoricAcid_Pos[0], sideLRBoricAcid_MotherLogical,
                              "BoricAcid_LR_MotherPV", logiInnerDetector, false, 0);
            new G4PVPlacement(nullptr, sideLRBoricAcid_Pos[1], sideLRBoricAcid_MotherLogical,
                              "BoricAcid_LR_MotherPV", logiInnerDetector, false, 1);

            new G4PVPlacement(nullptr, sideFBBoricAcid_Pos[0], sideFBBoricAcid_MotherLogical,
                              "BoricAcid_FB_MotherPV", logiInnerDetector, false, 0);
            new G4PVPlacement(nullptr, sideFBBoricAcid_Pos[1], sideFBBoricAcid_MotherLogical,
                              "BoricAcid_FB_MotherPV", logiInnerDetector, false, 1);

            new G4PVPlacement(nullptr, bottomBoricAcid_Pos, bottomBoricAcid_MotherLogical,
                              "BoricAcid_Bottom_MotherPV", logiInnerDetector, false, 0);
        }
    } else {
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), realPE_SideFB_BoxLogical,
                          "physPolyEthylene_SideFB", logiWorkArea, false, 0);
    }
    if (fNeutronMode) {
        // Neutron simulation section
        RockSolid->SetDeltaThetaAngle(180. * deg);
        logiRockCavity->SetSolid(
            new G4Sphere("RockCavitySolid", 0, NEUT_CavityRadius, 0, 360 * deg, 0, 180 * deg));
        logiRockCavityVis->SetVisibility(false);

        switch (whichNShieldingConf) {
            default:
                G4cerr << "Wrong neutron shielding configuration! I will fuck up this simulation!"
                       << G4endl;
                exit(-1);
                return;
                break;
            case kNS_RealConf:
                G4cout << "Real configuration has been selected." << G4endl;
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  logiRockCavity, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox +
                        G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                          boratedPE_thickness / 2.,
                                      0, 0) +
                        FB_Trans,
                    boratedPE_SideFB_BoxLogical, "physBoratedPE_SideFB", logiWorkArea, false, 0);
                new G4PVPlacement(
                    nullptr,
                    PosPbBox -
                        G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                          boratedPE_thickness / 2.,
                                      0, 0) +
                        FB_Trans,
                    boratedPE_SideFB_BoxLogical, "physBoratedPE_SideFB", logiWorkArea, false, 1);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox +
                        G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                          boratedPE_thickness + realPE_thickness1 / 2.,
                                      0, 0) +
                        FB_Trans,
                    realPE_SideFB_BoxLogical, "physPolyEthylene_SideFB", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox -
                        G4ThreeVector(pbbox_size + sideFBScint_MotherBox->GetZHalfLength() * 2. +
                                          boratedPE_thickness + realPE_thickness1 / 2.,
                                      0, 0) +
                        FB_Trans,
                    realPE_SideFB_BoxLogical, "physPolyEthylene_SideFB", logiWorkArea, false, 1);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox +
                        G4ThreeVector(0.,
                                      -sideLRScint_YPos -
                                          sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                          boratedPE_thickness / 2.,
                                      0) +
                        LR_Trans,
                    boratedPE_SideLR_BoxLogical, "physBoratedPE_SideLR", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox -
                        G4ThreeVector(0.,
                                      -sideLRScint_YPos -
                                          sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                          boratedPE_thickness / 2.,
                                      0) +
                        LR_Trans,
                    boratedPE_SideLR_BoxLogical, "physBoratedPE_SideLR", logiWorkArea, false, 1);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox -
                        G4ThreeVector(0.,
                                      -sideLRScint_YPos -
                                          sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                          boratedPE_thickness - realPE_thickness2 / 2.,
                                      0) +
                        G4ThreeVector(0, 0, -realPE_LRCuttingLength2 / 2.) + LR_Trans,
                    realPE_SideLR2_BoxLogical, "physPolyEthylene_SideLR2", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox +
                        G4ThreeVector(0.,
                                      -sideLRScint_YPos -
                                          sideLRScint_MotherBox->GetZHalfLength() * 2. -
                                          boratedPE_thickness - realPE_thickness1 / 2.,
                                      0) +
                        G4ThreeVector(0, 0, -realPE_LRCuttingLength1 / 2.) + LR_Trans,
                    realPE_SideLR1_BoxLogical, "physPolyEthylene_SideLR1", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox - G4ThreeVector(0, 0,
                                             (pbbox_Zinnerlength + pbbox_thick) / 2. +
                                                 boratedPE_thickness / 2. + boratedPE_BottomGap),
                    boratedPE_Bottom_BoxLogical, "physBoratedPE_Bottom", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox - G4ThreeVector(0, 0,
                                             (pbbox_Zinnerlength + pbbox_thick) / 2. +
                                                 boratedPE_thickness + boratedPE_BottomGap +
                                                 realPE_thickness1 / 2.),
                    realPE_Bottom_BoxLogical, "physPolyEthylene_Bottom", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox - G4ThreeVector(0, 0,
                                             (pbbox_Zinnerlength + pbbox_thick) / 2. +
                                                 boratedPE_thickness * 1.5 + boratedPE_BottomGap +
                                                 realPE_thickness1),
                    boratedPE_Bottom_BoxLogical, "physBoratedPE_Bottom", logiWorkArea, false, 1);

                new G4PVPlacement(
                    nullptr,
                    PosPbBox +
                        G4ThreeVector(0, 0, pbbox_zsize + sst_zsize * 2. + rubber_GapFromSSOVC),
                    realPbUnderRubberB4C_BoxLogical, "B4CRubber_PbUnder", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    G4ThreeVector(0, 0,
                                  physTopPbBox1->GetTranslation().z() + TopPbBox->GetZHalfLength() +
                                      topScint_MotherBox->GetZHalfLength() * 2. +
                                      rubber_thickness / 2.),
                    realPbTopRubberB4C_BoxLogical, "B4CRubber_PbTop", logiWorkArea, false, 0);

                new G4PVPlacement(
                    nullptr,
                    G4ThreeVector(0, 0,
                                  physTopPbBox1->GetTranslation().z() + TopPbBox->GetZHalfLength() +
                                      topScint_MotherBox->GetZHalfLength() * 2. + rubber_thickness +
                                      realPE_thickness1 / 2.),
                    realPE_Top_BoxLogical, "physPolyEthylene_Top", logiWorkArea, false, 0);

                new G4PVPlacement(nullptr, sideLRBoricAcid_Pos[0], sideLRBoricAcid_MotherLogical,
                                  "BoricAcid_LR_MotherPV", logiInnerDetector, false, 0);
                new G4PVPlacement(nullptr, sideLRBoricAcid_Pos[1], sideLRBoricAcid_MotherLogical,
                                  "BoricAcid_LR_MotherPV", logiInnerDetector, false, 1);

                new G4PVPlacement(nullptr, sideFBBoricAcid_Pos[0], sideFBBoricAcid_MotherLogical,
                                  "BoricAcid_FB_MotherPV", logiInnerDetector, false, 0);
                new G4PVPlacement(nullptr, sideFBBoricAcid_Pos[1], sideFBBoricAcid_MotherLogical,
                                  "BoricAcid_FB_MotherPV", logiInnerDetector, false, 1);

                new G4PVPlacement(nullptr, bottomBoricAcid_Pos, bottomBoricAcid_MotherLogical,
                                  "BoricAcid_Bottom_MotherPV", logiInnerDetector, false, 0);
                break;
            case kNS_Naked:
                G4cout << "Naked shielding configuration has been selected." << G4endl;
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  logiRockCavity, false, 0);
                break;
            case kNS_B4C:
                G4cout << "B4C rubber shielding configuration has been selected." << G4endl;
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  toyRubberB4C_BoxLogical, false, 0);

                new G4PVPlacement(nullptr, G4ThreeVector(), toyRubberB4C_BoxLogical,
                                  "physB4CRubber", logiRockCavity, false, 0);
                break;
            case kNS_B4CnPE:
                G4cout << "B4C rubber + PE 20cm shielding configuration has been selected."
                       << G4endl;
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  toyRubberB4C_BoxLogical, false, 0);

                new G4PVPlacement(nullptr, G4ThreeVector(), toyRubberB4C_BoxLogical,
                                  "physB4CRubber", toyPE_BoxLogical, false, 0);

                new G4PVPlacement(nullptr, G4ThreeVector(), toyPE_BoxLogical, "physPolyEthylene",
                                  logiRockCavity, false, 0);
                break;
            case kNS_PE20:
                G4cout << "PE 20cm shielding configuration has been selected." << G4endl;
                toyPE_thickness = 20 * cm;
                toyPE_Box->SetXHalfLength(workboxX + toyPE_thickness);
                toyPE_Box->SetYHalfLength(workboxY + toyPE_thickness);
                toyPE_Box->SetZHalfLength(workboxH + toyPE_thickness);
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  toyPE_BoxLogical, false, 0);

                new G4PVPlacement(nullptr, G4ThreeVector(), toyPE_BoxLogical, "physPolyEthylene",
                                  logiRockCavity, false, 0);
                break;
            case kNS_PE10:
                G4cout << "PE 10cm shielding configuration has been selected." << G4endl;
                toyPE_thickness = 10 * cm;
                toyPE_Box->SetXHalfLength(workboxX + toyPE_thickness);
                toyPE_Box->SetYHalfLength(workboxY + toyPE_thickness);
                toyPE_Box->SetZHalfLength(workboxH + toyPE_thickness);
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  toyPE_BoxLogical, false, 0);

                new G4PVPlacement(nullptr, G4ThreeVector(), toyPE_BoxLogical, "physPolyEthylene",
                                  logiRockCavity, false, 0);
                break;
            case kNS_B4CnBAcid:
                G4cout << "B4C rubber + BoricAcid shielding configuration has been selected."
                       << G4endl;
                new G4PVPlacement(nullptr, G4ThreeVector(), logiWorkArea, "physWorkArea",
                                  toyRubberB4C_BoxLogical, false, 0);

                new G4PVPlacement(nullptr, sideLRBoricAcid_Pos[0], sideLRBoricAcid_MotherLogical,
                                  "BoricAcid_LR_MotherPV", logiInnerDetector, false, 0);
                new G4PVPlacement(nullptr, sideLRBoricAcid_Pos[1], sideLRBoricAcid_MotherLogical,
                                  "BoricAcid_LR_MotherPV", logiInnerDetector, false, 1);

                new G4PVPlacement(nullptr, sideFBBoricAcid_Pos[0], sideFBBoricAcid_MotherLogical,
                                  "BoricAcid_FB_MotherPV", logiInnerDetector, false, 0);
                new G4PVPlacement(nullptr, sideFBBoricAcid_Pos[1], sideFBBoricAcid_MotherLogical,
                                  "BoricAcid_FB_MotherPV", logiInnerDetector, false, 1);

                new G4PVPlacement(nullptr, bottomBoricAcid_Pos, bottomBoricAcid_MotherLogical,
                                  "BoricAcid_Bottom_MotherPV", logiInnerDetector, false, 0);

                new G4PVPlacement(nullptr, PosInnerDet + topBoricAcid_Pos,
                                  topBoricAcid_MotherLogical, "BoricAcid_Top_MotherPV",
                                  logiWorkArea, false, 0);

                new G4PVPlacement(nullptr, G4ThreeVector(), toyRubberB4C_BoxLogical,
                                  "physB4CRubber", logiRockCavity, false, 0);
                break;
        }

        new G4PVPlacement(nullptr, G4ThreeVector(), logiRockCavity, "physRockCavity", logiHall,
                          false, 0);
    } else {
        new G4PVPlacement(tube_alignMtx, G4ThreeVector(), logiRockCavity, "physRockCavity",
                          logiRock, false, 0);
        new G4PVPlacement(nullptr,                // rotation
                          G4ThreeVector(0, 0, 0), // translation
                          logiRock,               // associated logical vol
                          "physRock",             // name
                          logiHall,               // parent
                          false,                  // no "Many"
                          0);                     // copy number
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -RockDiskThick_half), "physRockDisk",
                          logiRockDisk, physHall, false, 0, false);
        new G4PVPlacement(WA_alignMtx, G4ThreeVector(0, workboxH, 0), logiWorkArea, "physWorkArea",
                          logiRockCavity, false, 0);
    }

    //////////////////////////////
    // Set Attributes
    //////////////////////////////
    G4cout << "Set Geometry Attributes...\n";
    logiHallVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.1));
    logiHall->SetVisAttributes(logiHallVis);
    logiHall->SetVisAttributes(G4VisAttributes::GetInvisible);

    // llization
    G4VisAttributes *va_smallbox = new G4VisAttributes(brown);
    va_smallbox->SetForceSolid(true);
    lv_smallbox->SetVisAttributes(va_smallbox);

    G4VisAttributes *va_smallboxB = new G4VisAttributes(brown);
    va_smallboxB->SetForceSolid(true);
    lv_smallboxB->SetVisAttributes(va_smallboxB);

    G4VisAttributes *va_TopPlate = new G4VisAttributes(brown);
    va_TopPlate->SetForceSolid(true);
    lv_TopPlate->SetVisAttributes(va_TopPlate);

    G4VisAttributes *va_TopPlateB = new G4VisAttributes(brown);
    va_TopPlate->SetForceSolid(true);
    lv_TopPlateB->SetVisAttributes(va_TopPlateB);

    G4VisAttributes *va_copillar = new G4VisAttributes(brown);
    va_copillar->SetForceSolid(true);
    lv_copillar->SetVisAttributes(va_copillar);

    //	  IVC_Welding//////////

    G4VisAttributes *va_IVCwe = new G4VisAttributes(red);
    va_IVCwe->SetForceSolid(true);
    lv_IVCwe->SetVisAttributes(va_IVCwe);

    G4VisAttributes *va_OVCwe = new G4VisAttributes(red);
    va_OVCwe->SetForceSolid(true);
    lv_OVCwe->SetVisAttributes(va_OVCwe);

    logiTargetRoomVis = new G4VisAttributes(G4Colour(1, 1, 1, 0.8));
    logiTargetRoom->SetVisAttributes(logiTargetRoomVis);

    G4VisAttributes *logicMASSVis = new G4VisAttributes(orangel);
    logicMASSVis->SetVisibility(true);
    logicMASSVis->SetForceWireframe(true);
    logicMASS->SetVisAttributes(logicMASSVis);

    G4VisAttributes *lv_springVis = new G4VisAttributes(greenl);
    lv_springVis->SetVisibility(true);
    lv_springVis->SetForceSolid(true);
    lv_spring->SetVisAttributes(lv_springVis);

    logiRockVis = new G4VisAttributes(G4Colour(1, 0, 0, 0.5));
    logiRockVis->SetForceSolid(false);
    logiRockDiskVis = new G4VisAttributes(red);
    logiWorkAreaVis = new G4VisAttributes(white);
    logiWorkAreaVis->SetVisibility(false);
    logiInnerDetectorVis = new G4VisAttributes(G4Colour(1, 0, 0, 0.5));
    logiPbBoxVis         = new G4VisAttributes(grey);
    logiPbBoxVis         = new G4VisAttributes(greyl);
    logiSSOVCVis         = new G4VisAttributes(cyanl);
    logimuVis            = new G4VisAttributes(yellowl);
    logiCu4Vis           = new G4VisAttributes(bluel);
    logiCu3Vis           = new G4VisAttributes(green);
    logiCu2Vis           = new G4VisAttributes(cyan);
    logiCu1Vis           = new G4VisAttributes(brownl);
    logiSuperMSVis       = new G4VisAttributes(G4Colour(0.4, 0.2, 0.6, 0.8));
    logiSuperMS->SetVisAttributes(logiSuperMSVis);
    logiCuMCPVis = new G4VisAttributes(lblue);
    logiCuMCPVis->SetForceSolid(true);
    logiCuP1Vis       = new G4VisAttributes(lgreen);
    logiPbP2Vis       = new G4VisAttributes(grey);
    logiCuP3Vis       = new G4VisAttributes(red);
    lv_PhotonFrameVis = new G4VisAttributes(greenl);

    for (int nn = 0; nn < 6; nn++) {
        logiCMOVis[nn] = new G4VisAttributes(yellow);
    }
    logiSquareDiskUpVis = new G4VisAttributes(brown);

    logiPhotonDetVis = new G4VisAttributes(G4Colour(1, 1, 1, 0.3));
    logiGeWaferVis   = new G4VisAttributes(blue);
    logiVacDiskVis   = new G4VisAttributes(orange);

    if (flagInvisible == 999) {
        logiPhotonDetVis->SetVisibility(true);
        logiPhotonDetVis->SetForceSolid(true);
        logiGeWaferVis->SetVisibility(true);
        logiGeWaferVis->SetForceSolid(true);
        logiVacDiskVis->SetVisibility(true);
        logiPhotonDet->SetVisAttributes(logiPhotonDetVis);
        logiGeWafer->SetVisAttributes(logiGeWaferVis);
        logiVacDisk->SetVisAttributes(logiVacDiskVis);
    } else {
        logiVacDisk->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiPhotonDetVis->SetVisibility(true);
        logiPhotonDetVis->SetForceWireframe(true);
        logiPhotonDet->SetVisAttributes(logiPhotonDetVis);
    }
    if (flagInvisible) {
        logiRock->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiWorkArea->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiTopPbBox->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiPbBox->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiSSOVCOuter->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu4Outer->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu3Outer->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu2Outer->SetVisAttributes(G4VisAttributes::GetInvisible);
        logiCu1Outer->SetVisAttributes(G4VisAttributes::GetInvisible);

        if (flagInvisible == 999) {
            for (int nn = 0; nn < 5; nn++) {
                logiCMOCell[nn]->SetVisAttributes(G4VisAttributes::GetInvisible);
            }

            logiCuMCP->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiCuP1->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiPbP2->SetVisAttributes(G4VisAttributes::GetInvisible);
            logiCuP3->SetVisAttributes(G4VisAttributes::GetInvisible);
        }
    } else {
        logiCuP1Vis->SetForceSolid(true);
        logiPbP2Vis->SetForceSolid(true);
        logiCuP3Vis->SetForceSolid(true);
        for (int nn = 0; nn < 6; nn++) {
            logiCMOVis[nn]->SetVisibility(true);
            logiCMOVis[nn]->SetForceSolid(true);
            logiReflectorVis[nn]->SetVisibility(true);
            logiReflectorVis[nn]->SetForceSolid(true);
        }
        logiGeWaferVis->SetVisibility(true);
        logiGeWaferVis->SetForceSolid(true);

        logiRock->SetVisAttributes(logiRockVis);
        logiRockCavity->SetVisAttributes(logiRockCavityVis);
        logiRockDisk->SetVisAttributes(logiRockDiskVis);
        logiWorkArea->SetVisAttributes(logiWorkAreaVis);

        logiInnerDetector->SetVisAttributes(logiInnerDetectorVis);
        logiTopPbBox->SetVisAttributes(logiPbBoxVis);
        logiPbBox->SetVisAttributes(logiPbBoxVis);
        logiSSOVCVis->SetForceWireframe(true);
        logiSSOVCOuter->SetVisAttributes(logiSSOVCVis);

        logimuVis->SetForceWireframe(true);
        logimu->SetVisAttributes(logimuVis);

        logiCu4Vis->SetForceWireframe(true);
        logiCu4Outer->SetVisAttributes(logiCu4Vis);
        logiCu3Vis->SetForceWireframe(true);
        logiCu4Outer->SetVisAttributes(logiCu4Vis);
        logiCu3Vis->SetForceWireframe(true);
        logiCu3Outer->SetVisAttributes(logiCu3Vis);
        logiCu2Outer->SetVisAttributes(logiCu2Vis);
        logiCu1Outer->SetVisAttributes(logiCu1Vis);
        logiCuMCP->SetVisAttributes(logiCuMCPVis);
        logiCuP1->SetVisAttributes(logiCuP1Vis);
        logiPbP2->SetVisAttributes(logiPbP2Vis);
        logiCuP3->SetVisAttributes(logiCuP3Vis);
        lv_PhotonFrameVis->SetForceSolid(true);
        lv_PhotonFrame->SetVisAttributes(lv_PhotonFrameVis);

        for (int nn = 0; nn < 6; nn++) {
            logiCMOCell[nn]->SetVisAttributes(logiCMOVis[nn]);
            logiReflector[nn]->SetVisAttributes(logiReflectorVis[nn]);

            logiSquareDiskUpVis->SetForceSolid(true);
            logiSquareDiskUp[nn]->SetVisAttributes(logiSquareDiskUpVis);
        }
        logiGeWafer->SetVisAttributes(logiGeWaferVis);
    }

#if G4VERSION_NUMBER < 1000
    ConstructAMoREPilot_SDandField();
#endif

    /////////////////////////////
    // photocathode surfaces...
    /////////////////////////////

    new G4LogicalBorderSurface("MLCS_photocathode_logsurf",
                               physiPMT, // exiting glass into vac.
                               physiPMTVacu, Photocathode_opsurf);

    G4Region *PmtRegion = new G4Region("MLCS");
    PmtRegion->AddRootLogicalVolume(logicPMT);
    new CupPMTOpticalModel("MLCS_optical_model", physiPMT);

    /////////////////////////////////////////////////////////////////
    //  Set Region
    //

    G4Region *crystalsRegion = new G4Region("crystals");
    for (int nn = 0; nn < 6; nn++)
        crystalsRegion->AddRootLogicalVolume(logiCMOCell[nn]);
}

void AmoreDetectorConstruction::ConstructAMoREPilot_SDandField() {
    //////////////////////////////
    // --- TG sensitive detector
    //////////////////////////////

    CupScintSD *TGSD;
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    TGSD = new CupScintSD(SDname = "/CupDet/TGSD", 30);
    SDman->AddNewDetector(TGSD);
    for (int nn = 0; nn < 6; nn++)
        fPilot_logiCMOCell[nn]->SetSensitiveDetector(TGSD);
    fPilot_logiGOLDa->SetSensitiveDetector(TGSD);
    fPilot_logiGOLDb->SetSensitiveDetector(TGSD);

    CupVetoSD *MuonScintSD;
    MuonScintSD = new CupVetoSD("/CupDet/MuonVetoSD",
                                2 + 4 + 4); // Top 2, SideX 4, SideY 4
    SDman->AddNewDetector(MuonScintSD);
    fPilot_TopScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
    fPilot_TopScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
    fPilot_TopScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);

    fPilot_SideFBScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
    fPilot_SideFBScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
    fPilot_SideFBScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);

    fPilot_SideLRScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
    fPilot_SideLRScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
    fPilot_SideLRScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);

    CupPMTSD *pmtSD;
    pmtSD = new CupPMTSD("/cupdet/pmt/MLCS", 34); // number of pmts= 28
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD);
    fPilot_logicPMT->SetSensitiveDetector(pmtSD);
    fPilot_logicPMTVacu->SetSensitiveDetector(pmtSD);
}
