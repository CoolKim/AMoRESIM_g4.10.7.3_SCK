/**
 * @file AmoreDetectorStaticInfo.hh
 * @author BaseHardware
 * @date 190315
 * @brief Header files for static variable definitions
 * @details This file contains static variables for geometric dimenstions. If you want to define
 * dimensional variables which can be changed in run time, please set them in the data/settings*
 * files and use the variables with CupParam.
 * @see CupParam
 * */
#ifndef AmoreDetectorStaticInfo_HH
#define AmoreDetectorStaticInfo_HH 1

#include "G4Color.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "AmoreSim/AmoreDetectorStaticType.hh"

/// Static information for geometric dimension for AMoRE simulation
namespace AmoreDetectorStaticInfo {
	/// Tolerence for boolean solids, normal size
	constexpr G4double solidBooleanTol = 100 * um;
	/// Tolerence for boolean solids, small size
	constexpr G4double solidBooleanTolSmall = 5 * um;
	/// Tolerence for boolean solids, big size
	constexpr G4double solidBooleanTolBig = 1 * mm;
	/// Tolerence angle for boolean solids
	constexpr G4double solidBooleanTolAngle = 1. * deg;
	/// A conversion constant for inch to cm
	constexpr G4double AMoRE_unit_inch = 2.54 * cm;
	/// Not determined, for AMoRE-I crystal dimensions
	constexpr G4double toBeDetermined = 4.5 * cm;

	/// Color table for AMoRE simulation framework visualization attributes
	namespace ColorTable {
		const G4Colour white   = G4Colour::White();
		const G4Colour grey    = G4Colour::Gray();
		const G4Colour red     = G4Colour::Red();
		const G4Colour blue    = G4Colour::Blue();
		const G4Colour cyan    = G4Colour::Cyan();
		const G4Colour magenta = G4Colour::Magenta();
		const G4Colour yellow  = G4Colour::Yellow();
		const G4Colour green   = G4Colour::Green();
		const G4Colour brown   = G4Colour::Brown();
        const G4Colour orange(.75, .55, 0.0);
        const G4Colour lgrey(.85, .85, .85);
        const G4Colour lblue(0.0, 0.0, .75);
        const G4Colour lgreen(0.0, .75, 0.0);
        const G4Colour greyl(0.5, 0.5, 0.5, 0.3);
        const G4Colour cyanl(0.0, 1.0, 1.0, 0.5);
        const G4Colour orangel(.75, .55, 0.0, 0.5);
        const G4Colour bluel(0.0, 0.0, 1.0, 0.5);
        const G4Colour greenl(0.0, 1.0, 0.0, 0.3);
        const G4Colour brownl(0.7, 0.4, 0.1, 0.5);
        const G4Colour redl(1.0, .0, 0.3);
        const G4Colour whitel(1, 1, 1, 0.8);
        const G4Colour lgreenl(G4Colour(129 / 255., 193 / 255., 71 / 255., 0.9));
	} // namespace ColorTable

	/// Geometric dimension of AMoRE-Pilot run 7 environment. Defined for mimicking AMoRE-Pilot's
	/// neutron shielding in AMoRE-I geometry
	namespace AMoRE_Pilot_Run7Save {
		/// Thickness of PE for neutron shielding, thin part
		constexpr G4double polyEthylene_thickness1 = 10. * cm;
		/// Thickness of PE for neutron shielding, thick part
		constexpr G4double polyEthylene_thickness2 = 30. * cm;
		/// Cutting length 1 for LR part for PE to avoid overlapping with a gantry geometry
		constexpr G4double polyEthylene_LRCuttingLength1 = 60. * cm;
		/// Cutting length 2 for LR part for PE to avoid overlapping with a gantry geometry
		constexpr G4double polyEthylene_LRCuttingLength2 = 13. * cm;
		/// Thickness of borated PE for neutron shielding, all parts share the same value
		constexpr G4double boratedPE_thickness = 2.5 * cm;
		/// Gap between bottom outer surface of lead shield and borated PE shield at bottom
		constexpr G4double boratedPE_BottomGap = 19. * cm;

		/// Thickness of boric acid panel
		constexpr G4double thickness_BoricAcid = 3. * mm;
		/// Thickness of PE housing for boric acid
		constexpr G4double thickness_BoricAcidHousing = 1. * mm;
		/// Width of space at edge part of boric acid
		constexpr G4double thickness_BoricAcidSpacing = 10. * mm;

		/// Thickness of boron carbide rubber
		constexpr G4double boronCarbideRubber_thickness = 4. * mm;
		/// X and Y dimension of boron carbide rubber
		constexpr G4double boronCarbideRubber_UnderPb_XYsize = 100. * cm;
		/// Gap between top surface of SSOVC and B4C rubber panel
		constexpr G4double boronCarbideRubber_GapFromSSOVC = 10. * cm;
	} // namespace AMoRE_Pilot_Run7Save

	namespace AMoRE_I {
		struct ModuleSDInfo;

		/// Total number of tower in crystal array
		constexpr G4int totalNumOfTower = 4;
		/// Total number of modules in a single tower
		constexpr G4int maxModuleNumInTower = 6;
		/// Information list for crystals
		extern const CrystalModuleInfo crystalModuleInfoList[totalNumOfTower][maxModuleNumInTower];

		/// Total world size
		constexpr G4double bounding_size = 30. * m / 2.0;

		/// Radius of rock shell, it can be hemisphere or orb
		constexpr G4double RockRadius = 8.* m;
		/// Radius of cavity in the rock shell
		constexpr G4double RockCavityRadius = 5.0 * m;
		/// Radius of cavity in the rock shell in the neutron mode
		constexpr G4double NEUT_CavityRadius = 3.5 * m + 5 * mm;

		/// X size of work area in the cavern
		constexpr G4double workboxX = 3.7 * m / 2;
		/// Y size of work area in the cavern
		constexpr G4double workboxY = 3.6 * m / 2;
		/// Z size of work area in the cavern
		constexpr G4double workboxH = 4.5 * m / 2;

		/// Thickness of rock disk at the floor of laboratory
		constexpr G4double RockDiskThick_half = 1.0 * m;

		// OVC 
		constexpr G4double ss_radius      = 63.5 * cm / 2;
		constexpr G4double ss_thick       = 0.5 * cm;
		constexpr G4double ss_height_half = 153.8 * cm / 2; 
		constexpr G4double sstopbox_xsize = 74.8 * cm / 2.0;
		constexpr G4double sstopbox_ysize = 74.8 * cm / 2.0;
		constexpr G4double sst_zsize_half = 3. * cm / 2.;
		constexpr G4double ssb_zsize      = 2. * cm / 2;
		constexpr G4double plateGap       = 12. * cm; 

		// Layer4_ 50K-SHIELD Copper
		constexpr G4double cu4_radius       = 57.6 * cm / 2;
		constexpr G4double cu4_thick        = 0.3 * cm;
		constexpr G4double cu4_inner_height = 134.2 * cm / 2; // 136.9-1.0-1.7
		constexpr G4double cu4t_zsize       = 1.7 * cm / 2;
		constexpr G4double cu4b_zsize       = 1.0 * cm / 2;
		constexpr G4double cu4_gap_fromTop  = 120 * mm;

		// Layer3_Cu IVC
		constexpr G4double cu3_radius       = 51.4 * cm / 2;
		constexpr G4double cu3_thick        = 0.5 * cm;
		constexpr G4double cu3_inner_height = 116.18 * cm / 2; // 119.88 -1.8 - 1.9
		constexpr G4double cu3t_zsize       = 1.8 * cm / 2;
		constexpr G4double cu3b_zsize       = 1.9 * cm / 2;
		constexpr G4double cu3_gap_fromTop  = 120 * mm;

		// Layer2_Cu SHIELD-STILL
		constexpr G4double cu2_radius       = 45.5 * cm / 2;
		constexpr G4double cu2_thick        = 0.1 * cm;
		constexpr G4double cu2_inner_height = 100.2 * cm / 2;
		constexpr G4double cu2t_zsize       = 2.0 * cm / 2;
		constexpr G4double cu2b_zsize       = 0.1 * cm / 2;
		constexpr G4double cu2_gap_fromTop  = 120 * mm;

		// Layer1_Cu 50mK-SHIELD
		constexpr G4double cu1_radius       = 43.5 * cm / 2;
		constexpr G4double cu1_thick        = 0.1 * cm;
		constexpr G4double cu1_inner_height = 85.1 * cm / 2;
		constexpr G4double cu1t_zsize       = 2.0 * cm / 2;
		constexpr G4double cu1b_zsize       = 0.1 * cm / 2;
		constexpr G4double cu1_gap_fromTop  = 120 * mm;

		//  Pb Top Box
		constexpr G4double pbtopbox_size              = 150. * cm / 2.0;
		constexpr G4double pbtopbox_zsize             = 15. * cm / 2.0;
		constexpr G4double pbtopbox_housing_thickness = 3 * mm;
		constexpr G4double lifting_PbTopBox           = 22.0 * cm;

		//  Pb  Box
		constexpr G4double leadBoxThickness        = 20. * cm;
		constexpr G4double leadBoxHousingThickness = 0 * mm;
		constexpr G4double leadBoxXYinnerlength    = 80. * cm;
		constexpr G4double leadBoxZinnerlength     = 162.1 * cm;
		constexpr G4double verticalGapToMuVeto     = 885 * mm;

		// Cu Plate: Mixing Chamber SHIELD
		constexpr G4double copperMixingChamberPlateRadius   = 40.8 * cm / 2;
		constexpr G4double copperMixingChamberPlatThickness = 3. * cm / 2;

		// Plate1_Cu
		constexpr G4double copperPlate1ZSize_half = 0.5 * cm / 2;
		constexpr G4double plateGap2              = 1.0 * cm;

		// Plate2_Pb Shield
		constexpr G4double leadPlate2ZSize_half = 5. * cm / 2;

		// Plate3_Cu
		constexpr G4double cup3_zsize = 1. * cm / 2;

		// Polyethylene(+Borated) geometry
		constexpr G4double polyEthylene_commonThickness = 30. * cm;
		constexpr G4double polyEthyleneLR_cuttingLength = 13. * cm;
		constexpr G4double polyEthylene_topThickness    = 10. * cm;
		constexpr G4double boratedPE_thickness          = 3. * cm;
		constexpr G4double bottomNeutronShieldGap       = 19. * cm;

		constexpr G4double nShield_mufflerLR_SizeX = 1.5 * m;
		constexpr G4double nShield_mufflerLR_SizeY = 1.5 * m;
		constexpr G4double nShield_mufflerLR_DistX = 1.0 * m;
		constexpr G4double nShield_mufflerLR_DistY = 0.5 * m;

		constexpr G4double nShield_mufflerFB_SizeX = 1.2 * m;
		constexpr G4double nShield_mufflerFB_SizeY = 1.35 * m;
		// constexpr G4double nShield_mufflerFB_DistX = 1.0 * m;
		// constexpr G4double nShield_mufflerFB_DistY = 0.5 * m;

		// Boric acid panel geometry
		constexpr G4double thickness_BoricAcid        = 3. * mm;
		constexpr G4double thickness_BoricAcidHousing = 1. * mm;
		constexpr G4double thickness_BoricAcidSpacing = 10. * mm;

		// B4C rubber shield
		constexpr G4double boronCarbideRubber_thickness      = 4. * mm;
		constexpr G4double boronCarbideRubber_UnderPb_XYsize = 100. * cm;
		constexpr G4double boronCarbideRubber_GapFromSSOVC   = 10. * cm;

		// G10 material holders between SS and 50K-SHIELD Copper
		constexpr G4double Cuholder_R     = 2.540 * cm / 2; // outer radius
		constexpr G4double G10holderXY[5] = {162.5 * mm, 190 * mm, 190 * mm, 230 * mm, 260 * mm};

		// HBeam parameters
		constexpr G4double Bplate_size       = 300. / 2;
		constexpr G4double Bplate_thick_half = 40. / 2.;
		constexpr G4double Hbeam1_length     = 3460.1 / 2;
		constexpr G4double Hbeam2_length     = 800. / 2;

		/// Gantry geomtety: length of vertical beam
		constexpr G4double Horizontal_beamX_length = 3000. / 2;
		/// Gantry geomtety: length of vertical beam (Large?)
		constexpr G4double Horizontal_beamXLarge_length = 3000. / 2;
		/// Gantry geomtety: length of horizontal beam
		constexpr G4double Horizontal_beamY_length = 1300. / 2;
		/// Gantry geomtety: length of horizontal beam (Large?)
		constexpr G4double Horizontal_beamYLarge_length = 1300. / 2;

		constexpr double x_beam[4] = {1.0, -1.0, -1.0, 1.0};
		constexpr double y_beam[4] = {1.0, 1.0, -1.0, -1.0};

		constexpr double HoLayer[3] = {1000., 2925., 3480.1};
		constexpr double HoLine[2]  = {1, -1};

		constexpr double HoLayerY[3] = {1000., 2925., 3480.1};
		constexpr double HoLineY[2]  = {1, -1};

		/// X dimension of superconducting magnet shield (SMS) box
		constexpr G4double scShieldSizeX = 252 * mm; //248 * mm;
		/// Y dimension of superconducting magnet shield (SMS) box
		constexpr G4double scShieldSizeY = 252 * mm; //245 * mm;
		/// Z dimension of superconducting magnet shield (SMS) box
		constexpr G4double scShieldSizeZ = 371 * mm;

		//constexpr G4double type3_scShieldSizeX = 252 * mm;
		//constexpr G4double type3_scShieldSizeZ = 371 * mm;

		/// Thickness of lead shield at the surface of SMS geometry
		constexpr G4double scShieldLeadThickness = 2. * mm;

		/// Thickness of bottom panel of copper frame for SMS
		constexpr G4double scsCuFrameBottomThickness = 6 * mm;

		/// Cut position X of bottom panel of copper frame
		constexpr G4double scsCuFrameBottomCutPosX = 102.95 * mm;
		/// Cut position Y of bottom panel of copper frame
		constexpr G4double scsCuFrameBottomCutPosY = 101.45 * mm;
		/// Angle of cutting angle for the bottom panel (You should do the math a little bit to
		/// understand this values lol)
		constexpr G4double scsCuFrameBottomCutAngle = 45. * deg;

		constexpr G4double scsCuFrameBottomHoleRadius = 30 * mm;
		constexpr G4double scsCuFrameBottomHoleDistX  = 61.25 * mm;
		constexpr G4double scsCuFrameBottomHoleDistY  = 62. * mm;

		constexpr G4double scsCuFrameTopDiskRadius    = 204. * mm;
		constexpr G4double scsCuFrameTopDiskThickness = 6. * mm;

		constexpr G4double scsCuFrameTopBlock1SizeX     = 23.5 * mm;
		constexpr G4double scsCuFrameTopBlock1SizeY     = 90 * mm;
		constexpr G4double scsCuFrameTopBlock2SizeX     = 90 * mm;
		constexpr G4double scsCuFrameTopBlock2SizeY     = 22 * mm;
		constexpr G4double scsCuFrameTopLongBlock1SizeX = 7.5 * mm;
		constexpr G4double scsCuFrameTopLongBlock1SizeY = 123.5 * mm;
		constexpr G4double scsCuFrameTopLongBlock2SizeX = 122 * mm;
		constexpr G4double scsCuFrameTopLongBlock2SizeY = 7.5 * mm;
		constexpr G4double scsCuFrameTopBlocksThickness = 15 * mm;

		constexpr G4double scsCuFrameHoleAtTopBlockRadius = 4.5 * mm;
		constexpr G4double scsCuFrameHoleAtTopBlockLength = 15 * mm;
		constexpr G4double scsCuFrameHoleAtTopBlockDist   = 15 * mm;

		constexpr G4double scsCuFrameSkeletonThickness = 7.5 * mm;
		constexpr G4double scsCuFrameSkeletonWidth     = 15 * mm;
		constexpr G4double scsCuFrameSkeletonDistance  = 80. * mm;

		/// Thickness for upper panel of copper frame
		constexpr G4double upperPanelThick = 10. * mm;
		/// Size X for upper panel of copper frame
		constexpr G4double upperPanelSizeX = 196. * mm;
		/// Size Y for upper panel of copper frame
		constexpr G4double upperPanelSizeY = 196. * mm;

		constexpr G4double innerPartDepth = 5. * mm;
		constexpr G4double innerPartSizeX = 166. * mm;
		constexpr G4double innerPartSizeY = 151. * mm;

		constexpr G4double innerMajorHoleSizeX       = 56. * mm;
		constexpr G4double innerMajorHoleSizeY       = 54. * mm;
		constexpr G4double innerMajorHoleDistanceEX  = 9.7 * mm;
		constexpr G4double innerMajorHoleDistanceEY  = 5.5 * mm;
		constexpr G4double innerMinorHole1DistanceEC = 14.5 * mm;
		constexpr G4double innerMinorHole1SizeX      = 20. * mm;
		constexpr G4double innerMinorHole1SizeY      = 20. * mm;
		constexpr G4double innerMinorHole2SizeX      = 20. * mm;
		constexpr G4double innerMinorHole2SizeY      = 27. * mm;

		constexpr G4double edgeHoleRadius     = 5. * mm;
		constexpr G4double edgeHoleDistanceEX = 7. * mm;
		constexpr G4double edgeHoleDistanceEY = 7. * mm;
		constexpr G4double edgeHoleDistanceX  = 50. * mm;
		constexpr G4double edgeHoleDistanceY  = 50. * mm;

		constexpr G4double cornerCuttedAngle = 45. * deg;
		constexpr G4double innerCarvingTol   = 1. * mm;

		constexpr G4double type3_panelFrameThick = 5. * mm;
		constexpr G4double type3_subPanelThick   = 8. * mm;
		constexpr G4double type3_upperPanelThick = type3_panelFrameThick+type3_subPanelThick;
		constexpr G4double type3_innerPartDepth  = 8. * mm;
		constexpr G4double type3_innerPartSizeX  = 161. * mm;
		constexpr G4double type3_innerPartSizeY  = 131. * mm;

		constexpr G4double type3_innerMajorHoleSizeX = 58. * mm;
		constexpr G4double type3_innerMajorHoleSizeY = 50.5 * mm;
		constexpr G4double type3_innerMajorHoleDistanceX = 10. * mm;
		constexpr G4double type3_innerMajorHoleDistanceY = 7.5 * mm;

		constexpr G4double type3_innerMinorHole1SizeX = 15. * mm;
		constexpr G4double type3_innerMinorHole1SizeY = 15. * mm;
		constexpr G4double type3_innerMinorHole2SizeX = 20. * mm;
		constexpr G4double type3_innerMinorHole2SizeY = 28. * mm;
		constexpr G4double type3_innerMinorHoleDistance = 17.5 * mm;

		constexpr G4double type3_edgeHoleRadius = 4.5 * mm;
		constexpr G4double type3_edgeHoleDistanceEX = 8. * mm;
		constexpr G4double type3_edgeHoleDistanceEY = 8. * mm;
		constexpr G4double type3_edgeHoleDistanceX  = 45. * mm;
		constexpr G4double type3_edgeHoleDistanceY  = 45. * mm;

		// for Weight
		constexpr G4double weightForUpperPanelSizeY = 62.5 * mm;
		constexpr G4double weightForUpperPanelSizeX = 58. * mm;
		constexpr double weightForUpperPanelThick[4] = {6, 14.5, 12.5, 0};

		// for additional height
		constexpr G4double additionalHeightSizeX = 146. * mm;
		constexpr G4double additionalHeightSizeY = 131. * mm;
		constexpr G4double type3_additionalHeightThick = 7. * mm;

		// From Build_I_DetectorArrayBottomPanel
		constexpr G4double bottomPanelSizeX = 166. * mm;
		constexpr G4double bottomPanelSizeY = 151. * mm;
		constexpr G4double bottomPanelThick = 12.50 * mm;

		constexpr G4double bottomPanelFrameThick = 5 * mm;
		constexpr G4double bottomPanelFrameWidth = 15. * mm;
		constexpr G4double bottomPanelFrameBigWidth = 25. * mm;

		constexpr G4double topBoxThick      = 7.5 * mm;
		constexpr G4double topBoxSizeY      = 30. * mm;
		constexpr G4double topHoleSizeY = 20. * mm;
		constexpr G4double topHoleThick = topBoxThick - bottomPanelFrameThick;
		constexpr G4double topSmallDis = 13. * mm;

		constexpr G4double baseSmallHoleSizeX = 9. * mm;
		constexpr G4double baseBigHoleSizeX   = 19. * mm;
		constexpr G4double baseHoleSizeY = 7. * mm;
		constexpr G4double longFrameHoleSizeY = 9. * mm;
		constexpr G4double longFrameHoleSizeX = 25.5 * mm;
		constexpr G4double cutterSizeMinorX = 25. * mm;
		constexpr G4double cutterSizeMajorX = 47. * mm;

		// From Build_I_CopperFrameUpper & Bottom
		constexpr G4double copperFrameSizeX  = 77. * mm;
		constexpr G4double copperFrameSizeY  = 74. * mm;
		constexpr G4double copperFrameSizeZ  = 4.64 * mm;
		constexpr G4double copperFrameXWidth = 1 * cm;
		constexpr G4double copperFrameYWidth = 1.005 * cm;
		constexpr G4double copperFrameHoleDiameter = 5.5 * mm;
		constexpr G4double copperFrameHoleDistance = 5 * mm;

		constexpr G4double type3_copperFrameSizeX = 78. * mm;
		constexpr G4double type3_copperFrameSizeY = 73. * mm;
		constexpr G4double type3_copperFrameSizeZ = 2.5 * mm;
		constexpr G4double type3_copperFrameXWidth = 1 * cm;
		constexpr G4double type3_copperFrameYWidth = 1 * cm;

		constexpr G4double type3_copperFrameHole1SizeY = 7.5 * mm;
		constexpr G4double type3_copperFrameHole2SizeX = 14.5 * mm;
		constexpr G4double type3_copperFrameHole3SizeX = 8 * mm;
		constexpr G4double type3_frameHole1Dist = -10.75 * mm;
		constexpr G4double type3_frameHole2Dist = 19.75 * mm;

		constexpr G4double type3_copperFrameHoleWidth = 5 * mm;
		constexpr G4double type3_copperFrameHoleDiameter = 4.3 * mm;
		constexpr G4double type3_copperFrameHoleDistance = 5 * mm;

		constexpr G4double type3_smallBlockSizeX = 10 * mm;
		constexpr G4double type3_smallBlockSizeY = 10 * mm;
		constexpr G4double type3_smallBlock1SizeZ = 9.5 * mm;
		constexpr G4double type3_smallBlock2SizeZ = 3. * mm;
		constexpr G4double type3_smallBlock3SizeZ = 8.5 * mm;

		constexpr G4double type3_bottomSupportASizeX = 0.8 * cm;
		constexpr G4double type3_bottomSupportASizeY = type3_copperFrameSizeY;
		constexpr G4double type3_bottomSupportASizeZ = 3. * mm;
		constexpr G4double type3_bottomSupportBSizeX = 2.4 * cm;
		constexpr G4double type3_bottomSupportBSizeY = type3_copperFrameSizeY;
		constexpr G4double type3_bottomSupportBSizeZ = 3. * mm;
		constexpr G4double type3_bottomSupportADistX = (8.5 + 3) * mm;
		constexpr G4double type3_bottomSupportBDistX = (4+6+5.75) * mm;

		constexpr G4double smallBlock1SizeX = 10 * mm;
		constexpr G4double smallBlock1SizeY = 10 * mm;
		constexpr G4double smallBlock1SizeZ = 6 * mm;

		constexpr G4double smallBlock2SizeX = 1 * cm;
		constexpr G4double smallBlock2SizeY = 1 * cm;
		constexpr G4double smallBlock2SizeZ = 4 * mm;

		constexpr G4double smallBlock1DistX = -5 * mm;
		constexpr G4double smallBlock1DistY = -15 * mm;

		constexpr G4double smallBlock2DistX = 5 * mm;
		constexpr G4double smallBlock2DistY = 5 * mm;

		constexpr G4double bottomSupportASizeX = 1.5 * cm;
		constexpr G4double bottomSupportASizeY = copperFrameSizeY;
		constexpr G4double bottomSupportASizeZ = copperFrameSizeZ / 2.;
		constexpr G4double bottomSupportBSizeX = 1.9 * cm;
		constexpr G4double bottomSupportBSizeY = copperFrameSizeY;
		constexpr G4double bottomSupportBSizeZ = copperFrameSizeZ / 2.;
		constexpr G4double bottomSupportADistX = 4 * mm;
		constexpr G4double bottomSupportBDistX = 2 * mm;
		constexpr G4double bottomSupportDepth  = 1 * mm;

		// From Build_I_DetectorArray
		constexpr G4double arraySizeX             = 196. * mm;
		constexpr G4double arraySizeY             = 196. * mm;
		//constexpr G4double arrayHeight            = 357.90 * mm;
		constexpr G4double arrayHeight      = scShieldSizeZ-scsCuFrameBottomThickness-scShieldLeadThickness;
		constexpr G4double type1_detModuleMarginX = 0 * mm;
		constexpr G4double type1_detModuleMarginY = 4 * mm;
		constexpr G4double type2_detModuleMarginX = 0 * mm;
		constexpr G4double type2_detModuleMarginY = 3 * mm;

		constexpr G4double towerRodRadius    = 3.5 * mm;
		constexpr G4double towerRodTopRadius = 5 * mm;
		constexpr G4double towerRodTopLength = 7 * mm;

		constexpr G4double detModuleBoltM3SRadius  = 3.48 * mm;
		constexpr G4double detModuleBoltM3SLength  = 3.0 * mm;
		constexpr G4double dMBM3SDistanceXFromXEnd = 14.5 * mm;

		constexpr G4double type1_distanceBetweenTower = 8 * mm;
		constexpr G4double type2_distanceBetweenTower = 6 * mm;

		// From Build_I_PhotonDetector
		// PhotonFrame total size
		constexpr G4double photonFrameSizeX = 57. * mm;
		constexpr G4double photonFrameSizeY = 73.8 * mm;
		constexpr G4double photonFrameSizeZ = 6. * mm;

		constexpr G4double type3_photonFrameSizeX = 55 * mm;
		constexpr G4double type3_photonFrameSizeY = 66 * mm;
		constexpr G4double type3_photonFrameSizeZ = 8.5 * mm;
		constexpr G4double type3_photonFrameTopSpace = 1.0 * mm;

		// PhotonFrame Hole and Base
		constexpr G4double type1_photonFrameHoleRadius     = 22.5 * mm;
		constexpr G4double type1_photonFrameBaseWidth      = 1 * cm;
		constexpr G4double type1_photonFrameBaseSpaceWidth = 1.4 * cm;
		constexpr G4double type1_photonFrameBaseThick      = 2 * mm;

		constexpr G4double type2_photonFrameHoleRadius = 24.5 * mm;
		constexpr G4double type2_photonFrameRoofThick  = 0.2 * mm;
		constexpr G4double type2_photonFrameRoofWidth  = 0.5 * mm;
		constexpr G4double type2_photonFrameBaseThick  = 3.5 * mm;
		constexpr G4double type2_photonFrameBaseWidthX = 7.5 * mm;
		constexpr G4double type2_photonFrameBaseWidthY = 6.5 * mm;
		constexpr G4double type2_photonFrameHeartLeapMinorRadius = 5 * mm;
		constexpr G4double type2_photonFrameHeartAngle           = 20. * deg;

		constexpr G4double type3_photonFrameHoleRadius = 16. * mm;
		constexpr G4double type3_photonFrameHoleThick  = 2. * mm;
		constexpr G4double type3_photonFrameHoleBoxSizeX = 17.95 * mm;

		constexpr G4double type3_photonFrameBaseThick  = 6.5 * mm; //3.25 * mm;
		constexpr G4double type3_photonFrameBaseWidthX = 7.0 * mm;
		constexpr G4double type3_photonFrameBaseWidthY = 7.0 * mm;

		// PhotonFrame Peek
		constexpr G4double type2_photonFramePeekHoleRadius = 2.9 * mm;
		constexpr G4double type2_photonFramePeekHoleThick  = 2 * mm;

		constexpr G4double type2_photonFramePeekRdaius = 2.775 * mm;
		constexpr G4double type2_photonFramePeekThick  = 1.25 * mm;

		constexpr G4double type3_photonFramePeekHoleRadius     = 0.35 * mm;
		constexpr G4double type3_photonFramePeekHoleThick      = 1 * mm;
		constexpr G4double type3_photonFramePeekRdaius         = 1.5 * mm;
		constexpr G4double type3_photonFramePeekThick          = 2.5 * mm;
		constexpr G4double type3_photonFramePeekBoltRadius     = 2. * mm;
		constexpr G4double type3_photonFramePeekBoltHeadRadius = 2.5 * mm;

		// PhotonFrame Clamp
		constexpr G4double type2_photonFrameClampBox1SizeXY = 5.5 * mm;
		constexpr G4double type2_photonFrameClampBox2SizeX  = 2 * mm;
		constexpr G4double type2_photonFrameClampBox2SizeY  = 7. * mm;
		constexpr G4double type2_photonFrameClampBoltRadius = 1.6 * mm;
		constexpr G4double type2_photonFrameClampThick      = 2. * mm;

		constexpr G4double type3_photonFrameClampRadius = 6.5 / 2 * mm;
		constexpr G4double type3_photonFrameClampThick  = 1.2 * mm;
		constexpr G4double type3_photonFrameClampSizeX  = 2.5 * mm;
		constexpr G4double type3_photonFrameClampSizeZ  = 3.5 * mm;

		// PhotonFrame Germanium Wafer
		constexpr G4double germaniumWaferRadius   = 2.54 * cm;
		constexpr G4double germaniumWaferThick    = 0.003 * cm; // 300 um
		constexpr G4double germaniumWaferVacThick = 0.001 * cm; // 100 um
		constexpr G4double waferHolderWidth       = 5 * mm;
		constexpr G4double waferHolderLength      = 25 * mm;
		constexpr G4double waferHolderThick       = 0.5 * mm;

		constexpr G4double type1_germaniumWaferSlotDepth = 1 * mm;
		constexpr G4double type1_germaniumWaferSlotThick = 2 * mm;

		constexpr G4double type2_germaniumWaferSlotThick = 0.4 * mm;

		constexpr G4double type3_germaniumWaferRadius      = 2.54 * cm;
		constexpr G4double type3_germaniumWaferThick       = 0.0028 * cm; // 280 um;

		// 3 gold film above GeWafer
		constexpr G4double geWaferUpperGoldRadius         = 2.5 * mm;
		constexpr G4double geWaferUpperGoldThickness      = 0.0003 * mm;
		constexpr G4double geWaferUpperGoldRadialDistance = 6 * mm;

		constexpr G4double type3_geWaferUpperGoldRadius    = 3. * mm;
		constexpr G4double type3_geWaferUpperGoldThickness = 0.0003 * mm;

		// PhotonFrame Copper Bolts
		constexpr G4double type3_photonFrameCuBoltsRadius     = 3.0 * mm;
		//constexpr G4double type3_photonFrameCuBoltsHoleRadius = 3. * mm;
		//constexpr G4double type3_photonFrameCuBoltsThick      = 2 * mm;

		// PhotonFrame Brass Bolts
		constexpr G4double type3_photonPhononBoltsRadius     = 3.0 / 2. * mm;
		constexpr G4double type3_photonPhononBoltsThick      = 12 * mm;
		constexpr G4double type3_photonPhononBoltsHeadRadius = 4.0 / 2. * mm;
		constexpr G4double type3_photonPhononBoltsHeadThick  = 2. * mm;

		// From Build_I_SingleDetectorModule
		constexpr G4double type3_moduleSupportRodRadius = 4. * mm;

		constexpr G4double moduleSupportRodRadius = 3.5 * mm;

		constexpr G4double reflectorThick   = 64 * um;
		constexpr G4double reflectorLengthX = 5.5 * cm;
		constexpr G4double reflectorLengthY = 5.2 * cm;
		//constexpr G4double type3_reflectorLengthX = 5.7 * cm;
		//constexpr G4double type3_reflectorLengthY = 5.7 * cm; //5.2 * cm;

		// a gold film below crystal
		constexpr G4double crystalBottomGoldRadius    = 1 * cm;
		constexpr G4double crystalBottomGoldThickness = 0.0003 * mm;

		// From Build_I_Top/FB/LRScintillator
		constexpr G4double scint_Thick_half             = 50.000 * mm / 2.;
		constexpr G4double scint_muffler_Thick_half     = 30.000 * mm / 2.;
		constexpr G4double scint_bottom_Thick_half      = 30.000 * mm / 2.;
		constexpr G4double scint_FlatL_ratio            = 0.7;
		constexpr G4double scint_muffler_FlatL_ratio    = 1.0;
		constexpr G4double scint_reflector_thick        = 1 * mm;
		constexpr G4double scint_muffler_sideRefl_thick = 0.1 * mm;

		constexpr G4double topScint_trapZ_half     = 368.855 * mm / 2.;
		constexpr G4double topScint_trapX1_half    = 37.000 * mm / 2.;
		constexpr G4double topScint_trapX2_half    = 381.000 * mm / 2.;
		constexpr G4double topScint_boxHeight_half = 1725.000 * mm / 2.;
		constexpr G4double topScint_boxWidth_half  = 762.063 * mm / 2.;

		constexpr G4double sideFBScint_trapZ_half     = 250.000 * mm / 2.;
		constexpr G4double sideFBScint_trapX1_half    = 37.000 * mm / 2.;
		constexpr G4double sideFBScint_trapX2_half    = 275.000 * mm / 2.;
		constexpr G4double sideFBScint_boxHeight_half = 1700.000 * mm / 2.;
		constexpr G4double sideFBScint_boxWidth_half  = 550.000 * mm / 2.;

		constexpr G4double sideLRScint_trapZ_half     = 500.350 * mm / 2.;
		constexpr G4double sideLRScint_trapX1_half    = 46.699 * mm / 2.;
		constexpr G4double sideLRScint_trapX2_half    = 513.333 * mm / 2.;
		constexpr G4double sideLRScint_boxHeight_half = 1700.000 * mm / 2.;
		constexpr G4double sideLRScint_boxWidth_half  = 600.000 * mm / 2.;

		constexpr G4double mufflerFBScint_trapZ_half     = 570.000 * mm / 2.;
		constexpr G4double mufflerFBScint_trapX1_half    = 40.000 * mm / 2.;
		constexpr G4double mufflerFBScint_trapX2_half    = 600.000 * mm / 2.;
		constexpr G4double mufflerFBScint_boxHeight_half = 1030.000 * mm / 2.;
		constexpr G4double mufflerFBScint_boxWidth_half  = 600.000 * mm / 2.;
		constexpr G4double mufflerFBScint_YDistanceFromLead = 0.45 * m;

		constexpr G4double mufflerLRScint_trapZ_half     = 570.000 * mm / 2.;
		constexpr G4double mufflerLRScint_trapX1_half    = 40.000 * mm / 2.;
		constexpr G4double mufflerLRScint_trapX2_half    = 600.000 * mm / 2.;
		constexpr G4double mufflerLRScint_boxHeight_half = 1030.000 * mm / 2.;
		constexpr G4double mufflerLRScint_boxWidth_half  = 600.000 * mm / 2.;
		constexpr G4double mufflerLRScint_YDistanceFromLead = 0.45 * m;


		// From Build_I_MuonVetoPMT
		constexpr G4double muonVetoPMTBacksideThickness = 1 * mm;
		constexpr G4double muonVetoPMTTubeThickness     = 2 * mm;
		constexpr G4double muonVetoPMTWindowThickness   = 4 * mm;
		constexpr G4double muonVetoPMTSizeR             = 51.5 * mm / 2;
		constexpr G4double muonVetoPMTSizeZ             = 109 * mm;
		constexpr G4double muonVetoPMTTopSizeRmin       = 0;
		constexpr G4double muonVetoPMTVacuTopSizeRmin   = 0;
		constexpr G4double muonVetoPMTGreaseShieldThick = 5 * um;
		constexpr G4double muonVetoPMTGreaseThickness   = 2 * mm;

		// From GeometriesInDetModuleOf / GeometriesInOpticalDetOf
		// SSCB BOLTS: M4Top Bolts (Deprecated in AMoRE-I simulation geometry)
		constexpr G4double detModuleBoltM4TopRadius = 4.11 * mm;
		constexpr G4double detModuleBoltM4TopLength = 5.0 * mm;
		constexpr G4double detModuleBoltM4Radius    = 4.6 * mm;
		constexpr G4double detModuleBoltM4Length    = 4.0 * mm;

		// SSCB BOLTS: M3-6
		constexpr G4double detModuleBoltM3Radius  = 3.25 * mm;
		constexpr G4double detModuleBoltM3Length  = 3.0 * mm;
		constexpr G4double dMBM3DistanceXFromXEnd = 11.625 * mm;

		// SSCB BOLTS: M3-4
		constexpr G4double detModuleBoltM3FRadius     = 2.925 * mm;
		constexpr G4double detModuleBoltM3FLength     = 3.0 * mm;
		constexpr G4double detMBM3FAtSupportADistance = 0.6 * mm;
		constexpr G4double detMBM3FAtSupportBDistance = 2.5 * mm;

		// additional Bolts : AB
		// beside of pilar
		constexpr G4double detModuleBoltABRadius  = 3.7625 * mm;
		constexpr G4double detModuleBoltABLength  = 3.5 * mm;
		constexpr G4double dMBABDistanceXFromXEnd = 10.7525 * mm;

		constexpr G4double detModuleTopPinConnectorSizeX = 20. * mm;
		constexpr G4double detModuleTopPinConnectorSizeY = 4. * mm;
		constexpr G4double detModuleB1PinConnectorSizeX  = 8. * mm;
		constexpr G4double detModuleB1PinConnectorSizeY  = 4. * mm;
		constexpr G4double detModuleB2PinConnectorSizeX  = 4. * mm;
		constexpr G4double detModuleB2PinConnectorSizeY  = 4. * mm;
		constexpr G4double detModulePinConnectorSizeZ    = 1. * mm;

		// Adhesive for optical detector
		constexpr G4double detModuleFT2850Thickness = 0.13 * mm;

		constexpr G4double detModuleFT2850ASizeX = 5. * mm;
		constexpr G4double detModuleFT2850ASizeY = 8. * mm;

		constexpr G4double detModuleFT2850BSizeX = 5. * mm;
		constexpr G4double detModuleFT2850BSizeY = 6. * mm;

		constexpr G4double detModuleFT2850CSizeX = 20. * mm;
		constexpr G4double detModuleFT2850CSizeY = 8. * mm;

		constexpr G4double dMFT2850AYDistanceFromCenter = 8.75 * mm;
		constexpr G4double dMFT2850BYDistanceFromCenter = 11.25 * mm;

		constexpr G4double detModuleFT2850AForCrystalSizeX = 4. * mm;
		constexpr G4double detModuleFT2850AForCrystalSizeY = 4. * mm;

		constexpr G4double detModuleFT2850BForCrystalSizeX = 8. * mm;
		constexpr G4double detModuleFT2850BForCrystalSizeY = 5. * mm;

		constexpr G4double detModuleFT2850CForCrystalSizeX = 4. * mm;
		constexpr G4double detModuleFT2850CForCrystalSizeY = 4. * mm;

		constexpr G4double dMFT2850AForCrystalXDistanceFromCenter = -9.5 * mm;

		constexpr G4double dMFT2850BForCrystalXDistanceFromCenter = -6.5 * mm;
		constexpr G4double dMFT2850BForCrystalYDistanceFromCenter = 30. * mm;

		constexpr G4double dMFT2850CForCrystalXDistanceFromCenter = -6. * mm;
		constexpr G4double dMFT2850CForCrystalYDistanceFromCenter = -30. * mm;

		constexpr G4double photDetPCBThick           = 0.115 * mm;
		constexpr G4double photDetPCBBaseBoxSizeX    = 16.8 * mm;
		constexpr G4double photDetPCBBaseBoxSizeY    = 29.26 * mm;
		constexpr G4double photDetPCBTriangle0Length = 11.8 * mm;
		constexpr G4double photDetPCBTriangle1Length = 9.108 * mm;
		constexpr G4double photDetPCBTriangle2Length = 4.780 * mm;
		constexpr G4double photDetPCBOriginDistX     = 0 * mm;
		constexpr G4double photDetPCBOriginDistY     = -25 * mm;

		constexpr G4double type3_photDetPCBBaseBoxThick     = 5. * mm;
		constexpr G4double type3_photDetPCBBaseBoxSizeX     = 26.99 * mm;
		constexpr G4double type3_photDetPCBBaseBoxSizeY     = 64 * mm;
		constexpr G4double type3_photDetPCBBaseBoxSizeX1    = 5.53 * mm;
		constexpr G4double type3_photDetPCBBaseBoxSizeY1    = 20.48 * mm;
		constexpr G4double type3_photDetPCBBaseBoxSizeY2    = 15 * mm;
		constexpr G4double type3_photDetPCBBaseInnerThick   = 4. * mm;
		constexpr G4double type3_photDetPCBBaseInnerSizeX1  = 1. * mm;
		constexpr G4double type3_photDetPCBBaseInnerSizeY1  = 22.53 * mm;
		constexpr G4double type3_photDetPCBBaseInnerSizeY2  = 19. * mm;
		constexpr G4double type3_photDetPCBBaseInnerSizeY3  = 23.53 * mm;
		constexpr G4double type3_photDetPCBBaseInnerSizeY4  = 24.47 * mm;
		constexpr G4double type3_photDetPCBBaseInnerSizeY5  = 21.2 * mm;
		constexpr G4double type3_photDetPCBBaseSpaceSizeX   = 4.52 * mm;

		constexpr G4double photDetPCBHoleSizeX = 5.5 * mm;
		constexpr G4double photDetPCBHoleSizeY = 8.0 * mm;
		constexpr G4double photDetPCBHoleSizeZ = 10.0 * mm;

		constexpr G4double photDetPCBBoltsHeadHeight = 1 * mm;
		constexpr G4double photDetPCBBoltsHoleRadius = 1.25 * mm;
		constexpr G4double photDetPCBBoltsHoleHeight = 4.0 / 2. * mm;
		constexpr G4double photDetPCBBoltsHoleDistX  = 2.3 * mm;
		constexpr G4double photDetPCBBoltsHoleDistY  = 2.5 * mm;

		constexpr G4double photDetPCBSolderBox1SizeX = 1. * mm;
		constexpr G4double photDetPCBSolderBox1SizeY = 5.06 * mm;
		constexpr G4double photDetPCBSolderBox2SizeX = 4. * mm;
		constexpr G4double photDetPCBSolderBox2SizeY = 1.266 * mm;
		constexpr G4double photDetPCBEpoxyBoxSizeX   = 1.1248 * mm;
		constexpr G4double photDetPCBEpoxyBoxSizeY   = 1.1248 * mm;
		constexpr G4double photDetPCBEpoxyBoxSizeZ   = 0.46 * mm;

		constexpr G4double photDetPCBEpoxyBox1Dist = 15 * mm;
		constexpr G4double photDetPCBEpoxyBox2Dist = 18 * mm;

		constexpr G4double type3_crystalPCBWidth     = 19.8 * mm;
		constexpr G4double type3_crystalPCBLength    = 39. * mm;
		constexpr G4double type3_crystalPCBGrooveX   = 3.8 * mm;
		constexpr G4double type3_crystalPCBGrooveY   = 3.2 * mm;
		constexpr G4double type3_crystalPCBBigY1     = 32.8 * mm;
		constexpr G4double type3_crystalPCBBigY2     = 22.2 * mm;
		constexpr G4double type3_crystalPCBShortY    = 3. * mm;
		constexpr G4double type3_crystalPCBDownBoxX1 = 8. * mm;
		constexpr G4double type3_crystalPCBDownBoxX2 = 9.35 * mm;
		constexpr G4double type3_crystalPCBDownBoxY  = 6. * mm;

		constexpr G4double crystalPCBThick           = 0.115 * mm;
		constexpr G4double crystalPCBWidth           = 14.9 * mm;
		constexpr G4double crystalPCBUpBoxLength     = 11.2 * mm;
		constexpr G4double crystalPCBUpBoxXDistance  = 4.67 * mm;
		constexpr G4double crystalPCBUpBoxYDistance  = 23.4 * mm;
		constexpr G4double crystalPCBMiddleBoxLength = 23.2 * mm;
		constexpr G4double crystalPCBDownBoxWidth    = 6.25 * mm;
		constexpr G4double crystalPCBDownBoxLength   = 4.35 * mm;
		constexpr G4double crystalPCBOriginDistX     = -9.25 * mm;
		constexpr G4double crystalPCBOriginDistY     = 10.3 * mm;

		constexpr G4double crystalPCBSolderBox1SizeX = 1.4 * mm;
		constexpr G4double crystalPCBSolderBox1SizeY = 4.74 * mm;
		constexpr G4double crystalPCBSolderBox2SizeX = 3.68 * mm;
		constexpr G4double crystalPCBSolderBox2SizeY = 1.8 * mm;
		constexpr G4double crystalPCBSolderBox1XDist = 3.5 * mm;
		constexpr G4double crystalPCBSolderBox1YDist = 0. * mm;
		constexpr G4double crystalPCBSolderBox2XDist = -3.5 * mm;
		constexpr G4double crystalPCBSolderBox2YDist = 3. * mm;

		constexpr G4double crystalPCBEpoxyBox1SizeXY = 1.286 * mm;
		constexpr G4double crystalPCBEpoxyBox2SizeX  = 1.74 * mm;
		constexpr G4double crystalPCBEpoxyBox2SizeY  = 0.95 * mm;
		constexpr G4double crystalPCBEpoxyBox12SizeZ = 0.46 * mm;
		constexpr G4double crystalPCBEpoxyBox1XDist  = -4.1 * mm;
		constexpr G4double crystalPCBEpoxyBox1YDist  = -10. * mm;
		constexpr G4double crystalPCBEpoxyBox2XDist  = 1.5 * mm;
		constexpr G4double crystalPCBEpoxyBox2YDist  = -6. * mm;

		constexpr G4double crystalPeekThick  = 0.42 * mm;
		constexpr G4double crystalPeekRadius = 4.787 * mm;
		constexpr G4double crystalPeekDist1  = 15 * mm;
		constexpr G4double crystalPeekDist2  = 2 * mm;

		constexpr G4double crystalHeaterSpacing = 0.5 * mm;
		constexpr G4double crystalHeaterXDist   = 0. * mm;
		constexpr G4double crystalHeaterYDist   = 15.5 * mm;

		constexpr G4double crystalHeaterBoxSizeX       = 3.7 * mm;
		constexpr G4double crystalHeaterBoxSizeY       = 3.8 * mm;
		constexpr G4double crystalHeaterBoxSizeZ       = 0.38 * mm;
		constexpr G4double crystalHeaterElectrodeSizeX = 0.04 * mm;
		constexpr G4double crystalHeaterElectrodeSizeY = 2.695 * mm;
		constexpr G4double crystalHeaterElectrodeSizeZ = 100 * nm;

		constexpr G4double crystalHeaterElectrodeStartX  = -1.43 * mm;
		constexpr G4double crystalHeaterElectrodeSpacing = 0.26 * mm;
	} // namespace AMoRE_I

	namespace AMoRE_200 {
		// Outer Shield  ------------
		// Rock geometry
		//constexpr G4double cavern_nMode_radius = 6. * m;
		constexpr G4double cavern_nMode_radius = 12. * m;

		// Detector Room
		constexpr G4double room_dist_x = 3.17 * m;
		constexpr G4double room_dist_y = 4.86 * m;

		constexpr G4double top_plate_x = 2200 * mm;
		constexpr G4double top_plate_y = 2000 * mm;

		// Muon Veto (PS)
		constexpr G4double longPSlength  = 1680 * mm ; //1662 * mm; // real design 1670 mm ?? 
		constexpr G4double shortPSlength = 121 * cm;
		constexpr G4double profile_thickness = 10 * cm;
		constexpr G4double veto_frame_width = 25 * mm;
		constexpr G4double veto_frame_thickness = 1 * mm;
		constexpr G4double veto_frame_space = 5 * mm;

		constexpr G4double bottom_veto_housingX = 3720 * mm;
		constexpr G4double bottom_veto_housingY = 3375 * mm;

		// water tank
		constexpr G4double waterhousing_thickness = 50. * mm;
		constexpr G4double watertank_thickness = 700. * mm;
		constexpr G4double watertank_top_thickness = 800. * mm;
		constexpr G4double PMTroom_thickness = 250. * mm;

		// Hat Inner Size & H-beam
		constexpr G4double HatInnerX = 5000 * mm / 2.;
		constexpr G4double HatInnerY = 5000 * mm / 2.;
		constexpr G4double HatInnerZ = 2200 * mm;
		constexpr G4double HatHBeam_size = 150. * mm;
		constexpr G4double HatHBeam_thickness = 10. * mm;
		constexpr G4double HatHBeam_in_thickness = 7. * mm;
		constexpr G4double HatHBeam_heightL = 5000. * mm;
		constexpr G4double HatHBeam_heightS = 2050. * mm;
		constexpr G4double HatAlPlate_thickness = 6. * mm;

		// Detector H-beam
		constexpr G4double DetHbeam_size = 120. * mm;
		constexpr G4double DetHbeam_height = 1017. * mm;
		constexpr G4double DetHbeamS_height = 700. * mm;
		constexpr G4double DetHbeam_thickness = 7. * mm;

		// H-Beam 
		constexpr G4double HBeam_size = 300. *mm;
		constexpr G4double HBeam_thickness = 15. * mm;
		constexpr G4double HBeam_housingX = 13120 * mm;
		constexpr G4double HBeam_housingY = 8000 * mm;
		//constexpr G4double HBeam_housingZ = 4550 * mm;
		constexpr G4double HBeam_housingDist = 760 * mm;

		// Pit
		constexpr G4double pitBox_x = 9600 * mm / 2.;
		constexpr G4double pitBox_z = 6020 * mm;
		//constexpr G4double pitBox_z = 6900 * mm;
		//constexpr G4double pitBox_x = 9100 * mm / 2.;
		//constexpr G4double pitBox_z = 6120 * mm;

		// Real cavern (Handeok mine cavern) model
		constexpr G4double cavern_pizza_radius      = 18.4 * m; //8.4 m; // Original: 15.4 m
		constexpr G4double cavern_pizza_angle_real  = 60. * deg;
		constexpr G4double cavern_pizza_thickness   = 21. *m ; // 11. m;
		constexpr G4double cavern_subpizza_radius   = 2.6 * m; // Original: 5.6 m
		constexpr G4double cavern_subpizza_angle    = 60 * deg;
		constexpr G4double cavern_totalpizza_dist   = 2.5 * m; //1 m;
		constexpr G4double cavern_loaf_thickness    = 21. *m;  //11 m;
		constexpr G4double cavern_loaf_height       = 11.4 * m; 
		constexpr G4double cavern_pizza_angle_tol   = 0. * deg;

		// InnerDetector ------
		// Crystal Tower
		constexpr G4int nModuleInTower = 9;
		constexpr G4int maxNumOfTower  = 70;
		extern const AMoRE200CrystalModuleInfo crystalModuleInfoList[maxNumOfTower];

		constexpr G4double Module_gap         = 20. * mm / 2.;
		constexpr G4double Module_radius      = 74. * mm / 2.; 
		constexpr G4double Module_width       = 18. * mm / 2.;
		constexpr G4double Module_base_thick  = 2. * mm; 
		constexpr G4double Module_thick       = 2.5 * mm;
		constexpr G4double Tower_height       = 700. * mm; //683. * mm;
		constexpr G4double Array_radius       = 940. * mm / 2. + solidBooleanTol;
		constexpr G4double Array_height       = 700. * mm + solidBooleanTol; 

		// Heater
		constexpr G4double crystalHeaterSpacing = 0.5 * mm;
		constexpr G4double crystalHeaterXDist   = 0. * mm;
		constexpr G4double crystalHeaterYDist   = 15.5 * mm;

		constexpr G4double crystalHeaterBoxSizeX       = 3.7 * mm;
		constexpr G4double crystalHeaterBoxSizeY       = 3.8 * mm;
		constexpr G4double crystalHeaterBoxSizeZ       = 0.38 * mm;
		constexpr G4double crystalHeaterElectrodeSizeX = 0.04 * mm;
		constexpr G4double crystalHeaterElectrodeSizeY = 2.695 * mm;
		constexpr G4double crystalHeaterElectrodeSizeZ = 100 * nm;

		constexpr G4double crystalHeaterElectrodeStartX  = -1.43 * mm;
		constexpr G4double crystalHeaterElectrodeSpacing = 0.26 * mm;

		// Araldite
		constexpr G4double araldite_thick = 32 * um;

		// Reflector
		constexpr G4double reflector_thick   = 64 * um; // 64 micro-meter
		constexpr G4double reflector_gapz    = 1 * mm;
		constexpr G4double reflector_gapx    = 3 * mm;
		constexpr G4double reflector_gapy    = 1.7 * mm;

		// Bottom Module Frame
		constexpr G4double bottomframe_height= 7.5 * mm;
		constexpr G4double frame_hole_depth  = 4.78 * mm; // instead of frame_hole_r
		constexpr G4double frame_hole_r      = 32.22 * mm; // !!!!
		constexpr G4double leg_height        = 3.5 * mm;
		constexpr G4double bottomframe_thick = 1.5 * mm;
		constexpr G4double bottomadd_size    = 2.5 * mm;
		constexpr G4double bottomhole1_size  = 5. * mm;
		constexpr G4double bottomhole2_size  = 3.3 * mm;
		constexpr G4double bottomhole3_size  = 7.5 * mm;
		constexpr G4double bottomhole_w      = 4.  * mm;
		constexpr G4double bottomhole_l      = 2.6 * mm;

		// Top Module Frame
		constexpr G4double topframe_height = 9. * mm;
		constexpr G4double topleg_height = 7.5 * mm;
		constexpr G4double topadd_height = 1.5 * mm;
		constexpr G4double topadd_length = 3.8 * mm;
		constexpr G4double tophole0_size = 3. * mm;
		constexpr G4double tophole1_size = 5.3 * mm;
		constexpr G4double tophole2_size = 4.2 * mm;

		// Post
		constexpr G4double post_height = 48.7 * mm; // !!!!
		constexpr G4double post_top_h  = 2. * mm;
		constexpr G4double post_radius = 8.5/2 * mm;
		constexpr G4double post_top_r  = 7.9/2 * mm;
		constexpr G4double post_bolt_r = 4.2/2 * mm;

		// clamp
		constexpr G4double clamp_wafer_w        = 6. * mm;
		constexpr G4double clamp_wafer_l        = 12.5 * mm;
		constexpr G4double clamp_wafer_holesize = 4.1 * mm;

		constexpr G4double clamp_top_w        = 11. * mm;
		constexpr G4double clamp_top_l        = 12 * mm;
		constexpr G4double clamp_top_box1_w   = 2.2 * mm - solidBooleanTol;
		constexpr G4double clamp_top_box1_l   = 3. * mm;
		constexpr G4double clamp_top_box_h    = 3. * mm;
		constexpr G4double clamp_top_holesize = 7.9 * mm;
		constexpr G4double clamp_top_holepos  = 4.5 * mm;

		constexpr G4double clamp_bot_h     = 4. * mm;
		constexpr G4double clamp_thick     = 1. * mm;

		// Wafer
		constexpr G4double wafer_radius   = 2.54 * cm;
		constexpr G4double wafer_thick    = 0.03 * cm; // 300 um
		constexpr G4double wafer_vac_thick = 0.01 * cm; // 100 um

		// Gold Film
		constexpr G4double upperGold_radius   = 1.5 * mm;
		constexpr G4double upperGold_thick    = 0.0003 * mm;
		constexpr G4double bottomGold_radius  = 1.4 * cm; //1 * cm;
		constexpr G4double bottomGold_thick   = 0.0003 * mm;

		//PCB & sensor
		constexpr G4double pcb_thick = 115 * um;
		constexpr G4double pcbbox_size = 9 * mm;
		constexpr G4double pcbtub_radius = 9 * mm;
		constexpr G4double lightsensor_size = 4.5 * mm;
		constexpr G4double heatsensor_size = 7 * mm;

		// Light detector
		constexpr G4double light_box1_w    = 12. * mm;
		constexpr G4double light_box2_w    = 8. * mm;
		constexpr G4double light_box3_w    = 13. * mm;   
		constexpr G4double light_box1_l    = 66 * mm;// !!!!
		constexpr G4double light_box2_l    = 40. * mm;// !!!!
		constexpr G4double light_box3_l    = tophole1_size;
		constexpr G4double light_box_thick = 1.5 * mm;

		// Heat detector
		constexpr G4double heat_box1_w   = 8. * mm;
		constexpr G4double heat_box2_w   = 10. * mm;
		constexpr G4double heat_box1_l  = 48.43 * mm; //!!!!
		constexpr G4double heat_box2_l  = 62.89 * mm; //!!!!
		constexpr G4double heat_box1_thick = 2. * mm;
		constexpr G4double heat_box2_thick = 1. * mm;

		// Screws
		constexpr G4double M5_size      = 5;
		constexpr G4double M5_height    = 3;
		constexpr G4double M5_12_height = 12;
		constexpr G4double M5_15_height = 16;
		constexpr G4double M4_size      = 4;
		constexpr G4double M4_height    = 2;
		constexpr G4double M4_4_height  = 4;
		constexpr G4double M4_6_height  = 6;

		// CMO Module
		constexpr G4double CMOcell_r         = 5.0 * cm / 2.;
		constexpr G4double CMOcell_h         = 5. * cm / 2.;

		// CMOArray
		constexpr G4double ArrayTargetMass = 183 * kg; //210. * kg;
		constexpr G4double Tower_Ratio = 1.3;
		constexpr G4double RSpacing = 1 * mm;
		constexpr G4double ZSpacing = 0 * cm;

		// Copper Frame
		constexpr G4double Copper_Fx      = 7.3 * cm;
		constexpr G4double Copper_Fy      = 7.3 * cm;
		constexpr G4double CopperFz      = 4.64 * mm;
		constexpr G4double CopBoxAZ_half = 6 * mm / 2.;

		// Stycast
		constexpr G4double stycast_gap        = 0.2 * mm;
		constexpr G4double stycast_thick_half = 0.13 * mm / 2.;

		// Photon Detector: Ge Wafer + Vaccum disk
		constexpr G4double PhotonDet_r     = 2.54 * 2 * cm / 2.;
		constexpr G4double PhotonDet_zsize = 0.6 * cm / 2.;
		constexpr G4double GeWafer_r       = 2.54 * 2 * cm / 2.;
		constexpr G4double GeWafer_zsize   = 0.004 * cm / 2.;   // 300 nm + 100 nm Vac
		constexpr G4double VacDisk_zsize   = 0.00001 * cm / 2.; // 100 nm

		// SSOVC
		constexpr G4double sst_zsize_half = 5.* cm / 2.; //3. * cm / 2.;
		constexpr G4double sst_xsize_half = 1800. * mm / 2.; //2500. * mm / 2.;
		constexpr G4double sst_ysize_half = 1500. * mm / 2.; //2000. * mm / 2.;

		constexpr G4double OVCthick  = 0.5 * cm; 
		constexpr G4double OVCthickB = 15. * mm;
		constexpr G4double ovc_gap   = 30 * mm;
		constexpr G4double ss_radius = 1270. * mm /2.; 
		constexpr G4double ss_inner_height_half = 2890.* mm /2.;

		// Cu 1,2,3,4
		constexpr G4double cu4thick   = 4. * mm; 
		constexpr G4double cu4thickB  = 6. * mm;
		constexpr G4double cu4_radius = 1206. * mm /2.;
		constexpr G4double cu4_height = 2706. * mm /2.; 

		constexpr G4double cu3thick   = 6. * mm; 
		constexpr G4double cu3thickB  = 27.* mm;
		constexpr G4double cu3_radius = 1142. * mm / 2.;
		constexpr G4double cu3_height = 2449. * mm/ 2.; 

		constexpr G4double cu2thick  = 3 * mm; 
		constexpr G4double cu2thickB = 4 * mm;
		constexpr G4double cu2_radius = 1086. * mm / 2.; 
		constexpr G4double cu2_height = 2163. * mm / 2.;

		constexpr G4double cu1thick  = 3 * mm;
		constexpr G4double cu1thickB = 6 * mm;
		constexpr G4double cu1_radius = 1026. * mm / 2.; 
		constexpr G4double cu1_height = 1907. * mm / 2.; 

		constexpr G4double cumcp_radius = 1000. * mm / 2.;
		constexpr G4double inlead_radius = 980. * mm / 2.;

		constexpr G4double cu4to3diff = 124 * mm;
		constexpr G4double cu3to2diff = 127.5 * mm;
		constexpr G4double cu2to1diff = 124 * mm;

		// Cu, Pb Plates
		constexpr G4double plateGap     = 23. * cm; //12. * cm;
		constexpr G4double plateGap2    = 101.892 * mm ; //10. * cm; //5.0 * cm;
		constexpr G4double cumcp_height = 27. * mm / 2; //3. * cm / 2;
		constexpr G4double cup1_zsize   = 30. * mm / 2.; //2. * cm /2; //0.5 * cm / 2;
		constexpr G4double cup2_zsize   = 3. * mm / 2.;
		constexpr G4double cup3_zsize   = 25. * mm / 2.; //2. *cm /2; // 1. * cm / 2;
		constexpr G4double pbp1_zsize   = 8. * cm / 2;
		constexpr G4double pbp2_zsize   = 1. * cm / 2;

	} // namespace AMoRE_200
} // namespace AmoreDetectorStaticInfo

#endif
