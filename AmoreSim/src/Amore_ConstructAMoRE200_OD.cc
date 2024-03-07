//
//  Original by G. Horton-Smith 2004/12/02
//
//  Modified by E.J.Jeon 2007/06/14
//  Modified by Y.S.Yoon 2015/06/15
//  Updated by J.Seo 2024/02/15

#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreDetectorStaticInfo.hh"
#include "AmoreSim/AmoreEventAction.hh"
#include "CupSim/CupPMTSD.hh" // for "sensitive detector"
#include "CupSim/CupParam.hh"
#include "CupSim/CupScintSD.hh"
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4RotationMatrix.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include <fstream>
#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////
// declaration of "private" static utility functions that we
// don't need in class definition
static G4LogicalVolume *MakeHBeam1(G4Material *mat);
static G4LogicalVolume *MakePit(G4double pitBox_x, G4double pitBox_z, G4Material *mat, G4double HatBarrel_gap);

////////////////////////////////////////////////////////////////
/** ConstructXeDetector_OD() constructs the GenericLAND inner detector.
	It is a member of the AmoreDetectorConstruction class, so it has
	direct access to all the member variables in that class, including
	saved material pointers and such.
	*/
G4LogicalVolume *AmoreDetectorConstruction::ConstructAMoRE200_OD() {
	using namespace CLHEP;
	using namespace AmoreDetectorStaticInfo;
	using namespace AmoreDetectorStaticInfo::ColorTable;
	using namespace AmoreDetectorStaticInfo::AMoRE_200;

	// -- database
	CupParam &db(CupParam::GetDB());

	// -- Visualization Attributes
	G4VisAttributes *visRock = new G4VisAttributes(G4Colour(0.0, 0.8, 0.8, 0.3));
	G4VisAttributes *visCavern = new G4VisAttributes(orangel);
	G4VisAttributes *visAir = new G4VisAttributes(G4Colour(1,1,1,0.4));
	G4VisAttributes *leadVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.7, 0.3));
	G4VisAttributes *stainlessVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.9));
	G4VisAttributes *aluminiumVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.1));
	G4VisAttributes *ironVis = new G4VisAttributes(G4Colour(0,0,0,0.1));
	G4VisAttributes *CuShieldVisAttr = new G4VisAttributes(G4Colour(23/255., 189/255., 79/255., 0.8));
	G4VisAttributes *boricAcidVisAttr = new G4VisAttributes(G4Colour(0.4,0.4,0.4,0.2));
	G4VisAttributes *plasticScintVis = new G4VisAttributes(G4Colour(0.8, 1.0, 0.725, 0.4));
	G4VisAttributes *plasticVetoVis = new G4VisAttributes(G4Colour(1.0, 0.855, 0.725, 0.4));
	G4VisAttributes *shieldPE_VisAttr = new G4VisAttributes(G4Colour(129 / 255., 193 / 255., 71 / 255., 0.9));
	G4VisAttributes *shieldWaterTank_VisAttr = new G4VisAttributes(G4Colour(0, 0 , 1, 0.1));

	visAir->SetForceSolid(true);
	visRock->SetForceSolid(true);
	visCavern->SetForceSolid(true);
	leadVis->SetForceSolid(true);
	stainlessVis->SetForceSolid(true);
	aluminiumVis->SetForceSolid(true);
	CuShieldVisAttr->SetForceSolid(true);
	boricAcidVisAttr->SetForceSolid(true);
	plasticScintVis->SetForceSolid(true);
	plasticVetoVis->SetForceSolid(true);
	shieldPE_VisAttr->SetForceSolid(true);
	shieldWaterTank_VisAttr->SetForceSolid(true);


	////////////////////////////////////////////////////
	// Primitive values retrived from database
	////////////////////////////////////////////////////
	G4double rockshell_thickness       = db["rockshell_thickness"];
	G4double airbuffer_radius          = db["airbuffer_radius"];
	G4double airbuffer_height          = db["airbuffer_height"];

	G4double lead_shield_thickness     = db["lead_shield_thickness"];
	G4double lead_housing_thickness    = db["lead_housing_thickness"];

	G4double boricacid_thickness       = db["boricacid_thickness"];

	G4double plastic_scintillator_thickness = db["plastic_scintillator_thickness"];
	G4double plastic_veto_thickness    = db["plastic_veto_thickness"];
	G4double plastic_veto_width        = db["plastic_veto_width"];
	G4double al_plate_thickness        = db["al_plate_thickness"];

	G4double PE_shield_thickness       = db["PE_shield_thickness"];
	G4double thin_lead_shield_thickness= db["thin_lead_shield_thickness"];

	G4double nShield_hatGapFromLead    = db["nShield_hatGapFromLead"];
	G4double nShield_GapFromCeiling    = db["nShield_GapFromCeiling"];

	G4double cavern_sphere_radius      = db["cavern_sphere_radius"];
	G4double rock_shell_radius         = db["rock_shell_radius"];

	G4LogicalVolume *retvalLV = nullptr;

	////////////////////////////////////////////////////
	// Composite values calculated from primitive values (Temporary)
	////////////////////////////////////////////////////

	// Barrel shield
	const int nPStype = 2;
	G4double lead_shield_height    =  // Lead shield height
		airbuffer_height + thin_lead_shield_thickness + lead_shield_thickness; 
	G4double lead_shield_halfsize    =  // Lead shield  half x,y
		airbuffer_radius + thin_lead_shield_thickness + lead_shield_thickness; 

	G4double PS_housing_height =  // Barrel size
		lead_shield_height + lead_housing_thickness*2 + 
		boricacid_thickness + PE_shield_thickness;

	G4double PS_housing_halfsize =  // Barrel size
		longPSlength + profile_thickness*2.5 + plastic_veto_thickness;

	// Rock
	G4double rock_floor_thickness = pitBox_z/2.;

	// World
	G4double floor_spacing = rock_shell_radius;
	G4double world_size    = rock_shell_radius * sqrt(2);
	G4double world_size_Z  = rock_shell_radius * sqrt(2) + rock_floor_thickness + floor_spacing;


	////////////////////////////////////////////////////
	// Make the world
	////////////////////////////////////////////////////
	G4Box *worldSolid = new G4Box("World_solid", world_size, world_size, world_size_Z / 2.);
	G4LogicalVolume *logiWorld  = new G4LogicalVolume(worldSolid, _air, "logiWorld");
	world_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logiWorld, 
			"physWorld", nullptr, false, OverlapCheck);

	////////////////////////////////////////////////////
	// Build the rock geometry (Hemisphere and Floor)
	////////////////////////////////////////////////////
	G4Tubs *floorSolid = new G4Tubs("Rock_Floor", 0, rock_shell_radius, rock_floor_thickness, 0, 360 * deg);
	G4LogicalVolume *logiFloor  = new G4LogicalVolume(floorSolid, _rock, "logiFloor");

	G4Sphere *solidRock  = new G4Sphere("Rock_Solid", 0, rock_shell_radius, 0, 360 * deg, 0, 90 * deg);
	G4LogicalVolume *logiRock = new G4LogicalVolume(solidRock, _rock, "logiRock", 0, 0, 0);
	//G4LogicalVolume *logiVirtualRock   = new G4LogicalVolume(solidRock, _rock, "logiVirtualRock", 0, 0, 0);
	//G4Sphere *solidVirtualRock = new G4Sphere("VirtualRock_Solid", 0, 15 * m, 0, 360 * deg, 0, 90 * deg);
	//G4LogicalVolume *logiRock = new G4LogicalVolume(solidVirtualRock, _rock, "logiRock", 0, 0, 0);

	G4Box *RockBox       = new G4Box("RockBox", 
			PS_housing_halfsize - plastic_veto_thickness - profile_thickness + rockshell_thickness,
			PS_housing_halfsize - plastic_veto_thickness - profile_thickness + rockshell_thickness,
			PS_housing_height/2. + rockshell_thickness + ovc_gap/2. + sst_zsize_half);
	G4Box *TopRockBox    = new G4Box("TopRockBox",
			PS_housing_halfsize - plastic_veto_thickness - profile_thickness,
			PS_housing_halfsize - plastic_veto_thickness - profile_thickness,
			rockshell_thickness/2.);
	G4LogicalVolume *logiRockShell = new G4LogicalVolume(RockBox, _rock, "logiRockShell");
	G4LogicalVolume *logiTopRock   = new G4LogicalVolume(TopRockBox, _rock1, "logiTopRock");

	logiRock  -> SetVisAttributes(visRock);
	logiFloor -> SetVisAttributes(visRock);
	logiRockShell -> SetVisAttributes(visRock);
	logiTopRock   -> SetVisAttributes(visRock);
	//logiTopRock -> SetVisAttributes(G4VisAttributes::Invisible);

	// Virtural volume for rock gamma simulation
	G4Box *shieldHousingBox = new G4Box("shieldHousingBox",
			PS_housing_halfsize - plastic_veto_thickness - profile_thickness,
			PS_housing_halfsize - plastic_veto_thickness - profile_thickness,
			PS_housing_height/2.+ sst_zsize_half + ovc_gap/2.);
	G4LogicalVolume *shieldHousingLV = new G4LogicalVolume(shieldHousingBox, _air, "shieldHousingLV");
	shieldHousingLV->SetVisAttributes(visCavern);

	////////////////////////////////////////////
	// Make Shield Solid and Logical volumes 
	////////////////////////////////////////////
	G4Box *IDspaceBox   = new G4Box("IDspace_Box", 
			airbuffer_radius - boricacid_thickness, 
			airbuffer_radius - boricacid_thickness, 
			airbuffer_height / 2. - boricacid_thickness / 2.);

	// PE shield ----------
	G4Box *shieldPEBox   = new G4Box("IPEShield_Box", 
			lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
			lead_shield_halfsize + lead_housing_thickness + boricacid_thickness + PE_shield_thickness,
			PS_housing_height/2.);
	G4VSolid *shieldPESolid = new G4SubtractionSolid("IPEShield_Solid", shieldPEBox, IDspaceBox, 0,
			G4ThreeVector(0,0,shieldPEBox->GetZHalfLength() - IDspaceBox->GetZHalfLength()) );

	G4LogicalVolume *shieldPELV = new G4LogicalVolume(shieldPESolid, _polyethylene, "IPEShield_LV");
	shieldPELV -> SetVisAttributes(shieldPE_VisAttr);

	// Boric Acid shield -------
	G4Box *shieldBoricAcidBox = new G4Box("shieldBoricAcid_Box",
			lead_shield_halfsize + lead_housing_thickness + boricacid_thickness,
			lead_shield_halfsize + lead_housing_thickness + boricacid_thickness,
			(lead_shield_height + boricacid_thickness) / 2. + lead_housing_thickness);
	G4VSolid *shieldBoricAcidSolid = new G4SubtractionSolid("shieldBoricAcid_Solid", shieldBoricAcidBox, IDspaceBox, 0,
			G4ThreeVector(0, 0, shieldBoricAcidBox->GetZHalfLength() - IDspaceBox->GetZHalfLength()));

	G4LogicalVolume *shieldBoricAcidLV = new G4LogicalVolume(shieldBoricAcidSolid, _BoricAcidRubber, "shieldBoricAcid_LV");
	shieldBoricAcidLV -> SetVisAttributes(boricAcidVisAttr);

	// Lead shield ----------
	G4Box *leadShieldBox   = new G4Box("OuterVetoShield_Box", 
			lead_shield_halfsize, lead_shield_halfsize, lead_shield_height / 2.);
	G4VSolid *leadShieldSolid = new G4SubtractionSolid("LeadShield_Solid", leadShieldBox, IDspaceBox, 0,
			G4ThreeVector(0,0,leadShieldBox->GetZHalfLength() + lead_housing_thickness - IDspaceBox->GetZHalfLength()) );

	G4LogicalVolume *leadShieldLV = new G4LogicalVolume(leadShieldSolid, _lead, "OuterVetoShield_LV");
	leadShieldLV -> SetVisAttributes(leadVis);

	G4Box *leadHousingBox   = new G4Box("OuterVetoHousing_Box",
			lead_shield_halfsize + lead_housing_thickness,
			lead_shield_halfsize + lead_housing_thickness,
			lead_shield_height/2. + lead_housing_thickness);
	G4VSolid *leadHousingSolid = new G4SubtractionSolid("LeadHousing_Solid", leadHousingBox, IDspaceBox, 0,
			G4ThreeVector(0,0,leadHousingBox->GetZHalfLength() - IDspaceBox->GetZHalfLength()) );

	G4LogicalVolume *leadHousingLV = new G4LogicalVolume(leadHousingSolid, _stainless, "OuterVetoHousing_LV");
	leadHousingLV -> SetVisAttributes(stainlessVis);

	// Thin lead shield ------- 
	// Has been dicided to remove copper shield. Aug.2022.
	// Lead shield divided two part. 
	// (For low 214Bi contaminated lead(20 cm) and low 210Pb contaimnated lead(5 cm).) 
	G4Box *ThinLeadShieldBox = new G4Box("ThinLead_Box", 
			airbuffer_radius+thin_lead_shield_thickness, 
			airbuffer_radius+thin_lead_shield_thickness, 
			(airbuffer_height+thin_lead_shield_thickness)/2.);
	G4VSolid *ThinLeadShieldSolid = new G4SubtractionSolid("ThinLead_Solid", 
			ThinLeadShieldBox, IDspaceBox, 0,
			G4ThreeVector(0,0,
				ThinLeadShieldBox->GetZHalfLength() + lead_housing_thickness - IDspaceBox->GetZHalfLength()));
	G4LogicalVolume *ThinLeadShieldLV = new G4LogicalVolume(ThinLeadShieldSolid, _lead, "ThinLeadShield_LV");
	ThinLeadShieldLV -> SetVisAttributes(CuShieldVisAttr);

	// Boric Acid ----------
	G4Box *boricAcidBox   = new G4Box("BoricAcid_Box", 
			airbuffer_radius, airbuffer_radius, airbuffer_height / 2.);
	G4VSolid *boricAcidSolid = new G4SubtractionSolid("BoricAcid_Solid", boricAcidBox, IDspaceBox, 0,
			G4ThreeVector(0,0,
				boricAcidBox->GetZHalfLength() + lead_housing_thickness - IDspaceBox->GetZHalfLength()) );

	G4LogicalVolume *boricAcidLV = new G4LogicalVolume(boricAcidSolid, _BoricAcidRubber, "BoricAcid_LV");
	boricAcidLV -> SetVisAttributes(boricAcidVisAttr);

	////////////////////////////////////////////
	// Muon Veto (Plastic scintillator)
	////////////////////////////////////////////
	// veto housing
	G4Box *plasticVetoHousing1Box = new G4Box("PlasticVetoHousing1_Box",
			PS_housing_halfsize, PS_housing_halfsize, PS_housing_height/2.);
	G4VSolid *plasticVetoHousingSolid = new G4SubtractionSolid("PlasticVetoHousing_Solid", plasticVetoHousing1Box, IDspaceBox, 0,
			G4ThreeVector(0, 0, plasticVetoHousing1Box->GetZHalfLength() - IDspaceBox->GetZHalfLength()));
	G4LogicalVolume *plasticVetoHousing1LV = new G4LogicalVolume(plasticVetoHousingSolid, _air, "PlasticVetoHousing1_LV");
	plasticVetoHousing1LV->SetVisAttributes(G4VisAttributes::Invisible);

	G4Box *plasticVetoHousing2Box = new G4Box("PlasticVetoHousing2_Box",
			bottom_veto_housingX/2., bottom_veto_housingY/2., plastic_veto_thickness/2.);
	G4LogicalVolume *plasticVetoHousing2LV = new G4LogicalVolume(plasticVetoHousing2Box, _air, "PlasticVetoHousing2_LV");
	//plasticVetoHousing2LV->SetVisAttributes(G4VisAttributes::Invisible);


	////////////////////////////////////////////
	// Muon Veto (Water Cerenkov)
	////////////////////////////////////////////
	// HAT Air ----------
	G4Box *shieldHatAirBox = new G4Box("HatAir_Box", HatInnerX, HatInnerY, HatInnerZ/2.);
	G4Box *shieldHatAirBox1 = new G4Box("HatAir_Box1", 
			HatInnerX+waterhousing_thickness, HatInnerY+waterhousing_thickness, HatInnerZ/2.);

	// Hat Water Cerenkov Housing  -------------------
	G4Box *shieldWaterHousingBox = new G4Box("WaterHousing_Box",
			HatInnerX + waterhousing_thickness*2 + watertank_thickness,
			HatInnerY + waterhousing_thickness*2 + watertank_thickness,
			(HatInnerZ + waterhousing_thickness*2 + watertank_top_thickness+PMTroom_thickness)/2.);
	G4VSolid *HatWaterHousingSolid = new G4SubtractionSolid("HatWaterHousing_Solid",shieldWaterHousingBox, shieldHatAirBox,0,
			G4ThreeVector(0,0, -shieldWaterHousingBox->GetZHalfLength() + shieldHatAirBox->GetZHalfLength() ));
	G4LogicalVolume *shieldWaterHousingLV = new G4LogicalVolume(HatWaterHousingSolid,_stainless, "HatWaterHousing_LV");
	shieldWaterHousingLV -> SetVisAttributes(stainlessVis);

	// Air in Water Cerenkov Housing ----------------
	G4Box *shieldWaterAirBox = new G4Box("WCAir_Box",
			HatInnerX + waterhousing_thickness + watertank_thickness,
			HatInnerY + waterhousing_thickness + watertank_thickness,
			(HatInnerZ + watertank_top_thickness + PMTroom_thickness)/2.);
	G4VSolid *HatWaterAirSolid = new G4SubtractionSolid("WCAir_Solid", shieldWaterAirBox, shieldHatAirBox1,0,
			G4ThreeVector(0,0, -shieldWaterAirBox->GetZHalfLength() + shieldHatAirBox1->GetZHalfLength() ));
	G4LogicalVolume *shieldWaterTankAirLV = new G4LogicalVolume(HatWaterAirSolid, _air, "HatWaterTankAir_LV");
	shieldWaterTankAirLV->SetVisAttributes(G4VisAttributes::Invisible);

	// Hat Water Tank --------------------------------
	G4Box *shieldWaterTankBox = new G4Box("WaterTank_Box",
			HatInnerX + waterhousing_thickness + watertank_thickness,
			HatInnerY + waterhousing_thickness + watertank_thickness,
			(HatInnerZ + watertank_top_thickness)/2.);
	G4VSolid *HatWaterTankSolid = new G4SubtractionSolid("HatWaterTank_Solid", shieldWaterTankBox, shieldHatAirBox1,0,
			G4ThreeVector(0,0, -shieldWaterTankBox->GetZHalfLength() + shieldHatAirBox1->GetZHalfLength() ));
	G4LogicalVolume *shieldWaterTankLV = new G4LogicalVolume(HatWaterTankSolid, _water, "HatWaterTank_LV");
	shieldWaterTankLV -> SetVisAttributes(shieldWaterTank_VisAttr);

	// PMT ROOM -------------------------------------
	/*
	G4Box *shieldPMTroomBox = new G4Box("PMTroom_Box",
			HatInnerX + waterhousing_thickness + watertank_thickness,
			HatInnerY + waterhousing_thickness + watertank_thickness,
			PMTroom_thickness/2);
	G4LogicalVolume *shieldPMTroomLV = new G4LogicalVolume(shieldPMTroomBox,_air,"HatPMTroom_LV");
	*/

	// HAT H-beam ----------------
	G4Box *HatBeamHousingInBox = new G4Box("HatBeamHousing_Box",
			HatInnerX - HatHBeam_size, HatInnerY - HatHBeam_size, (HatInnerZ - HatHBeam_size)/2.);
	G4VSolid *HatBeamHousingSolid = new G4SubtractionSolid("HatBeamHousingSolid",
			shieldHatAirBox, HatBeamHousingInBox,0, G4ThreeVector(0,0,-HatHBeam_size/2.) );
	G4LogicalVolume *HatBeamHousingLV = new G4LogicalVolume(HatBeamHousingSolid, _air, "HatBeamHousing_LV");
	HatBeamHousingLV->SetVisAttributes(visAir);

	G4Box *HatBeamLong1Box = new G4Box("HatBeamLong1_Box",HatHBeam_size/2., HatHBeam_size/2., HatHBeam_heightL/2.);
	G4Box *HatBeamLong2Box = new G4Box("HatBeamLong2_Box",HatHBeam_size/2., HatHBeam_size/2., 
			HatHBeam_heightL/2.-HatHBeam_size);
	G4Box *HatBeamShortBox = new G4Box("HatBeamShort_Box",HatHBeam_size/2., HatHBeam_size/2., HatHBeam_heightS/2.);
	G4Box *HatBeamSpaceBox = new G4Box("HatBeamSpace_Box",
			HatHBeam_size/2.-HatHBeam_thickness, HatHBeam_size/2., HatHBeam_heightL);
	G4VSolid *HatBeamLong1Solid = new G4SubtractionSolid("HatBeamLong1_Solid0",
			HatBeamLong1Box, HatBeamSpaceBox, 0, {0, -HatHBeam_size/2.-HatHBeam_in_thickness/2.,0});
	HatBeamLong1Solid = new G4SubtractionSolid("HatBeamLong1_Solid",
			HatBeamLong1Solid, HatBeamSpaceBox, 0, {0, HatHBeam_size/2.+HatHBeam_in_thickness/2.,0});
	G4VSolid *HatBeamLong2Solid = new G4SubtractionSolid("HatBeamLong2_Solid0",
			HatBeamLong2Box, HatBeamSpaceBox, 0, {0, -HatHBeam_size/2.-HatHBeam_in_thickness/2.,0});
	HatBeamLong2Solid = new G4SubtractionSolid("HatBeamLong2_Solid",
			HatBeamLong2Solid, HatBeamSpaceBox, 0, {0, HatHBeam_size/2.+HatHBeam_in_thickness/2.,0});
	G4VSolid *HatBeamShortSolid = new G4SubtractionSolid("HatBeamShort_Solid0",
			HatBeamShortBox, HatBeamSpaceBox, 0, {0, -HatHBeam_size/2.-HatHBeam_in_thickness/2.,0});
	HatBeamShortSolid = new G4SubtractionSolid("HatBeamShort_Solid",
			HatBeamShortSolid, HatBeamSpaceBox, 0, {0, HatHBeam_size/2.+HatHBeam_in_thickness/2., 0});
	G4LogicalVolume *HatBeamLong1LV = new G4LogicalVolume(HatBeamLong1Solid, _iron2, "HatBeamLong1_LV");
	G4LogicalVolume *HatBeamLong2LV = new G4LogicalVolume(HatBeamLong2Solid, _iron2, "HatBeamLong2_LV");
	G4LogicalVolume *HatBeamShortLV = new G4LogicalVolume(HatBeamShortSolid, _iron2, "HatBeamShort_LV");
	HatBeamLong1LV->SetVisAttributes(ironVis);
	HatBeamLong2LV->SetVisAttributes(ironVis);
	HatBeamShortLV->SetVisAttributes(ironVis);

	// HAT Aluminium plate -------
	G4Box *HatAlPlateInBox = new G4Box("HatAlPlateIn_Box",
			HatInnerX - HatHBeam_size - HatAlPlate_thickness, 
			HatInnerY - HatHBeam_size - HatAlPlate_thickness,
			(HatInnerZ - HatHBeam_size - HatAlPlate_thickness)/2.);
	G4VSolid *HatAlPlateSolid = new G4SubtractionSolid("HatAlPlateSolid",
			HatBeamHousingInBox, HatAlPlateInBox, 0, G4ThreeVector(0,0,-HatAlPlate_thickness/2.) );

	G4LogicalVolume *HatAlPlateLV = new G4LogicalVolume(HatAlPlateSolid, _aluminium, "HatAlPlate_LV");
	HatAlPlateLV->SetVisAttributes(aluminiumVis);

	// HAT Boric Acid ------------ 
	G4Box *shieldHatSpaceBox = new G4Box("HatBoric_Box",
			HatInnerX - HatHBeam_size - HatAlPlate_thickness - boricacid_thickness,
			HatInnerY - HatHBeam_size - HatAlPlate_thickness - boricacid_thickness,
			(HatInnerZ - HatHBeam_size - HatAlPlate_thickness - boricacid_thickness)/2.);
	G4VSolid *HatBoricAcidSolid = new G4SubtractionSolid("HatBoricAcid_Solid", HatAlPlateInBox, shieldHatSpaceBox, 0,
			G4ThreeVector(0,0,-boricacid_thickness/2.));

	G4LogicalVolume *shieldHatBoricLV = new G4LogicalVolume(HatBoricAcidSolid, _BoricAcidRubber, "HatBoric_LV");
	shieldHatBoricLV -> SetVisAttributes(boricAcidVisAttr);

	// Detector supporting H-beam  
	G4Box *DetHbeamHousingBox = new G4Box("DetHbeamHousingBox",
			shieldHatSpaceBox->GetXHalfLength(), shieldHatSpaceBox->GetYHalfLength(), shieldHatSpaceBox->GetZHalfLength()-DetHbeam_size/2.);
	G4LogicalVolume *DetHbeamHousingLV = new G4LogicalVolume(DetHbeamHousingBox, _air, "DetHbeamHousing_LV");
	DetHbeamHousingLV->SetVisAttributes(G4VisAttributes::Invisible);

	G4Box *DetHbeamHBox = new G4Box("DetHbeamH_Box", DetHbeam_size/2., DetHbeam_size/2., shieldHatSpaceBox->GetYHalfLength());
	G4Box *DetHbeamVBox = new G4Box("DetHbeamV_Box", DetHbeam_size/2., DetHbeam_size/2., PS_housing_height/2.);
	G4Box *DetHbeamBox = new G4Box("DetHbeam_Box", DetHbeam_size/2., DetHbeam_size/2., DetHbeam_height/2.);
	G4Box *DetHbeamSBox = new G4Box("DetHbeamS_Box", DetHbeam_size/2., DetHbeam_size/2., DetHbeamS_height/2.);
	G4Box *DetHbeamMBox = new G4Box("DetHbeamM_Box", DetHbeam_size/2., DetHbeam_size/2., DetHbeamS_height/4.);
	G4Box *DetHbeamSpaceBox = new G4Box("DetHbeamSpace_Box", 
			DetHbeam_size/2.-DetHbeam_thickness/2., DetHbeam_size/2.-DetHbeam_thickness,shieldHatSpaceBox->GetYHalfLength()*2);

	G4VSolid *DetHbeamSolid = new G4SubtractionSolid("DetHbeam1_Solid", DetHbeamBox, DetHbeamSpaceBox, 0,
			G4ThreeVector(-DetHbeam_size/2.,0,0));
	DetHbeamSolid = new G4SubtractionSolid("DetHbeam_Solid", DetHbeamSolid, DetHbeamSpaceBox, 0,
			G4ThreeVector(DetHbeam_size/2.,0,0));
	G4LogicalVolume *DetHbeamLV = new G4LogicalVolume(DetHbeamSolid, _iron2, "DetHbeam_LV");
	DetHbeamLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamSSolid = new G4SubtractionSolid("DetHbeam2_Solid", DetHbeamSBox, DetHbeamSpaceBox, 0,
			G4ThreeVector(-DetHbeam_size/2.,0,0));
	DetHbeamSSolid = new G4SubtractionSolid("DetHbeamS_Solid", DetHbeamSSolid, DetHbeamSpaceBox, 0,
			G4ThreeVector(DetHbeam_size/2.,0,0));
	G4LogicalVolume *DetHbeamSLV = new G4LogicalVolume(DetHbeamSSolid, _iron2, "DetHbeamS_LV");
	DetHbeamSLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamMSolid = new G4SubtractionSolid("DetHbeam3_Solid", DetHbeamMBox, DetHbeamSpaceBox, 0,
			G4ThreeVector(-DetHbeam_size/2.,0,0));
	DetHbeamMSolid = new G4SubtractionSolid("DetHbeamM_Solid", DetHbeamMSolid, DetHbeamSpaceBox, 0,
			G4ThreeVector(DetHbeam_size/2.,0,0));
	G4LogicalVolume *DetHbeamMLV = new G4LogicalVolume(DetHbeamMSolid, _iron2, "DetHbeamM_LV");
	DetHbeamMLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamHSolid = new G4SubtractionSolid("DetHbeam4_Solid", DetHbeamHBox, DetHbeamSpaceBox, 0,
			G4ThreeVector(-DetHbeam_size/2.,0,0));
	DetHbeamHSolid = new G4SubtractionSolid("DetHbeamH_Solid", DetHbeamHSolid, DetHbeamSpaceBox, 0,
			G4ThreeVector(DetHbeam_size/2.,0,0));
	G4LogicalVolume *DetHbeamHLV = new G4LogicalVolume(DetHbeamHSolid, _iron2, "DetHbeamH_LV");
	DetHbeamHLV->SetVisAttributes(ironVis);

	G4VSolid *DetHbeamVSolid = new G4SubtractionSolid("DetHbeam5_Solid", DetHbeamVBox, DetHbeamSpaceBox, 0,
			G4ThreeVector(-DetHbeam_size/2.,0,0));
	DetHbeamVSolid = new G4SubtractionSolid("DetHbeamV_Solid", DetHbeamVSolid, DetHbeamSpaceBox, 0,
			G4ThreeVector(DetHbeam_size/2.,0,0));
	G4LogicalVolume *DetHbeamVLV = new G4LogicalVolume(DetHbeamVSolid, _iron2, "DetHbeamV_LV");
	DetHbeamVLV->SetVisAttributes(ironVis);

	//////////////////////////////////////////
	// Cavern and Rock geometry options
	//////////////////////////////////////////
	G4VSolid *cavern_solid = nullptr;
	G4LogicalVolume *logiCavern = nullptr;
	G4VPhysicalVolume *physCavern = nullptr;
	G4ThreeVector cavern_tlate = G4ThreeVector(0);

	G4Sphere *neutronmodeCavern = nullptr;
	G4Sphere *hemiSphereCavern = nullptr;

	G4Box *cavernLoafBox = nullptr;
	G4Tubs *cavernMainPizza = nullptr;
	G4Tubs *cavernSubPizza = nullptr;

	if (fNeutronMode) { /// for neutron simulation
		neutronmodeCavern = new G4Sphere("cavern_solid", 0, cavern_nMode_radius, 0., 360. * deg, 0, 180. * deg);
		cavern_solid = neutronmodeCavern;
		logiCavern   = new G4LogicalVolume(cavern_solid, _air, "logiCavern");
		logiCavern  -> SetVisAttributes(G4VisAttributes::Invisible);
		physCavern   = new G4PVPlacement(nullptr, {}, logiCavern, "physCavern", logiWorld, false, 0, OverlapCheck);
	} 	else if (fRockgammaMode) { /// for rock gammna simulation
		fRockPhysical  = new G4PVPlacement(NULL, {0,0,0}, logiRockShell, 
				"physRock", logiWorld, false, 0, OverlapCheck);
		new G4PVPlacement(NULL, {0, 0, RockBox->GetZHalfLength() - rockshell_thickness/2.}, logiTopRock,  
				"physTopRock", logiRockShell, false, 0, OverlapCheck);
		new G4PVPlacement(NULL, 
				{0, 0, RockBox->GetZHalfLength() - rockshell_thickness - shieldHousingBox->GetZHalfLength()}, shieldHousingLV,
				"shieldHousingPV", logiRockShell, false, 0, OverlapCheck);
	}	else { /// for muon and internal/external background simulation
		G4ThreeVector rockTlate = G4ThreeVector(0, 0, 0);
		G4ThreeVector rockFloorTlate = G4ThreeVector(0, 0, - rock_floor_thickness);

		//new G4PVPlacement(NULL, rockTlate, logiVirtualRock, "physVirtualRock", logiWorld, false, 0, OverlapCheck);
		fRockPhysical = new G4PVPlacement(NULL, rockTlate, logiRock, 
				"physRock", logiWorld, false, 0, OverlapCheck);
				//"physRock", logiVirtualRock, false, 0, OverlapCheck);
		fFloorPhysical = new G4PVPlacement(nullptr, rockFloorTlate, logiFloor, 
				"physFloor",logiWorld, false, 0, OverlapCheck);
		AmoreEventAction::SetPrimSkew(rockTlate);

		switch (whichCavernType) {
			case kCavern_Toy_Cylinder:{
				G4Exception(__PRETTY_FUNCTION__, "CAVERN", G4ExceptionSeverity::FatalException,
				"Cavern type of cylinder toy model hasn't been implemented in AMoRE-II.");
				return nullptr;
				break;
			}
			case kCavern_Toy_HemiSphere:{
				hemiSphereCavern = new G4Sphere("cavern_solid", 0, 
				cavern_sphere_radius, 0, 360 * deg, 0, 90 * deg);
				cavern_solid     = hemiSphereCavern;
				break;
			}
			case kCavern_RealModel:{
				G4double cavern_loaf_width = (cavern_subpizza_radius +
					(cavern_pizza_radius - cavern_subpizza_radius) * std::sin(cavern_pizza_angle_real / 2.)) * 2.;
				G4double cavern_loaf_totalheight = cavern_loaf_height;

				cavernLoafBox = new G4Box("CavernLoafBox", 
					cavern_loaf_thickness / 2., cavern_loaf_width / 2., cavern_loaf_totalheight / 2.);
				cavernMainPizza = new G4Tubs("CavernMainPizza", 0, 
					cavern_pizza_radius,
					cavern_pizza_thickness / 2., 0,
					cavern_pizza_angle_real + cavern_pizza_angle_tol * 2.);
				cavernSubPizza  = new G4Tubs("CavernSubPizza", 0, 
					cavern_subpizza_radius,
					cavern_pizza_thickness / 2., 0,
					cavern_pizza_angle_real + cavern_pizza_angle_tol * 2.);

				G4RotationMatrix *pizzaAlignRotMtx = new G4RotationMatrix();
				pizzaAlignRotMtx->rotateY(90 * deg);
				pizzaAlignRotMtx->rotateZ(cavern_subpizza_angle / 2. + cavern_pizza_angle_tol);

				G4RotationMatrix *pizzaPart1RotMtx = new G4RotationMatrix();
				pizzaPart1RotMtx->rotateZ(cavern_subpizza_angle + cavern_pizza_angle_tol);

				G4RotationMatrix *pizzaPart2RotMtx = new G4RotationMatrix();
				pizzaPart2RotMtx->rotateZ(-cavern_subpizza_angle - cavern_pizza_angle_tol);

				G4ThreeVector pizzaTlateInLoaf = G4ThreeVector(0, 0, 
					-cavern_loaf_totalheight / 2. - cavern_totalpizza_dist);

				G4UnionSolid *pizzaUnionStage1 = new G4UnionSolid("PizzaUnion1_Solid", 
					cavernMainPizza, cavernSubPizza, pizzaPart1RotMtx,
					G4ThreeVector(cavern_pizza_radius - cavern_subpizza_radius, 0, 0)
						.rotateZ(cavern_pizza_angle_tol));
				G4UnionSolid *pizzaUnionStage2 = new G4UnionSolid("PizzaUnion2_Solid", 
					pizzaUnionStage1, cavernSubPizza, pizzaPart2RotMtx,
					G4ThreeVector(cavern_pizza_radius - cavern_subpizza_radius, 0, 0)
						.rotateZ(cavern_pizza_angle_real + cavern_pizza_angle_tol));
				G4SubtractionSolid *pizzaUnionStage3 = new G4SubtractionSolid("PizzaUnion3_Solid",
					pizzaUnionStage2, cavernLoafBox,pizzaAlignRotMtx,{0,0,0});

				G4UnionSolid *cavernFinalSolid = new G4UnionSolid("CavernFinalSolid", 
					cavernLoafBox, pizzaUnionStage3,pizzaAlignRotMtx, pizzaTlateInLoaf);

				cavern_solid = cavernFinalSolid;

				cavern_tlate = G4ThreeVector(0, 0, cavern_loaf_totalheight / 2.);
				break;
			}
			default:{
				G4Exception(__PRETTY_FUNCTION__, "CAVERN",
				G4ExceptionSeverity::FatalErrorInArgument, "Cavern type is wrong.");
				return nullptr;
			}
		}
		logiCavern = new G4LogicalVolume(cavern_solid, _air, "logiCavern");
		logiCavern->SetVisAttributes(visCavern);
		physCavern = new G4PVPlacement(nullptr, cavern_tlate, 
				logiCavern, "physCavern",logiRock, false, OverlapCheck);
		fCavernPhysical = physCavern;
	}

	//////////////////////////////////////////
	// Detector root shield volumes positioning
	//////////////////////////////////////////

	// Muon veto housing ------------
	G4PVPlacement *vetoHousing1PV = nullptr;
	G4PVPlacement *vetoHousing2PV = nullptr;
	if(fNeutronMode){
		vetoHousing1PV = new G4PVPlacement(nullptr, 
				{0, 0, - plasticVetoHousing1Box->GetZHalfLength() - nShield_GapFromCeiling},
				plasticVetoHousing1LV, "PlasticVetoHousing_PV", logiCavern, false, 0, 0);
		vetoHousing2PV = new G4PVPlacement(nullptr,
				{0, 0, -plasticVetoHousing1Box->GetZHalfLength()*2 - nShield_GapFromCeiling - HBeam_size 
				- plasticVetoHousing2Box->GetZHalfLength()},
				plasticVetoHousing2LV, "PlasticVetoHousing2_PV", logiCavern, false, 0, 0);
	} else if(!fRockgammaMode){
		vetoHousing1PV = new G4PVPlacement(nullptr, 
				{0, 0, -plasticVetoHousing1Box->GetZHalfLength()},
				plasticVetoHousing1LV, "PlasticVetoHousing_PV", logiCavern, false, 0, 0);
		vetoHousing2PV = new G4PVPlacement(nullptr,
				{0, 0, - HBeam_size - plasticVetoHousing2Box->GetZHalfLength()},
				plasticVetoHousing2LV, "PlasticVetoHousing2_PV", logiFloor, false, 0, 0);
	}

	// PE shield --------------------
	G4PVPlacement *shieldPEPV;
	if(fRockgammaMode){
		shieldPEPV = new G4PVPlacement(nullptr, {0, 0, -sst_zsize_half-ovc_gap/2.},
				shieldPELV, "PEShield_PV", shieldHousingLV, false, 0, OverlapCheck);
	} else{
		shieldPEPV = new G4PVPlacement(nullptr, {0, 0, 0},
				shieldPELV, "PEShield_PV", plasticVetoHousing1LV, false, 0, OverlapCheck);
	}

	// BoricAcid Shield -------------
	G4PVPlacement *shieldBoricAcidPV = new G4PVPlacement(nullptr,
			G4ThreeVector(0, 0, PE_shield_thickness / 2.),
			shieldBoricAcidLV, "BoricAcidShield_PV", shieldPELV, false, 0, OverlapCheck);

	// Lead shield -----------------
	f200_VetoMaterialPhysical = new G4PVPlacement(nullptr, 
			{0, 0, boricacid_thickness / 2.},
			leadHousingLV, "LeadShield_PV", shieldBoricAcidLV, false, 0, OverlapCheck);

	G4VPhysicalVolume *BufferMother = new G4PVPlacement(nullptr,	{}, 
			leadShieldLV, "LeadShield_PV", leadHousingLV, false, 0, OverlapCheck);

	// Thin Lead shield -----------
	new G4PVPlacement(nullptr,{0,0,lead_shield_thickness/2.},
			ThinLeadShieldLV, "ThinLeadShield_PV", leadShieldLV, false, 0, OverlapCheck);

	// Boric Acid ----------
	new G4PVPlacement(nullptr,{0, 0, thin_lead_shield_thickness/2.}, 
			boricAcidLV, "InnerBoricAcid_PV", ThinLeadShieldLV, false, 0, OverlapCheck);



	// Upper part positioning --------------------------
	G4PVPlacement *shieldWaterHousingPV = nullptr;
	G4PVPlacement *shieldWaterTankAirPV = nullptr;
	G4PVPlacement *shieldWaterTankPV = nullptr;
	G4PVPlacement *HatAlPlatePV = nullptr;
	G4PVPlacement *shieldHatBoricPV = nullptr;
	G4PVPlacement *HatBeamHousingPV = nullptr;
	G4PVPlacement *DetHbeamHousingPV = nullptr;
	G4PVPlacement *DetHbeamHPV[4] = {nullptr};
	G4PVPlacement *DetHbeamVPV[8] = {nullptr};

	if(!fRockgammaMode) { // Upper part shields.
		// Hat Water cerenkov detector ----------
		shieldWaterHousingPV = new G4PVPlacement( nullptr,
				{0,0,shieldWaterHousingBox->GetZHalfLength() + nShield_hatGapFromLead + solidBooleanTol},
				shieldWaterHousingLV,"HatWaterHousing_PV",logiCavern, false, 0, 0);
		shieldWaterTankAirPV = new G4PVPlacement( nullptr,
				{0,0,0},
				shieldWaterTankAirLV, "HatWaterTankAir_PV", shieldWaterHousingLV, false, 0, 0);
		shieldWaterTankPV = new G4PVPlacement( nullptr,
				{0,0,
				-shieldWaterHousingBox->GetZHalfLength()+waterhousing_thickness+shieldWaterTankBox->GetZHalfLength()},
				shieldWaterTankLV, "HatWaterTank_PV", shieldWaterTankAirLV, false, 0, OverlapCheck);
				/*
		new G4PVPlacement( nullptr,	{0,0,
				shieldWaterHousingBox->GetZHalfLength()-waterhousing_thickness - shieldPMTroomBox->GetZHalfLength()},
				shieldPMTroomLV, "HatPMTroom_PV", shieldWaterHousingLV, false, 0, OverlapCheck);
				*/

		// Hat Aluminium plate -------------------
		HatAlPlatePV = new G4PVPlacement( nullptr,
				{0, 0, HatBeamHousingInBox->GetZHalfLength() + nShield_hatGapFromLead},
				HatAlPlateLV, "HatAlPlate_PV", logiCavern, false, 0, 0);

		// Hat Boric Acid ------------------------
		shieldHatBoricPV = new G4PVPlacement( nullptr,
				//{0,0, shieldHatAirBox->GetZHalfLength() + nShield_hatGapFromLead},
				{0,0, HatAlPlateInBox->GetZHalfLength() + nShield_hatGapFromLead},
				shieldHatBoricLV, "HatBoricAcid_PV", logiCavern, false, 0, 0);

		if(fEnable_Gantry){
			// Hat H-beam ----------------------------
			HatBeamHousingPV = new G4PVPlacement( nullptr, 
					{0,0, shieldHatAirBox->GetZHalfLength() + nShield_hatGapFromLead},
					HatBeamHousingLV, "HatBeamHousing_PV", logiCavern, false, 0, 0);


			double beamdist_y = (HatInnerY*2-HatHBeam_size)/5.;
			double beamdist_x = (HatInnerX*2-HatHBeam_size)/5.;
			for(int ih = 0; ih < 6; ih++){
				new G4PVPlacement( nullptr, {-HatInnerX+HatHBeam_size/2., -HatInnerY+HatHBeam_size/2.+ ih*beamdist_y, 
						HatInnerZ/2.-HatHBeam_size - HatHBeam_heightS/2.},
						HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
				new G4PVPlacement( nullptr, {HatInnerX-HatHBeam_size/2., -HatInnerY+HatHBeam_size/2.+ ih*beamdist_y, 
						HatInnerZ/2.-HatHBeam_size - HatHBeam_heightS/2.},
						HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
				if( 0 <ih && ih< 5){
					new G4PVPlacement( nullptr, {-HatInnerX+HatHBeam_size/2.+ih*beamdist_x, -HatInnerY+HatHBeam_size/2., 
							HatInnerZ/2.-HatHBeam_size - HatHBeam_heightS/2.},
							HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
					new G4PVPlacement( nullptr, {-HatInnerX+HatHBeam_size/2.+ih*beamdist_x, HatInnerY-HatHBeam_size/2., 
							HatInnerZ/2.-HatHBeam_size - HatHBeam_heightS/2.},
							HatBeamShortLV, "HatBeamShort_PV", HatBeamHousingLV, false, 0, OverlapCheck);
				}
			}

			G4RotationMatrix *hatbeamRotMtx = new G4RotationMatrix();
			hatbeamRotMtx->rotateY(90*deg);
			new G4PVPlacement(G4Transform3D(*hatbeamRotMtx,G4ThreeVector(
							0,-HatInnerY+HatHBeam_size/2.,HatInnerZ/2.-HatHBeam_size/2.)),
					HatBeamLong1LV, "HatBeamLong1_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*hatbeamRotMtx,G4ThreeVector(
							0,HatInnerY-HatHBeam_size/2.,HatInnerZ/2.-HatHBeam_size/2.)),
					HatBeamLong1LV, "HatBeamLong1_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			hatbeamRotMtx->rotateZ(90*deg);
			new G4PVPlacement(G4Transform3D(*hatbeamRotMtx,G4ThreeVector(
							-HatInnerX+HatHBeam_size/2.,0,HatInnerZ/2.-HatHBeam_size/2.)),
					HatBeamLong2LV, "HatBeamLong2_PV", HatBeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*hatbeamRotMtx,G4ThreeVector(
							HatInnerX-HatHBeam_size/2.,0,HatInnerZ/2.-HatHBeam_size/2.)),
					HatBeamLong2LV, "HatBeamLong2_PV", HatBeamHousingLV, false, 0, OverlapCheck);

			// Detector supporting H-beam
			DetHbeamHousingPV = new G4PVPlacement( nullptr, 
					{0, 0, shieldHatSpaceBox->GetZHalfLength() + nShield_hatGapFromLead + DetHbeam_size/2.},
					DetHbeamHousingLV, "DetHbeamHousing_PV", logiCavern, false, 0, 0);

			G4RotationMatrix *detbeamRotMtx = new G4RotationMatrix();
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							- DetHbeam_height/2.-DetHbeam_size/2., - DetHbeam_height/2.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							- DetHbeam_height/2.-DetHbeam_size/2., DetHbeam_height/2.-DetHbeam_size/2., 
							- DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							DetHbeam_height/2.+DetHbeam_size/2., - DetHbeam_height/2.+DetHbeam_size/2., 
							- DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							DetHbeam_height/2.+DetHbeam_size/2., DetHbeam_height/2.-DetHbeam_size/2., 
							- DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							- DetHbeamS_height/4.-DetHbeam_size/2., - DetHbeamS_height/4.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height/2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							- DetHbeamS_height/4.-DetHbeam_size/2., DetHbeamS_height/4.-DetHbeam_size/2.,
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height/2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							DetHbeamS_height/4.+DetHbeam_size/2.,	- DetHbeamS_height/4.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height/2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							DetHbeamS_height/4.+DetHbeam_size/2., DetHbeamS_height/4.-DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height/2.)),
					DetHbeamSLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

			DetHbeamVPV[0] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.-DetHbeam_size/2., 
							-DetHbeamHousingBox->GetYHalfLength()+DetHbeam_size/2.,
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[1] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.+DetHbeam_size/2., 
							DetHbeamHousingBox->GetYHalfLength()-DetHbeam_size/2., 
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[2] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							+DetHbeam_height/2.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetYHalfLength()+DetHbeam_size/2.,
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[3] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							+DetHbeam_height/2.-DetHbeam_size/2., 
							DetHbeamHousingBox->GetYHalfLength()-DetHbeam_size/2., 
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[4] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.-DetHbeam_size/2., 
							-DetHbeamHousingBox->GetYHalfLength()+DetHbeam_size/2.-DetHbeamS_height,
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[5] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.+DetHbeam_size/2., 
							DetHbeamHousingBox->GetYHalfLength()-DetHbeam_size/2.+DetHbeamS_height, 
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[6] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							+DetHbeam_height/2.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetYHalfLength()+DetHbeam_size/2.-DetHbeamS_height,
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);
			DetHbeamVPV[7] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							+DetHbeam_height/2.-DetHbeam_size/2., 
							DetHbeamHousingBox->GetYHalfLength()-DetHbeam_size/2.+DetHbeamS_height, 
							-DetHbeamVBox->GetZHalfLength())),
					DetHbeamVLV, "HatBeamV_PV", logiCavern, false, 0, 0);

			detbeamRotMtx->rotateX(90*deg);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.-DetHbeam_size/2., 0, 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.+ DetHbeam_size/2.)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							DetHbeam_height/2.+DetHbeam_size/2., 0, 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.+ DetHbeam_size/2.)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeamS_height/4.-DetHbeam_size/2., 0, 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height-DetHbeam_size/2.)),
					DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							DetHbeamS_height/4.+DetHbeam_size/2., 0, 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height-DetHbeam_size/2.)),
					DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

			DetHbeamHPV[0] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.-DetHbeam_size/2., 0, +nShield_hatGapFromLead + DetHbeam_size/2.)),
					DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, 0);
			DetHbeamHPV[1] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							-DetHbeam_height/2.+DetHbeam_size/2., 0, +nShield_hatGapFromLead + DetHbeam_size/2.)),
					DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, 0);
			DetHbeamHPV[2] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							+DetHbeam_height/2.+DetHbeam_size/2., 0, +nShield_hatGapFromLead + DetHbeam_size/2.)),
					DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, 0);
			DetHbeamHPV[3] = new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							+DetHbeam_height/2.-DetHbeam_size/2., 0, +nShield_hatGapFromLead + DetHbeam_size/2.)),
					DetHbeamHLV, "HatBeamH_PV", logiCavern, false, 0, 0);

			detbeamRotMtx->rotateZ(90*deg);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							0, -DetHbeam_height/2.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.+ DetHbeam_size/2.)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							0, DetHbeam_height/2.-DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.+ DetHbeam_size/2.)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							0, DetHbeamS_height/4.-DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.+ DetHbeam_size/2.)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							0, -DetHbeamS_height/4.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.+ DetHbeam_size/2.)),
					DetHbeamLV, "HatBeamLong_PV", DetHbeamHousingLV, false, 0, OverlapCheck);

			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							0, -DetHbeamS_height/4.+DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height-DetHbeam_size/2.)),
					DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
			new G4PVPlacement(G4Transform3D(*detbeamRotMtx,G4ThreeVector(
							0, DetHbeamS_height/4.-DetHbeam_size/2., 
							-DetHbeamHousingBox->GetZHalfLength()+DetHbeam_height/2.+solidBooleanTol
							+DetHbeam_height/2.-DetHbeamS_height-DetHbeam_size/2.)),
					DetHbeamMLV, "HatBeamShort_PV", DetHbeamHousingLV, false, 0, OverlapCheck);
		}
	}

	//////////////////////////////////////////
	// reposition the chambers to fit the cavern type
	//////////////////////////////////////////
	if(!fNeutronMode && !fRockgammaMode){
		G4double HatBottomPosZ = plasticVetoHousing1Box->GetZHalfLength()*2 + nShield_hatGapFromLead + nShield_GapFromCeiling;

		G4ThreeVector topTlate = G4ThreeVector(0,0,
				HatBottomPosZ + shieldWaterHousingBox->GetZHalfLength());
		G4ThreeVector topBeamTlate = G4ThreeVector(0,0,
				HatBottomPosZ + shieldHatAirBox->GetZHalfLength());
		G4ThreeVector topAlPlateTlate = G4ThreeVector(0,0,
				HatBottomPosZ + HatBeamHousingInBox->GetZHalfLength());
		G4ThreeVector topBoricTlate = G4ThreeVector(0, 0, 
				HatBottomPosZ + HatAlPlateInBox->GetZHalfLength());
		G4ThreeVector detHbeamTlate = G4ThreeVector(0, 0, 
				HatBottomPosZ + shieldHatSpaceBox->GetZHalfLength() + DetHbeam_size/2.);

		G4ThreeVector barrelTlate = G4ThreeVector(0,0, plasticVetoHousing1Box->GetZHalfLength());
		G4ThreeVector bottomTlate = G4ThreeVector(0, 0, rock_floor_thickness - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());

		G4ThreeVector cavernTlate = G4ThreeVector(0);

		if(whichCavernType == kCavern_RealModel){
			HatBottomPosZ = -cavernLoafBox->GetZHalfLength() + plasticVetoHousing1Box->GetZHalfLength()*2
				+ nShield_hatGapFromLead + nShield_GapFromCeiling;

			topTlate = G4ThreeVector(room_dist_x, -room_dist_y,
					HatBottomPosZ+ shieldWaterHousingBox->GetZHalfLength());
			topBeamTlate = G4ThreeVector(room_dist_x, -room_dist_y,
					HatBottomPosZ + shieldHatAirBox->GetZHalfLength());
			topAlPlateTlate = G4ThreeVector(room_dist_x, -room_dist_y,
					HatBottomPosZ + HatBeamHousingInBox->GetZHalfLength());
			topBoricTlate = G4ThreeVector(room_dist_x, -room_dist_y, 
					HatBottomPosZ + HatAlPlateInBox->GetZHalfLength());
			detHbeamTlate = G4ThreeVector(room_dist_x, -room_dist_y,
					HatBottomPosZ + shieldHatSpaceBox->GetZHalfLength() + DetHbeam_size/2.);

			barrelTlate = G4ThreeVector(room_dist_x, -room_dist_y, 
					-cavernLoafBox->GetZHalfLength()+plasticVetoHousing1Box->GetZHalfLength());
			//bottomTlate = G4ThreeVector(room_dist_x, -room_dist_y,
			bottomTlate = G4ThreeVector(0,0,
					rock_floor_thickness - HBeam_size - plasticVetoHousing2Box->GetZHalfLength());

			//cavernTlate = G4ThreeVector(-room_dist_x, room_dist_y,rock_floor_thickness*2);
			cavernTlate = G4ThreeVector(-room_dist_x, room_dist_y, cavernLoafBox->GetZHalfLength());
		}
		vetoHousing1PV->SetTranslation(barrelTlate);
		vetoHousing2PV->SetTranslation(bottomTlate);
		shieldWaterHousingPV->SetTranslation(topTlate);
		HatAlPlatePV->SetTranslation(topAlPlateTlate);
		shieldHatBoricPV->SetTranslation(topBoricTlate);
		physCavern->SetTranslation(cavernTlate);

		if(fEnable_Gantry){
			G4ThreeVector detHbeamHPos[4] = {DetHbeamHPV[0]->GetTranslation(), DetHbeamHPV[1]->GetTranslation(),
				DetHbeamHPV[2]->GetTranslation(), DetHbeamHPV[3]->GetTranslation()};
			G4ThreeVector detHbeamVPos[8] = {DetHbeamVPV[0]->GetTranslation(), DetHbeamVPV[1]->GetTranslation(),
				DetHbeamVPV[2]->GetTranslation(), DetHbeamVPV[3]->GetTranslation(),
				DetHbeamVPV[4]->GetTranslation(), DetHbeamVPV[5]->GetTranslation(),
				DetHbeamVPV[6]->GetTranslation(), DetHbeamVPV[7]->GetTranslation()};
			HatBeamHousingPV->SetTranslation(topBeamTlate);
			DetHbeamHousingPV->SetTranslation(detHbeamTlate);
			for(int iH = 0; iH < 8; iH++){
				if(iH < 4) {
					DetHbeamHPV[iH]->SetTranslation(detHbeamHPos[iH]+detHbeamTlate+G4ThreeVector(
								0,0,-DetHbeamHousingBox->GetZHalfLength()-nShield_hatGapFromLead - DetHbeam_size ));
				}
				DetHbeamVPV[iH]->SetTranslation(detHbeamVPos[iH]+detHbeamTlate+G4ThreeVector(
							0, 0, -DetHbeamHousingBox->GetZHalfLength()-DetHbeam_size-HBeam_size));
				if(OverlapCheck) DetHbeamVPV[iH]->CheckOverlaps();
			}
			if(OverlapCheck){
				HatBeamHousingPV->CheckOverlaps();
				DetHbeamHousingPV->CheckOverlaps();
			}
		}
		if(OverlapCheck){
			vetoHousing1PV->CheckOverlaps();
			vetoHousing2PV->CheckOverlaps();
			shieldWaterHousingPV->CheckOverlaps();
			HatAlPlatePV->CheckOverlaps();
			shieldHatBoricPV->CheckOverlaps();
		}
	}

	//////////////////////////////////////////
	// H-beam and Pit option
	//////////////////////////////////////////
	G4double dist[7] = {2150,380,3800,565,465,1680,1680};
	G4double len[4] = {1192,2208,2218,882};
	// WC supporting H-beam housing----------
	G4Box *HbeamHousingOut = new G4Box("hbeamhousingOut_Box", 
			HBeam_housingX/2., HBeam_housingY/2., (PS_housing_height+HBeam_size+nShield_GapFromCeiling)/2.);
	G4Box *HbeamHousingIn = new G4Box("hbeamhousingIn_Box",
			HBeam_housingX/2. - HBeam_size,
			HBeam_housingY/2. - HBeam_size,
			(PS_housing_height+nShield_GapFromCeiling)/2.);
	//G4Box *HbeamHousingCut = new G4Box("hbeamHousingCut", top_plate_x/2., top_plate_y/2., PS_housing_height);
	G4Box *HbeamHousingCut = new G4Box("hbeamHousingCut", sst_xsize_half, sst_ysize_half, PS_housing_height);
	G4VSolid *HbeamHousing = new G4SubtractionSolid("hbeamhousing_Solid", HbeamHousingOut, HbeamHousingIn, 0, 
			//G4ThreeVector(0,0, - HbeamHousingOut->GetZHalfLength() - HBeam_size));
					 G4ThreeVector(0,0, - HbeamHousingOut->GetZHalfLength() +HbeamHousingIn->GetZHalfLength()));
	HbeamHousing = new G4SubtractionSolid("hbeamhousing_Solid", HbeamHousing, HbeamHousingCut, 0, 
			G4ThreeVector( HBeam_housingDist, 0, 0));
	G4LogicalVolume *HbeamHousingLV = new G4LogicalVolume(HbeamHousing, _air, "HbeamHousing_LV");
	HbeamHousingLV->SetVisAttributes(G4VisAttributes::Invisible);

	G4LogicalVolume *HbeamBotLV = nullptr;
	G4LogicalVolume *PitLV = nullptr;

	if(fEnable_Gantry){
		// H-beam Bottom ----------
		HbeamBotLV    = MakeHBeam1(_iron1);
		HbeamBotLV->SetVisAttributes(ironVis);

		// Pit ----------
		PitLV    = MakePit(pitBox_x, pitBox_z,_rebar, HBeam_size);
		PitLV->SetVisAttributes(ironVis);

		// WC supporting H-beam -------------------
		G4Box *beamCut = new G4Box("beamCut", HBeam_size, HBeam_size/2. - HBeam_thickness, HBeam_housingX);
		G4RotationMatrix *beamRotMtx = new G4RotationMatrix();

		// H-beam vertical
		G4Box *HbeamLong = new G4Box("HbeamLong", HBeam_size/2., HBeam_size/2., PS_housing_height/2.);
		G4VSolid *HbeamV = new G4SubtractionSolid("HbeamV1", HbeamLong, beamCut, 0, 
				{HBeam_size + HBeam_thickness, 0, 0}); 
		HbeamV = new G4SubtractionSolid("HbeamV", HbeamV, beamCut, 0, 
				{-HBeam_size - HBeam_thickness, 0, 0});
		G4LogicalVolume *HbeamV_LV = new G4LogicalVolume(HbeamV, _iron2, "HbeamV_LV");
		HbeamV_LV->SetVisAttributes(ironVis);

		// H-beam horizontal
		beamRotMtx->rotateX(90*deg);
		G4Box *HbeamLong2 = new G4Box("HbeamLong2", HBeam_size/2., HBeam_housingY/2., HBeam_size/2.);
		G4VSolid *HbeamH = new G4SubtractionSolid("HbeamH1", HbeamLong2, beamCut, 
				G4Transform3D(*beamRotMtx, {HBeam_size + HBeam_thickness,0, 0}));
		HbeamH = new G4SubtractionSolid("HbeamH", HbeamH, beamCut, 
				G4Transform3D(*beamRotMtx, {-HBeam_size + -HBeam_thickness,0, 0}));
		G4LogicalVolume *HbeamH_LV = new G4LogicalVolume(HbeamH, _iron2, "HbeamH_LV");
		HbeamH_LV->SetVisAttributes(ironVis);

		// Small H-beams
		G4Box *HbeamShort[7];
		G4VSolid *HbeamSH[7];
		G4LogicalVolume *HbeamSH_LV[7];
		beamRotMtx->rotateZ(90*deg);
		for(int i = 0; i < 7; i++){
			HbeamShort[i] = new G4Box(Form("HbeamShort%d",i), dist[i]/2., HBeam_size/2., HBeam_size/2.);
			HbeamSH[i] = new G4SubtractionSolid(Form("HbeamSH1%d",i), HbeamShort[i], beamCut, 
					G4Transform3D(*beamRotMtx, {0, HBeam_size + HBeam_thickness, 0}));
			HbeamSH[i] = new G4SubtractionSolid(Form("HbeamSH1%d",i), HbeamSH[i], beamCut, 
					G4Transform3D(*beamRotMtx, {0, -HBeam_size - HBeam_thickness, 0}));
			HbeamSH_LV[i] = new G4LogicalVolume(HbeamSH[i], _iron2, Form("HbeamSH%d_LV",i));
			HbeamSH_LV[i]->SetVisAttributes(ironVis);
		}


		////////////////////////////////////////////////////////////////
		//// H-beam, pit positioning 
		////////////////////////////////////////////////////////////////
		// WC supporting H-beam
		G4ThreeVector vbeamPos = G4ThreeVector(HBeam_housingX/2.-HBeam_size/2., HBeam_housingY/2.-HBeam_size/2., -HBeam_size/2.);
		G4ThreeVector hbeamPos = G4ThreeVector(HBeam_housingX/2.-HBeam_size/2., 0, HbeamHousingOut->GetZHalfLength()-HBeam_size/2.);

		new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
		vbeamPos[1] *= -1;
		new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
		new G4PVPlacement(nullptr, hbeamPos, HbeamH_LV, "HbeamH_PV", HbeamHousingLV, false, 0, OverlapCheck); 
		for(int i = 0; i < 7; i++){
			vbeamPos[0] -= dist[i] + HBeam_size;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
			vbeamPos[1] *= -1;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
			hbeamPos[0] -= dist[i] + HBeam_size;
			new G4PVPlacement(nullptr, hbeamPos, HbeamH_LV, "HbeamH_PV", HbeamHousingLV, false, 0, OverlapCheck); 
		}
		for(int i = 0; i < 4; i++){
			vbeamPos[1] -= len[i] + HBeam_size;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
			vbeamPos[0] *= -1;
			new G4PVPlacement(nullptr, vbeamPos, HbeamV_LV, "HbeamV_PV", HbeamHousingLV, false, 0, OverlapCheck);
		}

		G4ThreeVector sbeamPos = G4ThreeVector(HBeam_housingX/2., HBeam_housingY/2. - HBeam_size/2., 
				HbeamHousingOut->GetZHalfLength() - HBeam_size/2.);
		for(int i = 0; i < 7; i++){
			sbeamPos[0] -= dist[i]/2.+HBeam_size;
			new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
			for(int j = 0; j < 4; j++){
				if((i==1 || i == 2 || i==3) && (j==0 || j == 3)){
					sbeamPos[1] -= len[j]/2. + HBeam_size/2.;
					new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
					sbeamPos[1] -= len[j]/2. + HBeam_size/2.;
					new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
				}
				else if (i==2 && j==1){
					sbeamPos[1] -= len[j] + HBeam_size;
				}
				else{
					sbeamPos[1] -= len[j] + HBeam_size;
					new G4PVPlacement(nullptr, sbeamPos, HbeamSH_LV[i], "HbeamSH_PV", HbeamHousingLV, false, 0, OverlapCheck);
				}
			}
			sbeamPos[1] = HBeam_housingY/2. - HBeam_size/2.;
			sbeamPos[0] -= dist[i]/2.;
		}

		//////////////////////////////////////////////////////////////
		G4PVPlacement *HbeamBotPV = nullptr;
		G4PVPlacement *pitPV = nullptr;
		G4PVPlacement *HbeamHousingPV = nullptr;

		G4ThreeVector HbeamTlate = G4ThreeVector( 0, 0, 
				-plasticVetoHousing1Box->GetZHalfLength()*2 - nShield_GapFromCeiling - HBeam_size/2.);
		G4ThreeVector pitTlate = G4ThreeVector(-pitBox_x/2., 0, 
				-plasticVetoHousing1Box->GetZHalfLength()*2 -nShield_GapFromCeiling - pitBox_z/2.);
		G4ThreeVector HbeamHousingTlate = G4ThreeVector(-HBeam_housingDist,0,
				-plasticVetoHousing1Box->GetZHalfLength()*2 - nShield_GapFromCeiling + HbeamHousingOut->GetZHalfLength());

		if(fNeutronMode){
			HbeamBotPV = new G4PVPlacement(nullptr, HbeamTlate, HbeamBotLV,"HbeamBot_PV",logiCavern,false,0, OverlapCheck);
			pitPV = new G4PVPlacement(nullptr, pitTlate, PitLV,"pit_PV",logiCavern,false,0, OverlapCheck);
			HbeamHousingPV = new G4PVPlacement (nullptr, HbeamHousingTlate, HbeamHousingLV, 
					"HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
		}
		else if(!fRockgammaMode){
			if(whichCavernType == kCavern_RealModel){
				//HbeamTlate = G4ThreeVector( room_dist_x, -room_dist_y, rock_floor_thickness-HBeam_size/2.);
				//pitTlate = G4ThreeVector( -pitBox_x/2. + room_dist_x, -room_dist_y, 0);
				HbeamTlate = G4ThreeVector( 0, 0, rock_floor_thickness-HBeam_size/2.);
				pitTlate = G4ThreeVector( -pitBox_x/2., 0, 0);
				HbeamHousingTlate = G4ThreeVector( room_dist_x - HBeam_housingDist, -room_dist_y, 
						-cavernLoafBox->GetZHalfLength() + HbeamHousingOut->GetZHalfLength());
			}
			else if(whichCavernType == kCavern_Toy_HemiSphere){
				HbeamTlate = G4ThreeVector(0,0, rock_floor_thickness- HBeam_size/2.);
				pitTlate = G4ThreeVector(-pitBox_x/2., 0, rock_floor_thickness - pitBox_z/2.);
				HbeamHousingTlate = G4ThreeVector( - HBeam_housingDist, 0, HbeamHousingOut->GetZHalfLength());
			}
			HbeamBotPV = new G4PVPlacement(nullptr, HbeamTlate, HbeamBotLV,"HbeamBot_PV",logiFloor,false,0, OverlapCheck);
			pitPV = new G4PVPlacement(nullptr, pitTlate, PitLV,"pit_PV",logiFloor,false,0, OverlapCheck);
			HbeamHousingPV = new G4PVPlacement (nullptr, HbeamHousingTlate, HbeamHousingLV, 
					"HbeamHousing_PV", logiCavern, false, 0, OverlapCheck);
		}

		if (fDbgMsgOn) {
			cout << " ======================== H-beams =========" << endl;
			cout << "  ( length: mm, weight: kg, global coordinate)" << endl;
			cout << " Bottom H-beam " << endl;
			cout << "     mass      : " << HbeamBotLV->GetMass(true,false)/kg;
			if(fFloorPhysical!=nullptr){
				cout << "     coordinate: " << HbeamBotPV->GetTranslation() + fFloorPhysical->GetTranslation() << endl;
			} else cout << "     coordinate: " << HbeamBotPV->GetTranslation() << endl;
			cout << " Pit " << endl;
			cout << "     mass      : " << PitLV->GetMass(true,false)/kg;
			if(fFloorPhysical!=nullptr){
				cout << "     coordinate: " << pitPV->GetTranslation() + fFloorPhysical->GetTranslation() << endl;
			} else cout << "     coordinate: " << pitPV->GetTranslation() << endl;
			cout << " H-beam housing 1" << endl;
			cout << "     dimension                     : " << HbeamHousingOut->GetXHalfLength()*2 << " x " 
				<< HbeamHousingOut->GetYHalfLength()*2 << " x " 
				<< HbeamHousingOut->GetZHalfLength()*2 << endl;
			cout << "     coordinaate                   : " 
				<< HbeamHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
			cout << "     mass (HbeamV, HbeamH, HbeamSH): " << HbeamV_LV->GetMass(true,false)/kg << ", "
				<< HbeamH_LV->GetMass(true,false)/kg << ", "
				<< HbeamSH_LV[0]->GetMass(true,false)/kg << endl;
			cout << " H-beam housing 2" << endl;
			cout << "     dimension                 : " << shieldHatAirBox->GetXHalfLength()*2 << " x " 
				<< shieldHatAirBox->GetYHalfLength()*2 << " x " 
				<< shieldHatAirBox->GetZHalfLength()*2 << endl;
			cout << "     coordinaate               : " << HatBeamHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
			cout << "     mass (Long1, Long2, Short): " << HatBeamLong1LV->GetMass(true,false)/kg << ", "
				<< HatBeamLong2LV->GetMass(true,false)/kg << ", " 
				<< HatBeamShortLV->GetMass(true,false)/kg << endl;
			cout << " H-beam housing 3" << endl;
			cout << "     dimension                               : " << DetHbeamHousingBox->GetXHalfLength()*2 << " x " 
				<< DetHbeamHousingBox->GetYHalfLength()*2 << " x " 
				<< DetHbeamHousingBox->GetZHalfLength()*2 << endl;
			cout << "     coordinaate                             : " << 
				DetHbeamHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
			cout << "     mass (Long, Short, Horizontal, Vertical): " << DetHbeamLV->GetMass(true,false)/kg << ", "
				<< DetHbeamSLV->GetMass(true,false)/kg << ", " 
				<< DetHbeamHLV->GetMass(true,false)/kg << ", "
				<< DetHbeamVLV->GetMass(true,false)/kg << endl;
		}
	}
	else{ // Without H-beams
		G4PVPlacement *HbeamHousingPV = nullptr;
		if(fRockgammaMode){ HbeamHousingPV = nullptr;	}
		else {
			G4ThreeVector HbeamHousingTlate = G4ThreeVector(-HBeam_housingDist,0,
					-plasticVetoHousing1Box->GetZHalfLength()*2 - nShield_GapFromCeiling + HbeamHousingOut->GetZHalfLength());
			HbeamHousingPV = new G4PVPlacement (nullptr, HbeamHousingTlate, HbeamHousingLV, 
					"HbeamHousing_PV", logiCavern, false, 0, 0);

			if(!fNeutronMode){
				if(whichCavernType == kCavern_RealModel){
					HbeamHousingTlate = G4ThreeVector( room_dist_x - HBeam_housingDist, -room_dist_y, 
							-cavernLoafBox->GetZHalfLength() + HbeamHousingOut->GetZHalfLength());
				}
				else if(whichCavernType == kCavern_Toy_HemiSphere){
					HbeamHousingTlate = G4ThreeVector( - HBeam_housingDist, 0, HbeamHousingOut->GetZHalfLength());
				}
			}
			HbeamHousingPV->SetTranslation(HbeamHousingTlate);
			if(OverlapCheck) HbeamHousingPV->CheckOverlaps();
		}
	}

	//////////////////////////////////////////
	// Additional PE
	//////////////////////////////////////////
	G4Box *addPEBox1 = nullptr;
	G4Box *addPEBox2 = nullptr;
	G4LogicalVolume *addPE1LV = nullptr;
	G4LogicalVolume *addPE2LV = nullptr;

	if(fAdditionalPE)	{
		G4ThreeVector addPEPos = G4ThreeVector( 
				HBeam_housingX/2. - HBeam_size - dist[0]/2.,
				HBeam_housingY/2. - HBeam_size*3/2 - len[0]/2.,
				HbeamHousingOut->GetZHalfLength()-HBeam_size/2.);

		for(int i = 1; i < 5; i ++){
			addPEPos[0] -= dist[i-1]/2. + HBeam_size + dist[i]/2.;
			for(int j = 0; j < 4; j ++){
				if(j == 0 || j == 3){
					addPEBox1 = new G4Box("addPEBox", 
							dist[i]/2.+ HBeam_size/2.-HBeam_thickness,
							len[j]/4. + HBeam_size/4.- HBeam_thickness - 5*solidBooleanTolBig,
							HBeam_size/2.- HBeam_thickness);

					addPE1LV = new G4LogicalVolume(addPEBox1, _polyethylene, "additionalPE_LV" );
					addPE1LV->SetVisAttributes(shieldPE_VisAttr);

					if(j == 3){
						addPEPos[1] -= (len[j] - HBeam_size*3)/4.;
					}
					addPEPos[1] -= (len[j]-HBeam_size)/4.;
					new G4PVPlacement(nullptr, addPEPos, addPE1LV, "additionalPE1_PV", HbeamHousingLV, false, 0, OverlapCheck);
					addPEPos[1] -= (len[j] + HBeam_size*3)/4.;
				}
				else{
					addPEBox2 = new G4Box("addPEBox", 
							dist[i]/2.+ HBeam_size/2.-HBeam_thickness,
							len[j]/2. + HBeam_size/2.- HBeam_thickness - solidBooleanTol,
							HBeam_size/2.- HBeam_thickness);

					addPE2LV = new G4LogicalVolume(addPEBox2, _polyethylene, "additionalPE_LV" );
					addPE2LV->SetVisAttributes(shieldPE_VisAttr);

					addPEPos[1] -= len[j]/2.;
					if( i!=2 && i!=3 ) {
						new G4PVPlacement(nullptr, addPEPos, addPE2LV, "additionalPE2_PV", HbeamHousingLV, false, 0, OverlapCheck);
					}
					addPEPos[1] -= len[j]/2. + HBeam_size;
				}
			}
			addPEPos[1] = HBeam_housingY/2. - HBeam_size*3/2 - len[0]/2.;
		}
	}


	//////////////////////////////////////////////////////
	if (fDbgMsgOn) {
		cout << "\n ================================ Outer shields ======" << endl;
		cout << " (length: mm, weight: kg, global coordinate)\n" << endl;
		cout << " World dimension(box): " 
			<< worldSolid->GetXHalfLength()*2 << " x " 
			<< worldSolid->GetYHalfLength()*2 << " x " 
			<< worldSolid->GetZHalfLength()*2 << endl;
		if(fRockgammaMode){
			cout << " Rock for rockgamma simulation" << endl;
			cout << "     shell mass  : " << logiRockShell->GetMass(true,false)/kg << endl;
			cout << "     top mass    : " << logiTopRock->GetMass(true,false)/kg << endl;
			cout << " PE (box) " << endl; //..............................
			cout << "     mass      : " << shieldPELV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " 
				<< shieldPEBox->GetXHalfLength()*2 << " x " 
				<< shieldPEBox->GetYHalfLength()*2 << " x " 
				<< shieldPEBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< shieldPEBox->GetXHalfLength()-shieldBoricAcidBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< fRockPhysical->GetTranslation() + shieldPEPV->GetTranslation() << endl;
		}
		else if(fNeutronMode){
			cout << " Cavern" << endl;
			cout << "     dimension (r): " << neutronmodeCavern->GetOuterRadius() << endl; 
			cout << "     coordinate   : " << physCavern->GetTranslation() << endl;
		}
		else{
			cout << " Rock floor (tube)" << endl;
			cout << "     mass             : " << logiFloor->GetMass(true,false)/kg << endl;
			cout << "     dimension (r x h): " 
				<< floorSolid->GetOuterRadius() << " x " << floorSolid->GetZHalfLength()*2 << endl;
			cout << "     coordinate       : " << fFloorPhysical->GetTranslation() << endl;
			cout << " Rock body (hemisphere)" << endl;
			cout << "     mass         : " << logiRock->GetMass(true,false)/kg << endl;
			cout << "     dimension (r): " << solidRock->GetOuterRadius() << endl;
			cout << "     coordinate   : " << fRockPhysical->GetTranslation() << endl;
			cout << " Cavern ( '" << AmoreDetectorConstruction::GetCavernTypeName(whichCavernType) 
				<< "' option selected )" << endl;
			if(whichCavernType==kCavern_RealModel){
				cout << "     dimension (box part): " 
					<< cavernLoafBox->GetXHalfLength()*2 << " x " 
					<< cavernLoafBox->GetYHalfLength()*2 << " x " 
					<< cavernLoafBox->GetZHalfLength()*2 << endl;
			}
			else if(whichCavernType==kCavern_Toy_HemiSphere)	cout << "     dimension (r): " << hemiSphereCavern->GetOuterRadius() << endl;
			cout << "     coordinate   : " << fRockPhysical->GetTranslation() + physCavern->GetTranslation() << endl;
		}

		if(!fRockgammaMode){
			/*
			cout << " ============================ Muon Veto ======" << endl;
			cout << " PS veto housing " << endl; //.......................
			cout << "     dimension : " 
				<< plasticVetoHousing1Box->GetXHalfLength()*2 << " x " 
				<< plasticVetoHousing1Box->GetYHalfLength()*2 << " x " 
				<< plasticVetoHousing1Box->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< plasticVetoHousing1Box->GetXHalfLength() - shieldPEBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() << endl;
			cout << " PS supporter (Aluminium profile) mass: \n"
				<< "     " << plasticVetoSupporterV1LV->GetMass(true,false)/kg << "(" 
				<< plasticVetoSupporterV1Box->GetXHalfLength()*2. << "x" 
				<< plasticVetoSupporterV1Box->GetYHalfLength()*2. << "x" 
				<< plasticVetoSupporterV1Box->GetZHalfLength()*2. << ") \n" 
				<< "     " << plasticVetoSupporterV2LV->GetMass(true,false)/kg << "(" 
				<< plasticVetoSupporterV2Box->GetXHalfLength()*2. << "x" 
				<< plasticVetoSupporterV2Box->GetYHalfLength()*2. << "x" 
				<< plasticVetoSupporterV2Box->GetZHalfLength()*2. << ") \n" 
				<< "     " << plasticVetoSupporterH1LV->GetMass(true,false)/kg << "(" 
				<< plasticVetoSupporterH1Box->GetXHalfLength()*2. << "x" 
				<< plasticVetoSupporterH1Box->GetYHalfLength()*2. << "x" 
				<< plasticVetoSupporterH1Box->GetZHalfLength()*2. << ") \n" 
				<< "     " << plasticVetoSupporterH2LV->GetMass(true,false)/kg << "(" 
				<< plasticVetoSupporterH2Box->GetXHalfLength()*2. << "x" 
				<< plasticVetoSupporterH2Box->GetYHalfLength()*2. << "x" 
				<< plasticVetoSupporterH2Box->GetZHalfLength()*2. << ")" << endl;
			cout << " PS veto stainless flame" << endl;
			cout << "     mass       : " << plasticVetoLV[0]->GetMass(true,false)/kg << endl;
			cout << "     dimension  : " 
				<< plasticVetoBox[0]->GetXHalfLength()*2 << " x "
				<< plasticVetoBox[0]->GetYHalfLength()*2 << " x " 
				<< plasticVetoBox[0]->GetZHalfLength()*2 << endl;
			cout << "     coordinate : " <<  physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() 
				+ GetPhysicalVolumeByName("PlasticVeto0_PV")->GetTranslation() << endl;
			cout << " PS veto Aluminium plate " << endl;
			cout << "     mass      : " << aluminiumHolderLV[0]->GetMass(true,false)/kg << endl;
			cout << "     dimension : " <<  plasticScintHolderBox[0]->GetXHalfLength()*2 << " x "
				<<  plasticScintHolderBox[0]->GetYHalfLength()*2 << " x "
				<<  plasticScintHolderBox[0]->GetZHalfLength()*2 << endl;
			cout << " Plastic scintillator" << endl;
			cout << "     dimension: " 
				<< plasticScintBox[0]->GetXHalfLength()*2 << " x " 
				<< plasticScintBox[0]->GetYHalfLength()*2 << " x " 
				<< plasticScintBox[0]->GetZHalfLength()*2 << endl;
			cout << "     mass     : " << plasticScintOLV[0]->GetMass(true,false)/kg << endl;
			cout << " PS veto detector ( # of veto: " << f200_VetoTotCNum
				<< " =>  barrel: " << nVetoZ*8 << ", bottom: " << nVetoB*2 << ")" << endl;
			cout << " ============================= Bottom shields ======" << endl;
			cout << " PE (box) " << endl; //..............................
			cout << "     mass      : " << shieldPELV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " 
				<< shieldPEBox->GetXHalfLength()*2 << " x " 
				<< shieldPEBox->GetYHalfLength()*2 << " x " 
				<< shieldPEBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< shieldPEBox->GetXHalfLength()-shieldBoricAcidBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() 
				+ shieldPEPV->GetTranslation() << endl;
			cout << " BoricAcid " << endl; //..............................
			cout << "     mass      : " << shieldBoricAcidLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " 
				<< shieldBoricAcidBox->GetXHalfLength()*2 << " x "
				<< shieldBoricAcidBox->GetYHalfLength()*2 << " x " 
				<< shieldBoricAcidBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< shieldBoricAcidBox->GetXHalfLength() - leadShieldBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() 
				+ shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() << endl;
			cout << " Lead shield " << endl; //..............................
			cout << "     mass      : " << leadShieldLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " 
				<< leadShieldBox->GetXHalfLength()*2 << " x "
				<< leadShieldBox->GetYHalfLength()*2 << " x " 
				<< leadShieldBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< leadShieldBox->GetXHalfLength() - ThinLeadShieldBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() 
				+ shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() 
				+ f200_VetoMaterialPhysical->GetTranslation() + BufferMother->GetTranslation() << endl;
			cout << " Thin lead shield " << endl; //..............................
			cout << "     mass      : " << ThinLeadShieldLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " 
				<< ThinLeadShieldBox->GetXHalfLength()*2 << " x "
				<< ThinLeadShieldBox->GetYHalfLength()*2 << " x " 
				<< ThinLeadShieldBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< ThinLeadShieldBox->GetXHalfLength() - boricAcidBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() 
				+ shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() 
				+ f200_VetoMaterialPhysical->GetTranslation() + BufferMother->GetTranslation() 
				+ GetPhysicalVolumeByName("ThinLeadShield_PV")->GetTranslation() << endl;
			cout << " Inner boric acid " << endl; //..............................
			cout << "     mass      : " << boricAcidLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " 
				<< boricAcidBox->GetXHalfLength()*2 << " x "
				<< boricAcidBox->GetYHalfLength()*2 << " x " 
				<< boricAcidBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " 
				<< boricAcidBox->GetXHalfLength() - IDspaceBox->GetXHalfLength() << endl;
			cout << "     coordinate: " 
				<< physCavern->GetTranslation() + vetoHousing1PV->GetTranslation() 
				+ shieldPEPV->GetTranslation() + shieldBoricAcidPV->GetTranslation() 
				+ f200_VetoMaterialPhysical->GetTranslation() + BufferMother->GetTranslation() 
				+ GetPhysicalVolumeByName("ThinLeadShield_PV")->GetTranslation()
				+ GetPhysicalVolumeByName("InnerBoricAcid_PV")->GetTranslation() << endl;
			cout << " Inner detector space " << endl; //..............................
			cout << "     dimension : " 
				<< IDspaceBox->GetXHalfLength()*2 << " x "
				<< IDspaceBox->GetYHalfLength()*2 << " x " 
				<< IDspaceBox->GetZHalfLength()*2 << endl;
			cout << " ============================= Upper shields ======" << endl;
			if(fAdditionalPE){
				cout << " Additional PE shield" << endl;
				cout << "     mass        : " << addPE1LV->GetMass(true,false)/kg 
					<< ", " << addPE2LV->GetMass(true,false)/kg << endl;
				cout << "     dimension(1): " << addPEBox1->GetXHalfLength()*2 << " x " 
					<< addPEBox1->GetYHalfLength()*2 << " x " 
					<< addPEBox1->GetZHalfLength()*2 << endl;
				cout << "     dimension(2): " << addPEBox2->GetXHalfLength()*2 << " x " 
					<< addPEBox2->GetYHalfLength()*2 << " x " 
					<< addPEBox2->GetZHalfLength()*2 << endl;
				cout << "     coordinate  : " << GetPhysicalVolumeByName("additionalPE1_PV")->GetTranslation()
					+ physCavern->GetTranslation() + GetPhysicalVolumeByName("HbeamHousing_PV")->GetTranslation() << ", " 
					<< GetPhysicalVolumeByName("additionalPE2_PV")->GetTranslation()
					+ physCavern->GetTranslation() + GetPhysicalVolumeByName("HbeamHousing_PV")->GetTranslation() << endl;
			}
			*/
			cout << " WC veto housing " << endl;
			cout << "     mass              : " << shieldWaterHousingLV->GetMass(true,false)/kg << endl;
			cout << "     dimension (out)   : " 
				<< shieldWaterHousingBox->GetXHalfLength()*2 << " x " 
				<< shieldWaterHousingBox->GetYHalfLength()*2 << " x " 
				<< shieldWaterHousingBox->GetZHalfLength()*2 << endl;
			cout << "     dimension (inner) : " 
				<< shieldHatAirBox->GetXHalfLength()*2 << " x " 
				<< shieldHatAirBox->GetYHalfLength()*2 << " x " 
				<< shieldHatAirBox->GetZHalfLength()*2 << endl;
			cout << "     thickness         : " << waterhousing_thickness << endl;
			cout << "     coordinate        : " << shieldWaterHousingPV->GetTranslation() 
				+ physCavern->GetTranslation() << endl;
			cout << " WC water" << endl;
			cout << "     mass      : " << shieldWaterTankLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " << shieldWaterTankBox->GetXHalfLength()*2 << " x "
				<< shieldWaterTankBox->GetYHalfLength()*2 << " x "
				<< shieldWaterTankBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " << watertank_thickness << "(top thickness: " << watertank_top_thickness << ")" << endl;
			cout << "     coordinate: " << shieldWaterTankPV->GetTranslation() 
				+ shieldWaterHousingPV->GetTranslation() + physCavern->GetTranslation() << endl;
			cout << " Hat Aluminium plate " << endl;
			cout << "     mass      : " << HatAlPlateLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " << HatBeamHousingInBox->GetXHalfLength()*2 << " x " 
				<< HatBeamHousingInBox->GetYHalfLength()*2 << " x " 
				<< HatBeamHousingInBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " << HatAlPlate_thickness << endl; 
			cout << "     coordinate: " << HatAlPlatePV->GetTranslation() + physCavern->GetTranslation() << endl;
			cout << " Hat boric acid " << endl;
			cout << "     mass      : " << shieldHatBoricLV->GetMass(true,false)/kg << endl;
			cout << "     dimension : " << HatAlPlateInBox->GetXHalfLength()*2 << " x " 
				<< HatAlPlateInBox->GetYHalfLength()*2 << " x " 
				<< HatAlPlateInBox->GetZHalfLength()*2 << endl;
			cout << "     thickness : " << boricacid_thickness << endl; 
			cout << "     coordinate: " << shieldHatBoricPV->GetTranslation() + physCavern->GetTranslation() << endl; 
		}
	}

	cout << " end " << endl;
	if(fRockgammaMode){
		retvalLV = shieldHousingLV;
	}
	else {
		retvalLV  = logiCavern;
	}
	return retvalLV;

}

////////////////////////////////////////////////////////////////
/// H-beam shape is not H shape but square pillar.
static G4LogicalVolume *MakeHBeam1(G4Material *mat)
{
	G4double beam1Boxsize_x = 8800/2.;
	G4double beam1Boxsize_y = 4975/2.;
	G4double beam1Boxsize_z = 300;
	G4double ssize_x        = 300/2.;
	G4double ssize_y        = 675/2.;
	G4double ssize_z        = 600;
	G4double lsize_x        = 3600/2.;
	G4double lsize_y        = 3375/2.;
	G4double lsize_z        = 600;

	G4ThreeVector lposHole;
	G4ThreeVector rposHole;
	G4SubtractionSolid *HbeamBot1Solid[20];
	G4SubtractionSolid *HbeamBot2Solid[20];

	G4Box *Hbeam1Box       = new G4Box("Hbeam1_Box",beam1Boxsize_x,beam1Boxsize_y,beam1Boxsize_z/2.);
	G4Box *smallsquarehole = new G4Box("smallsquarehole",ssize_x,ssize_y,ssize_z/2.);
	G4Box *largesquarehole = new G4Box("largesquarehole",lsize_x,lsize_y,lsize_z/2.);

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<5; j++)
		{
			lposHole = G4ThreeVector(
					-beam1Boxsize_x + 400 + ssize_x + ssize_x*2*i + 200*i,
					-beam1Boxsize_y + 400 + ssize_y + ssize_y*2*j + 200*j,
					0);
			if(j==0 && i==0)
			{
				HbeamBot1Solid[j+i*5] = new G4SubtractionSolid("HbeamBot_Solid",
						Hbeam1Box,smallsquarehole,0,lposHole);
			}else
			{
				HbeamBot1Solid[j+i*5] = new G4SubtractionSolid("HbeamBot_Solid",
						HbeamBot1Solid[j-1+i*5],smallsquarehole,0,lposHole);
			}
		}
	}
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<5; j++)
		{
			rposHole = G4ThreeVector(
					beam1Boxsize_x - 400 - ssize_x - ssize_x*2*i - 200*i,
					-beam1Boxsize_y + 400 + ssize_y + ssize_y*2*j + 200*j,
					0);
			if(j==0 && i==0)
			{
				HbeamBot2Solid[j+i*5] = new G4SubtractionSolid("HbeamBot_Solid",
						HbeamBot1Solid[19],smallsquarehole,0,rposHole);
			}else{
				HbeamBot2Solid[j+i*5] = new G4SubtractionSolid("HbeamBot_Solid",
						HbeamBot2Solid[j-1+i*5],smallsquarehole,0,rposHole);
			}
		}
	}
	G4SubtractionSolid *HbeamBotSolid = new G4SubtractionSolid("HbeamBot_Solid",
			HbeamBot2Solid[19],largesquarehole,0,G4ThreeVector(0,0,0));
	G4LogicalVolume *HBeam1LV         = new G4LogicalVolume(HbeamBotSolid,mat,"HbeamBot_LV");

	return HBeam1LV;
}

static G4LogicalVolume *MakePit(G4double pitBox_x, G4double pitBox_z, G4Material *mat, G4double HatBarrel_gap)
{
	//G4double pitBox_x = 9600/2.;
	//G4double pitBox_z = 6900;
	G4double pitBox_y = 3600/2.;
	G4double inBox1_y = (3600-600)/2;
	G4double inBox1_x = 4200/2;
	G4double inBox2_x = 6000/2; 
	G4double inBox3_y = 5000/2.;
	G4double beam1Boxsize_x = 8800/2.;
	G4double beam1Boxsize_y = 4975/2.;
	//G4double beam1Boxsize_z = 300;
	G4double beam1Boxsize_z = 460;

	G4Box *Hbeam1Box    = new G4Box("Hbeam1_Box",beam1Boxsize_x,beam1Boxsize_y,beam1Boxsize_z/2.);
	G4Box *pitMotherBox = new G4Box("pitMother_Box",pitBox_x,pitBox_y,pitBox_z/2.);
	G4Box *pitInBox1    = new G4Box("pitInner1_Box",inBox1_x,inBox1_y,pitBox_z/2.);
	G4Box *pitInBox2    = new G4Box("pitInner2_Box",inBox2_x,inBox1_y,pitBox_z/2.);
	G4Box *pitInBox3    = new G4Box("pitInner2_Box",inBox2_x,inBox3_y,pitBox_z/2.);

	G4VSolid *pit1Solid    = new G4SubtractionSolid("pitSolid",pitMotherBox,Hbeam1Box,0,
			G4ThreeVector(pitBox_x-2400,0,pitBox_z/2.-beam1Boxsize_z/2.));
	G4VSolid *pit2Solid    = new G4SubtractionSolid("pitSolid",pit1Solid,pitInBox1,0,
			G4ThreeVector(pitBox_x-inBox1_x-HatBarrel_gap,0.,HatBarrel_gap));
	G4VSolid *pit3Solid    = new G4SubtractionSolid("pitSolid",pit2Solid,pitInBox2,0,
			G4ThreeVector(-pitBox_x+inBox2_x+HatBarrel_gap,0,2600+HatBarrel_gap));
	G4VSolid *pitSolid    = new G4SubtractionSolid("pitSolid",pit3Solid,pitInBox3,0,
			G4ThreeVector(-inBox2_x,0.,-4300));

	G4LogicalVolume *pitLV                 = new G4LogicalVolume(pitSolid,mat,"pit_LV");
	return pitLV;
}
