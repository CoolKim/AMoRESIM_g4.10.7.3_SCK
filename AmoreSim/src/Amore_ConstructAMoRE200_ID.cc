//
//  Original by G. Horton-Smith 2004/12/02
//
//  Modified by E.J.Jeon 2007/06/14
//  Modified by Y.S.Yoon 2015/06/15
//  Updated by J.Seo 2024/02/15

#include "globals.hh"

#include "AmoreSim/AmoreCMOParameterisation.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreDetectorStaticInfo.hh"
#include "CupSim/CupPMTOpticalModel.hh" // for same PMT optical model as main sim
#include "CupSim/CupParam.hh"
#include "CupSim/Cup_PMT_LogicalVolume.hh" // for making PMT assemblies

#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Torus.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Version.hh"
#include "G4VisAttributes.hh"

#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"
//#include "CLHEP/Matrix/Matrix.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "G4Region.hh"
#include "G4Types.hh"

#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//using namespace CLHEP;
using namespace std;
using namespace AmoreDetectorStaticInfo;
using namespace AmoreDetectorStaticInfo::AMoRE_200;
using namespace AmoreDetectorStaticInfo::ColorTable;
using AMoRE200CrystalModuleInfo = AmoreDetectorStaticInfo::AMoRE200CrystalModuleInfo;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreDetectorConstruction::ConstructAMoRE200_ID(G4LogicalVolume *aWorkAreaLV) {

	G4int flagInvisible = 0;
	G4int flagOneCell   = 0;

	CupParam &db(CupParam::GetDB());

	////////////////////////////////////////////////////
	///// Primitive values retrived from database
	////////////////////////////////////////////////////
	G4double airbuffer_radius    = db["airbuffer_radius"];
	G4double airbuffer_height    = db["airbuffer_height"];

	G4double thin_lead_shield_thickness  = db["thin_lead_shield_thickness"];
	G4double lead_shield_thickness  = db["lead_shield_thickness"];
	G4double lead_housing_thickness = db["lead_housing_thickness"];

	G4double boricacid_thickness     = db["boricacid_thickness"];
	G4double PE_shield_thickness     = db["PE_shield_thickness"];
	G4double nShield_GapFromCeiling  = db["nShield_GapFromCeiling"];

	G4double Cu1GapFromBot = db["Cu1Gap"];
	G4double Cu2GapFromBot = db["Cu2Gap"];
	G4double Cu3GapFromBot = db["Cu3Gap"];
	G4double Cu4GapFromBot = db["Cu4Gap"];
	G4double TRGapFromCuP6 = db["TRGap"];

	///////////////////////////////////////////////////////
	///// composite variables
	///////////////////////////////////////////////////////

	/// Work Area
	G4LogicalVolume *logiWorkArea = aWorkAreaLV;
	G4Material *WA_mat = logiWorkArea->GetMaterial();

	///////////////////////////////////////////////////////
	// Build Source Housing: teflon tube
	// diameter 12 mm and thickness 1 mm 
	///////////////////////////////////////////////////////
	G4double housing_thickness = 1 * mm;
	G4double housing_radius = 12 / 2. * mm;
	G4double source_radius = 1 / 2. * mm;
	G4double source_xpos = ss_radius + housing_radius + housing_thickness + solidBooleanTol;
	G4double source_height[2] = {ss_inner_height_half*5/10, ss_inner_height_half*3/10};

	G4Torus *SHousing1 = new G4Torus("SHousing1", 
			housing_radius-housing_thickness, housing_radius, source_xpos, 0, 360. * degree);
	G4LogicalVolume *logiSourceHousing = new G4LogicalVolume(SHousing1, _teflon, "SourceHousing1LV");

	//////////////////////////////////////////////////////
	// Build Source Volume: Thorium wire
	// diameter 1 mm
	//////////////////////////////////////////////////////
	G4Torus *Source1 = new G4Torus("Source1", 0, source_radius,	source_xpos, 0, 360. * degree);
	//G4LogicalVolume *logiSourceVolume = new G4LogicalVolume(Source1, _ThWire, "SourceVolume1LV");
	G4LogicalVolume *logiSourceVolume = new G4LogicalVolume(Source1, _air, "SourceVolume1LV");


	///////////////////////////////////////////////////////
	// Build Layer5_OVC Stainless Steel (SS)
	///////////////////////////////////////////////////////
	G4ThreeVector RealModel_shield_topPos =  G4ThreeVector(room_dist_x, -room_dist_y, 
			-cavern_loaf_height/2. + 
			airbuffer_height + lead_shield_thickness + 
			lead_housing_thickness*2 + boricacid_thickness + PE_shield_thickness);

	G4Tubs *SSCylinder = new G4Tubs("SSCylinder", 0, 
			ss_radius, ss_inner_height_half, 0, 360. * deg);
	G4VSolid *SSTop      = new G4Box("SSTop", 
			sst_xsize_half, sst_ysize_half, sst_zsize_half);
	G4VSolid *SSOVC0     = new G4UnionSolid("SSOVC0", SSCylinder, SSTop, nullptr,
			{0, 0, ss_inner_height_half + sst_zsize_half});
	G4LogicalVolume *logiSSOVC  = new G4LogicalVolume(SSOVC0, _stainless, "logiSSOVC");

	///////////////////////////////////////////////////////
	// For Radon Air 
	///////////////////////////////////////////////////////
	G4Box *radonAirBox = new G4Box("RadonAir_Box",
			airbuffer_radius - boricacid_thickness, airbuffer_radius - boricacid_thickness, 
			(airbuffer_height - boricacid_thickness) / 2.);
	G4VSolid *radonAirSolid = new G4SubtractionSolid("RadonAir_Solid", radonAirBox, SSCylinder, 0,
			G4ThreeVector(0,0,
				radonAirBox->GetZHalfLength() - SSCylinder->GetZHalfLength()) );
	G4LogicalVolume *radonAirLV = new G4LogicalVolume(radonAirSolid, _air, "RadonAir_LV");


	if (fRockgammaMode){
		G4double workarea_posz= static_cast<G4Box *>(logiWorkArea->GetSolid())->GetZHalfLength();
		f200_OVCPhysical = new G4PVPlacement( nullptr, 
				{0, 0, -workarea_posz + PE_shield_thickness + lead_shield_thickness 
				+ boricacid_thickness + ovc_gap + SSCylinder->GetZHalfLength()}, 
				logiSSOVC, "physSSOVC", logiWorkArea, false, 0,OverlapCheck);
	}
	else if (fNeutronMode) {
		f200_OVCPhysical = new G4PVPlacement( nullptr, 
				{0, 0, - SSCylinder->GetZHalfLength() - nShield_GapFromCeiling + ovc_gap}, 
				logiSSOVC, "physSSOVC", logiWorkArea, false, 0,OverlapCheck);
	} 
	else{
		G4double ovcPosZ = boricacid_thickness + thin_lead_shield_thickness + lead_shield_thickness + 
			lead_housing_thickness*2 + boricacid_thickness + PE_shield_thickness;
		G4ThreeVector ovcTlate = G4ThreeVector(0,0, 
				SSCylinder->GetZHalfLength() + ovcPosZ + ovc_gap);
		G4ThreeVector airTlate = G4ThreeVector(0,0,
				radonAirBox->GetZHalfLength() + ovcPosZ);
		if(whichCavernType==kCavern_RealModel){
			ovcTlate = RealModel_shield_topPos + 
				G4ThreeVector(0, 0, - SSCylinder->GetZHalfLength() + ovc_gap + nShield_GapFromCeiling);
			airTlate = RealModel_shield_topPos + G4ThreeVector(0,0, 
					- radonAirBox->GetZHalfLength() + nShield_GapFromCeiling + ovc_gap);
		}

		f200_OVCPhysical = new G4PVPlacement( nullptr, ovcTlate, logiSSOVC, 
				"physSSOVC", logiWorkArea, false, 0,OverlapCheck);

		// Radon Air -----------
		new G4PVPlacement(nullptr, airTlate,
				radonAirLV, "RadonAir_PV", logiWorkArea, false, 0, OverlapCheck);

		// Calibration source
		//G4ThreeVector SourcePosition = ovcTlate - G4ThreeVector(0,0,source_height[0]);
		G4ThreeVector SourcePosition = G4ThreeVector(0,0,-source_height[0]);
		new G4PVPlacement ( nullptr, SourcePosition, 
				logiSourceHousing, "physSourceHousing1", radonAirLV, false, 0, OverlapCheck);
		new G4PVPlacement( nullptr, SourcePosition,
				logiSourceVolume, "physSourceVolume1", radonAirLV, false, 0, OverlapCheck);

		SourcePosition = G4ThreeVector(0,0, -source_height[1]);
		new G4PVPlacement ( nullptr, SourcePosition, 
				logiSourceHousing, "physSourceHousing2", radonAirLV, false, 0, OverlapCheck);
		new G4PVPlacement( nullptr, SourcePosition,
				logiSourceVolume, "physSourceVolume2", radonAirLV, false, 0, OverlapCheck);
	}

	G4Tubs *SSOVC1      = new G4Tubs("SSOVC_InSpaceSol", 0, ss_radius - OVCthick, 
			//ss_inner_height_half - OVCthick/2. - sst_zsize_half, 0, 360 * deg);
			ss_inner_height_half - OVCthickB/2., 0, 360 * deg);

	G4LogicalVolume *logiSSOVCIn = new G4LogicalVolume(SSOVC1, WA_mat, "SSOVC_InnerSpaceLV");
	//new G4PVPlacement(0, G4ThreeVector(0, 0, OVCthick/2.), 
	new G4PVPlacement(0, G4ThreeVector(0, 0, OVCthickB/2.), 
			logiSSOVCIn, "physSSOVCIn", logiSSOVC, false, 0, OverlapCheck);

	///////////////////////////////////////////////////////
	// Layer4_ 50K-SHIELD Copper
	///////////////////////////////////////////////////////
	//G4double cu4zthick = cu4thick/2.;
	G4Tubs *Cu4Cylinder = new G4Tubs("Cu4Cylinder", 0, cu4_radius, cu4_height, 0, 360. * deg);
	G4LogicalVolume *logiCu4     = new G4LogicalVolume(Cu4Cylinder, _copper, "logiCu4", 0, 0, 0);
	new G4PVPlacement(0, // no rotation
			G4ThreeVector(0, 0, -SSOVC1->GetZHalfLength() + cu4_height + Cu4GapFromBot),
			logiCu4, "physCu4", logiSSOVCIn, false, 0,OverlapCheck);

	G4Tubs *Cu4CylIn  = new G4Tubs("Cu4CylinderInSol", 0, cu4_radius - cu4thick,
			//cu4_height - cu4zthick, 0, 360 * deg);
			cu4_height - cu4thickB/2., 0, 360 * deg);
	G4LogicalVolume *logiCu4In = new G4LogicalVolume(Cu4CylIn, WA_mat, "Cu4CylinderInLV");
	//new G4PVPlacement(0, G4ThreeVector(0, 0, cu4thick/2.), 
	new G4PVPlacement(0, G4ThreeVector(0, 0, cu4thickB/2.), 
			logiCu4In, "physCu4In", logiCu4, false, 0, OverlapCheck);


	///////////////////////////////////////////////////////
	// Layer3_Cu IVC
	///////////////////////////////////////////////////////
	//G4double cu3zthick = cu3thick/2.;
	G4Tubs *Cu3Cylinder = new G4Tubs("Cu3Cylinder", 0, cu3_radius, cu3_height, 0, 360. * deg);
	G4LogicalVolume *logiCu3     = new G4LogicalVolume(Cu3Cylinder, _copper, "logiCu3", 0, 0, 0);
	G4Tubs *Cu3CylIn    = new G4Tubs("Cu3CylinderInSol", 0, 
			cu3_radius - cu3thick, cu3_height - cu3thickB/2., 0, 360 * deg);
	G4LogicalVolume *logiCu3In   = new G4LogicalVolume(Cu3CylIn, WA_mat, "Cu3CylinderInLV");
	new G4PVPlacement(0, // no rotation
			G4ThreeVector(0, 0, -Cu4CylIn->GetZHalfLength() + cu3_height + Cu3GapFromBot),
			logiCu3, "physCu3", logiCu4In, false,0, OverlapCheck);
	//new G4PVPlacement(0, G4ThreeVector(0, 0, cu3thick/2.), 
	new G4PVPlacement(0, G4ThreeVector(0, 0, cu3thickB/2.), 
			logiCu3In, "physCu3In", logiCu3, false, 0, OverlapCheck);


	///////////////////////////////////////////////////////
	// Layer2_Cu3 SHIELD-STILL
	///////////////////////////////////////////////////////
	//G4double cu2zthick = cu2thick/2.;
	G4Tubs *Cu2Cylinder = new G4Tubs("Cu2Cylinder", 0, cu2_radius, cu2_height, 0, 360. * deg);
	G4LogicalVolume *logiCu2     = new G4LogicalVolume(Cu2Cylinder, _copper, "logiCu2", 0, 0, 0);
	G4Tubs *Cu2CylIn    = new G4Tubs("Cu2CylinderInSol", 0, 
			//cu2_radius - cu2thick, cu2_height - cu2zthick, 0, 360 * deg);
			cu2_radius - cu2thick, cu2_height - cu2thickB/2., 0, 360 * deg);
	G4LogicalVolume *logiCu2In   = new G4LogicalVolume(Cu2CylIn, WA_mat, "Cu2CylinderInLV");
	new G4PVPlacement(0, // no rotation
			G4ThreeVector(0, 0, -Cu3CylIn->GetZHalfLength() + cu2_height + Cu2GapFromBot),
			logiCu2, "physCu2", logiCu3In, false, 0,OverlapCheck);
	//new G4PVPlacement(0, G4ThreeVector(0, 0, cu2thick/2.), 
	new G4PVPlacement(0, G4ThreeVector(0, 0, cu2thickB/2.), 
			logiCu2In, "physCu2In", logiCu2, false, 0, OverlapCheck);

	///////////////////////////////////////////////////////
	// Layer1_Cu4 50mk-SHIELD
	///////////////////////////////////////////////////////
	//G4double cu1zthick = cu1thick/2.;
	G4Tubs *Cu1Cylinder = new G4Tubs("Cu1Cylinder", 0, cu1_radius, cu1_height, 0, 360. * deg);
	G4LogicalVolume *logiCu1     = new G4LogicalVolume(Cu1Cylinder, _copper, "logiCu1", 0, 0, 0);
	G4Tubs *Cu1CylIn    = new G4Tubs("Cu1CylInSol", 0, 
			//cu1_radius - cu1thick, cu1_height - cu1zthick, 0, 360 * deg);
			cu1_radius - cu1thick, cu1_height - cu1thickB/2., 0, 360 * deg);
	G4LogicalVolume *logiCu1In   = new G4LogicalVolume(Cu1CylIn, WA_mat, "Cu1CylinderInLV");
	new G4PVPlacement(0, // no rotation
			G4ThreeVector(0, 0, -Cu2CylIn->GetZHalfLength() + cu1_height + Cu1GapFromBot),
			logiCu1, "physCu1", logiCu2In, false, 0,OverlapCheck);
	//new G4PVPlacement(0, G4ThreeVector(0, 0, cu1thick/2.), 
	new G4PVPlacement(0, G4ThreeVector(0, 0, cu1thickB/2.), 
			logiCu1In, "physCu1In", logiCu1, false, 0, OverlapCheck);

	///////////////////////////////////////////////////////
	// Build Mixing Chamber Cu, Pb Plates
	///////////////////////////////////////////////////////
	/// Cu Plates
	G4double CenterCuMCPZ  = cu1_height - cumcp_height - plateGap;
	G4ThreeVector PosCuMCP = G4ThreeVector(0., 0., CenterCuMCPZ);

	G4double CenterCuP1Z  = CenterCuMCPZ - cumcp_height - cup1_zsize - plateGap2;
	G4ThreeVector PosCuP1 = G4ThreeVector(0., 0., CenterCuP1Z);

	G4double CenterPbP1Z  = CenterCuP1Z - cup1_zsize - pbp1_zsize;
	G4ThreeVector PosPbP1 = G4ThreeVector(0., 0., CenterPbP1Z);

	G4double CenterCuP2Z  = CenterPbP1Z - pbp1_zsize - cup2_zsize;
	G4ThreeVector PosCuP2 = G4ThreeVector(0., 0., CenterCuP2Z);

	G4double CenterPbP2Z  = CenterCuP2Z - cup2_zsize - pbp1_zsize;
	G4ThreeVector PosPbP2 = G4ThreeVector(0., 0., CenterPbP2Z);

	G4double CenterCuP3Z  = CenterPbP2Z - pbp1_zsize - cup3_zsize;
	G4ThreeVector PosCuP3 = G4ThreeVector(0., 0., CenterCuP3Z);

	G4double CenterPbP3Z  = CenterCuP3Z - cup3_zsize - pbp1_zsize;
	G4ThreeVector PosPbP3 = G4ThreeVector(0., 0., CenterPbP3Z);

	G4double CenterCuP4Z  = CenterPbP3Z - pbp1_zsize - cup2_zsize;
	G4ThreeVector PosCuP4 = G4ThreeVector(0., 0., CenterCuP4Z);

	G4double CenterPbP4Z  = CenterCuP4Z - cup2_zsize - pbp2_zsize; 
	G4ThreeVector PosPbP4 = G4ThreeVector(0., 0., CenterPbP4Z);

	G4double CenterCuP5Z  = CenterPbP4Z - pbp2_zsize - cup1_zsize;
	G4ThreeVector PosCuP5 = G4ThreeVector(0., 0., CenterCuP5Z);

	///////////////////////////////////////////////////////
	///// Build Cu and Lead Plates
	///////////////////////////////////////////////////////
	G4Tubs *CuMCPlate = new G4Tubs("CuMCPlate", 0, cumcp_radius, cumcp_height, 0, 360. * deg);
	G4Tubs *CuPlate1  = new G4Tubs("CuPlate1", 0, inlead_radius, cup1_zsize, 0, 360. * deg);
	G4Tubs *CuPlate2  = new G4Tubs("CuPlate2", 0, inlead_radius, cup2_zsize, 0, 360. * deg);
	G4Tubs *CuPlate3  = new G4Tubs("CuPlate3", 0, inlead_radius, cup3_zsize, 0, 360. * deg);
	G4Tubs *PbPlate1  = new G4Tubs("PbPlate1", 0, inlead_radius, pbp1_zsize, 0, 360. * deg);
	G4Tubs *PbPlate2  = new G4Tubs("PbPlate2", 0, inlead_radius, pbp2_zsize, 0, 360. * deg);
	G4LogicalVolume *logiCuMCP = new G4LogicalVolume(CuMCPlate, _copper, "logiCuMCP", 0, 0, 0);
	G4LogicalVolume *logiCuP1  = new G4LogicalVolume(CuPlate1, _copper, "logiCuP1", 0, 0, 0);
	G4LogicalVolume *logiCuP2  = new G4LogicalVolume(CuPlate2, _copper, "logiCuP2", 0, 0, 0);
	G4LogicalVolume *logiCuP3  = new G4LogicalVolume(CuPlate3, _copper, "logiCuP3", 0, 0, 0);
	G4LogicalVolume *logiPbP1  = new G4LogicalVolume(PbPlate1, _lead, "logiPbP1", 0, 0, 0);
	G4LogicalVolume *logiPbP2  = new G4LogicalVolume(PbPlate2, _lead, "logiPbP2", 0, 0, 0);


	new G4PVPlacement(0, // no rotation
			PosCuMCP, logiCuMCP, "physCuMCP", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosCuP1, logiCuP1, "physCuP1", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosPbP1, logiPbP1, "physPbP1", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosCuP2, logiCuP2, "physCuP2", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosPbP2, logiPbP1, "physPbP2", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosCuP3, logiCuP3, "physCuP3", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosPbP3, logiPbP1, "physPbP3", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosCuP4, logiCuP2, "physCuP4", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosPbP4, logiPbP2, "physPbP4", logiCu1In, false, 0,OverlapCheck);
	new G4PVPlacement(0, // no rotation
			PosCuP5, logiCuP1, "physCuP5", logiCu1In, false, 0,OverlapCheck);

	///////////////////////////////////////////////////////
	// SC shield
	///////////////////////////////////////////////////////
	G4double SC_radius = 991./2.;
	G4double SC_height = 1200./2.; //840./2.;
	G4double SC_thick = 1.;
	G4double SC_overlap = 150.;
	G4ThreeVector PosSC = G4ThreeVector(0., 0., CenterCuP5Z - cup1_zsize - SC_height + SC_overlap);
	G4Tubs *SCshield_out = new G4Tubs("SCshield_out", 0, SC_radius, SC_height, 0, 360. * deg);
	G4Tubs *SCshield_in = new G4Tubs("SCshield_in", 0, SC_radius-SC_thick, SC_height-SC_thick/2., 0, 360. * deg); 
	G4VSolid *SCshield = new G4SubtractionSolid("SCshield", SCshield_out, SCshield_in, 0, G4ThreeVector(0, 0, SC_thick/2. + solidBooleanTol) );
	G4LogicalVolume *logiSCshield = new G4LogicalVolume(SCshield, _lead, "logiSCshield", 0, 0, 0);
	new G4PVPlacement(0, PosSC, logiSCshield, "physSCshield", logiCu1In, false, 0, OverlapCheck);


	///////////////////////////////////////////////////////
	// CMO Target area
	///////////////////////////////////////////////////////
	G4double TR_radius = inlead_radius;
	G4double TR_height = 800./2.; // 900/2.; // This value should be changed !!!!!

	G4ThreeVector PosTR = G4ThreeVector(0., 0., CenterCuP5Z - TR_height - cup1_zsize - TRGapFromCuP6);
	G4Tubs *TargetRoom  = new G4Tubs("TargetRoom", 0, TR_radius, TR_height, 0, 360. * deg);
	G4LogicalVolume *logiTargetRoom = new G4LogicalVolume(TargetRoom, _vacuum, "logiTargetRoom", 0, 0, 0);

	new G4PVPlacement(0, // no rotation
			PosTR, logiTargetRoom, "physTargetRoom", logiCu1In, false, 0,OverlapCheck);

	//////////////////////////////////////////////////////
	/// CRYSTAL ARRAY
	//////////////////////////////////////////////////////
	G4Tubs *CrystalArray = new G4Tubs("CrystalArray", 0, Array_radius, Array_height/2., 0, 360. * deg);
	G4LogicalVolume *logiCrystalArray  = new G4LogicalVolume(CrystalArray, _vacuum, "CrystalArrayLV");
	new G4PVPlacement(0, {0,0,TR_height - Array_height/2.}, logiCrystalArray, "physCrystalArray", logiTargetRoom, false, 0, OverlapCheck);

	//////////////////////////////////////////////////////
	/// CRYSTAL TOWER
	//////////////////////////////////////////////////////
	G4ThreeVector PosTower = G4ThreeVector(0., 0., Array_height/2. - Tower_height/2.);
	switch (whichAMoRE200Phase) {
		case kPhase1:{
									 f200_TotTowerNum     = 10;
									 f200_TotCrystalNum   = f200_TotTowerNum * nModuleInTower;
									 break;
								 }
		case kPhase2:{
									 f200_TotTowerNum   = 47; // for 178kg Crystal, 70 for proposed tower design;
									 //f200_TotTowerNum   = 48; // for 178kg Crystal, 70 for proposed tower design;
									 f200_TotCrystalNum = f200_TotTowerNum * nModuleInTower;
									 break;
								 }
		default:
								 G4Exception("AmoreDetectorConstruction::ConstructAMoRE200_ID", "", 
										 G4ExceptionSeverity::FatalException, "AMoRE200 Phase type was not selected.");
								 break;
	}


	//////////////////////////////////////////////////////
	/// CRYSTAL MODULES
	//////////////////////////////////////////////////////
	for(int itower = 0; itower < f200_TotTowerNum; itower++){
		const AMoRE200CrystalModuleInfo &nowInfo = crystalModuleInfoList[itower];
		G4double cell_h = nowInfo.fCrystalHeight;
		G4double tower_x = nowInfo.fTowerX;
		G4double tower_y = nowInfo.fTowerY;
		G4double module_h = cell_h + bottomframe_height + topframe_height + solidBooleanTol;

		G4ThreeVector modulePos = G4ThreeVector(tower_x, tower_y, -Tower_height/2 + module_h/2.);
		for(int imodule = 0; imodule < nModuleInTower; imodule++){
			G4LogicalVolume *logiCrystalModule = MakeModule(_vacuum, _Li2MoO4, _vm2000, _copper, _copper3, _teflon2, _SiWafer, _gold, itower, imodule);
			new G4PVPlacement(0, modulePos, logiCrystalModule, ("physCrystalModule" + to_string(itower) + "_" + to_string(imodule)).c_str(), logiCrystalArray, false, 0, OverlapCheck);
			modulePos[2] += module_h;
		}
	}

	//f200_logiCrystalCell = G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalCellLV",false);


	//////////////////////////////////////////////////////
	G4ThreeVector globalPos_OVC = f200_OVCPhysical->GetTranslation();
	G4ThreeVector globalPos_Cu1In = f200_OVCPhysical->GetTranslation() + 
		GetPhysicalVolumeByName("physSSOVCIn")->GetTranslation() + 
		GetPhysicalVolumeByName("physCu4")->GetTranslation() + GetPhysicalVolumeByName("physCu4In")->GetTranslation() + 
		GetPhysicalVolumeByName("physCu3")->GetTranslation() + GetPhysicalVolumeByName("physCu3In")->GetTranslation() + 
		GetPhysicalVolumeByName("physCu2")->GetTranslation() + GetPhysicalVolumeByName("physCu2In")->GetTranslation() + 
		GetPhysicalVolumeByName("physCu1")->GetTranslation() + GetPhysicalVolumeByName("physCu1In")->GetTranslation();
	G4ThreeVector globalPos_Source1 = {0,0,0};
	G4ThreeVector globalPos_Source2 = {0,0,0};
	if(!fNeutronMode){
		globalPos_Source1 = GetPhysicalVolumeByName("physSourceHousing1")->GetTranslation();
		globalPos_Source2 = GetPhysicalVolumeByName("physSourceHousing2")->GetTranslation();
	}
	if(!fRockgammaMode) {
		globalPos_Cu1In += GetPhysicalVolumeByName("physCavern")->GetTranslation();
		globalPos_OVC  += GetPhysicalVolumeByName("physCavern")->GetTranslation();
		if(!fNeutronMode){
			globalPos_Source1 += GetPhysicalVolumeByName("physCavern")->GetTranslation();
			globalPos_Source2 += GetPhysicalVolumeByName("physCavern")->GetTranslation();
		}
	}
	G4ThreeVector globalPos_PbP1 = globalPos_Cu1In + GetPhysicalVolumeByName("physPbP1")->GetTranslation();
	G4ThreeVector globalPos_PbP2 = globalPos_Cu1In + GetPhysicalVolumeByName("physPbP2")->GetTranslation();
	G4ThreeVector globalPos_PbP3 = globalPos_Cu1In + GetPhysicalVolumeByName("physPbP3")->GetTranslation();
	G4ThreeVector globalPos_PbP4 = globalPos_Cu1In + GetPhysicalVolumeByName("physPbP4")->GetTranslation();
	G4ThreeVector globalPos_TR = globalPos_Cu1In + GetPhysicalVolumeByName("physTargetRoom")->GetTranslation();

	if (fDbgMsgOn) {
		cout << "\n ========================== Inner Volumes ============" << endl;
		cout << " (length: mm, weight: kg, global coordinate) \n" << endl;
		cout << " OVC " << endl;
		cout << "     mass            : " << logiSSOVC->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh) : " << SSCylinder->GetOuterRadius() << " x " << SSCylinder->GetZHalfLength()*2 << endl;
		cout << "     thickness       : " << OVCthick << "(" << OVCthickB << ")"<< endl;
		cout << "     coordinate      : " << globalPos_OVC<< endl;
		cout << " 50K Cu Chamber " << endl;
		cout << "     mass            : " << logiCu4->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh) : " << Cu4Cylinder->GetOuterRadius() << " x " << Cu4Cylinder->GetZHalfLength()*2 << endl;
		cout << "     thickness (bot) : " << cu4thick << "(" << cu4thickB << ")" << endl;
		cout << "     coordinate      : " << GetPhysicalVolumeByName("physCu4")->GetTranslation()
			+ globalPos_OVC + GetPhysicalVolumeByName("physSSOVCIn")->GetTranslation() << endl;
		cout << " IVC " << endl;
		cout << "     mass            : " << logiCu3->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh) : " << Cu3Cylinder->GetOuterRadius() << " x " << Cu3Cylinder->GetZHalfLength()*2 << endl;
		cout << "     thickness (bot) : " << cu3thick << "(" << cu3thickB << ")" << endl;
		cout << "     coordinate      : " << GetPhysicalVolumeByName("physCu3")->GetTranslation()
			+ globalPos_OVC + GetPhysicalVolumeByName("physSSOVCIn")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu4")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu4In")->GetTranslation()<< endl;
		cout << " 1K Cu Chamber " << endl;
		cout << "     mass            : " << logiCu2->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh) : " << Cu2Cylinder->GetOuterRadius() << " x " << Cu2Cylinder->GetZHalfLength()*2 << endl;
		cout << "     thickness (bot) : " << cu2thick << "(" << cu2thickB << ")" << endl;
		cout << "     coordinate      : " << GetPhysicalVolumeByName("physCu2")->GetTranslation()
			+ globalPos_OVC + GetPhysicalVolumeByName("physSSOVCIn")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu4")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu4In")->GetTranslation()
			+ GetPhysicalVolumeByName("physCu3")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu3In")->GetTranslation() << endl;
		cout << " 50mK Cu Chamber " << endl;
		cout << "     mass            : " << logiCu1->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh) : " << Cu1Cylinder->GetOuterRadius() << " x " << Cu1Cylinder->GetZHalfLength()*2 << endl;
		cout << "     thickness (bot) : " << cu1thick << "(" << cu1thickB << ")" << endl;
		cout << "     coordinate      : " << GetPhysicalVolumeByName("physCu1")->GetTranslation()
			+ globalPos_OVC + GetPhysicalVolumeByName("physSSOVCIn")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu4")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu4In")->GetTranslation()
			+ GetPhysicalVolumeByName("physCu3")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu3In")->GetTranslation()
			+ GetPhysicalVolumeByName("physCu2")->GetTranslation() 
			+ GetPhysicalVolumeByName("physCu2In")->GetTranslation() << endl;
		cout << "    Cu1In           : " << globalPos_Cu1In << endl;
		cout << "=============== Inner plate shields ==================" << endl;
		cout << " Mix Plate " << endl;
		cout << "     mass           : " << logiCuMCP->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << CuMCPlate->GetOuterRadius() << " x " << CuMCPlate->GetZHalfLength()*2 << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physCuMCP")->GetTranslation()
			+ globalPos_Cu1In << endl;
		cout << " Cu plate(top,bot) " << endl;
		cout << "     mass           : " << logiCuP1->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << CuPlate1->GetOuterRadius() << " x " << CuPlate1->GetZHalfLength()*2 << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physCuP1")->GetTranslation()
			+ globalPos_Cu1In << ", " << GetPhysicalVolumeByName("physCuP5")->GetTranslation()+globalPos_Cu1In << endl;
		cout << " Pb plate(thick) " << endl;
		cout << "     mass           : " << logiPbP1->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << PbPlate1->GetOuterRadius() << " x " << PbPlate1->GetZHalfLength()*2 << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physPbP1")->GetTranslation()
			+ globalPos_Cu1In << ", " << GetPhysicalVolumeByName("physPbP2")->GetTranslation()+globalPos_Cu1In
			<< ", " << GetPhysicalVolumeByName("physPbP3")->GetTranslation()+globalPos_Cu1In << endl;
		cout << " Cu plate(thin) " << endl;
		cout << "     mass           : " << logiCuP2->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << CuPlate2->GetOuterRadius() << " x " << CuPlate2->GetZHalfLength()*2 << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physCuP2")->GetTranslation()
			+ globalPos_Cu1In << ", " << GetPhysicalVolumeByName("physCuP4")->GetTranslation()+globalPos_Cu1In << endl;
		cout << " Pb plate(thin) " << endl;
		cout << "     mass           : " << logiPbP2->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << PbPlate2->GetOuterRadius() << " x " << PbPlate2->GetZHalfLength()*2 << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physPbP4")->GetTranslation()
			+ globalPos_Cu1In << endl;
		cout << " Cu plate(additional) " << endl;
		cout << "     mass           : " << logiCuP3->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << CuPlate3->GetOuterRadius() << " x " << CuPlate3->GetZHalfLength()*2 << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physCuP3")->GetTranslation()
			+ globalPos_Cu1In << endl;
		cout << " SC shield " << endl;
		cout << "     mass           : " << logiSCshield->GetMass(true,false)/kg << endl;
		cout << "     dimension (rxh): " << SCshield_out->GetOuterRadius() << " x " << SCshield_out->GetZHalfLength()*2 << endl;
		cout << "     thickness      : " << SC_thick << endl;
		cout << "     coordinate     : " << GetPhysicalVolumeByName("physSCshield")->GetTranslation() + globalPos_Cu1In << endl;
		cout << " ================== Calibration system ==========" << endl;
		cout << " Source housing" << endl;
		cout << "     mass               : " << logiSourceHousing->GetMass(true,false)/kg << endl;
		cout << "     dimension (r x dia): " << SHousing1->GetRtor() << " x " << SHousing1->GetRmax()*2 << endl;
		cout << "     tube thickness     : " << SHousing1->GetRmax()-SHousing1->GetRmin() << endl;
		cout << " Source wire" << endl;
		cout << "     mass                 : " << 
			G4LogicalVolumeStore::GetInstance()->GetVolume("SourceVolume1LV",false)->GetMass(true,false)/kg << endl;
		cout << "     dimension (r x thick): " << Source1->GetRtor() << " x " << Source1->GetRmax()*2 << endl;
		cout << "     coordinate           : " << globalPos_Source1 << ", " << globalPos_Source2 << endl;
		cout << " ================== Crystal modules  ==========" << endl;
		cout << " The number of towers                 : " << f200_TotTowerNum << endl;
		cout << " The number of crystals in the tower  : " << nModuleInTower << endl;
		cout << " Total number of crystals             : " << f200_TotTowerNum * nModuleInTower << endl;
		cout << "     mass (lmo6, lmo5)    : " << f200_logiCrystalCell[0]->GetMass(true,false)/kg << ", " 
			<< f200_logiCrystalCell[9]->GetMass(true,false)/kg << endl;
		double totcellmass = 0;
		for (int i = 0; i < f200_TotTowerNum * nModuleInTower; i++){
			totcellmass += f200_logiCrystalCell[i]->GetMass(true,false)/kg;
		}
		cout << "     mass (total crystal) : " << totcellmass << endl;
		cout << " Reflector " << endl;
		cout << "     mass (side,bottom)           : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("ReflectorLV",false)->GetMass(true,false)/kg << ", " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("ReflectorBLV",false)->GetMass(true,false)/kg << endl;
		cout << "     thickness                    : " << reflector_thick << endl;
		cout << "     local coordinate(side,bottom): " << GetPhysicalVolumeByName("physReflector0")->GetTranslation() << ", " 
			<< GetPhysicalVolumeByName("physReflectorBottom0")->GetTranslation() << endl;
		cout << " Bottom frame " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BottomModuleLV",false)->GetMass(true,false)/kg << endl;
		cout << "     thickness       : " << Module_thick << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physBottomModule0")->GetTranslation() << endl;
		cout << " Posts " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("PostLV",false)->GetMass(true,false)/kg << endl;
		cout << "     thickness       : " << post_radius << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physPost0")->GetTranslation() << endl;
		cout << " Top frame " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("TopModuleLV",false)->GetMass(true,false)/kg << endl;
		cout << "     thickness       : " << Module_thick << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physTopModule0")->GetTranslation() << endl;
		cout << " Photon frame " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("LightDetectorLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physLightDetector0")->GetTranslation() << endl;
		cout << " Phonon frame " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("HeatDetectorLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physHeatDetector0")->GetTranslation() << endl;
		cout << " Crystal module frame total mass: " << 
			G4LogicalVolumeStore::GetInstance()->GetVolume("BottomModuleLV",false)->GetMass(true,false)/kg 
			+ G4LogicalVolumeStore::GetInstance()->GetVolume("TopModuleLV",false)->GetMass(true,false)/kg
			+ 3 * G4LogicalVolumeStore::GetInstance()->GetVolume("PostLV",false)->GetMass(true,false)/kg
			+ G4LogicalVolumeStore::GetInstance()->GetVolume("LightDetectorLV",false)->GetMass(true,false)/kg
			+ G4LogicalVolumeStore::GetInstance()->GetVolume("HeatDetectorLV",false)->GetMass(true,false)/kg << endl;
		cout << " Clamp (top)" << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("ClampTopLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physClampTop0")->GetTranslation() << endl;
		cout << " Clamp (wafer) "<< endl; 
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("ClampWaferLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physClampWafer0")->GetTranslation() << endl;
		cout << " Clamp (bottomp)" << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("ClampBottomLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physClampBottom0")->GetTranslation() << endl;
		cout << " Teflon sheet " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("TeflonSheetLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physTeflonSheet0")->GetTranslation() << endl;
		cout << " Wafer " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("WaferLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physWafer0")->GetTranslation() << endl;
		cout << " Upper gold film " << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("UpperGoldFilmLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physUpperGold0")->GetTranslation() << endl;
		cout << " Lower gold film" << endl;
		cout << "     mass            : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BottomGoldFilmLV",false)->GetMass(true,false)/kg << endl;
		cout << "     local coordinate: " << GetPhysicalVolumeByName("physBottomGold0")->GetTranslation() << endl;
		cout << " Bolts mass " << endl;
		cout << "     M5 post : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BoltsPostLV",false)->GetMass(true,false)/kg << endl;
		cout << "     M4 long : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BoltsClampLV",false)->GetMass(true,false)/kg << endl;
		cout << "     M4 short: " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BoltsBodyLV",false)->GetMass(true,false)/kg << endl;
		cout << "     M4 thin : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BoltsLightLV",false)->GetMass(true,false)/kg << endl;
		cout << "     M4 head : " 
			<< G4LogicalVolumeStore::GetInstance()->GetVolume("BoltsHeatLV",false)->GetMass(true,false)/kg << endl;

		cout << " ...... Global Position ...... " << endl;
		cout << "     crystal Array : " << globalPos_TR + GetPhysicalVolumeByName("physCrystalArray")->GetTranslation() << endl;
		cout << "     crystal module0_0: " << globalPos_TR + GetPhysicalVolumeByName("physCrystalArray")->GetTranslation() + GetPhysicalVolumeByName("physCrystalModule0_0")->GetTranslation() << endl;
		cout << "     clamp top0    : " << globalPos_TR + GetPhysicalVolumeByName("physCrystalArray")->GetTranslation() + GetPhysicalVolumeByName("physCrystalModule0_0")->GetTranslation() + GetPhysicalVolumeByName("physClampTop0")->GetTranslation() << endl;

/*
		cout << " post global position------" <<endl;
		for(int i = 0; i < f200_TotTowerNum; i++){
			for(int j = 0; j < nModuleInTower; j++){
				G4int detid = i * 9 + j;
				cout << "		physPost"<< detid << "_" << 0 <<": " << globalPos_TR + GetPhysicalVolumeByName("physCrystalArray")->GetTranslation() +
					GetPhysicalVolumeByName(("physCrystalModule" + to_string(i) + "_" + to_string(j)).c_str())->GetTranslation() +
					GetPhysicalVolumeByName(("physPost"+to_string(detid)+"_0").c_str())->GetTranslation() + G4ThreeVector(3,0,0) << endl;
				cout << "		physPost"<< detid << "_" << 1 <<": " << globalPos_TR + GetPhysicalVolumeByName("physCrystalArray")->GetTranslation() +
					GetPhysicalVolumeByName(("physCrystalModule" + to_string(i) + "_" + to_string(j)).c_str())->GetTranslation() +
					GetPhysicalVolumeByName(("physPost"+to_string(detid)+"_1").c_str())->GetTranslation() + G4ThreeVector(3,0,0) << endl;
				cout << "		physPost"<< detid << "_" << 2 <<": " << globalPos_TR + GetPhysicalVolumeByName("physCrystalArray")->GetTranslation() +
					GetPhysicalVolumeByName(("physCrystalModule" + to_string(i) + "_" + to_string(j)).c_str())->GetTranslation() +
					GetPhysicalVolumeByName(("physPost"+to_string(detid)+"_2").c_str())->GetTranslation() + G4ThreeVector(3,0,0) << endl;
			}
		}
*/
		cout << " ===================================================" << endl;
	}

	//////////////////////////////
	// Set Attributes
	//////////////////////////////
	G4VisAttributes *logiTargetRoomVis = new G4VisAttributes(G4Colour(0,0,0,0.2));
	logiTargetRoom->SetVisAttributes(logiTargetRoomVis);
	logiTargetRoom->SetVisAttributes(G4VisAttributes::Invisible);
	

	G4VisAttributes *logiSSOVCVis   = new G4VisAttributes(cyanl);
	G4VisAttributes *logiCu4Vis     = new G4VisAttributes(bluel);
	G4VisAttributes *logiCu3Vis     = new G4VisAttributes(greenl);
	G4VisAttributes *logiCu2Vis     = new G4VisAttributes(cyanl);
	G4VisAttributes *logiCu1Vis     = new G4VisAttributes(brownl);
	G4VisAttributes *logiCuMCPVis   = new G4VisAttributes(lblue);
	G4VisAttributes *logiCuP1Vis    = new G4VisAttributes(lgreen);
	G4VisAttributes *logiPbP2Vis    = new G4VisAttributes(lgrey);
	G4VisAttributes *logiSourceHousingVis = new G4VisAttributes(whitel);
	G4VisAttributes *logiSourceVis  = new G4VisAttributes(redl);

	if (flagInvisible || flagOneCell) {
		logiSSOVC->SetVisAttributes(G4VisAttributes::Invisible);
		logiCu4->SetVisAttributes(G4VisAttributes::Invisible);
		logiCu3->SetVisAttributes(G4VisAttributes::Invisible);
		logiCu2->SetVisAttributes(G4VisAttributes::Invisible);
		logiCu1->SetVisAttributes(G4VisAttributes::Invisible);

		if (flagInvisible == 999) {
			logiCuMCP->SetVisAttributes(G4VisAttributes::Invisible);
			logiCuP1->SetVisAttributes(G4VisAttributes::Invisible);
			logiCuP2->SetVisAttributes(G4VisAttributes::Invisible);
			logiCuP3->SetVisAttributes(G4VisAttributes::Invisible);
			logiPbP1->SetVisAttributes(G4VisAttributes::Invisible);
			logiPbP2->SetVisAttributes(G4VisAttributes::Invisible);
		}
	} 
	else {
		logiTargetRoomVis->SetForceSolid(true);
		logiSSOVCVis->SetForceSolid(true);
		logiCu4Vis->SetForceSolid(true);
		logiCu3Vis->SetForceSolid(true);
		logiCu2Vis->SetForceSolid(true);
		logiCu1Vis->SetForceSolid(true);
		logiCuMCPVis->SetForceSolid(true);
		logiCuP1Vis->SetForceSolid(true);
		logiPbP2Vis->SetForceSolid(true);
		//logiSSOVCVis->SetForceSolid(true);
		//logiSourceVis->SetForceSolid(true);
		//logiSourceHousingVis->SetForceSolid(true);

		logiSSOVC->SetVisAttributes(logiSSOVCVis);
		//logiSSOVCIn->SetVisAttributes(logiTargetRoomVis);
		logiCu4->SetVisAttributes(logiCu4Vis);
		//logiCu4In->SetVisAttributes(logiTargetRoomVis);
		logiCu3->SetVisAttributes(logiCu3Vis);
		//logiCu3In->SetVisAttributes(logiTargetRoomVis);
		logiCu2->SetVisAttributes(logiCu2Vis);
		//logiCu2In->SetVisAttributes(logiTargetRoomVis);
		logiCu1->SetVisAttributes(logiCu1Vis);
		//logiCu1In->SetVisAttributes(logiTargetRoomVis);
		logiCuMCP->SetVisAttributes(logiCuMCPVis);
		logiCuP1->SetVisAttributes(logiCuP1Vis);
		logiCuP2->SetVisAttributes(logiCuP1Vis);
		logiCuP3->SetVisAttributes(logiCuP1Vis);
		logiPbP1->SetVisAttributes(logiPbP2Vis);
		logiPbP2->SetVisAttributes(logiPbP2Vis);
		//logiSCshield->SetVisAttributes(logiPbP2Vis);
		//logiSourceHousing->SetVisAttributes(logiSourceHousingVis);
		//logiSourceVolume->SetVisAttributes(logiSourceVis);
	}

#if G4VERSION_NUMBER < 1000
	ConstructAMoRE200_SDandField();
#endif

	/////////////////////////////////////////////////////////////////
	//  Set Region
	G4Region *crystalsRegion = new G4Region("crystals");
	//crystalsRegion->AddRootLogicalVolume(f200_logiCrystalCell);
	for(int icry = 0; icry < f200_TotCrystalNum; icry++){
		crystalsRegion->AddRootLogicalVolume(f200_logiCrystalCell[icry]);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Crystal module making function
///////////////////////////////////////////////////////////////////////////////
G4LogicalVolume *AmoreDetectorConstruction::MakeModule(G4Material *towerMat, G4Material *crystalMat, G4Material *reflectorMat, G4Material *frameMat, G4Material *frameMat1, G4Material *clampMat, G4Material *waferMat, G4Material *filmMat, G4int TowerNum, G4int ModuleNum){

	G4int CrystalIndex = TowerNum * nModuleInTower + ModuleNum;

	const AMoRE200CrystalModuleInfo &nowInfo = crystalModuleInfoList[TowerNum];
	G4double cell_h = nowInfo.fCrystalHeight;
	G4double cell_r = nowInfo.fCrystalRadius;
	G4double module_h = cell_h + bottomframe_height + topframe_height + solidBooleanTol;
	G4double module_r = cell_r + 12; 

	///_______________________________________________
	/// Crystal Module
	///-----------------------------------------------
	G4Tubs *CrystalModule = new G4Tubs("CrystalModule", 0, module_r + 0.5, module_h/2, 0, 360. * deg);
	G4LogicalVolume *resultLV = new G4LogicalVolume(CrystalModule, towerMat, "CrystalModuleLV");

	///_______________________________________________
	/// Crystal Cell 
	///-----------------------------------------------
	G4Tubs *CrystalCell = new G4Tubs("CrystalCell", 0, cell_r, cell_h/2, 0, 360. * deg);
	//G4LogicalVolume *logiCrystalCell = new G4LogicalVolume(CrystalCell, crystalMat, ("CrystalCellLV_" + to_string(CrystalIndex)).c_str());
	G4LogicalVolume *logiCrystalCell = new G4LogicalVolume(CrystalCell, crystalMat, "CrystalCellLV");
	G4VisAttributes *logiCrystalVis = new G4VisAttributes(yellow);   
	logiCrystalVis->SetVisibility(true);
	logiCrystalVis->SetForceSolid(true);
	logiCrystalCell->SetVisAttributes(logiCrystalVis);
	f200_logiCrystalCell[CrystalIndex] = logiCrystalCell;


	///_______________________________________________
	/// Heater
	///-----------------------------------------------
	G4Box *crystalHeaterBox = new G4Box("HeaterForCrystalBox", 
			crystalHeaterBoxSizeX / 2., crystalHeaterBoxSizeY / 2., crystalHeaterBoxSizeZ / 2.); 
	G4Box *crystalHeaterElectrodeBox = new G4Box("HeaterForCrystalElectrodeBox", 
			crystalHeaterElectrodeSizeX / 2., crystalHeaterElectrodeSizeY / 2., crystalHeaterElectrodeSizeZ / 2.); 
	G4LogicalVolume *logiHeater = new G4LogicalVolume(crystalHeaterBox, _SiWafer, "HeaterForCrystal_LV");
	G4LogicalVolume *logiHeaterElectrode = new G4LogicalVolume( 
			crystalHeaterElectrodeBox, _AuPd, "HeaterForCrystal_Electrode_LV");
	G4VisAttributes *logiHeaterVis = new G4VisAttributes(orange);
	logiHeaterVis->SetForceSolid(true);
	logiHeater->SetVisAttributes(logiHeaterVis);
	G4VisAttributes *logiHeaterElecVis = new G4VisAttributes(yellow);
	logiHeaterElectrode->SetVisAttributes(logiHeaterElecVis);

	G4ThreeVector nowSolidPosition = {crystalHeaterElectrodeStartX, 0,
		-crystalHeaterBoxSizeZ / 2. + crystalHeaterElectrodeSizeZ / 2.}; 
	for (G4int i = 0; i < 12; i++) 
	{    
		new G4PVPlacement(nullptr, nowSolidPosition, logiHeaterElectrode,
				"HeaterForCrystal_Electrode_PV", logiHeater, false, 0, OverlapCheck);
		nowSolidPosition[0] += crystalHeaterElectrodeSpacing;
	}  


	///_______________________________________________
	/// Araldite
	///-----------------------------------------------
	G4Box *AralditeBox = new G4Box("AralditeBox", crystalHeaterBoxSizeX / 2., crystalHeaterBoxSizeY/2., araldite_thick/2.);
	G4Tubs *AlspacerTub = new G4Tubs("Alspacer", 0, araldite_thick/2, crystalHeaterBoxSizeX/2., 0, 360.*deg);
	G4LogicalVolume *logiAraldite = new G4LogicalVolume(AralditeBox, _araldite, "AralditeLV");
	G4LogicalVolume *logiAlspacer = new G4LogicalVolume(AlspacerTub, _aluminium, "AlspacerLV");

	G4VisAttributes *logiAlspacerVis = new G4VisAttributes(grey);
	logiAlspacerVis->SetForceSolid(true);
	logiAlspacer->SetVisAttributes(logiAlspacerVis);
	G4VisAttributes *logiEpoxiVis = new G4VisAttributes(orangel);
	logiEpoxiVis->SetForceSolid(true);
	logiAraldite->SetVisAttributes(logiEpoxiVis);

	G4RotationMatrix *alspRotMtx = new G4RotationMatrix();
	alspRotMtx->rotateY(90.*deg);
	new G4PVPlacement(G4Transform3D(*alspRotMtx, {0, AralditeBox->GetXHalfLength()/2.,0}), logiAlspacer, "physAraldite" ,logiAraldite, false, 0, OverlapCheck);
	new G4PVPlacement(G4Transform3D(*alspRotMtx, {0, -AralditeBox->GetXHalfLength()/2.,0}), logiAlspacer, "physAraldite" ,logiAraldite, false, 0, OverlapCheck);



	///_______________________________________________
	/// Reflector 
	///-----------------------------------------------
	//G4double reflector_r = cell_r + reflector_thick;
	G4double reflector_rb = cell_r + 2.5;
	G4double reflector_r = cell_r + 3;
	G4double reflector_h = cell_h + 3.6;
	if ( cell_r > 53./2 ) reflector_h = cell_h + 3.2;
	G4Box *ReflectorHoleL = new G4Box("ReflectorHoleL", clamp_top_l, clamp_top_l*0.6, clamp_top_l/2.);
	G4Box *ReflectorHoleM = new G4Box("ReflectorHoleM", clamp_wafer_w, clamp_wafer_w*0.6, clamp_wafer_w/3.);
	G4Box *ReflectorHoleS = new G4Box("ReflectorHoleS", clamp_bot_h*0.5, clamp_bot_h*1.2, clamp_bot_h);
	G4Box *ReflectorHoleH = new G4Box("ReflectorHoleH", bottomhole_w/2., bottomhole_w/2, bottomhole_w);
	G4VSolid *Reflector = new G4Tubs("Reflector_tmp", reflector_r-reflector_thick, reflector_r, reflector_h/2, 0, 360. * deg);
	G4VSolid *Reflector_B = new G4Tubs("Reflector_Btmp", bottomGold_radius, reflector_rb, reflector_thick/2., 0, 360. * deg);

	G4RotationMatrix *ReflectorholeSRotMtx = new G4RotationMatrix();
	G4ThreeVector sholePos = G4ThreeVector(reflector_rb, 0, 0);
	G4ThreeVector lholePos = G4ThreeVector(reflector_r, 0, reflector_h/2.);
	G4ThreeVector hholePos = G4ThreeVector(0, -reflector_rb/2, 0);

	for (int i = 0; i < 24; i++){
		if (i==1 || i==9 || i ==17) { 
			Reflector = new G4SubtractionSolid("Reflector", Reflector, ReflectorHoleL, ReflectorholeSRotMtx, lholePos);
		}
		if (i==6 || i==14 || i==22) {
			Reflector = new G4SubtractionSolid("Reflector", Reflector, ReflectorHoleM, ReflectorholeSRotMtx, lholePos);
			Reflector_B = new G4SubtractionSolid("Reflector_B", Reflector_B, ReflectorHoleS, ReflectorholeSRotMtx, sholePos);
			ReflectorholeSRotMtx->rotateZ(60.*deg);
		}
		if (i==18){
			Reflector_B = new G4SubtractionSolid("Reflector_B", Reflector_B, ReflectorHoleH, nullptr, hholePos);
		}
		sholePos.rotateZ(360/24. * deg);
		lholePos.rotateZ(360/24. * deg);
	}

	G4LogicalVolume *logiReflector = new G4LogicalVolume(Reflector, reflectorMat, "ReflectorLV");
	G4LogicalVolume *logiReflector_B = new G4LogicalVolume(Reflector_B, reflectorMat, "ReflectorBLV");
	G4VisAttributes *logiReflectorVis = new G4VisAttributes(G4Colour(0.75, 0.0, 0.75, 0.3));
	logiReflectorVis->SetForceSolid(true);
	logiReflector->SetVisAttributes(logiReflectorVis);
	logiReflector_B->SetVisAttributes(logiReflectorVis);

	///_______________________________________________
	/// Bottom Module Frame
	///-----------------------------------------------
	G4VSolid *BottomModule     = new G4Tubs("BottomModule", module_r - Module_width, module_r, Module_thick/2., 0, 360. * deg);
	G4VSolid *BottomModuleBase = new G4Tubs("BottomModuleBase", module_r - Module_width - Module_base_thick, module_r - Module_width, bottomframe_thick/2., 0, 360. * deg);
	G4Tubs *BottomHole1 = new G4Tubs("BottomHole1", 0, bottomhole1_size/2., Module_thick * 10, 0, 360. *  deg); 
	G4Tubs *BottomHole2 = new G4Tubs("BottomHole2", 0, bottomhole2_size/2., Module_thick * 10, 0, 360. *  deg); 
	G4Box *BottomHole3  = new G4Box("BottomHole3", bottomhole3_size/2, bottomhole3_size/2, Module_thick);
	G4Box *BottomHole4  = new G4Box("BottomHole4", bottomhole_w/2., bottomhole_l/2, Module_thick);
	G4Tubs *BottomLeg1   = new G4Tubs("BottomLeg1", module_r - Module_width, module_r, leg_height/2., -360./24/2 * deg, 360./24 * deg);
	G4Tubs *BottomLeg2   = new G4Tubs("BottomLeg2", module_r - Module_width, module_r, leg_height/2., 360./24*7.5 * deg, 360./24 * deg);
	G4Tubs *BottomLeg3   = new G4Tubs("BottomLeg3", module_r - Module_width, module_r, leg_height/2., 360./24*15.5 * deg, 360./24 * deg);
	G4VSolid *BottomLeg = new G4UnionSolid("BottomLeg", BottomLeg1, BottomLeg2, nullptr, {0,0,0});
	BottomLeg = new G4UnionSolid("BottomLeg", BottomLeg, BottomLeg3, nullptr, {0,0,0});
	G4VSolid *BottomAdditional = new G4Tubs("BottomAdditional", 0, bottomadd_size, bottomframe_thick/2., 0, 360. * deg);

	G4ThreeVector holePos  = G4ThreeVector(module_r-frame_hole_depth, 0, 0);
	G4ThreeVector hole3Pos = G4ThreeVector(module_r, 0, 0);
	G4ThreeVector addPos   = G4ThreeVector(module_r - Module_width - Module_base_thick+0.5, 0, -Module_thick/2. + bottomframe_thick/2.);
	G4ThreeVector legPos   = G4ThreeVector(0, 0, -Module_thick/2.-leg_height/2.);
	hole3Pos.rotateZ(360/24. * 18 * deg);

	BottomModule = new G4UnionSolid("BottomModule", BottomModule, BottomModuleBase, nullptr, {0,0,-Module_thick/2 + bottomframe_thick/2.});
	BottomModule = new G4UnionSolid("BottomModule", BottomModule, BottomLeg, nullptr, legPos);

	G4RotationMatrix *Bottomhole4RotMtx = new G4RotationMatrix();

	for(int i = 0; i < 24; i++){
		if(i!=4 && i!=18 && i!=19 ) BottomModule = new G4SubtractionSolid("BottomModule", BottomModule, BottomHole1, nullptr, holePos);
		else if( i!=18 ) BottomModule = new G4SubtractionSolid("BottomModule", BottomModule, BottomHole2, nullptr, holePos);
		else BottomModule = new G4SubtractionSolid("BottomModule", BottomModule, BottomHole3, nullptr, hole3Pos);

		holePos.rotateZ(360/24. * deg);
		if (i==6 || i==14) { // || i==22) {
			BottomModule = new G4UnionSolid("BottomModule", BottomModule, BottomAdditional, nullptr, addPos);
			BottomModule = new G4SubtractionSolid("BottomModule", BottomModule, BottomHole4, Bottomhole4RotMtx, addPos);
			Bottomhole4RotMtx->rotateZ(60.*deg);
		}
		addPos.rotateZ(360/24. * deg);
		}
		G4LogicalVolume *logiBottomModule = new G4LogicalVolume(BottomModule, frameMat1, "BottomModuleLV");
		G4VisAttributes *logiFrameVis = new G4VisAttributes(brown);
		logiFrameVis->SetForceSolid(true);
		logiBottomModule->SetVisAttributes(logiFrameVis);

		///_______________________________________________
		/// Top Module Frame
		///-----------------------------------------------
		G4VSolid *TopModule = new G4Tubs("TopModule", module_r - Module_width, module_r, Module_thick/2., 0, 360. *deg);

		G4Tubs *TopLeg1 = new G4Tubs("TopLeg1",module_r - Module_width, module_r, topleg_height/2., -360./24/2 * deg, 360./24 * deg);
		G4Tubs *TopLeg2 = new G4Tubs("TopLeg2",module_r - Module_width, module_r, topleg_height/2., 360./24*7.5 * deg, 360./24 * deg);
		G4Tubs *TopLeg3 = new G4Tubs("TopLeg3",module_r - Module_width, module_r, topleg_height/2., 360./24*15.5 * deg, 360./24 * deg);
		G4VSolid *TopLeg = new G4UnionSolid("TopLeg", TopLeg1, TopLeg2, nullptr, {0,0,0});
		TopLeg = new G4UnionSolid("TolLeg", TopLeg, TopLeg3, nullptr, {0,0,0});

		G4Tubs *TopAdd1 = new G4Tubs("TopAdd1", module_r - Module_width, module_r, topadd_height/2., 360./24*5.5*deg, 360./24*deg);
		G4Tubs *TopAdd2 = new G4Tubs("TopAdd2", module_r - Module_width, module_r, topadd_height/2., 360./24*13.5*deg, 360./24*deg);
		G4Tubs *TopAdd3 = new G4Tubs("TopAdd3", module_r - Module_width, module_r, topadd_height/2., 360./24*21.5*deg, 360./24*deg);
		G4Tubs *TopAdd4 = new G4Tubs("TopAdd4", module_r - Module_width - topadd_length, module_r - Module_width, topadd_height/2., 360./24*5.75*deg, 360./24/2*deg);
		G4Tubs *TopAdd5 = new G4Tubs("TopAdd5", module_r - Module_width - topadd_length, module_r - Module_width, topadd_height/2., 360./24*13.75*deg, 360./24/2*deg);
		G4Tubs *TopAdd6 = new G4Tubs("TopAdd6", module_r - Module_width - topadd_length, module_r - Module_width, topadd_height/2., 360./24*21.75*deg, 360./24/2*deg);
		G4VSolid *TopAdd = new G4UnionSolid("TopAdd", TopAdd1, TopAdd2, nullptr, {0,0,0});
		TopAdd = new G4UnionSolid("TopAdd", TopAdd, TopAdd3, nullptr, {0,0,0});
		TopAdd = new G4UnionSolid("TopAdd", TopAdd, TopAdd4, nullptr, {0,0,0});
		TopAdd = new G4UnionSolid("TopAdd", TopAdd, TopAdd5, nullptr, {0,0,0});
		TopAdd = new G4UnionSolid("TopAdd", TopAdd, TopAdd6, nullptr, {0,0,0});

		G4Tubs *TopAddTub = new G4Tubs("TopAddTub", 0, tophole0_size/2., topadd_height/2., 0, 360. * deg);
		G4Tubs *TopHole1 = new G4Tubs("TopHole1", 0, tophole1_size/2, Module_thick * 10, 0, 360. * deg);
		G4Tubs *TopHole2 = new G4Tubs("TopHole2", 0, tophole2_size/2, Module_thick * 10, 0, 360. * deg);

		addPos  = G4ThreeVector(module_r - Module_width - topadd_length, 0, Module_thick/2. + topadd_height/2.);
		legPos  = G4ThreeVector(0, 0, Module_thick/2. + topleg_height/2.);
		holePos = G4ThreeVector(module_r-frame_hole_depth, 0, 0);

		TopModule = new G4UnionSolid("TopModule", TopModule, TopLeg, nullptr, legPos);
		TopModule = new G4UnionSolid("TopModule", TopModule, TopAdd, nullptr, {0, 0, Module_thick/2. + topadd_height/2.});
		for(int i = 0; i < 24; i++){
			if(i==0 || i==8 || i ==16){ // post position
				TopModule = new G4SubtractionSolid("TopModule", TopModule, TopHole2, nullptr, holePos);
			}	
			else if(i==5 || i==6 || i==14 || i==19 || i==22){ // smallest hole
				if(i==6 || i==14 || i==22) {
					TopModule = new G4UnionSolid("TopModule", TopModule, TopAddTub, nullptr, addPos);
				}
				TopModule = new G4SubtractionSolid("TopModule", TopModule, BottomHole2, nullptr, holePos);
			} 
			else if(i==18) { 
				TopModule = new G4SubtractionSolid("TopModule", TopModule, BottomHole3, nullptr, hole3Pos);
			} 
			else TopModule = new G4SubtractionSolid("TopModule", TopModule, TopHole1, nullptr, holePos);

			holePos.rotateZ(360/24. * deg);
			addPos.rotateZ(360/24. * deg);
		}
		G4LogicalVolume *logiTopModule = new G4LogicalVolume(TopModule, frameMat, "TopModuleLV");
		logiTopModule->SetVisAttributes(logiFrameVis);

		///_______________________________________________
		/// POSTs
		///-----------------------------------------------
		G4double post_h = cell_h - 0.9;
		G4VSolid *Post  = new G4Tubs("Post", post_bolt_r, post_radius-solidBooleanTol*2, (post_h-post_top_h)/2., 0, 360. * deg); 
		G4Tubs *PostTop = new G4Tubs("PostTop1", post_bolt_r, post_top_r, post_top_h/2., 0, 360. * deg);
		Post = new G4UnionSolid("Post", Post, PostTop, nullptr, {0,0,(post_h-post_top_h)/2+post_top_h/2});
		G4LogicalVolume *logiPost = new G4LogicalVolume(Post, frameMat, "PostLV");
		logiPost->SetVisAttributes(logiFrameVis);

		///_______________________________________________
		// Top module teflon sheets
		///-----------------------------------------------
		G4double teflonSheet_thick = 0.15 * mm;
		G4Tubs *teflon1 = new G4Tubs("teflon1", module_r - Module_width - topadd_length - teflonSheet_thick, module_r - Module_width - 1 * mm , topadd_height/2. + teflonSheet_thick, 360./24*5.725*deg, 360./24/1.8*deg);
		G4Tubs *teflon2 = new G4Tubs("teflon2", module_r - Module_width - topadd_length - teflonSheet_thick, module_r - Module_width - 1 * mm , topadd_height/2. + teflonSheet_thick, 360./24*13.725*deg, 360./24/1.8*deg);
		G4Tubs *teflon3 = new G4Tubs("teflon3", module_r - Module_width - topadd_length - teflonSheet_thick, module_r - Module_width - 1 * mm , topadd_height/2. + teflonSheet_thick, 360./24*21.725*deg, 360./24/1.8*deg);
		G4Tubs *teflon1_1 = new G4Tubs("teflon1_1", module_r - Module_width - topadd_length - teflonSheet_thick , module_r - Module_width, topadd_height/2. + solidBooleanTol, 360./24*5.74*deg, 360./24/1.9*deg);
		G4Tubs *teflon2_1 = new G4Tubs("teflon2_1", module_r - Module_width - topadd_length - teflonSheet_thick , module_r - Module_width, topadd_height/2. + solidBooleanTol, 360./24*13.74*deg, 360./24/1.9*deg);
		G4Tubs *teflon3_1 = new G4Tubs("teflon3_1", module_r - Module_width - topadd_length - teflonSheet_thick , module_r - Module_width, topadd_height/2. + solidBooleanTol, 360./24*21.74*deg, 360./24/1.9*deg);
		G4VSolid *teflon4 = new G4Tubs("teflon4", 0, tophole0_size/2. + teflonSheet_thick, topadd_height/2. + teflonSheet_thick, 0, 360. * deg);
		G4Tubs *teflon4_1 = new G4Tubs("teflon4_1", 0, tophole0_size/2. + solidBooleanTol, topadd_height/2. + solidBooleanTol, 0, 360. * deg);
		teflon4 = new G4SubtractionSolid("teflon4", teflon4, teflon4_1, nullptr, {0,0,0}); 

		G4VSolid *TeflonSheet = new G4UnionSolid("TeflonSheet", teflon1, teflon2, nullptr, {0,0,0});
		TeflonSheet = new G4UnionSolid("TeflonSheet", TeflonSheet, teflon3, nullptr, {0,0,0});
		G4VSolid *TeflonSheet_1 = new G4UnionSolid("TeflonSheet_1", teflon1_1, teflon2_1, nullptr,{0,0,0});
		TeflonSheet_1 = new G4UnionSolid("TeflonSheet_1", TeflonSheet_1, teflon3_1, nullptr, {0,0,0});

		addPos  = G4ThreeVector(module_r - Module_width - topadd_length, 0, 0);
		addPos.rotateZ(360/24. * 6 * deg);
		for(int i = 0; i < 3; i ++){
			TeflonSheet = new G4UnionSolid("TeflonSheet", TeflonSheet, teflon4, nullptr, addPos);
			addPos.rotateZ(120. * deg);
		}
		TeflonSheet = new G4SubtractionSolid("TeflonSheet", TeflonSheet, TeflonSheet_1, nullptr, {0,0,0});

		G4LogicalVolume *logiTeflonSheet = new G4LogicalVolume(TeflonSheet, clampMat, "TeflonSheetLV");
		G4VisAttributes *logiClampVis = new G4VisAttributes(white);
		logiClampVis->SetForceSolid(true);
		logiTeflonSheet->SetVisAttributes(logiClampVis);

		///_______________________________________________
		/// Clamps
		///-----------------------------------------------
		G4VSolid *ClampWafer   = new G4Box("ClampWafer", clamp_wafer_w/2., clamp_wafer_l/2, post_top_h/2.);
		G4Tubs *ClampWaferHole = new G4Tubs("ClampWaferHole", 0, clamp_wafer_holesize/2, post_top_h, 0, 360. * deg);
		ClampWafer = new G4SubtractionSolid("ClampWafer", ClampWafer, ClampWaferHole, nullptr, {0, 12.5/2-3, 0});
		G4LogicalVolume *logiClampWafer = new G4LogicalVolume(ClampWafer, clampMat, "ClampWaferLV");
		logiClampWafer->SetVisAttributes(logiClampVis);

		G4VSolid *ClampTop  = new G4Box("ClampTop", clamp_top_w/2, clamp_top_l/2., post_top_h/2.);
		G4Box *ClampTopBox1 = new G4Box("ClampTopBox1", clamp_top_box1_w/2, clamp_top_box1_l/2., clamp_top_box_h/2.);
		G4Box *ClampTopBox2 = new G4Box("ClampTopBox2", clamp_top_box_h/2., clamp_top_box_h/2., post_top_h/2.);
		G4Tubs *ClampTopHole = new G4Tubs("ClampTopHole", 0, clamp_top_holesize/2, post_top_h, 0, 360.*deg);
		ClampTop = new G4UnionSolid("ClampTop", ClampTop, ClampTopBox1, nullptr, {-clamp_top_w/2+clamp_top_box1_w/2-AmoreDetectorStaticInfo::solidBooleanTol*2, 0, post_top_h/2.+clamp_top_box_h/2.});
		ClampTop = new G4UnionSolid("ClampTop", ClampTop, ClampTopBox2, nullptr, {-clamp_top_w/2-clamp_top_box_h/2-AmoreDetectorStaticInfo::solidBooleanTol*2., 0, post_top_h+clamp_thick+AmoreDetectorStaticInfo::solidBooleanTol});
		ClampTop = new G4SubtractionSolid("ClampTop", ClampTop, ClampTopHole, nullptr, {clamp_top_w/2-clamp_top_holepos,0,0});
		G4LogicalVolume *logiClampTop = new G4LogicalVolume(ClampTop, clampMat, "ClampTopLV");
		logiClampTop->SetVisAttributes(logiClampVis);

		G4VSolid *ClampBottom = new G4Box("ClampBottom", bottomhole_w/2., bottomhole_l/2.-solidBooleanTol, clamp_bot_h/2.);
		G4VSolid *ClampBottom2 = new G4Box("ClampBottom2", bottomhole_w/2., bottomhole_l/2.-solidBooleanTol, (clamp_bot_h/2.+clamp_thick/2.)/2.);
		G4VSolid *ClampBot_tub = new G4Tubs("ClampBot_tub", 0, bottomhole_w/2., clamp_thick/2., 180.* deg, 360.*deg);
		ClampBottom = new G4UnionSolid("ClampBottom", ClampBottom, ClampBot_tub, nullptr, {0, -bottomhole_l/2, 0});
		ClampBottom2 = new G4UnionSolid("ClampBottom2", ClampBottom2, ClampBot_tub, nullptr, {0, -bottomhole_l/2, -(clamp_bot_h/2.-clamp_thick/2.)/2.});
		G4LogicalVolume *logiClampBottom = new G4LogicalVolume(ClampBottom, clampMat, "ClampBottomLV");
		G4LogicalVolume *logiClampBottom2 = new G4LogicalVolume(ClampBottom2, clampMat, "ClampBottomLV");
		logiClampBottom->SetVisAttributes(logiClampVis);
		logiClampBottom2->SetVisAttributes(logiClampVis);

		///_______________________________________________
		/// Wafer
		///-----------------------------------------------
		G4double wafer_radius_adj = wafer_radius;
		if (cell_r > 50./2) wafer_radius_adj = 60.0 / 2.;

		G4Tubs *Wafer = new G4Tubs("Wafer", 0, wafer_radius_adj, wafer_thick/2., 0, 360. *deg);
		G4LogicalVolume *logiWafer = new G4LogicalVolume(Wafer, waferMat, "WaferLV");
		G4VisAttributes *logiWaferVis = new G4VisAttributes(bluel);
		logiWaferVis->SetForceSolid(true);
		logiWafer->SetVisAttributes(logiWaferVis);

		///_______________________________________________
		/// Gold film
		///_______________________________________________
		G4Tubs *UpperGoldFilm = new G4Tubs("UpperGoldFilm", 0, upperGold_radius, upperGold_thick/2., 0, 360. * deg);
		G4LogicalVolume *logiUpperGoldFilm = new G4LogicalVolume(UpperGoldFilm, filmMat, "UpperGoldFilmLV");

		G4Tubs *BottomGoldFilm = new G4Tubs("BottomGoldFilm", 0, bottomGold_radius, bottomGold_thick/2., 0, 360. * deg);
		G4LogicalVolume *logiBottomGoldFilm = new G4LogicalVolume(BottomGoldFilm, filmMat, "BottomGoldFilmLV");

		G4VisAttributes *logiGoldFilmVis = new G4VisAttributes(G4Colour(1, 0.87, 0));
		logiGoldFilmVis->SetForceSolid(true);
		logiUpperGoldFilm->SetVisAttributes(logiGoldFilmVis);
		logiBottomGoldFilm->SetVisAttributes(logiGoldFilmVis);

		///_______________________________________________
		/// light detector
		///_______________________________________________
		G4double ldet_box1_l = cell_r * 2 + 16;
		G4double ldet_box2_l = cell_r * 2 - 10;
		G4Box *LightBox1 = new G4Box("LightBox1", light_box1_w/2., ldet_box1_l/2., light_box_thick/2.);
		G4Box *LightBox2 = new G4Box("LightBox2", light_box2_w/2., ldet_box2_l/2., light_box_thick/2.);
		G4Box *LightBox3 = new G4Box("LightBox3", light_box3_w/2., light_box3_l/2., light_box_thick);
		G4VSolid *LightDetector = new G4SubtractionSolid("LightDetector", LightBox1, TopHole1, nullptr, {0, 0, 0} );
		LightDetector = new G4SubtractionSolid("LightDetector", LightDetector, LightBox3, nullptr, {0, -ldet_box1_l/2., 0});
		LightDetector = new G4UnionSolid("LightDetector", LightDetector, LightBox2, nullptr, {light_box1_w/2+light_box2_w/2,-ldet_box1_l/2.+ldet_box2_l/2., 0});
		G4LogicalVolume *logiLightDetector = new G4LogicalVolume(LightDetector, frameMat, "LightDetectorLV");
		logiLightDetector->SetVisAttributes(logiFrameVis);

		///_______________________________________________
		/// PCB & Solder for light detector
		///_______________________________________________
		G4Tubs *PCBtubs = new G4Tubs("PCBtubs", tophole1_size/2., pcbtub_radius, pcb_thick/2., 0, 210.*deg);
		G4LogicalVolume *logiPCBlight = new G4LogicalVolume(PCBtubs, _copper2, "PCBlightLV");
		G4VisAttributes *logiPCBVis = new G4VisAttributes(orange);
		logiPCBVis->SetForceSolid(true);
		logiPCBlight->SetVisAttributes(logiPCBVis);

		G4LogicalVolume *logiSolderLight = new G4LogicalVolume(PCBtubs, _solder, "SolderLightLV");
		G4VisAttributes *logiSolderVis = new G4VisAttributes(greyl);
		logiSolderVis->SetForceSolid(true);
		logiSolderLight->SetVisAttributes(logiSolderVis);


		///_______________________________________________
		/// Stycast for light sensor
		///_______________________________________________
		G4Box *StycastLight = new G4Box("StycastLight", lightsensor_size/2., lightsensor_size/2., pcb_thick/2.);
		G4LogicalVolume *logiStycastLight = new G4LogicalVolume(StycastLight, _stycast, "StycastLightLV");
		logiStycastLight->SetVisAttributes(logiEpoxiVis);


		///_______________________________________________
		/// heat detector
		///_______________________________________________
		G4double hdet_box1_l = cell_r * 2 - 10;
		G4double hdet_box2_l = cell_r * 2 + 10;
		G4Box *HeatBox1 = new G4Box("HeatBox1", heat_box1_w/2., hdet_box1_l/2., heat_box1_thick/2.);
		G4Box *HeatBox2 = new G4Box("HeatBox2", heat_box2_w/2., hdet_box2_l/2., heat_box2_thick/2.);
		G4VSolid *HeatDetector = new G4SubtractionSolid("HeatDetector", HeatBox1, HeatBox2, nullptr, {0, -hdet_box2_l/3. -hdet_box1_l/2., heat_box1_thick/2.-heat_box2_thick/2.});
		HeatDetector = new G4UnionSolid("HeatDetector", HeatDetector, HeatBox2, nullptr, {heat_box1_w/2+heat_box2_w/2, -hdet_box1_l/2.+hdet_box2_l/2, -heat_box1_thick/2.+heat_box2_thick/2.});
		G4LogicalVolume *logiHeatDetector = new G4LogicalVolume(HeatDetector, frameMat1, "HeatDetectorLV");
		logiHeatDetector->SetVisAttributes(logiFrameVis);


		///_______________________________________________
		/// PCB for heat detector 
		///_______________________________________________
		G4Box *PCBbox = new G4Box("PCBbox", pcbbox_size/2., pcbbox_size/2., pcb_thick/2.);
		G4LogicalVolume *logiPCBheat = new G4LogicalVolume(PCBbox, _copper2, "PCBheatLV"); 
		logiPCBheat->SetVisAttributes(logiPCBVis);
		G4LogicalVolume *logiSolderHeat = new G4LogicalVolume(PCBbox, _solder, "SolderHeatLV");
		logiSolderHeat->SetVisAttributes(logiSolderVis);
		

		///_______________________________________________
		/// Stycast for heat sensor
		///_______________________________________________
		G4Box *StycastHeat = new G4Box("StycastHeat", heatsensor_size/2., heatsensor_size/2., pcb_thick/2.);
		G4LogicalVolume *logiStycastHeat = new G4LogicalVolume(StycastHeat, _stycast, "StycastHeatLV");
		logiStycastHeat->SetVisAttributes(logiEpoxiVis);


		///_______________________________________________
		/// Screws
		///_______________________________________________
		G4Tubs *Screw_M5_head = new G4Tubs("M5_head",0, M5_size/2., M5_height/2., 0, 360. * deg);
		G4Tubs *Screw_M5_12 = new G4Tubs("M5_12", 0, post_bolt_r, M5_12_height/2., 0, 360. * deg);
		//G4Tubs *Screw_M5_15 = new G4Tubs("M5_15", 0, post_bolt_r, M5_15_height/2., 0, 360. * deg);
		G4Tubs *Screw_M4_head = new G4Tubs("M4_head", 0, M4_size/2., M4_height/2., 0, 360. * deg);
		G4Tubs *Screw_M4_6 = new G4Tubs("M4_6", 0, tophole0_size/2., M4_6_height/2., 0, 360. * deg);
		G4Tubs *Screw_M4_4 = new G4Tubs("M4_4", 0, tophole0_size/2., M4_4_height/2., 0, 360. * deg);

		G4VSolid *BoltsPost = new G4UnionSolid("BoltsPost", Screw_M5_head, Screw_M5_12, nullptr, {0,0,-M5_height/2.-M5_12_height/2.});
		G4LogicalVolume *logiBoltsPost = new G4LogicalVolume(BoltsPost, frameMat, "BoltsPostLV");
		G4VisAttributes *logiBoltsVis = new G4VisAttributes(greyl);
		logiBoltsVis->SetForceSolid(true);
		logiBoltsPost->SetVisAttributes(logiBoltsVis);

		G4VSolid *BoltsClamp = new G4UnionSolid("BoltsClamp", Screw_M4_head, Screw_M4_6, nullptr, {0, 0, -M4_height/2. -M4_6_height/2.});
		G4LogicalVolume *logiBoltsClamp = new G4LogicalVolume(BoltsClamp, frameMat, "BoltsClampLV");
		logiBoltsClamp->SetVisAttributes(logiBoltsVis);

		G4VSolid *BoltsLight = new G4SubtractionSolid("BoltsLight", Screw_M4_head, Screw_M5_head, nullptr, {0,0, -M4_height/2.});
		G4LogicalVolume *logiBoltsLight = new G4LogicalVolume(BoltsLight, frameMat, "BoltsLightLV");
		G4LogicalVolume *logiBoltsHeat  = new G4LogicalVolume(Screw_M4_head, frameMat, "BoltsHeatLV");
		G4LogicalVolume *logiBoltsBody  = new G4LogicalVolume(Screw_M4_4, frameMat, "BoltsBodyLV");
		logiBoltsLight->SetVisAttributes(logiBoltsVis);
		logiBoltsHeat->SetVisAttributes(logiBoltsVis);
		logiBoltsBody->SetVisAttributes(logiBoltsVis);


		///_______________________________________________
		/// POSITIONING of the module componants
		///-----------------------------------------------
		G4int CrystalOffset = 1.;
		//G4int CrystalIndex = TowerNum * nModuleInTower + ModuleNum;
		G4ThreeVector cellPos = G4ThreeVector(0, 0, module_h/2 - topleg_height - Module_thick - cell_h/2. + CrystalOffset); 

		/// Crystal .................................
		new G4PVPlacement(0, cellPos, logiCrystalCell, ("physCrystalCell" + to_string(TowerNum) + "_" + to_string(ModuleNum)).c_str(), resultLV, false, CrystalIndex, OverlapCheck);
		

		/// Araldite ................................
		G4ThreeVector epoxyPos = cellPos;
		epoxyPos[2] += -cell_h/2. -araldite_thick/2.- bottomGold_thick;
		epoxyPos[1] += -reflector_rb/2;
		new G4PVPlacement(0, epoxyPos, logiAraldite, ("physAraldite" +to_string(TowerNum) + "_" + to_string(ModuleNum)).c_str(), resultLV, false, 0, OverlapCheck);


		/// Heater ..................................
		epoxyPos[2] += -araldite_thick/2. - crystalHeaterBoxSizeZ/2.;
		new G4PVPlacement(0, epoxyPos, logiHeater, ("physHeater" + to_string(TowerNum) + "_" + to_string(ModuleNum)).c_str(), resultLV, false, 0, OverlapCheck);


		/// Reflector ...............................
		G4ThreeVector reflectorPos = cellPos;
		reflectorPos[2] += -cell_h/2. + reflector_h/2. - reflector_gapz + solidBooleanTol; 
		new G4PVPlacement(0, reflectorPos, logiReflector, ("physReflector" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		reflectorPos[2] += -reflector_h/2. - reflector_thick/2.;
		new G4PVPlacement(0, reflectorPos, logiReflector_B, ("physReflectorBottom" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// top Frame...............................
		cellPos[2] += Module_thick/2. + cell_h/2. - CrystalOffset;
		G4ThreeVector postPos = cellPos;
		G4ThreeVector teflonSheetPos = cellPos;
		G4ThreeVector clamptopPos = cellPos;
		G4ThreeVector clampPos = cellPos + G4ThreeVector(0, 0, topadd_height + Module_thick);
		G4ThreeVector boltsPos = cellPos;
		G4ThreeVector boltsLightPos = cellPos;
		new G4PVPlacement(0, cellPos, logiTopModule, ("physTopModule" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Bottom Frame .............................
		//reflectorPos[2] += -Module_thick + CrystalOffset;
		reflectorPos[2] += -Module_base_thick + CrystalOffset;
		G4ThreeVector clampbotPos = reflectorPos;
		G4ThreeVector boltsBottomPos = reflectorPos;
		G4ThreeVector boltsHeatPos = reflectorPos;
		new G4PVPlacement(0, reflectorPos, logiBottomModule, ("physBottomModule" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Wafer ..................................
		G4ThreeVector waferPos = clampPos;
		waferPos[2] += -post_top_h/2. + wafer_thick;
		new G4PVPlacement(0, waferPos, logiWafer, ("physWafer"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Upper Gold Film ........................
		waferPos[0] = light_box1_w/2.;
		waferPos[2] += wafer_thick/2. + upperGold_thick/2.;
		new G4PVPlacement(0, waferPos, logiUpperGoldFilm, ("physUpperGold"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		waferPos.rotateZ(120.*deg);
		new G4PVPlacement(0, waferPos, logiUpperGoldFilm, ("physUpperGold"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		waferPos.rotateZ(120.*deg);
		new G4PVPlacement(0, waferPos, logiUpperGoldFilm, ("physUpperGold"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// LightDetector ...........................
		waferPos[0] = light_box1_w/2.;
		waferPos[1] = 0;
		waferPos[2] += upperGold_thick/2. + light_box_thick/2.;
		G4RotationMatrix *ldetRotMtx = new G4RotationMatrix();
		ldetRotMtx->rotateZ(-15.*deg);
		new G4PVPlacement(G4Transform3D(*ldetRotMtx, waferPos), logiLightDetector, ("physLightDetector"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		//new G4PVPlacement(G4Transform3D(*ldetRotMtx, waferPos), logiLightDetector, "physLightDetector", resultLV, false, 0, OverlapCheck);

		/// PCB solder for light ............................
		G4ThreeVector PCBpos = waferPos;
		G4RotationMatrix *pcbRotMtx = new G4RotationMatrix();
		pcbRotMtx->rotateZ(215*deg);
		PCBpos[2] += light_box_thick/2. + pcb_thick/2.;
		new G4PVPlacement(G4Transform3D(*pcbRotMtx, PCBpos), logiSolderLight, ("physSolderLight"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		/// PCB for light ............................
		PCBpos[2] += pcb_thick;
		new G4PVPlacement(G4Transform3D(*pcbRotMtx,PCBpos), logiPCBlight, ("physPCBlight"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Stycast for light ............................
		PCBpos[0] -= tophole1_size/3;
		PCBpos[1] += tophole1_size*5/4;
		PCBpos[2] += pcb_thick;
		new G4PVPlacement(G4Transform3D(*ldetRotMtx, PCBpos), logiStycastLight, ("physStycastLight"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Bottom Gold Film .........................
		G4ThreeVector goldPos = cellPos;
		goldPos[2] -= bottomGold_thick/2. + cell_h - CrystalOffset + Module_thick/2.;
		new G4PVPlacement(0, goldPos, logiBottomGoldFilm, ("physBottomGold"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		/// Heat Detector ............................
		reflectorPos[0] += heat_box1_w/2.; 
		reflectorPos[1] += -module_r + hdet_box1_l/2. + bottomhole3_size*3/4.;
		reflectorPos[2] -= Module_thick/2.;
		new G4PVPlacement(G4Transform3D(*ldetRotMtx, reflectorPos), logiHeatDetector, ("physHeatDetector"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// PCB solder for heat .............................
		PCBpos = reflectorPos;
		PCBpos[2] += -Module_thick/2. - pcb_thick/2.;
		new G4PVPlacement(G4Transform3D(*ldetRotMtx, PCBpos), logiSolderHeat, ("physSolderHeat"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		/// PCB for heat .............................
		PCBpos[2] -= pcb_thick;
		new G4PVPlacement(G4Transform3D(*ldetRotMtx, PCBpos), logiPCBheat, ("physPCBheat"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Stycast for heat .............................
		PCBpos[1] -= pcbbox_size/2. + heatsensor_size/2.;
		PCBpos[2] -= pcb_thick;
		new G4PVPlacement(G4Transform3D(*ldetRotMtx, PCBpos), logiStycastHeat, ("physStycastHeat"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Teflon Sheet for top module ...............
		teflonSheetPos[2] += Module_thick/2. + topadd_height/2.;
		new G4PVPlacement(0, teflonSheetPos, logiTeflonSheet, ("physTeflonSheet"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// POST, clamp, bolts...........................
		/// Post ///
		postPos[0] += module_r-frame_hole_depth;
		postPos[2] -= Module_thick/2. + post_top_h + (post_h-post_top_h)/2.;
		postPos.rotateZ(360./24*deg);

		/// Clamp top ///
		clamptopPos[0] += -clamp_top_w/2.+4.5+module_r-frame_hole_depth;
		clamptopPos[2] += -Module_thick/2. - post_top_h/2.;
		G4RotationMatrix *clamptopMtx = new G4RotationMatrix();
		clamptopPos.rotateZ(360./24*deg);
		clamptopMtx->rotateZ(360./24*deg);

		/// Clamp wafer ///
		clampPos[1] += -clamp_wafer_l/2 + 3 + module_r-frame_hole_depth;
		clampPos[2] += wafer_thick * 2;
		G4RotationMatrix *clampRotMtx = new G4RotationMatrix();

		/// Clamp bottom ///
		G4ThreeVector clampbot2Pos = clampbotPos;
		clampbotPos[1] += module_r - Module_width - Module_base_thick+0.5;
		clampbotPos[2] += -Module_thick/2. + bottomframe_thick + clamp_thick/2.;
		clampbot2Pos[1] += module_r - Module_width - Module_base_thick+0.5;
		clampbot2Pos[2] += -Module_thick/2. + bottomframe_thick + (clamp_bot_h/2. + clamp_thick/2.)/2.;
		G4RotationMatrix *clampbotMtx = new G4RotationMatrix();

		/// Bolts post top ///
		boltsPos[0] += module_r-frame_hole_depth;
		boltsPos[2] += M5_height/2. + Module_thick/2.;
		boltsPos.rotateZ(360/24. * deg);

		/// Bolts post bottom ///
		boltsBottomPos[0] += module_r-frame_hole_depth;
		boltsBottomPos[2] -= M5_height/2. + Module_thick/2.;
		boltsBottomPos.rotateZ(360/24. * deg);
		G4RotationMatrix *boltsRotMtx = new G4RotationMatrix();
		boltsRotMtx->rotateX(180.*deg);

		/// Bolts clamp ///
		G4ThreeVector boltsClampPos = G4ThreeVector(0,0,0);
		boltsClampPos     = clampPos;
		boltsClampPos[1]  = module_r-frame_hole_depth;
		boltsClampPos[2] += post_top_h/2. + M4_height/2.;

		for(int j = 0; j < 3; j++){
			new G4PVPlacement(0, postPos, logiPost, ("physPost"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			postPos.rotateZ(120.*deg);

			new G4PVPlacement(G4Transform3D(*clamptopMtx, clamptopPos), logiClampTop, ("physClampTop"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			clamptopPos.rotateZ(120.*deg);
			clamptopMtx->rotateZ(120.*deg);

			new G4PVPlacement(G4Transform3D(*clampRotMtx, clampPos), logiClampWafer, ("physClampWafer"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			clampPos.rotateZ(120.*deg);
			clampRotMtx->rotateZ(120.*deg);

			if(j < 2) {
				new G4PVPlacement(G4Transform3D(*clampbotMtx, clampbotPos), logiClampBottom, ("physClampBottom"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			} 
			else new G4PVPlacement(G4Transform3D(*clampbotMtx, clampbot2Pos), logiClampBottom2, ("physClampBottom"+to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			clampbotPos.rotateZ(120.*deg);
			clampbot2Pos.rotateZ(120.*deg);
			clampbotMtx->rotateZ(120.*deg);

			new G4PVPlacement(0, boltsPos, logiBoltsPost, ("physBolts_Post" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			boltsPos.rotateZ(120. * deg);

			new G4PVPlacement(G4Transform3D(*boltsRotMtx, boltsBottomPos), logiBoltsPost, ("physBolts_Post" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			boltsBottomPos.rotateZ(120. * deg);

			new G4PVPlacement(0, boltsClampPos,logiBoltsClamp, ("physBolts_Clamp" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
			boltsClampPos.rotateZ(120.*deg);
		}

		/// Bolts Light detector.................
		boltsLightPos[0] += module_r-frame_hole_depth;
		boltsLightPos[2] += -M4_4_height/2. + Module_thick/2.;
		boltsLightPos.rotateZ(360/24.* 4 * deg);
		new G4PVPlacement(0, boltsLightPos, logiBoltsBody, ("physBolts_Body" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		boltsLightPos.rotateZ(360/24.* 15 * deg);
		new G4PVPlacement(0, boltsLightPos, logiBoltsBody, ("physBolts_Body" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		boltsLightPos[2] = waferPos[2] + light_box_thick/2. + (M4_height-M5_height/2.)/2.;
		boltsLightPos.rotateZ(360/24.* 9 * deg);
		new G4PVPlacement(0, boltsLightPos, logiBoltsLight, ("physBolts_Light" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		boltsLightPos.rotateZ(360/24.* 15 * deg);
		new G4PVPlacement(0, boltsLightPos, logiBoltsLight, ("physBolts_Light" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		/// Bolts Heat detector............................
		boltsHeatPos[0] += module_r-frame_hole_depth;
		boltsHeatPos[2] += -Module_thick/2. + M4_4_height/2.; 
		boltsHeatPos.rotateZ(360/24.*3*deg);
		new G4PVPlacement(0, boltsHeatPos, logiBoltsBody, ("physBolts_Body" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		boltsHeatPos.rotateZ(360/24.*16*deg);
		new G4PVPlacement(0, boltsHeatPos, logiBoltsBody, ("physBolts_Body" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		boltsHeatPos[2] -= M4_4_height/2. + heat_box2_thick + M4_height/2.;
		boltsHeatPos.rotateZ(360/24.* 8 * deg);
		new G4PVPlacement(0, boltsHeatPos, logiBoltsHeat, ("physBolts_Heat" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);
		boltsHeatPos.rotateZ(360/24.* 16 * deg);
		new G4PVPlacement(0, boltsHeatPos, logiBoltsHeat, ("physBolts_Heat" + to_string(CrystalIndex)).c_str(), resultLV, false, 0, OverlapCheck);

		return resultLV;
	}

	/////////////////////////////////////////////////////////////////////////
	// Crystal information from Hyejin at Oct 2021.
	// Crystal Type name, Crystal radius, crystal height, tower position x, y
	const AMoRE200CrystalModuleInfo AmoreDetectorStaticInfo::AMoRE_200::crystalModuleInfoList[maxNumOfTower] = {

		{"LMO6_0", 6.0 * cm / 2, 6.0 * cm,   Module_radius+Module_gap, 0},
		{"LMO5_0", 5.0 * cm / 2, 5.0 * cm, -(Module_radius+Module_gap),0},

		{"LMO5_1_0", 5.0 * cm / 2, 5.0 * cm,  3*(Module_radius+Module_gap),  0},
		{"LMO5_1_1", 5.0 * cm / 2, 5.0 * cm,  2*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO5_1_2", 5.0 * cm / 2, 5.0 * cm,  0,                            -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO5_1_3", 5.0 * cm / 2, 5.0 * cm, -2*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO5_1_4", 5.0 * cm / 2, 5.0 * cm, -3*(Module_radius+Module_gap),  0},
		{"LMO5_1_5", 5.0 * cm / 2, 5.0 * cm, -2*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO5_1_6", 5.0 * cm / 2, 5.0 * cm,  0,                             2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO5_1_7", 5.0 * cm / 2, 5.0 * cm,  2*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},

		{"LMO6_2_0" , 5.0 * cm / 2, 5.0 * cm,  5*(Module_radius+Module_gap),  0},
		{"LMO6_2_1" , 5.0 * cm / 2, 5.0 * cm,  4*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_2" , 5.0 * cm / 2, 5.0 * cm,  3*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_3" , 6.0 * cm / 2, 6.0 * cm,     Module_radius+Module_gap,  -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_4" , 6.0 * cm / 2, 6.0 * cm,   -(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_5" , 5.0 * cm / 2, 5.0 * cm, -3*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_6" , 5.0 * cm / 2, 5.0 * cm, -4*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_7" , 5.0 * cm / 2, 5.0 * cm, -5*(Module_radius+Module_gap),  0},
		{"LMO6_2_8" , 5.0 * cm / 2, 5.0 * cm, -4*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_9" , 5.0 * cm / 2, 5.0 * cm, -3*(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_10", 6.0 * cm / 2, 6.0 * cm,   -(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_11", 6.0 * cm / 2, 6.0 * cm,     Module_radius+Module_gap,   4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_12", 5.0 * cm / 2, 5.0 * cm,  3*(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_2_13", 5.0 * cm / 2, 5.0 * cm,  4*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},

		{"LMO6_3_0" , 6.0 * cm / 2, 6.0 * cm,  7*(Module_radius+Module_gap),  0},
		{"LMO6_3_1" , 6.0 * cm / 2, 6.0 * cm,  6*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_2" , 6.0 * cm / 2, 6.0 * cm,  5*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_3" , 6.0 * cm / 2, 6.0 * cm,  4*(Module_radius+Module_gap), -6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_4" , 6.0 * cm / 2, 6.0 * cm,  2*(Module_radius+Module_gap), -6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_5" , 6.0 * cm / 2, 6.0 * cm,  0,                            -6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_6" , 6.0 * cm / 2, 6.0 * cm, -2*(Module_radius+Module_gap), -6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_7" , 6.0 * cm / 2, 6.0 * cm, -4*(Module_radius+Module_gap), -6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_8" , 6.0 * cm / 2, 6.0 * cm, -5*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_9" , 6.0 * cm / 2, 6.0 * cm, -6*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_10", 6.0 * cm / 2, 6.0 * cm, -7*(Module_radius+Module_gap),  0},
		{"LMO6_3_11", 6.0 * cm / 2, 6.0 * cm, -6*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_12", 6.0 * cm / 2, 6.0 * cm, -5*(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_13", 6.0 * cm / 2, 6.0 * cm, -4*(Module_radius+Module_gap),  6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_14", 6.0 * cm / 2, 6.0 * cm, -2*(Module_radius+Module_gap),  6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_15", 6.0 * cm / 2, 6.0 * cm,  0,                             6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_16", 6.0 * cm / 2, 6.0 * cm,  2*(Module_radius+Module_gap),  6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_17", 6.0 * cm / 2, 6.0 * cm,  4*(Module_radius+Module_gap),  6*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_18", 6.0 * cm / 2, 6.0 * cm,  5*(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_3_19", 6.0 * cm / 2, 6.0 * cm,  6*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},

		{"LMO6_4_0" , 6.0 * cm / 2, 6.0 * cm,  7*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_4_1", 6.0 * cm / 2, 6.0 * cm, -7*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		{"LMO6_4_2", 6.0 * cm / 2, 6.0 * cm,   -(Module_radius+Module_gap),  8*(Module_radius+Module_gap)*sin(60*deg)}
		//{"LMO6_4_3", 6.0 * cm / 2, 6.0 * cm,     Module_radius+Module_gap ,  8*(Module_radius+Module_gap)*sin(60*deg)}

		/* for 70 module
		   {"LMO6_4_0" , 6.0 * cm / 2, 6.0 * cm,  9*(Module_radius+Module_gap),  0},
		   {"LMO6_4_1" , 6.0 * cm / 2, 6.0 * cm,  8*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_2" , 6.0 * cm / 2, 6.0 * cm,  7*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_3" , 6.0 * cm / 2, 6.0 * cm,  6*(Module_radius+Module_gap), -6*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_4" , 6.0 * cm / 2, 6.0 * cm,  5*(Module_radius+Module_gap), -8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_5" , 6.0 * cm / 2, 6.0 * cm,  3*(Module_radius+Module_gap), -8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_6" , 6.0 * cm / 2, 6.0 * cm,     Module_radius+Module_gap , -8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_7" , 6.0 * cm / 2, 6.0 * cm,   -(Module_radius+Module_gap), -8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_8" , 6.0 * cm / 2, 6.0 * cm, -3*(Module_radius+Module_gap), -8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_9" , 6.0 * cm / 2, 6.0 * cm, -5*(Module_radius+Module_gap), -8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_10", 6.0 * cm / 2, 6.0 * cm, -6*(Module_radius+Module_gap), -6*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_11", 6.0 * cm / 2, 6.0 * cm, -7*(Module_radius+Module_gap), -4*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_12", 6.0 * cm / 2, 6.0 * cm, -8*(Module_radius+Module_gap), -2*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_13", 6.0 * cm / 2, 6.0 * cm, -9*(Module_radius+Module_gap),  0},
		   {"LMO6_4_14", 6.0 * cm / 2, 6.0 * cm, -8*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_15", 6.0 * cm / 2, 6.0 * cm, -7*(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_16", 6.0 * cm / 2, 6.0 * cm, -6*(Module_radius+Module_gap),  6*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_17", 6.0 * cm / 2, 6.0 * cm, -5*(Module_radius+Module_gap),  8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_18", 6.0 * cm / 2, 6.0 * cm, -3*(Module_radius+Module_gap),  8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_19", 6.0 * cm / 2, 6.0 * cm,   -(Module_radius+Module_gap),  8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_20", 6.0 * cm / 2, 6.0 * cm,     Module_radius+Module_gap ,  8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_21", 6.0 * cm / 2, 6.0 * cm,  3*(Module_radius+Module_gap),  8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_22", 6.0 * cm / 2, 6.0 * cm,  5*(Module_radius+Module_gap),  8*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_23", 6.0 * cm / 2, 6.0 * cm,  6*(Module_radius+Module_gap),  6*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_24", 6.0 * cm / 2, 6.0 * cm,  7*(Module_radius+Module_gap),  4*(Module_radius+Module_gap)*sin(60*deg)},
		   {"LMO6_4_25", 6.0 * cm / 2, 6.0 * cm,  8*(Module_radius+Module_gap),  2*(Module_radius+Module_gap)*sin(60*deg)}
		   */
	};
