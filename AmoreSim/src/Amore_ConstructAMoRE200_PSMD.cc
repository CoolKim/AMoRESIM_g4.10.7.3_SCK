// Written by Jeewon Seo (2024.02.13.)

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreDetectorStaticInfo.hh"

#include "CupSim/CupParam.hh"
#include "CupSim/CupVetoSD.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4VisAttributes.hh"

// #include "TString.h"

#include <fstream>
#include <sstream>

void AmoreDetectorConstruction::ConstructAMoRE200_PSMD()
{
    CupParam &db ( CupParam::GetDB() );
	using namespace CLHEP;
	using namespace std;
	using namespace AmoreDetectorStaticInfo;
	using namespace AmoreDetectorStaticInfo::AMoRE_200;


	/////////////////////////////////////////////
	// --- Getting the physical volume of the veto housing
	/////////////////////////////////////////////
    G4VPhysicalVolume* VetoHousingSide_PV = GetPhysicalVolumeByName("PlasticVetoHousing_PV"); 
    G4VPhysicalVolume* VetoHousingBottom_PV = GetPhysicalVolumeByName("PlasticVetoHousing2_PV");


	/////////////////////////////////////////////
    // --- Retrieve various geometry parameters
	/////////////////////////////////////////////
    //G4double rockshell_thickness       = db["rockshell_thickness"];
	//G4double airbuffer_radius          = db["airbuffer_radius"];
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

	//G4double nShield_hatGapFromLead    = db["nShield_hatGapFromLead"];
	//G4double nShield_GapFromCeiling    = db["nShield_GapFromCeiling"];

	//G4double cavern_sphere_radius      = db["cavern_sphere_radius"];
	//G4double rock_shell_radius         = db["rock_shell_radius"];

	// G4double lead_shield_halfsize    =  // Lead shield  half x,y
		// airbuffer_radius + thin_lead_shield_thickness + lead_shield_thickness; 

	// G4double PS_housing_halfsize =  // Barrel size
		// longPSlength + profile_thickness*2.5 + plastic_veto_thickness;

	G4double PS_housing_height =  airbuffer_height + thin_lead_shield_thickness + lead_shield_thickness 
		+ lead_housing_thickness*2 + boricacid_thickness + PE_shield_thickness;


	/////////////////////////////////////////////
	// --- Visualization attributes
	/////////////////////////////////////////////
	G4VisAttributes *plasticScintVis = new G4VisAttributes(G4Colour(0.8, 1.0, 0.725, 0.4));
	plasticScintVis->SetForceSolid(true);
	G4VisAttributes *aluminiumVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.1));
	aluminiumVis->SetForceSolid(true);
	G4VisAttributes *stainlessVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.7, 0.9));
	stainlessVis->SetForceSolid(true);


	/////////////////////////////////////////////
	// --- Read the coordinates of the PSs from the file
	/////////////////////////////////////////////
	ifstream wherePS;
	const char *basic_fn = "pspos_amore200.dat";
	if (getenv("AmoreDATA") != NULL){
		wherePS.open((G4String(getenv("AmoreDATA")) + "/" + G4String(basic_fn)).c_str());
	}
	// print error message on failure of file open
	if (wherePS.fail()) {
		G4cerr << "Error, " << basic_fn << " could not be opened.\n";
		if (getenv("AmoreDATA") == NULL) {
			G4cerr << 
					"AmoreDATA environment variable is not set, so I was looking for"
					<< basic_fn << " in the current directory." 
			<< G4endl;
		} else {
			G4cerr << 
					"I was looking for it in the AmoreDATA directory, "
					<< getenv("AmoreDATA") 
			<< G4endl;
		}
		G4Exception(" ", " ", JustWarning, "Error, ps coordinates file could not be opened.\n");
	}
	// read max number of ps
	int maxPSNo;
	wherePS >> maxPSNo;

	f200_TotPSNum = maxPSNo;

	/////////////////////////////////////////////
	// --- PS sensitive detector
	/////////////////////////////////////////////
	/*
	G4SDManager *fSDman = G4SDManager::GetSDMpointer();
	CupVetoSD *PS_SD = new CupVetoSD("/cupdet/MuonVetoSD", maxPSNo);
	fSDman->AddNewDetector(PS_SD);
	*/
	// CupVetoSD *PSO_SD = new CupVetoSD("/cupdet/MuonVetoSD/PSO", maxPSNo);
	// // CupVetoSD *PSI_SD = new CupVetoSD("/cupdet/MuonVetoSD/PSI", maxPSNo);
	// fSDman->AddNewDetector(PSO_SD);
	// fSDman->AddNewDetector(PSI_SD);


    ////////////////////////////////////////////
	// --- Muon Veto (Plastic scintillator)
	////////////////////////////////////////////

    const int nPStype = 2;
	G4double PSlength[nPStype] = {longPSlength/2., shortPSlength/2.};

	G4Box *plasticVetoBox[nPStype];
	G4Box *plasticVetoCut[nPStype];
	G4Box *plasticVetoAirBox[nPStype];
	G4Box *plasticScintHolderBox[nPStype];
	G4Box *plasticScintBox[nPStype];

	G4VSolid *plasticVetoSolid1[nPStype];
	G4VSolid *plasticVetoSolid[nPStype];

	G4LogicalVolume *aluminiumHolderLV[nPStype];
	G4LogicalVolume *plasticScintOLV[nPStype];
	G4LogicalVolume *plasticScintILV[nPStype];
	G4LogicalVolume *plasticVetoAirLV[nPStype];
	G4LogicalVolume *plasticVetoLV[nPStype];

	for(G4int PStype = 0; PStype < nPStype; PStype++){
		// plastic scintillator
		// plasticScintBox[PStype] = new G4Box(Form("PlasticScint_Box%d",PStype), 
		plasticScintBox[PStype] = new G4Box((G4String("PlasticScint_Box") + std::to_string(PStype)), 
				plastic_scintillator_thickness/2.,
				PSlength[PStype] ,//- veto_frame_space,
				plastic_veto_width/2. - veto_frame_thickness);// - veto_frame_space);

		plasticScintOLV[PStype] = new G4LogicalVolume(plasticScintBox[PStype], _vinylt, "PlasticScintO_LV");
		// plasticScintOLV[PStype] -> SetSensitiveDetector(PS_SD);
		// plasticScintOLV[PStype] -> SetSensitiveDetector(PSO_SD);
		plasticScintOLV[PStype] -> SetVisAttributes(plasticScintVis);

		plasticScintILV[PStype] = new G4LogicalVolume(plasticScintBox[PStype], _vinylt, "PlasticScintI_LV");
		// plasticScintILV[PStype] -> SetSensitiveDetector(PSI_SD);
		// plasticScintILV[PStype] -> SetSensitiveDetector(PS_SD);
		plasticScintILV[PStype] -> SetVisAttributes(plasticScintVis);

		// aluminium plate
		//plasticScintHolderBox[PStype] = new G4Box(Form("AluminiumHolder_Box%d",PStype),
		plasticScintHolderBox[PStype] = new G4Box(G4String("AluminiumHolder_Box") + std::to_string(PStype),
				al_plate_thickness/2., 
				PSlength[PStype],
				plastic_veto_width/2. - veto_frame_thickness);

		aluminiumHolderLV[PStype] = new G4LogicalVolume(plasticScintHolderBox[PStype], _aluminium, "AluminiumHolder_LV");
		aluminiumHolderLV[PStype]->SetVisAttributes(aluminiumVis);

		// air volume
		//plasticVetoAirBox[PStype] = new G4Box(Form("PlasticVetoAir_Box%d",PStype),
		plasticVetoAirBox[PStype] = new G4Box(G4String("PlasticVetoAir_Box") + std::to_string(PStype),
				plastic_veto_thickness/2. - veto_frame_thickness,
				PSlength[PStype],
				plastic_veto_width/2. - veto_frame_thickness);
		plasticVetoAirLV[PStype] = new G4LogicalVolume(plasticVetoAirBox[PStype], _air, "PlasticVetoAir_LV");

		// muon veto box
		plasticVetoBox[PStype] = new G4Box(G4String("PlasticVeto_Box") + std::to_string(PStype),
				plastic_veto_thickness/2., 
				PSlength[PStype] + veto_frame_thickness,// + veto_frame_space,
				plastic_veto_width/2. );// + veto_frame_space);
		//plasticVetoCut[PStype] = new G4Box(Form("PlasticVetoCut%d",PStype),
		plasticVetoCut[PStype] = new G4Box(G4String("PlasticVetoCut") + std::to_string(PStype),
				plastic_veto_thickness/2.,
				PSlength[PStype] + veto_frame_thickness - veto_frame_width,
				plastic_veto_width/2. - veto_frame_width);

		// plasticVetoSolid1[PStype] = new G4SubtractionSolid(Form("PlasticVetoSolid1%d",PStype),
		plasticVetoSolid1[PStype] = new G4SubtractionSolid(G4String("PlasticVetoSolid1") + std::to_string(PStype),
				plasticVetoBox[PStype], plasticVetoCut[PStype],0,
				{plastic_veto_thickness - veto_frame_thickness/2. + solidBooleanTol,0,0});
		// plasticVetoSolid[PStype] = new G4SubtractionSolid(Form("PlasticVetoSolid%d",PStype),
		plasticVetoSolid[PStype] = new G4SubtractionSolid(G4String("PlasticVetoSolid") + std::to_string(PStype),
				plasticVetoSolid1[PStype], plasticVetoCut[PStype],0,
				{-plastic_veto_thickness + veto_frame_thickness/2. - solidBooleanTol,0,0});

		plasticVetoLV[PStype] = new G4LogicalVolume(plasticVetoSolid[PStype], _stainless, "PlasticVetoEnvelope_LV");
		plasticVetoLV[PStype]->SetVisAttributes(stainlessVis);

		// place the inner components of the Veto Box

		new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,0.), plasticVetoAirLV[PStype], 
			"PlasticVetoAir_PV", plasticVetoLV[PStype], false, 0, OverlapCheck);

		new G4PVPlacement(nullptr, G4ThreeVector(
			plastic_veto_thickness/2. - veto_frame_thickness - al_plate_thickness/2.,0,0), 
			aluminiumHolderLV[PStype], "AluminiumHolder_PV", plasticVetoAirLV[PStype], false, 0, OverlapCheck);

		new G4PVPlacement(nullptr, G4ThreeVector(
			-plastic_veto_thickness/2. + veto_frame_thickness + al_plate_thickness/2.,0,0), 
			aluminiumHolderLV[PStype], "AluminiumHolder_PV", plasticVetoAirLV[PStype], false, 0, OverlapCheck);

		new G4PVPlacement(nullptr, G4ThreeVector(
			plastic_veto_thickness/2. - veto_frame_thickness - al_plate_thickness - plastic_scintillator_thickness/2.,0, 0), 
			plasticScintOLV[PStype], "PlasticScintO_PV", plasticVetoAirLV[PStype], false, 0, OverlapCheck);

		new G4PVPlacement(nullptr, G4ThreeVector(
			-plastic_veto_thickness/2. + veto_frame_thickness + al_plate_thickness + plastic_scintillator_thickness/2.,0, 0),
			plasticScintILV[PStype], "PlasticScintI_PV", plasticVetoAirLV[PStype], false, 0, OverlapCheck);
	}

	G4Region *PSRegion = new G4Region("PS_muonveto");
	for (int i = 0; i < nPStype; i++){
		PSRegion->AddRootLogicalVolume(plasticScintOLV[i]);
		PSRegion->AddRootLogicalVolume(plasticScintILV[i]);
	}


    ////////////////////////////////////////////
	// --- position the plastic scintillators
	////////////////////////////////////////////
	G4int region, ps_id;
	G4double cood_x, cood_y, cood_z, ps_size;

	G4RotationMatrix *psRotMtx = new G4RotationMatrix();
	psRotMtx->rotateZ(-90*deg);

	G4double grad = 3;
	G4RotationMatrix *vpsRotMtx = new G4RotationMatrix();
	vpsRotMtx->rotateY(90*deg);
	vpsRotMtx->rotateX((90+grad)*deg);

	G4RotationMatrix *bpsRotMtx = new G4RotationMatrix();
	bpsRotMtx->rotateY(90*deg);

	while (wherePS.good()) {
		// get a line from the file
		char linebuffer[128];
		wherePS.getline(linebuffer, sizeof(linebuffer) - 1);
		if (wherePS.fail()) break;

		// skip blank lines and lines beginning with '#'
		if (linebuffer[0] == '#' || linebuffer[0] == '\0') continue;

		// put the line in an istrstream for convenient parsing
		istringstream lineStream(linebuffer);

		// parse out region, coordinates,
		region = ps_id = -1;
		cood_x = cood_y = cood_z = ps_size = 0.0;
		lineStream >> region >> ps_id >> cood_x >> cood_y >> cood_z >> ps_size;

		// check for bad data
		if (lineStream.fail() || region < 0 || (cood_x == 0. && cood_y == 0. && cood_z == 0.) || ps_size == 0.) {
			G4cerr << "BAD DATA in PS position file:  line=\"" << linebuffer << "\"\n";
			G4cerr.flush();
			continue;
		}

		// a possible DAQ numbering
		int ith = ps_id;

		// name this PS
		char PSname[64];
		sprintf(PSname, "physMuonVetoPS_Envelope%d", ith);

		// position the PS
		G4ThreeVector pos(cood_x, cood_y, cood_z);

		// place the PS
		if (region < 9){ // side PS
			new G4PVPlacement(G4Transform3D(*psRotMtx, pos), PSname, plasticVetoLV[0], VetoHousingSide_PV, false, ith, OverlapCheck);
		} else if (region < 11){ // bottom PS
			new G4PVPlacement(G4Transform3D(*bpsRotMtx, pos), PSname, plasticVetoLV[0], VetoHousingBottom_PV, false, ith, OverlapCheck);
		} else if (ps_size == longPSlength) { // vertical side PS
			new G4PVPlacement(G4Transform3D(*vpsRotMtx, pos), PSname, plasticVetoLV[0], VetoHousingSide_PV, false, ith, OverlapCheck);
		} else{ // vertical side PS
			new G4PVPlacement(G4Transform3D(*vpsRotMtx, pos), PSname, plasticVetoLV[1], VetoHousingSide_PV, false, ith, OverlapCheck);
		}

		if (ith == 12 || ith == 3*12 || ith == 5*12+17 || ith == 7*12+17) {
			psRotMtx->rotateZ(-90*deg);
		}
		if (ith == 132 || ith == 135 || ith == 138 || ith == 141){
			vpsRotMtx->rotateZ(90*deg);
		}
	}
	wherePS.close();

	f200_logiVetoPSO = *plasticScintOLV;
	f200_logiVetoPSI = *plasticScintILV;


	/////////////////////////////////////////////
	// --- PS's supporter (aluminium profile)
	/////////////////////////////////////////////
	// veto supporter (aluminium profile)
	G4Box *plasticVetoSupporterV1Box = new G4Box("PlasticVetoSupporterV1_Box",
			profile_thickness/2., profile_thickness/2., PS_housing_height/2.);
	G4Box *plasticVetoSupporterV2Box = new G4Box("PlasticVetoSupporterV2_Box",
			(profile_thickness - plastic_veto_thickness)/2., 
			(profile_thickness - plastic_veto_thickness)/2., 
			(PS_housing_height-(profile_thickness - plastic_veto_thickness))/2.);
	G4Box *plasticVetoSupporterH1Box = new G4Box("PlasticVetoSupporterH1_Box",
			profile_thickness/2., shortPSlength/2., profile_thickness/2.);
	G4Box *plasticVetoSupporterH2Box = new G4Box("PlasticVetoSupporterH2_Box",
			(profile_thickness - plastic_veto_thickness)/2.,
			shortPSlength/4.,
			(profile_thickness - plastic_veto_thickness)/2.);

	G4LogicalVolume *plasticVetoSupporterV1LV = new G4LogicalVolume(plasticVetoSupporterV1Box, 
			_alprofile, "PlasticVetoSupporterV1_LV");
	G4LogicalVolume *plasticVetoSupporterV2LV = new G4LogicalVolume(plasticVetoSupporterV2Box, 
			_alprofile, "PlasticVetoSupporterV2_LV");
	G4LogicalVolume *plasticVetoSupporterH1LV = new G4LogicalVolume(plasticVetoSupporterH1Box, 
			_alprofile, "PlasticVetoSupporterH1_LV");
	G4LogicalVolume *plasticVetoSupporterH2LV = new G4LogicalVolume(plasticVetoSupporterH2Box, 
			_alprofile, "PlasticVetoSupporterH2_LV");
	plasticVetoSupporterV1LV->SetVisAttributes(aluminiumVis);
	plasticVetoSupporterV2LV->SetVisAttributes(aluminiumVis);
	plasticVetoSupporterH1LV->SetVisAttributes(aluminiumVis);
	plasticVetoSupporterH2LV->SetVisAttributes(aluminiumVis);

	/////////////////////////////////////////////
	// --- Read the coordinates of the PS's supporter
	/////////////////////////////////////////////
	ifstream whereProf;
	const char *prof_fn = "profpos_amore200.dat";
	if (getenv("AmoreDATA") != NULL){
		whereProf.open((G4String(getenv("AmoreDATA")) + "/" + G4String(prof_fn)).c_str());
	}
	// print error message on failure of file open
	if (whereProf.fail()) {
		G4cerr << "Error, " << prof_fn << " could not be opened.\n";
		if (getenv("AmoreDATA") == NULL) {
			G4cerr << 
					"AmoreDATA environment variable is not set, so I was looking for"
					<< prof_fn << " in the current directory." 
			<< G4endl;
		} else {
			G4cerr << 
					"I was looking for it in the AmoreDATA directory, "
					<< getenv("AmoreDATA") 
			<< G4endl;
		}
		G4Exception(" ", " ", JustWarning, "Error, profile coordinates file could not be opened.\n");
	}

	G4int prof_id;
	G4double prof_size;
	G4RotationMatrix *profRotMtx = new G4RotationMatrix();
	profRotMtx->rotateZ(-90*deg);

	while (whereProf.good()) {
		// get a line from the file
		char linebuffer[128];
		whereProf.getline(linebuffer, sizeof(linebuffer) - 1);
		if (whereProf.fail()) break;

		// skip blank lines and lines beginning with '#'
		if (linebuffer[0] == '#' || linebuffer[0] == '\0') continue;

		// put the line in an istrstream for convenient parsing
		istringstream lineStream(linebuffer);

		// parse out region, coordinates,
		region = prof_id = -1;
		cood_x = cood_y = cood_z = prof_size = 0.0;
		lineStream >> region >> prof_id >> cood_x >> cood_y >> cood_z >> prof_size;

		// check for bad data
		if (lineStream.fail() || region < 0 || (cood_x == 0. && cood_y == 0. && cood_z == 0.) || prof_size == 0.) {
			G4cerr << "BAD DATA in profile position file:  line=\"" << linebuffer << "\"\n";
			G4cerr.flush();
			continue;
		}

		// a possible DAQ numbering
		int ith = prof_id;

		// name this PS
		char profname[64];
		sprintf(profname, "physPSSupporter%d", ith);

		// position the PS
		G4ThreeVector pos(cood_x, cood_y, cood_z);

		if (prof_size==1){
			new G4PVPlacement(nullptr, pos, profname, plasticVetoSupporterV1LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		} else if (prof_size==2){
			new G4PVPlacement(nullptr, pos, profname, plasticVetoSupporterV2LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		} else if (prof_size==21){
			new G4PVPlacement(G4Transform3D(*profRotMtx, pos), profname, plasticVetoSupporterH1LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		} else if (prof_size==22){
			new G4PVPlacement(G4Transform3D(*profRotMtx, pos), profname, plasticVetoSupporterH2LV, VetoHousingSide_PV, false, 0, OverlapCheck);
		}

		if (ith == 2 || ith == 10 || ith == 16 || ith == 24) {
			profRotMtx->rotateZ(-90*deg);
		}
	}
	whereProf.close();

	/////////////////////////////////////////////
	// --- Print out the geometry information
	/////////////////////////////////////////////
/*
	G4cout << " ======      Muon Veto      ======" << G4endl;
	G4cout << " PS supporter (Aluminium profile)" << G4endl;
	G4cout << "     mass: " << plasticVetoSupporterV1LV->GetMass(true,false)/kg << "(" 
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
		<< plasticVetoSupporterH2Box->GetZHalfLength()*2. << ")" << G4endl;

	G4cout << " PS veto stainless flame" << G4endl;
	G4cout << "     mass       : " << plasticVetoLV[0]->GetMass(true,false)/kg << "("
		<< plasticVetoBox[0]->GetXHalfLength()*2 << " x "
		<< plasticVetoBox[0]->GetYHalfLength()*2 << " x "
		<< plasticVetoBox[0]->GetZHalfLength()*2 << ") \n"
		<< "     coordinate : " <<  VetoHousingSide_PV->GetTranslation() 
			+ GetPhysicalVolumeByName("physMuonVetoPS0")->GetTranslation() << G4endl;

	G4cout << " PS veto Aluminium plate " << G4endl;
	G4cout << "     mass       : " << aluminiumHolderLV[0]->GetMass(true,false)/kg << "("
		<< plasticScintHolderBox[0]->GetXHalfLength()*2 << "x"
		<< plasticScintHolderBox[0]->GetYHalfLength()*2 << "x"
		<< plasticScintHolderBox[0]->GetZHalfLength()*2 << ")" << G4endl;
	
	G4cout << " Plastic scintillator " << G4endl;
	G4cout << "     mass       : " << plasticScintOLV[0]->GetMass(true,false)/kg << "("
		<< plasticScintBox[0]->GetXHalfLength()*2 << "x"
		<< plasticScintBox[0]->GetYHalfLength()*2 << "x"
		<< plasticScintBox[0]->GetZHalfLength()*2 << ")" << G4endl;
		*/
	
	G4cout << " The total number of PS veto detector = " << maxPSNo << G4endl; 

}