#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "AmoreSim/AmoreDetectorStaticInfo.hh"

#include "CupSim/CupPMTSD.hh" // for "sensitive detector"
#include "CupSim/CupParam.hh"
#include "CupSim/Cup_PMT_LogicalVolume.hh"     // for making PMT assemblies
#include "CupSim/CupTorusStack.hh"
#include "CupSim/CupPMTOpticalModel.hh"
#include "CupSim/CupDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Region.hh"
#include "G4OpticalSurface.hh"
#include "G4ThreeVector.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"

#include <fstream>
#include <sstream>

////////////////////////////////////////////////////////////////
// declaration of "private" static utility functions that we
// don't need in class definition
/*
static void MakeID_PMT_Support(Cup_PMT_LogicalVolume *that,
			       G4Material *SupportMat,
			       G4Material *ExteriorMat
			       );
                   */

void AmoreDetectorConstruction::ConstructAMoRE200_WCMD(){
    // -- database --------------------------------------
	using namespace CLHEP;
	using namespace std;
    using namespace AmoreDetectorStaticInfo::AMoRE_200;
    CupParam &db ( CupParam::GetDB() );

	/////////////////////////////////////////////
	// --- Getting the physical volume of the Water cerenkov muon detector
	/////////////////////////////////////////////
    G4VPhysicalVolume* WCTankAir_PV = GetPhysicalVolumeByName("HatWaterTankAir_PV");
    G4VPhysicalVolume* WCTank_PV = GetPhysicalVolumeByName("HatWaterTank_PV");

    G4Material* wc_air_material = WCTankAir_PV->GetLogicalVolume()->GetMaterial();
    G4Material* wc_material = WCTank_PV->GetLogicalVolume()->GetMaterial();

    G4double Air_hh = (HatInnerZ + watertank_top_thickness + PMTroom_thickness)/2.;
    G4double WT_hh = (HatInnerZ + watertank_top_thickness)/2.;

    
	/////////////////////////////////////////////
	// --- Read the coordinates of the PMTs from the file
	/////////////////////////////////////////////
	ifstream wherePMT;
	const char *basic_fn = "pmtpos_amore200.dat";
	if (getenv("AmoreDATA") != NULL){
		wherePMT.open((G4String(getenv("AmoreDATA")) + "/" + G4String(basic_fn)).c_str());
	}
	// print error message on failure of file open
	if (wherePMT.fail()) {
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
		G4Exception(" ", " ", JustWarning, "Error, pmt coordinates file could not be opened.\n");
	}
	// read max number of pmt
	int maxPMTNo;
	wherePMT >> maxPMTNo;

	/////////////////////////////////////////////
	// --- PMT sensitive detector
	/////////////////////////////////////////////
	G4SDManager *fSDman = G4SDManager::GetSDMpointer();
    CupPMTSD * pmtSD = new CupPMTSD("/cupdet/pmt/inner", maxPMTNo, 0, 10);
    fSDman->AddNewDetector(pmtSD);

    /* // test for photon simulation using Cup_PMT_LogicalVolume
    Cup_PMT_LogicalVolume* _logiInnerPMT10
        = new Cup_10inch_LogicalVolume
        ( "InnerPMT",
        //_mineralOil,
        _air,
        _glass,
        Photocathode_opsurf,
        PMT_Vac,
        _stainless,  // dynode material
      //( db["omit_pmt_masks"] != 0.0 ?
      //NULL :       // no mask
      //_blackAcryl // physical mask on tubes to block non-sensitive areas
      //),
        NULL,
        pmtSD,  // sensitive detector hook
        whichPmtStyle);
  MakeID_PMT_Support(_logiInnerPMT10,
		     _air,     // support material
		     _air);
    */



    /////////////////////////////////////////////
    // --- Make PMT assembly
    /////////////////////////////////////////////
    // R7081 size used 
    const G4int n_edge = 6;
    const G4String vname = "InnerPMT";
    G4double outer_z_edge[n_edge + 1] = {96.5, 37.0, 0, -37.0, -68.5, -92.5, -203.5};
    G4double outer_rho_edge[n_edge + 1] = {0, 108.5, 126.5, 108.5, 58.7, 42.25, 42.25};
    G4double outer_z_o[n_edge] = {-40.7, 0.0, 0.0, 83.3, -92.5, -203.5};
    G4double r_dynode = 27.5;
    G4double z_dynode = -30.0;
    G4double d_wall = 3.0;
    G4double pmt_r = 130.0;
    G4double pmt_h = outer_z_edge[0]-outer_z_edge[n_edge]+d_wall*2;
    G4double pmt_cathode_h = outer_z_edge[0]+d_wall;
    G4double pmt_body_h = pmt_h - pmt_cathode_h; 

    // --- Make PMT Solid --------------------------------
    G4Tubs* body_envelope_solid = new G4Tubs("body_envelope_solid", 0, pmt_r, pmt_body_h/2., 0, 2*M_PI);
    G4Tubs* cathode_envelope_solid = new G4Tubs("cathode_envelope_solid", 0, pmt_r, pmt_cathode_h/2., 0, 2*M_PI);

    CupTorusStack *body_solid = new CupTorusStack(vname + "_body_solid");
    body_solid->SetAllParameters(n_edge, outer_z_edge, outer_rho_edge, outer_z_o);

    CupTorusStack *body1_solid = new CupTorusStack(vname + "_body1_solid");
    CupTorusStack *body2_solid = new CupTorusStack(vname + "_body2_solid");
    CupTorusStack *inner1_solid = new CupTorusStack(vname + "_inner1_solid");
    CupTorusStack *inner2_solid = new CupTorusStack(vname + "_inner2_solid");

    G4double z_lowest_dynode = z_dynode;
    {
        G4double *dscratch       = new G4double[4 * (n_edge + 1)];
        G4double *inner_z_edge   = dscratch;
        G4double *inner_rho_edge = dscratch + n_edge + 1;
        G4ThreeVector norm;
        G4int iedge_equator = -1;
        // calculate inner surface edges, check dynode position, and find equator
        inner_z_edge[0]   = outer_z_edge[0] - d_wall;
        inner_rho_edge[0] = 0.0;
        for (int i = 1; i < n_edge; i++) {
            norm =
                body_solid->SurfaceNormal(G4ThreeVector(0.0, outer_rho_edge[i], outer_z_edge[i]));
            inner_z_edge[i]   = outer_z_edge[i] - d_wall * norm.z();
            inner_rho_edge[i] = outer_rho_edge[i] - d_wall * norm.y();
            if (inner_rho_edge[i] > r_dynode && inner_z_edge[i] < z_lowest_dynode)
                z_lowest_dynode = inner_z_edge[i];
            if (outer_z_edge[i] == 0.0 || inner_z_edge[i] == 0.0) iedge_equator = i;
        }
        inner_z_edge[n_edge]   = outer_z_edge[n_edge] + d_wall;
        inner_rho_edge[n_edge] = outer_rho_edge[n_edge] - d_wall;
        // one final check on dynode allowed position
        if (inner_rho_edge[n_edge] > r_dynode && inner_z_edge[n_edge] < z_lowest_dynode)
            z_lowest_dynode = inner_z_edge[n_edge];
        // sanity check equator index
        if (iedge_equator < 0) {
            iedge_equator = (1 + n_edge) / 2;
            G4cerr << "CupSim/Cup_PMT_LogicalVolume::ConstructPMT_UsingTorusStack: "
                "Warning, pathological PMT shape, equator edge not found!"
            << G4endl;
        }
        // sanity check on dynode height
        if (z_dynode > inner_z_edge[iedge_equator]) {
            z_dynode = inner_z_edge[iedge_equator];
            G4cerr << "CupSim/Cup_PMT_LogicalVolume::ConstructPMT_UsingTorusStack: "
                "Warning, dynode higher than equator, dynode truncated!"
            << G4endl;
        }
        // set body surfaces
        body1_solid->SetAllParameters(iedge_equator, outer_z_edge, outer_rho_edge, outer_z_o);
        body2_solid->SetAllParameters(n_edge - iedge_equator, outer_z_edge + iedge_equator,
            outer_rho_edge + iedge_equator, outer_z_o + iedge_equator);
        // set inner surfaces
        inner1_solid->SetAllParameters(iedge_equator, inner_z_edge, inner_rho_edge, outer_z_o);
        inner2_solid->SetAllParameters(n_edge - iedge_equator, inner_z_edge + iedge_equator,
            inner_rho_edge + iedge_equator, outer_z_o + iedge_equator);
        // CupTorusStack keeps its own copy of edges, so we can delete our workspace
        delete[] dscratch;
    }

    // --- Dynode volume
    G4double hh_dynode = (z_dynode - z_lowest_dynode) / 2.0;
    G4Tubs *dynode_solid = new G4Tubs(vname + "_dynode_solid", 0.0, 
        r_dynode, hh_dynode, 0., 2. * M_PI); 


    /////////////////////////////////////////////
    // --- Make Logical Volume
    /////////////////////////////////////////////
    G4LogicalVolume *body_envelope_LV = 
        new G4LogicalVolume(body_envelope_solid, wc_air_material, vname + "_body_LV");
    body_envelope_LV->SetVisAttributes(G4VisAttributes::Invisible);
    G4LogicalVolume *cathode_envelope_LV = 
        new G4LogicalVolume(cathode_envelope_solid, wc_material, vname + "_cathode_LV");
    cathode_envelope_LV->SetVisAttributes(G4VisAttributes::Invisible);

    G4LogicalVolume *body1_log = 
        new G4LogicalVolume(body1_solid, _glass, vname + "_body1_log");
    body1_log->SetSensitiveDetector(pmtSD);
    body1_log->SetVisAttributes(G4VisAttributes::Invisible);
    G4LogicalVolume *body2_log = 
        new G4LogicalVolume(body2_solid, _glass, vname + "_body2_log");
    body2_log->SetVisAttributes(G4VisAttributes::Invisible);

    G4LogicalVolume *inner1_log =
        new G4LogicalVolume(inner1_solid, PMT_Vac, vname + "_inner1_log");
    inner1_log->SetSensitiveDetector(pmtSD);
    inner1_log->SetVisAttributes(G4Color(0.7, 0.5, 0.3, 0.27)); // orange

    G4LogicalVolume *inner2_log =
        new G4LogicalVolume(inner2_solid, PMT_Vac, vname + "_inner2_log");
    inner2_log->SetVisAttributes(G4Color(0.6, 0.7, 0.8, 0.67)); // light blue

    G4LogicalVolume *dynode_log =
        new G4LogicalVolume(dynode_solid, _stainless, vname + "_dynode_log");
    dynode_log->SetVisAttributes(G4Color(0.5, 0.5, 0.5, 1.0)); // medium gray


    ////////////////////////////////////////////////////////////////
    // Attach optical surfaces to borders
    ////
    G4OpticalSurface *our_Mirror_opsurf = new G4OpticalSurface("Mirror_opsurf");
    our_Mirror_opsurf->SetFinish(polishedfrontpainted); // needed for mirror
    our_Mirror_opsurf->SetModel(glisur);
    our_Mirror_opsurf->SetType(dielectric_metal);
    our_Mirror_opsurf->SetPolish(0.999); // a guess -- FIXME
    G4MaterialPropertiesTable *propMirror = NULL;
    G4Material *matMirror                 = G4Material::GetMaterial("PMT_Mirror");
    if (matMirror) propMirror = matMirror->GetMaterialPropertiesTable();
    if (propMirror == NULL) {
        G4cerr << "Warning: setting PMT mirror reflectivity to 0.9999 because no PMT_Mirror "
            "material properties defined"
        << G4endl;
        propMirror = new G4MaterialPropertiesTable();
        propMirror->AddProperty("REFLECTIVITY", new G4MaterialPropertyVector());
        propMirror->AddEntry("REFLECTIVITY", twopi * hbarc / (800.0e-9 * m), 0.9999);
        propMirror->AddEntry("REFLECTIVITY", twopi * hbarc / (200.0e-9 * m), 0.9999);
    }
    our_Mirror_opsurf->SetMaterialPropertiesTable(propMirror);


    ////////////////////////////////////////////
	// --- position the PMTs
	////////////////////////////////////////////
	G4int region, pmt_id;
	G4double cood_x, cood_y, cood_z, pmt_size;
    G4RotationMatrix* _pmtRotMtx = new G4RotationMatrix();
    _pmtRotMtx->rotateX(M_PI);

    G4PVPlacement *body1_phys = new G4PVPlacement(
        _pmtRotMtx, G4ThreeVector(0,0,pmt_cathode_h/2.), body1_log,  
        vname + "_body1_phys", cathode_envelope_LV , false, 0, OverlapCheck);

    G4PVPlacement *body2_phys = new G4PVPlacement(
        _pmtRotMtx, G4ThreeVector(0,0, -pmt_body_h/2.), body2_log, 
        vname + "_body2_phys", body_envelope_LV, false, 0, OverlapCheck);

    G4PVPlacement *inner1_phys = new G4PVPlacement( 0, G4ThreeVector(0,0,0), 
        vname + "_inner1_phys", inner1_log, body1_phys, false, 0, OverlapCheck);

    G4PVPlacement *inner2_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0),
        vname + "_inner2_phys", inner2_log, body2_phys, false, 0, OverlapCheck);

    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, z_dynode - hh_dynode), 
        vname + "_dynode_phys", dynode_log, inner2_phys, false, 0, OverlapCheck); 

    new G4LogicalBorderSurface(vname + "_photocathode_logsurf1", inner1_phys, body1_phys,
        Photocathode_opsurf);
    new G4LogicalBorderSurface(vname + "_photocathode_logsurf2", body1_phys, inner1_phys,
        Photocathode_opsurf);
    new G4LogicalBorderSurface(vname + "_mirror_logsurf1", inner2_phys, body2_phys,
        our_Mirror_opsurf);
    new G4LogicalBorderSurface(vname + "_mirror_logsurf2", body2_phys, inner2_phys,
        our_Mirror_opsurf);


    /////////////////////////////////////////////
    // --- FastSimulationModel
    /////////////////////////////////////////////
    G4Region *PmtRegion = new G4Region(vname);
    // PmtRegion->AddRootLogicalVolume(body_envelope_LV);
    // PmtRegion->AddRootLogicalVolume(cathode_envelope_LV);
    PmtRegion->AddRootLogicalVolume(body1_log);
    new CupPMTOpticalModel(vname + "_optical_model1", body1_phys);
    PmtRegion->AddRootLogicalVolume(body2_log);
    new CupPMTOpticalModel(vname + "_optical_model2", body2_phys);



    while (wherePMT.good()) {
        // get a line from the file
        char linebuffer[128];
        wherePMT.getline( linebuffer, sizeof(linebuffer)-1);
        if ( wherePMT.fail()) break;

        // skip blank lines and lines beginning with '#'
        if (linebuffer[0] == '#' || linebuffer[0] == '\0') continue;

        // put the line in an istrstream for convenient parsing
        istringstream lineStream(linebuffer);

        // parse out region, coordinates,
        region = pmt_id = -1;
        cood_x = cood_y = cood_z = pmt_size = -1;
        lineStream >> region >> pmt_id >> cood_x >> cood_y >> cood_z >> pmt_size;

        // check for bad data
        if (lineStream.fail() || region < 0 || pmt_id < 0 || (cood_x == -1 && cood_y == -1 && cood_z == -1) || pmt_size == -1) {
            G4cerr << "BAD DATA in PMT file:   line=\"" << linebuffer << "\"\n";
            G4cerr.flush();
            continue;
        }

        // DAQ numbering 
        int ith = pmt_id;

/* // test for Photon simulation using Cup_PMT_LogicalVolume
        char PMTname[64];
        sprintf(PMTname,"physInnerPMT%d",ith);
        G4ThreeVector pmtpos( cood_x, cood_y, cood_z+pmt_body_h/2.+Air_hh-PMTroom_thickness);
        if ( db["omit_id_pmts"] == 0.0 ) {
            new G4PVPlacement(_pmtRotMtx,
			    pmtpos,
			    PMTname,
			    _logiInnerPMT10,
			    WCTankAir_PV,    // physical parent
			    false,
			    ith);
        }
        */

        // PMT name
        char PMTname1[64];
        char PMTname2[64];

        sprintf(PMTname1, "physInnerPMT%d_cathode", ith);
        sprintf(PMTname2, "physInnerPMT%d_body", ith);
        // sprintf(PMTname1, "physWCVetoPMT%d_cathode", ith);
        // sprintf(PMTname2, "physWCVetoPMT%d_body", ith);

        // position the PMTs
        G4ThreeVector pmtpos1(cood_x, cood_y, cood_z-pmt_cathode_h/2.+WT_hh);
        G4ThreeVector pmtpos2(cood_x, cood_y, cood_z+pmt_body_h/2.+Air_hh-PMTroom_thickness);

        new G4PVPlacement(0, pmtpos1, PMTname1, 
            cathode_envelope_LV, WCTank_PV, false, ith, OverlapCheck); 
        new G4PVPlacement(0, pmtpos2, PMTname2, 
            body_envelope_LV, WCTankAir_PV, false, ith, OverlapCheck);
    }
    wherePMT.close();


    f200_logiPMTbody = body1_log; 
    f200_logiPMTinner = inner1_log;

    G4cout << " The total number of PMT for Water cerenkov detector = " << maxPMTNo << G4endl; 

}

/* // test for photon simulation using Cup_PMT_LogicalVolume
////////////////////////////////////////////////////////////////
// definition of "private" static utility functions that we
// don't need in class definition
static void MakeID_PMT_Support(Cup_PMT_LogicalVolume *that,
			       G4Material *SupportMat,
			       G4Material *ExteriorMat)
{
  // ... this should build the PMT support geometry
  // ... empty function for now
}
*/