/*
   AmoreDetectorConstruction
   */

#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "AmoreSim/AmoreDetectorStaticInfo.hh"
#include "CupSim/CupDetectorMessenger.hh"
#include "CupSim/CupInputDataReader.hh"

#include "AmoreSim/AmoreDetectorMessenger.hh"
#include "CupSim/CupParam.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Version.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include "G4ios.hh"

#include <map>

using namespace std;
//using namespace CLHEP;

AmoreDetectorConstruction::eDetGeometry AmoreDetectorConstruction::whichDetGeometry =
    kNumDetGeometries;
AmoreDetectorConstruction::eVetoGeometry AmoreDetectorConstruction::whichVetoGeometry =
    kNumVetoGeometries;
AmoreDetectorConstruction::ePhaseAMoRE200 AmoreDetectorConstruction::whichAMoRE200Phase =
		kNumAMoRE200Phase;
AmoreDetectorConstruction::eNShieldConf AmoreDetectorConstruction::whichNShieldingConf =
    kNumNSConfs;

AmoreDetectorConstruction::eCavernType AmoreDetectorConstruction::whichCavernType = kNumCavernTypes;

AmoreDetectorConstruction::AmoreDetectorConstruction() : CupDetectorConstruction() {
    whichDetector  = kDetector_AmoreDetector;
    AmoreMessenger = new AmoreDetectorMessenger(this);

    // Default values for geometry settings

    fEnable_OriginalGeom  = true;
    fEnable_Scintillator  = true;
    fEnable_Gantry        = true;
    fEnable_InnerDetector = true;
    fEnable_Innermost     = true;
    fEnable_NeutronShield = true;

    fDbgMsgOn      = true;
		OverlapCheck  = true;
    fNeutronMode   = false;
		fRockgammaMode = false;

		// For AMoreII
		fAdditionalPE      = true;
    fNeutShieldingConf = -1;

    // AMoRE-Pilot LV and PV for SD information & Geometry settings
    fPilot_Enable_TargetRoom = true;
    fPilot_Enable_Mumetal    = false;

    fPilot_TopScint_BoxLogical    = nullptr;
    fPilot_SideFBScint_BoxLogical = nullptr;
    fPilot_SideLRScint_BoxLogical = nullptr;

    fPilot_TopScint_FlatTrapLogical    = nullptr;
    fPilot_SideFBScint_FlatTrapLogical = nullptr;
    fPilot_SideLRScint_FlatTrapLogical = nullptr;

    fPilot_SideFBScint_PMTTrapLogical = nullptr;
    fPilot_TopScint_PMTTrapLogical    = nullptr;
    fPilot_SideLRScint_PMTTrapLogical = nullptr;

    fPilot_logiCMOCell  = nullptr;
    fPilot_logicPMT     = nullptr;
    fPilot_logicPMTVacu = nullptr;
    fPilot_logiGOLDa    = nullptr;
    fPilot_logiGOLDb    = nullptr;

    // AMoRE-I LV and PV for SD information & Geometry settings
    fI_Enable_SuperConductingShield = true;
    fI_Enable_CrystalArray          = true;

    fI_TopScint_BoxLogical       = nullptr;
    fI_SideFBScint_BoxLogical    = nullptr;
    fI_MufflerFBScint_BoxLogical = nullptr;
    fI_SideLRScint_BoxLogical    = nullptr;

    fI_TopScint_FlatTrapLogical             = nullptr;
    fI_SideFBScint_FlatTrapLogical          = nullptr;
    fI_MufflerFBScint_FlatTrapLogical       = nullptr;
    fI_SideLRScint_FlatTrapLogical          = nullptr;

    fI_TopScint_PMTTrapLogical       = nullptr;
    fI_SideFBScint_PMTTrapLogical    = nullptr;
    fI_MufflerFBScint_PMTTrapLogical = nullptr;
    fI_SideLRScint_PMTTrapLogical    = nullptr;

    fI_DetectorModuleRegion = nullptr;
		fI_crystalsRegion = nullptr;

    fI_logicPMT     = nullptr;
    fI_logicPMTVacu = nullptr;

    fI_CopperDefault_VisAttr = nullptr;
    fI_LeadDefault_VisAttr   = nullptr;
    fI_GoldDefault_VisAttr   = nullptr;
    fI_CopperFrame_VisAttr   = nullptr;
    fI_OpticalFrame_VisAttr  = nullptr;
    fI_GeWafer_VisAttr       = nullptr;
    fI_LeadShield_VisAttr    = nullptr;
    fI_Reflector_VisAttr     = nullptr;
    fI_Crystal_VisAttr       = nullptr;
		fI_BrassBolt_VisAttr     = nullptr;

    // AMoRE-II LV and PV for SD information
		f200_TotTowerNum          = 0;
		f200_TotCrystalNum        = 0;
		f200_VetoTotCNum          = 0;
		f200_HatVetoTotCNum       = 0;
		f200_TotalVetoCNum        = 0;
		//f200_logiCrystalCell      = 0;
		for(int i = 0; i < 1000; i++){
			f200_logiCrystalCell[i] = nullptr;
		}
		f200_logiVetoPSO          = nullptr;
		f200_logiVetoPSI          = nullptr;
		f200_logiHatPSO           = nullptr;
		f200_logiHatPSI           = nullptr;
    f200_physGeWafer          = nullptr;
    f200_physVacDisk          = nullptr;
		f200_HatVetoMaterialPV    = nullptr;
		f200_PSO_PV               = nullptr;
		f200_PSI_PV               = nullptr;
		f200_HatPSO_PV            = nullptr;
		f200_HatPSI_PV            = nullptr;
    f200_VetoActiveMaterialPV = nullptr;
    f200_FloorPEPhysical      = nullptr;
    f200_physGeWafer          = nullptr;
    f200_physVacDisk          = nullptr;
		f200_HatVetoMaterialPV    = nullptr;
    f200_FloorPEPhysical      = nullptr;
    f200_CeilingPEPhysical    = nullptr;
    f200_RealPEPhysical       = nullptr;
    f200_VetoMaterialPhysical = nullptr;
    f200_AirBufferPhysical    = nullptr;
}

AmoreDetectorConstruction::~AmoreDetectorConstruction() { delete fPilot_logiCMOCell; }

G4bool
    AmoreDetectorConstruction::CheckSanity_CrystalModuleInfoArray(const CrystalModuleInfo *aTarget,
                                                                  size_t aSize) const {
    size_t i;
    map<string, const CrystalModuleInfo *> informationArray;
    for (i = 0; i < aSize; i++) {
        const CrystalModuleInfo *nowInfo = aTarget + i;

        auto findRes = informationArray.find(nowInfo->fName);
        if (findRes == informationArray.end())
            informationArray[nowInfo->fName] = nowInfo;
        else {
            if (!nowInfo->operator=(*(findRes->second))) {
                G4Exception(__PRETTY_FUNCTION__, "AMORE_SANITY", G4ExceptionSeverity::JustWarning,
                            ("Sanity check for a given CrystalModuleInfo array has been failed. "
                             "Please check values or data of the array. Name: " +
                             nowInfo->fName)
                                .c_str());
                return false;
            }
        }
    }
    return true;
}

G4VPhysicalVolume *AmoreDetectorConstruction::Construct() {
    // make materials if needed
    if (!materials_built) {
        ConstructMaterials();
    }

    // delete the old detector if we are constructing a new one
    if (world_phys) {
        delete world_phys;
        world_phys = NULL;
    }

    // construct the new detector
    switch (whichDetector) {
        case kDetector_AmoreDetector:
            ConstructAmoreDetector();
            break;
        default:
            CupDetectorConstruction::Construct();
            break;
    }

    return world_phys;
}
// end of AmoreDetectorConstruction::Construct()

// ----------------------------------------------------------------
G4String AmoreDetectorConstruction::GetDetectorTypeName(int i) {
    if (i < kNumGenericDetectors)
        return G4String("generic_") + CupDetectorConstruction::GetDetectorTypeName(i);
    switch (i) {
        case kDetector_AmoreDetector:
            return "AmoreDetector";
        default:
            return "detector-unknown";
    }
}

// ----------------------------------------------------------------
G4String AmoreDetectorConstruction::GetDetGeometryTypeName(eDetGeometry i) {
    switch (i) {
        case kDetector_AMoRE200:
            return "amore200";
        case kDetector_AMoRE10:
            return "amore10";
        case kDetector_AMoREPilot:
            return "amorepilot";
        case kDetector_AMoREPilotRUN5:
            return "amorepilotRUN5";
        case kDetector_AMoRE_I:
            return "amoreI";
        case kDetector_MyDetector:
            return "mydetector";
        default:
            return "detector-unknown";
    }
}

// ----------------------------------------------------------------
G4String AmoreDetectorConstruction::GetVetoGeometryTypeName(eVetoGeometry i) {
    switch (i) {
        case kVeto_200323:
            return "200323";
				case kVeto_Beam:
						return "WithBeam";
        default:
            return "unknown";
    }
}

// ----------------------------------------------------------------
G4String AmoreDetectorConstruction::GetAMoRE200PhaseName(ePhaseAMoRE200 i) {
    switch (i) {
        case kPhase1:
            return "AMoRE200_Phase1";
				case kPhase2:
						return "AMoRE200_Phase2";
        default:
            return "unknown";
    }
}

G4String AmoreDetectorConstruction::GetNShieldConfTypeName(eNShieldConf i) {
    switch (i) {
        case kNS_Naked:
            return "Naked";
        case kNS_B4C:
            return "B4C";
        case kNS_B4CnPE:
            return "B4C+PE";
        case kNS_B4CnBAcid:
            return "B4C+BoricAcid";
        case kNS_PE20:
            return "PE20cm";
        case kNS_PE10:
            return "PE10cm";
        case kNS_RealConf:
            return "RealConf";
        default:
            return "unknown";
    }
}
// ----------------------------------------------------------------
G4String AmoreDetectorConstruction::GetCavernTypeName(eCavernType i) {
    switch (i) {
        case kCavern_Toy_HemiSphere:
            return "ToyHemiSphere";
        case kCavern_Toy_Cylinder:
            return "ToyCylinder";
        case kCavern_RealModel:
            return "Real_ARF";
        default:
            return "unknown";
    }
}
// ----------------------------------------------------------------
void AmoreDetectorConstruction::ConstructMaterials() {
    using namespace CLHEP;
    G4String name, symbol;
    G4int nisotope, iz, n, nelements, natoms, ncomponents;
    G4double abundance, density, a;

    CupDetectorConstruction::ConstructMaterials();

    // ... insert any additional material definitions here
    G4NistManager *nist = G4NistManager::Instance();
    G4Isotope *Mo98 =
        new G4Isotope(name = "Molybdenum98", iz = 42, n = 98, a = 97.9054073 * g / mole);
    G4Element *_elementMo98 =
        new G4Element(name = "enriched Molybdenum", symbol = "Mo", nisotope = 1);
    _elementMo98->AddIsotope(Mo98, abundance = 100. * perCent);

    _elementB = nist->FindOrBuildElement("B");
		_elementW = nist->FindOrBuildElement("W");
		_elementTh = nist->FindOrBuildElement("Th");

		name      = "ThO2";
		density   = 10.0 * g / cm3;
		nelements = 2;
		
		_ThoriumDioxide = new G4Material(name, density, nelements);
		_ThoriumDioxide->AddElement(_elementTh, 1);
		_ThoriumDioxide->AddElement(_elementO, 2);

		name	=	"ThoriatedTungstenWire";
		density = 19.36 * g / cm3;

		_ThWire = new G4Material(name, density, ncomponents = 2, G4State::kStateSolid);
		_ThWire->AddElement(_elementW, 0.99);
		_ThWire->AddMaterial(_ThoriumDioxide, 0.01);

    name      = "MuMetal";
    density   = 7.87 * g / cm3;
    nelements = 3;

    _mumetal = new G4Material(name, density, nelements);
    _mumetal->AddElement(_elementNi, 0.80);
    _mumetal->AddElement(_elementFe, 0.15);
    _mumetal->AddElement(_elementMo98, 0.05);

    density   = 1.78 * g / cm3;
    nelements = 3;

    _stycast = new G4Material("Stycast", density, nelements);
    _stycast->AddElement(_elementC, natoms = 18);
    _stycast->AddElement(_elementH, natoms = 18);
    _stycast->AddElement(_elementO, natoms = 5);

    _vinylt = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    _vinylt->GetIonisation()->SetBirksConstant(1.50e-2 * (g / (MeV * cm2)) / _vinylt->GetDensity());

    //////////////////////////////////
    // G10 Material
    name      = "g10material";
    density   = 1.700 * g / cm3;
    nelements = 4;

    g10material = new G4Material(name, density, nelements);
    g10material->AddElement(_elementSi, natoms = 1);
    g10material->AddElement(_elementO, natoms = 2);
    g10material->AddElement(_elementC, natoms = 3);
    g10material->AddElement(_elementH, natoms = 3);

    //////////////////////////////////
    // Sillicon rubber
    //////////////////////////////////
    name      = "Silicon rubber";
    density   = 1.64 * g / cm3;
    nelements = 4;

    _SiRubber = new G4Material(name, density, nelements);
    _SiRubber->AddElement(_elementSi, natoms = 1);
    _SiRubber->AddElement(_elementO, natoms = 1);
    _SiRubber->AddElement(_elementH, natoms = 6);
    _SiRubber->AddElement(_elementC, natoms = 2);

    //////////////////////////////////
    // B4C
    //////////////////////////////////
    name      = "Boron Carbide";
    density   = 2.52 * g / cm3;
    nelements = 2;

    _B4C = new G4Material(name, density, nelements);
    _B4C->AddElement(_elementB, natoms = 4);
    _B4C->AddElement(_elementC, natoms = 1);

    //////////////////////////////////
    // Silicon rubber + 24% B4C
    //////////////////////////////////
    name    = "24% B4C + SiRubber";
    density = 1.64 * g / cm3;

    _B4CRubber24perCent = new G4Material(name, density, ncomponents = 2, G4State::kStateSolid);
    _B4CRubber24perCent->AddMaterial(_B4C, 24. * perCent);
    _B4CRubber24perCent->AddMaterial(_SiRubber, (100. - 24.) * perCent);

    //////////////////////////////////
    // Boric acid powder(H3BO3)
    //////////////////////////////////
    name    = "Boric acid powder";
    density = 1.435 * g / cm3;
    density *= 0.9; // for powder

    _BoricAcidPowder = new G4Material(name, density, ncomponents = 3, G4State::kStateSolid);
    _BoricAcidPowder->AddElement(_elementH, natoms = 3);
    _BoricAcidPowder->AddElement(_elementB, natoms = 1);
    _BoricAcidPowder->AddElement(_elementO, natoms = 3);

		/////////////////////////////////////////////
		// Silicon rubber + Boric acid powder(H3BO3)
		/////////////////////////////////////////////
		name	=	"BoricAcidRubber";
		density = 5.4 * g / 4 / cm3;

		_BoricAcidRubber = new G4Material(name, density, ncomponents = 2, G4State::kStateSolid);
		_BoricAcidRubber->AddMaterial(_BoricAcidPowder, 2.8/5.4 * 100 * perCent);
		_BoricAcidRubber->AddMaterial(_SiRubber, (100.- 2.8/5.4*100) * perCent);

    //////////////////////////////////
    // Pure boron
    //////////////////////////////////
    name    = "Pure boron";
    density = 2.52 * g / cm3;

    _PureBoron = new G4Material(name, density, ncomponents = 1, G4State::kStateSolid);
    _PureBoron->AddElement(_elementO, natoms = 1);

    //////////////////////////////////
    // Borated PE (5% boron content by weight)
    //////////////////////////////////
    name                           = "Borated PE";
    density                        = 1.01 * g / cm3;
    G4double boron_ratio_in_weight = 5. * perCent;
    G4double boron_ratio           = density * boron_ratio_in_weight / _PureBoron->GetDensity();

    _BoratedPE_5perCent = new G4Material(name, density, ncomponents = 2, G4State::kStateSolid);
    _BoratedPE_5perCent->AddMaterial(_PureBoron, boron_ratio);
    _BoratedPE_5perCent->AddMaterial(_polyethylene, 1 - boron_ratio);

    // --- Rock  SiO2 ---------------
    name      = "Rock1";
    density   = 2.7 * g / cm3;
    nelements = 2;

    _rock1 = new G4Material(name, density, nelements);
    _rock1->AddElement(_elementSi, natoms = 1);
    _rock1->AddElement(_elementO, natoms = 2);

		// ___ Iron for beam1
		name      = "Iron1";
		density   = 0.7086 * g / cm3;
		nelements = 1;

		_iron1    = new G4Material(name, density, nelements);
		_iron1->AddElement(_elementFe,natoms = 1);

		// ___ Iron for beam2
		name      = "Iron2";
		density   = 0.41 * g / cm3;
		nelements = 1;

		_iron2    = new G4Material(name, density, nelements);
		_iron2->AddElement(_elementFe,natoms = 1);

		// ___ Rebar Fe
		name      = "Rebar";
		density   = 0.3173 * g/ cm3;

		_rebar    = new G4Material(name, density, nelements);
		_rebar->AddElement(_elementFe, natoms = 1);

    name      = "Aluminium";
    density   = 2.7 * g / cm3;
    nelements = 1;

    _aluminium = new G4Material(name, density, nelements);
    _aluminium->AddElement(_elementAl, natoms = 1);

    name      = "AlProfile";
    density   = 2.7 / 2. * g / cm3;
    nelements = 1;

    _alprofile = new G4Material(name, density, nelements);
    _alprofile->AddElement(_elementAl, natoms = 1);

    name      = "Teflon2";
    density   = 2.2 * g / cm3;
    nelements = 2;
    _teflon2  = new G4Material(name, density, nelements);
    _teflon2->AddElement(_elementC, natoms = 2);
    _teflon2->AddElement(_elementF, natoms = 4);

    name      = "Copper3";
    density   = 8.96 * g / cm3;
    nelements = 1;

    _copper3 = new G4Material(name, density, nelements);
    _copper3->AddElement(_elementCu, natoms = 1);

    name      = "PEEK";
    density   = 1.32 * g / cm3;
    nelements = 3;

    _peek1 = new G4Material(name, density, nelements);
    _peek1->AddElement(_elementC, natoms = 19);
    _peek1->AddElement(_elementH, natoms = 12);
    _peek1->AddElement(_elementO, natoms = 3);

    name      = "PEEK2";
    density   = 1.32 * g / cm3;
    nelements = 3;

    _peek2 = new G4Material(name, density, nelements);
    _peek2->AddElement(_elementC, natoms = 19);
    _peek2->AddElement(_elementH, natoms = 12);
    _peek2->AddElement(_elementO, natoms = 3);

    G4Element *Pb = new G4Element(name = "Lead", symbol = "Pb", iz = 82, a = 207.19 * g / mole);
    G4Element *Sn = new G4Element(name = "Tin", symbol = "Sn", iz = 50, a = 118.71 * g / mole);
    G4Element *Pd = new G4Element(name = "Palladium", symbol = "Pd", iz = 46, a = 106.42 * g / mole);

    name      = "solder";
    density   = 8.5 * g / cm3;
    nelements = 2;

    _solder = new G4Material(name, density, nelements);
    _solder->AddElement(Sn, natoms = 60);
    _solder->AddElement(Pb, natoms = 40);

		G4Element *Zn = new G4Element(name = "Zinc", symbol = "Zn", iz = 30, a = 65.38 * g / mole);
		name    = "brass";
		density = 8.73 * g / cm3;
		nelements = 2;
		_brass = new G4Material(name, density, nelements);
		_brass->AddElement(_elementCu, natoms = 66);
		_brass->AddElement(Zn, natoms = 34);

    name      = "SiWafer";
    density   = 2.33 * g / cm3;
    nelements = 1;

    _SiWafer = new G4Material(name, density, nelements);
    _SiWafer->AddElement(_elementSi, natoms = 1);

    name      = "AuPd";
    density   = 16 * g / cm3;
    nelements = 2;

    _AuPd = new G4Material(name, density, nelements);
    _AuPd->AddElement(_elementAu, 60 * perCent);
    _AuPd->AddElement(Pd, 40 * perCent);

		// ------ Silver ------
		_elementAg = new G4Element(name = "Silver", symbol = "Ag", iz = 47, a = 107.87 * g / mole);
		name      = "Silver";
		density   = 10.49 * g / cm3;
		nelements = 1;

		_silver = new G4Material(name, density, nelements);
		_silver->AddElement(_elementAg, natoms = 1);

    _elementPb = new G4Element(name = "Lead", symbol = "Pb", iz = 82, a = 207.19 * g / mole);
    name       = "Lead";
    density    = 11.35 * g / cm3;
    nelements  = 1;

    G4Isotope *Pb210 = new G4Isotope(name = "Lead", iz = 82, n = 210, a = 207.19 * g / mole);
    _elementPb210    = new G4Element(name = "enriched Lead", symbol = "Pb", nisotope = 1);
    _elementPb210->AddIsotope(Pb210, abundance = 100. * perCent);

		// Added by Jeewon for 96.5% enriched Molybdenum LMO crystal
		//G4Isotope *Mo92 = new G4Isotope(name = "Molybdenum92", iz = 42, n = 92, a = 91.906811 * g / mole);
		//G4Isotope *Mo94 = new G4Isotope(name = "Molybdenum94", iz = 42, n = 94, a = 94.9050883 * g / mole);
		//G4Isotope *Mo95 = new G4Isotope(name = "Molybdenum95", iz = 42, n = 95, a = 94.9058421 * g / mole);
		//G4Isotope *Mo96 = new G4Isotope(name = "Molybdenum96", iz = 42, n = 96, a = 95.9046795 * g / mole);
		//G4Isotope *Mo97 = new G4Isotope(name = "Molybdenum97", iz = 42, n = 97, a = 96.9060215 * g / mole);
    G4Isotope *Mo100 = new G4Isotope(name = "Molybdenum100", iz = 42, n = 100, a = 99.907477 * g / mole);
    G4Element *_elementMo = new G4Element(name = "enriched Molybdenum", symbol = "Mo", nisotope = 1);
    _elementMo->AddIsotope(Mo100, abundance = 100 * perCent); 
		/*
    _elementMo->AddIsotope(Mo100, abundance = 96.5 * perCent); 
		_elementMo->AddIsotope(Mo98, abundance = 3.35 * perCent);
		_elementMo->AddIsotope(Mo97, abundance = 0.028 * perCent);
		_elementMo->AddIsotope(Mo96, abundance = 0.031 * perCent);
		_elementMo->AddIsotope(Mo95, abundance = 0.027 * perCent);
		_elementMo->AddIsotope(Mo94, abundance = 0.014 * perCent);
		_elementMo->AddIsotope(Mo92, abundance = 0.021 * perCent);
		*/

    name      = "PbMoO4";
    density   = 6.95 * g / cm3;
    nelements = 3;
    PbMoO4    = new G4Material(name, density, nelements);
    PbMoO4->AddElement(_elementPb210, natoms = 1);
    PbMoO4->AddElement(_elementMo, natoms = 1);
    PbMoO4->AddElement(_elementO, natoms = 4);
    PbMoO4->GetIonisation()->SetBirksConstant(0.117 * mm / MeV);

    _elementLi = nist->FindOrBuildElement("Li");

    name      = "Li2MoO4";
    //density   = 2.66 * g / cm3;
		density = 3.03 * g / cm3;
		//density = 4.33 * g / cm3;
    nelements = 3;
    _Li2MoO4  = new G4Material(name, density, nelements);
    _Li2MoO4->AddElement(_elementLi, natoms = 2);
    _Li2MoO4->AddElement(_elementMo, natoms = 1);
    _Li2MoO4->AddElement(_elementO, natoms = 4);
    _Li2MoO4->GetIonisation()->SetBirksConstant(0.117 * mm / MeV);



    const char *basic_fname = "amore_materials.dat";
    ifstream amore_mat;

    if (getenv("AmoreDATA") != NULL)
        amore_mat.open((G4String(getenv("AmoreDATA")) + "/" + G4String(basic_fname)).c_str());
    else
        amore_mat.open(("data/" + G4String(basic_fname)).c_str());

    int error_cnt = CupInputDataReader::ReadMaterials(amore_mat);
    if (error_cnt)
        G4cerr << "Error was occured in MaterialPropertiesTable generation step, Error cnt.: "
               << error_cnt << G4endl;
}
// == Construct Geometry for main detector ======================================
// uses parameters from database or file
void AmoreDetectorConstruction::ConstructAmoreDetector() {
    using namespace std;
    // -- database
    CupParam &db(CupParam::GetDB());
    // add spherical-geometry-specific parameters to parameter list

    switch (whichDetGeometry) {
        case kDetector_AMoRE200:
            if (getenv("AmoreDATA") != NULL)
                db.ReadFile((G4String(getenv("AmoreDATA")) + "/settings_amore200.dat").c_str());
            else
                db.ReadFile("data/settings_amore200.dat");
            ConstructAMoRE200();
            break;
        case kDetector_AMoREPilot:
            ConstructAMoREPilot();
            break;
        case kDetector_AMoREPilotRUN5:
            ConstructAMoREPilotRUN5();
            break;
        case kDetector_AMoRE_I: {
            using namespace AmoreDetectorStaticInfo::AMoRE_I;
            if (getenv("AmoreDATA") != NULL)
                db.ReadFile((G4String(getenv("AmoreDATA")) + "/settings_amoreI.dat").c_str());
            else
                db.ReadFile("data/settings_amoreI.dat");

            size_t totalArraySize = maxModuleNumInTower * totalNumOfTower;

            const CrystalModuleInfo **crystalIArray = new const CrystalModuleInfo *[totalArraySize];
            for (size_t i = 0; i < totalNumOfTower; i++) {
                for (size_t j = 0; j < maxModuleNumInTower; j++) {
                    crystalIArray[i * maxModuleNumInTower + j] = &crystalModuleInfoList[i][j];
                }
            }
            CheckSanity_CrystalModuleInfoArray(*crystalIArray, totalArraySize);
            delete[] crystalIArray;

            ConstructAMoRE_I();
            break;
        }
        case kDetector_AMoRE10:
            ConstructAMoRE10();
            break;
        case kDetector_MyDetector:
            ConstructMyDetector();
            break;
        default:
            G4cerr << "ERROR: INVALID VALUE for AmoreDetectorConstruction.whichDetGeometry"
                   << G4endl;
            // fall through and make a GenericLAND
    }
}
void AmoreDetectorConstruction::ConstructSDandField() {
    switch (whichDetGeometry) {
        case kDetector_AMoRE200: {
            ConstructAMoRE200_SDandField();
        } break;
        case kDetector_AMoREPilot:
            ConstructAMoREPilot_SDandField(); 
            break;
        case kDetector_AMoREPilotRUN5:
            ConstructAMoREPilotRUN5_SDandField(); 
            break;
        case kDetector_AMoRE_I: {
            Construct_I_SDandField();
        } break;
        case kDetector_AMoRE10:
            break;
        case kDetector_MyDetector:
            break;
        default:
            G4cerr << "ERROR: INVALID VALUE for AmoreDetectorConstruction.whichDetGeometry"
                   << G4endl;
            // fall through and make a GenericLAND
    }
}
