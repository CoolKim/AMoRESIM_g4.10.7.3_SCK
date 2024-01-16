/**@file AmoreDetectorConstruction.hh
 * @brief A source code file of AmoreDetectorConstruction main class for AMoRE simulation
 * @detail Original by Glenn Horton-Smith, December 2004.  (Based on earlier work first written Dec
 1999.) Some more materials and their classification by Dario Motta, Jan 2005.
 */

/**
 * @brief A class for initialization and controlling of AMoRE-specific geometries
 * @details This class provides various methods for constructing geometries for AMoRE simulation,
 * judging methods of incident into border and so on.
 * @author AMoRE Collaboration
 * */

#ifndef AmoreDetectorConstruction_HH
#define AmoreDetectorConstruction_HH 1
#include <set>

#include "G4GeometryTolerance.hh"
#include "G4NavigationHistory.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"

#include "AmoreSim/AmoreDetectorStaticType.hh"
#include "AmoreSim/AmoreModuleHit.hh"
#include "CupSim/CupDetectorConstruction.hh"
#include "G4Version.hh"

class G4SurfaceProperty;

class AmoreDetectorConstruction : public CupDetectorConstruction {
	private:
		G4bool fEnable_OriginalGeom;
		G4bool fEnable_Scintillator;
		G4bool fEnable_Gantry;
		G4bool fEnable_InnerDetector;
		G4bool fEnable_Innermost;
		G4bool fEnable_NeutronShield;
		//G4bool fEnable_AdditionalPEShield;

		G4int fNeutShieldingConf;
		G4bool fNeutronMode;
		G4bool fRockgammaMode;
		G4bool fAdditionalPE;
		G4bool fDbgMsgOn;
		G4bool OverlapCheck;

		std::set<AmoreModuleSDInfo> fModuleSDInfos;

	protected:
		using CrystalModuleInfo = AmoreDetectorStaticInfo::CrystalModuleInfo;

		G4VPhysicalVolume *fCavernPhysical;
		G4VPhysicalVolume *fRockPhysical;
		G4VPhysicalVolume *fFloorPhysical;

		// AMoRE-Pilot LV and PV for SD information & Geometry settings
		G4bool fPilot_Enable_TargetRoom;
		G4bool fPilot_Enable_Mumetal;

		G4LogicalVolume *fPilot_TopScint_BoxLogical;
		G4LogicalVolume *fPilot_SideFBScint_BoxLogical;
		G4LogicalVolume *fPilot_SideLRScint_BoxLogical;

		G4LogicalVolume *fPilot_TopScint_FlatTrapLogical;
		G4LogicalVolume *fPilot_SideFBScint_FlatTrapLogical;
		G4LogicalVolume *fPilot_SideLRScint_FlatTrapLogical;

		G4LogicalVolume *fPilot_SideFBScint_PMTTrapLogical;
		G4LogicalVolume *fPilot_TopScint_PMTTrapLogical;
		G4LogicalVolume *fPilot_SideLRScint_PMTTrapLogical;

		G4LogicalVolume **fPilot_logiCMOCell;
		G4LogicalVolume *fPilot_logicPMT;
		G4LogicalVolume *fPilot_logicPMTVacu;
		G4LogicalVolume *fPilot_logiGOLDa;
		G4LogicalVolume *fPilot_logiGOLDb;

		// AMoRE-I LV and PV for SD information & Geometry settings
		G4bool fI_Enable_SuperConductingShield;
		G4bool fI_Enable_CrystalArray;

		G4LogicalVolume *fI_TopScint_BoxLogical;
		G4LogicalVolume *fI_SideFBScint_BoxLogical;
		G4LogicalVolume *fI_MufflerFBScint_BoxLogical;
		G4LogicalVolume *fI_MufflerLRScint_BoxLogical;
		G4LogicalVolume *fI_SideLRScint_BoxLogical;

		G4LogicalVolume *fI_TopScint_FlatTrapLogical;
		G4LogicalVolume *fI_SideFBScint_FlatTrapLogical;
		G4LogicalVolume *fI_MufflerFBScint_FlatTrapLogical;
		G4LogicalVolume *fI_MufflerLRScint_FlatTrapLogical;
		G4LogicalVolume *fI_SideLRScint_FlatTrapLogical;

		G4LogicalVolume *fI_TopScint_PMTTrapLogical;
		G4LogicalVolume *fI_SideFBScint_PMTTrapLogical;
		G4LogicalVolume *fI_MufflerFBScint_PMTTrapLogical = nullptr;
		G4LogicalVolume *fI_MufflerLRScint_PMTTrapLogical = nullptr;
		G4LogicalVolume *fI_SideLRScint_PMTTrapLogical;

		G4Region *fI_DetectorModuleRegion;
		G4Region *fI_crystalsRegion;

		std::set<G4LogicalVolume *> fI_CrystalLVs;
		std::set<G4LogicalVolume *> fI_GeWaferLVs;
		std::set<G4LogicalVolume *> fI_GeWaferGoldFilmLVs;
		std::set<G4LogicalVolume *> fI_CrystalGoldFilmLVs;

		G4LogicalVolume *fI_logicPMT;
		G4LogicalVolume *fI_logicPMTVacu;

		G4VisAttributes *fI_CopperDefault_VisAttr;
		G4VisAttributes *fI_LeadDefault_VisAttr;
		G4VisAttributes *fI_GoldDefault_VisAttr;
		G4VisAttributes *fI_Vikuiti_VisAttr;
		G4VisAttributes *fI_CopperFrame_VisAttr;
		G4VisAttributes *fI_OpticalFrame_VisAttr;
		G4VisAttributes *fI_GeWafer_VisAttr;
		G4VisAttributes *fI_LeadShield_VisAttr;
		G4VisAttributes *fI_Reflector_VisAttr;
		G4VisAttributes *fI_Crystal_VisAttr;
		G4VisAttributes *fI_Bolt_VisAttr;
		G4VisAttributes *fI_FT2850_VisAttr;
		G4VisAttributes *fI_PEEK_VisAttr;
		G4VisAttributes *fI_Clamp_VisAttr;
		G4VisAttributes *fI_ClampBolt_VisAttr;
		G4VisAttributes *fI_PCB_VisAttr;
		G4VisAttributes *fI_Epoxy_VisAttr;
		G4VisAttributes *fI_PE_VisAttr;
		G4VisAttributes *fI_BoratedPE_VisAttr;
		G4VisAttributes *fI_BoricAcid_VisAttr;
		G4VisAttributes *fI_BoricAcidHousing_VisAttr;
		G4VisAttributes *fI_B4CRubber_VisAttr;
		G4VisAttributes *fI_Invisible_VisAttr;
		G4VisAttributes *fI_BrassBolt_VisAttr;

		// AMoRE-II LV and PV for SD information
		G4int f200_VetoTotCNum, f200_HatVetoTotCNum, f200_TotalVetoCNum;
		G4int f200_TotCrystalNum, f200_TotTowerNum;
		G4LogicalVolume *f200_logiVetoPSO;
		G4LogicalVolume *f200_logiVetoPSI;
		G4LogicalVolume *f200_logiHatPSO;
		G4LogicalVolume *f200_logiHatPSI;
		//G4LogicalVolume *f200_logiCrystalCell;
		G4LogicalVolume *f200_logiCrystalCell[1000];
		G4VPhysicalVolume *f200_physGeWafer;
		G4VPhysicalVolume *f200_physVacDisk;
		G4VPhysicalVolume *f200_HatVetoMaterialPV;
		G4VPhysicalVolume *f200_PSO_PV;
		G4VPhysicalVolume *f200_PSI_PV;
		G4VPhysicalVolume *f200_HatPSO_PV;
		G4VPhysicalVolume *f200_HatPSI_PV;
		G4VPhysicalVolume *f200_VetoActiveMaterialPV;
		G4VPhysicalVolume *f200_FloorPEPhysical;
		G4VPhysicalVolume *f200_CeilingPEPhysical;
		G4VPhysicalVolume *f200_RealPEPhysical;
		G4VPhysicalVolume *f200_VetoMaterialPhysical;
		G4VPhysicalVolume *f200_AirBufferPhysical;
		G4VPhysicalVolume *f200_OVCPhysical;

	public:
		typedef enum {
			kDetector_AmoreDetector = CupDetectorConstruction::kNumGenericDetectors,
			kNumDetectors
		} eDetector;

		typedef enum {
			kDetector_AMoREPilot = 0,
			kDetector_AMoREPilotRUN5,
			kDetector_AMoRE_I,
			kDetector_AMoRE10,
			kDetector_AMoRE200,
			kDetector_MyDetector,
			kNumDetGeometries
		} eDetGeometry;

		typedef enum{
			kPhase1 = 0,
			kPhase2,
			kNumAMoRE200Phase
		} ePhaseAMoRE200;

		typedef enum {
			kVeto_200323 = 0,
			kVeto_Beam,
			kNumVetoGeometries
		} eVetoGeometry;

		typedef enum {
			kCavern_Toy_HemiSphere = 0,
			kCavern_Toy_Cylinder,
			kCavern_RealModel,
			kNumCavernTypes
		} eCavernType;

		typedef enum {
			kNS_Naked = 0,
			kNS_B4C,
			kNS_B4CnPE,
			kNS_B4CnBAcid,
			kNS_PE20,
			kNS_PE10,
			kNS_RealConf,
			kNumNSConfs
		} eNShieldConf;

		AmoreDetectorConstruction();          // constructor
		virtual ~AmoreDetectorConstruction(); // destructor

		virtual G4VPhysicalVolume *Construct(); // make the volumes, return ptr to world
		virtual void ConstructSDandField();

		virtual int GetNumDetectorTypes() { return kNumDetectors; }
		virtual G4String GetDetectorTypeName(int i);

		static eDetGeometry GetNumDetGeometryTypes() { return kNumDetGeometries; }
		virtual G4String GetDetGeometryTypeName(eDetGeometry i);
		eDetGeometry GetWhichDetGeometry(void) { return whichDetGeometry; }
		static inline eDetGeometry GetDetGeometryType() { return whichDetGeometry; }
		virtual void SetWhichDetGeometry(eDetGeometry w) { whichDetGeometry = w; }

		static inline int GetNumVetoGeometryTypes() { return kNumVetoGeometries; };
		virtual G4String GetVetoGeometryTypeName(eVetoGeometry i);
		static inline eVetoGeometry GetVetoGeometryType() { return whichVetoGeometry; }
		inline void SetWhichVetoGeometry(eVetoGeometry w) { whichVetoGeometry = w; }

		static inline int GetAMoRE200PhaseNumber() { return kNumAMoRE200Phase; };
		virtual G4String GetAMoRE200PhaseName(ePhaseAMoRE200 i);
		static inline ePhaseAMoRE200 GetAMoRE200PhaseType() { return whichAMoRE200Phase; }
		inline void SetWhichAMoRE200Phase(ePhaseAMoRE200 w) { whichAMoRE200Phase = w;}

		static inline int GetNumNShieldConfTypes() { return kNumNSConfs; };
		virtual G4String GetNShieldConfTypeName(eNShieldConf i);
		static inline eNShieldConf GetNShieldConfType() { return whichNShieldingConf; }
		inline void SetNShieldConfType(eNShieldConf w) { whichNShieldingConf = w; }

		static inline int GetNumCavernType() { return kNumCavernTypes; }
		virtual G4String GetCavernTypeName(eCavernType i);
		static inline eCavernType GetCavernType() { return whichCavernType; }
		virtual void SetWhichCavernType(eCavernType w) { whichCavernType = w; }

		inline void SetEnableOriginalGeometry(G4bool a) { fEnable_OriginalGeom = a; }
		inline void SetEnableScintillator(G4bool a) { fEnable_Scintillator = a; }
		inline void SetEnableGantry(G4bool a) { fEnable_Gantry = a; }
		inline void SetEnableInnerDetector(G4bool a) { fEnable_InnerDetector = a; }
		inline void SetEnableInnermost(G4bool a) { fEnable_Innermost = a; }
		inline void SetEnableRealConf(G4bool a) { fEnable_NeutronShield = a; }
		inline void SetNeutronMode(G4bool a) { fNeutronMode = a; }
		inline void SetRockgammaMode(G4bool a) { fRockgammaMode = a; }
		inline void SetAdditionalPE(G4bool a) { fAdditionalPE = a; }
		inline void SetOverlapCheck(G4bool a) { OverlapCheck = a; }
		inline void SetDebugMessage(G4bool a) { fDbgMsgOn = a; }

		inline G4bool GetEnableOriginalGeometry() const { return fEnable_OriginalGeom; }
		inline G4bool GetEnableScintillator() const { return fEnable_Scintillator; }
		inline G4bool GetEnableGantry() const { return fEnable_Gantry; }
		inline G4bool GetEnableInnerDetector() const { return fEnable_InnerDetector; }
		inline G4bool GetEnableInnermost() const { return fEnable_Innermost; }
		inline G4bool GetEnableNeutronShield() const { return fEnable_NeutronShield; }
		//inline G4bool GetEnableAdditionalPEShield() const { return fEnable_AdditionalPEShield; }
		inline G4bool GetNeutronMode() const { return fNeutronMode; }
		inline G4bool GetRockgammaMode() const { return fRockgammaMode; }
		inline G4bool GetAdditionalPE() const { return fAdditionalPE; }
		inline G4bool GetOverlapCheck() const { return OverlapCheck; }
		inline G4bool GetDebugMessage() const { return fDbgMsgOn; } 

		// For common uses
		inline bool JudgeBorderIncident(const G4Step *aStep, const G4VPhysicalVolume *const *aTargetPV,
				G4int aNumOfTarget = 1) const;
		static inline bool DoesNavHistContainPV(const G4NavigationHistory *aNavHist,
				const G4VPhysicalVolume *const *aTargetPV,
				G4int aNumOfTarget = 1);
		inline bool Judge_CavernBorder(const G4Step *aStep) const;

		const std::set<AmoreModuleSDInfo> &GetModuleSDInfoList() const { return fModuleSDInfos; }

		// For AMoRE Pilot (can be moved to common section in the future)
		inline void Set_Pilot_EnableTargetRoom(G4bool a) { fPilot_Enable_TargetRoom = a; }
		inline void Set_Pilot_EnableMumetal(G4bool a) { fPilot_Enable_Mumetal = a; }
		inline void Set_PilotRUN5_EnableTargetRoom(G4bool a) { fPilot_Enable_TargetRoom = a; }
		inline void Set_PilotRUN5_EnableMumetal(G4bool a) { fPilot_Enable_Mumetal = a; }

		inline G4bool Get_Pilot_EnableTargetRoom() const { return fPilot_Enable_TargetRoom; }
		inline G4bool Get_Pilot_EnableMumetal() const { return fPilot_Enable_Mumetal; }
		inline G4bool Get_PilotRUN5_EnableTargetRoom() const { return fPilot_Enable_TargetRoom; }
		inline G4bool Get_PilotRUN5_EnableMumetal() const { return fPilot_Enable_Mumetal; }

		inline G4VPhysicalVolume *GetCavernPV() const { return fCavernPhysical; }
		inline G4VPhysicalVolume *GetRockPV() const { return fRockPhysical; }
		inline G4VPhysicalVolume *GetFloorPV() const { return fFloorPhysical; }
		inline G4VPhysicalVolume *GetFloorPEPV() const { return f200_FloorPEPhysical; }
		inline G4VPhysicalVolume *GetCeilingPEPV() const { return f200_CeilingPEPhysical; }
		inline G4VPhysicalVolume *GetRealPEPV() const { return f200_RealPEPhysical; }

		// For AMoRE I (can be moved to common section in the future)
		inline void Set_I_EnableSuperConductingShield(G4bool a) { fI_Enable_SuperConductingShield = a; }
		inline void Set_I_EnableCrystalArray(G4bool a) { fI_Enable_CrystalArray = a; }

		inline G4bool Get_I_EnableSuperConductingShield() const {
			return fI_Enable_SuperConductingShield;
		}
		inline G4bool Get_I_EnableCrystalArray() const { return fI_Enable_CrystalArray; }

		// For AMoRE 200 (can be moved to common section in the future)
		inline bool Judge_200_PEBorderForRealCT(const G4Step *aStep) const;
		inline bool Judge_200_PEBorderForHSCT(const G4Step *aStep) const;
		inline bool Judge_200_VABorder(const G4Step *aStep) const;
		inline bool Judge_200_VMBorder(const G4Step *aStep) const;
		inline bool Judge_200_ABBorder(const G4Step *aStep) const;
		inline bool Judge_200_OVCBorder(const G4Step *aStep) const;

		//protected:
		G4Material *_mumetal, *_stycast, *_vinylt, *g10material, *PbMoO4;
		G4Material *_B4CRubber24perCent;
		G4Material *_SiRubber;
		G4Material *_B4C;
		G4Material *_BoricAcidPowder;
		G4Material *_BoricAcidRubber;
		G4Material *_BoratedPE_5perCent;
		G4Material *_PureBoron;
		G4Material *_aluminium;
		G4Material *_alprofile;
		G4Material *_teflon1;
		G4Material *_teflon2;
		G4Material *_copper3;
		G4Material *_peek1;
		G4Material *_peek2;
		G4Material *_solder;
		G4Material *_brass;
		G4Material *_SiWafer;
		G4Material *_AuPd;
		G4Material *_Li2MoO4;
		G4Material *_silver;
		G4Material *_rock1;
		G4Material *_ThWire;
		G4Material *_ThoriumDioxide;
		G4Material *_iron1;
		G4Material *_iron2;
		G4Material *_rebar;

		G4Element *_elementB;
		G4Element *_elementPb210;
		G4Element *_elementLi;
		G4Element *_elementAg;
		G4Element *_elementW;
		G4Element *_elementTh;

		//pilot RUN5 elements
		//G4Element *Sn, *Pb, *Pd;
		//G4Material *_copper2, *mylar, *_siwafer, *_aupd, *peek, *solder;


	protected:
		G4bool CheckSanity_CrystalModuleInfoArray(const CrystalModuleInfo *aTarget, size_t aSize) const;

		void ConstructMaterials();
		void ConstructAmoreDetector();
		void ConstructAMoRE10();  ///< make the AMoRE10 detector

		void ConstructAMoRE200(); ///< make the AMoRE200 detector
		void ConstructAMoRE200_ID(G4LogicalVolume *aWorkAreaLV);       ///< make the AMoRE200 inner detector (cryostat)
		G4LogicalVolume *MakeModule(G4Material *towerMat, G4Material *crystalMat, G4Material *reflectorMat, 
				G4Material *frameMat, G4Material *frameMat1, G4Material *clampMat, G4Material *waferMat, 
				G4Material *filmMat, G4int TowerNum, G4int ModuleNum);
		G4LogicalVolume *MakeTower_phase2(G4Material *towerMat, G4Material *crystalMat, G4Material *reflectorMat, 
				G4Material *frameMat, G4Material *clampMat, G4Material *waferMat, G4Material *filmMat, G4int TowerNum);
		G4LogicalVolume *ConstructAMoRE200_OD(); ///< make the AMoRE200 outer detector (water tank)

		void ConstructMyDetector();              ///< make the MyDetector

		void ConstructAMoREPilot(); ///< make the AMoRE-Pilot detector RUN6 setting
		void ConstructAMoREPilotRUN5(); ///< make the AMoRE-Pilot detector RUN6 setting

		void ConstructAMoRE_I();    ///< make the AMoRE10 detector

		void ConstructAMoRE10_SDandField(); // make the AMoRE200 detector
		void ConstructAMoRE200_SDandField();
		void ConstructMyDetector_SDandField(); // make the MyDetector

		void ConstructAMoREPilot_SDandField(); // make the AMoRE-Pilot detector
		void ConstructAMoREPilotRUN5_SDandField(); // make the AMoRE-Pilot detector

		// For AMoRE I (can be moved to common section in the future)
		G4LogicalVolume *Build_I_DetectorArrayUpperPanel(G4bool aRealistic, G4Material *aFrameMaterial,
				G4Material *aSpaceMaterial);
		G4LogicalVolume *Build_I_DetectorArrayBottomPanel(G4bool aRealistic, G4Material *aFrameMaterial,
				G4Material *aSpaceMaterial);
		G4LogicalVolume *Build_I_PhotonDetector(G4bool aRealistic,G4int aType, G4Material *aFrameMaterial,
				G4Material *aSpaceMaterial,
				AmoreModuleSDInfo *aMSDForThisModule);
		G4LogicalVolume *Build_I_CopperFrameUpper(G4int aType);
		G4LogicalVolume *Build_I_CopperFrameBottom(G4int aType);
		G4LogicalVolume *Build_I_SingleDetectorModule(G4bool aRealistic, G4int aType, const CrystalModuleInfo &aCrystal,
				G4Material *aSpaceMaterial,
				AmoreModuleSDInfo *aMSDForThisModule,
				G4bool aIsFloor);
		G4LogicalVolume *Build_I_DetectorArray(G4int aType, G4Material *aArrayFrameMaterial,
				G4Material *aSpaceMaterial);
		G4LogicalVolume *Build_I_SCMagnetShieldAt(G4LogicalVolume *aMotherLV, G4ThreeVector aTlate,
				G4int aType, G4Material *aFrameMaterial,
				G4Material *aSpaceMaterial,
				G4Material *aShieldMaterial);

		G4LogicalVolume *Build_I_MuonVetoPMT(G4Material *aShieldMaterial, G4Material *aGreaseMaterial);
		G4LogicalVolume *Build_I_TopScintillator(G4Material *aScintMaterial,
				G4Material *aReflectorMaterial,
				G4SurfaceProperty *aSurfaceProperty);
		G4LogicalVolume *Build_I_Side_LRScintillator(G4Material *aScintMaterial,
				G4Material *aReflectorMaterial,
				G4SurfaceProperty *aSurfaceProperty);
		G4LogicalVolume *Build_I_Side_FBScintillator(G4Material *aScintMaterial,
				G4Material *aReflectorMaterial,
				G4SurfaceProperty *aSurfaceProperty);
		G4LogicalVolume *Build_I_Muffler_FBScintillator(G4Material *aScintMaterial,
				G4Material *aReflectorMaterial,
				G4Material *aSideReflectorMaterial,
				G4SurfaceProperty *aSurfaceProperty,
				G4SurfaceProperty *aSideReflectorSurfProp);
		G4LogicalVolume *Build_I_Muffler_LRScintillator(G4Material *aScintMaterial,
				G4Material *aReflectorMaterial,
				G4Material *aSideReflectorMaterial,
				G4SurfaceProperty *aSurfaceProperty,
				G4SurfaceProperty *aSideReflectorSurfProp);

		void Place_I_PCBs(G4bool aRealistic, G4LogicalVolume *aWhere, G4int aType, 
				G4ThreeVector aTlate);
		void Place_I_GeometriesInDetModuleOf(G4LogicalVolume *aWhere, G4int aType,
				const CrystalModuleInfo &aCrystal,
				G4bool aPlacePartsForJustBelow);
		void Place_I_GeometriesInRooftopForTowerOf(G4LogicalVolume *aWhere, G4int aType,
				G4ThreeVector aOriginOfTower);
		void Place_I_RealNeutronShieldingAt(G4LogicalVolume *aWhere, G4ThreeVector aTlate, G4int aType);

		void Construct_I_SDandField(); // make the AMoRE-Pilot detector

		G4VSolid *BuildCopperFrameSolid(G4double, G4double, G4double, G4double);
		G4VSolid *BuildPhotonFrameSolid(G4double, G4double);

		static eDetGeometry whichDetGeometry;
		static eVetoGeometry whichVetoGeometry;
		static eCavernType whichCavernType;
		static eNShieldConf whichNShieldingConf;
		static ePhaseAMoRE200 whichAMoRE200Phase;

		class AmoreDetectorMessenger *AmoreMessenger;
};

bool AmoreDetectorConstruction::DoesNavHistContainPV(const G4NavigationHistory *aNavHist,
		const G4VPhysicalVolume *const *aTargetPV,
		G4int aNumOfTarget) {
	G4int tHavHistDepth = aNavHist->GetDepth();
	for (int i = tHavHistDepth; i >= 0; i--) {
		for (int j = 0; j < aNumOfTarget; j++) {
			if (aNavHist->GetVolume(i) == aTargetPV[j]) return true;
		}
	}
	return false;
}

bool AmoreDetectorConstruction::JudgeBorderIncident(const G4Step *aStep,
		const G4VPhysicalVolume *const *aTargetPV,
		G4int aNumOfTarget) const {
	G4bool tPreContainsTargetPV = false;
	G4StepPoint *tPreStepPt     = aStep->GetPreStepPoint();
	if (tPreStepPt != nullptr)
		tPreContainsTargetPV =
			DoesNavHistContainPV(tPreStepPt->GetTouchable()->GetHistory(), aTargetPV, aNumOfTarget);

	G4bool tPostContainsTargetPV = false;
	G4StepPoint *tPostStepPt     = aStep->GetPostStepPoint();
	if (tPostStepPt != nullptr)
		tPostContainsTargetPV = DoesNavHistContainPV(tPostStepPt->GetTouchable()->GetHistory(),
				aTargetPV, aNumOfTarget);

	return (!tPreContainsTargetPV && tPostContainsTargetPV);
}

// For AMoRE-II simulation (don't use those method)

bool AmoreDetectorConstruction::Judge_CavernBorder(const G4Step *aStep) const {
	return JudgeBorderIncident(aStep, &fCavernPhysical, 1);
}

bool AmoreDetectorConstruction::Judge_200_PEBorderForRealCT(const G4Step *aStep) const {
	return JudgeBorderIncident(aStep, &f200_RealPEPhysical, 1);
}

bool AmoreDetectorConstruction::Judge_200_PEBorderForHSCT(const G4Step *aStep) const {
	const G4VPhysicalVolume *lTargetList[2] = {f200_CeilingPEPhysical, f200_FloorPEPhysical};
	return JudgeBorderIncident(aStep, lTargetList, 2);
}

bool AmoreDetectorConstruction::Judge_200_VABorder(const G4Step *aStep) const {
	return JudgeBorderIncident(aStep, &f200_VetoActiveMaterialPV, 1);
}

bool AmoreDetectorConstruction::Judge_200_VMBorder(const G4Step *aStep) const {
	return JudgeBorderIncident(aStep, &f200_VetoMaterialPhysical, 1);
}

bool AmoreDetectorConstruction::Judge_200_ABBorder(const G4Step *aStep) const {
	return JudgeBorderIncident(aStep, &f200_AirBufferPhysical, 1);
}
bool AmoreDetectorConstruction::Judge_200_OVCBorder(const G4Step *aStep) const {
	return JudgeBorderIncident(aStep, &f200_OVCPhysical, 1);
}
#endif
