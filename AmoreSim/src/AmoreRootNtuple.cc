// AmoreRootNtuple.cc
// 29-Sep-2000 Bill Seligman
// This is an example class for booking histograms of hit information
// in GEANT4.
//
// 2015-05-26 Modified by E.J.Jeon for CUP detector
//
//

#include <sstream>
#include <string>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "AmoreSim/AmoreModuleSD.hh"
#include "AmoreSim/AmoreRootNtuple.hh"
#include "AmoreSim/AmoreRootNtupleMessenger.hh"
#include "AmoreSim/AmoreScintSD.hh"
#include "AmoreSim/AmoreScintillation.hh"
#include "AmoreSim/AmoreTrackInformation.hh"
#include "CupSim/CupScintHit.hh"
#include "CupSim/CupParam.hh"
#include "CupSim/CupPrimaryGeneratorAction.hh"
#include "CupSim/CupScintSD.hh"
#include "CupSim/CupScintillation.hh"
#include "CupSim/CupVertexGen.hh"
#include "CupSim/CupVetoHit.hh"
#include "CupSim/CupVetoSD.hh"

// Include files for ROOT.
#include "Rtypes.h"

// Include files for the G4 classes
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4IonTable.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4VHitsCollection.hh"
#include "G4VProcess.hh"
#include "globals.hh"
using namespace std;

TROOT theROOT("AmoreSim/Amoresim", "AMoRE Geant4 simulation output tree");

AmoreRootNtuple::AmoreRootNtuple()
    : CupRootNtuple(), fRecordedEvt(0), fRecordWithCut(false), fRecordPrimary(false),
      myAmoreNtupleMessenger(nullptr), fEvtInfos(nullptr), fPrimAtCB(nullptr), fPrimAtOVC(nullptr),
			fOutputForPrim(nullptr) {
    fModuleArray           = nullptr;
    EndTrackList           = new std::vector<TTrack *>;
    myAmoreNtupleMessenger = new AmoreRootNtupleMessenger(this);
    fEvtInfo_VolumeTbl     = new std::map<std::string, int>;
    fEvtInfo_VolumeTbl->clear();
}

AmoreRootNtuple::~AmoreRootNtuple() {
    CloseFile();
    ClearET();

    delete EndTrackList;
    delete myAmoreNtupleMessenger;
    delete fEvtInfo_VolumeTbl;
}

void AmoreRootNtuple::OpenFile(const G4String filename, G4bool outputmode) {
    if (fRecordPrimary) {
        fOutputForPrim = new TFile((filename + "_prim" + ".root").c_str(), "RECREATE",
                                   "Output file for primrary generation");

        fValuesForPrim[0] = fValuesForPrim[1] = fValuesForPrim[2] = fValuesForPrim[3] =
            fValuesForPrim[4] = fValuesForPrim[5] = fValuesForPrim[6] = fValuesForPrim[7] =
                fValuesForPrim[8] = fValuesForPrim[9] = fValuesForPrim[10] = fValuesForPrim[11] =
                    0.;

        fEvtInfo_EvtID        = 0;
        fEvtInfo_EdepOV[0]       = 0;
        fEvtInfo_EdepOV[1]       = 0;
        fEvtInfo_HittedCMONum = 0;
        fEvtInfo_InciAtCB     = 0;
        fEvtInfo_InciAtOVC     = 0;
    }

    CupRootNtuple::OpenFile(filename, outputmode);
}

void AmoreRootNtuple::CreateTree() {
    CupRootNtuple::CreateTree();

    // Sensitive Detector for crystal detectors
    AmoreScintSD *tgsd;
    tgsd = (AmoreScintSD *)(CupRootNtuple::sdman->FindSensitiveDetector("/CupDet/TGSD"));
    AmoreModuleSD *moduleSD;
    moduleSD = (AmoreModuleSD *)(CupRootNtuple::sdman->FindSensitiveDetector("/CupDet/MDSD"));
    if (tgsd != nullptr) {
        fROOTOutputTree->Branch("TGSD", &(Ctgsd), 256000, 2);
    } else if (moduleSD != nullptr) {
        const AmoreDetectorConstruction *theDetCons =
            static_cast<const AmoreDetectorConstruction *>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        const std::set<AmoreModuleSDInfo> &theSDInfoList = theDetCons->GetModuleSDInfoList();
        std::vector<std::string> nameList(theSDInfoList.size());
        for (const auto &nowSDInfo : theSDInfoList) {
            nameList[nowSDInfo.fModuleID] = nowSDInfo.fModuleName;
        }
        fModuleArray = new DetectorArray_Amore(
            "AmoreModuleSD",
            "A DetectorArray object of detector module array for simulation of AMoRE", nameList);
        fROOTOutputTree->Branch("MDSD", &fModuleArray, 512000, 2);
    }

    fROOTOutputTree->Branch("EndTrack", &EndTrackList);

    if (fRecordPrimary) {
        fOutputForPrim->cd();
        fEvtInfos = new TTree("EvtInfos", "Event information for primary records");
        fEvtInfos->SetDirectory(fOutputForPrim);

        fEvtInfos->Branch("EvtID", &fEvtInfo_EvtID, "EvtID/I");
        fEvtInfos->Branch("HittedCMONum", &fEvtInfo_HittedCMONum, "HittedCMONum/I");
        fEvtInfos->Branch("VolTbl", fEvtInfo_VolumeTbl);
        fEvtInfos->Branch("InciAtCB", &fEvtInfo_InciAtCB, "InciAtCB/I");

        fPrimAtCB = new TNtupleD("PrimAtCB", "Primaries at Cavern Border", fgcVarListForPrim);
        fPrimAtCB->SetDirectory(fOutputForPrim);

        switch (AmoreDetectorConstruction::GetDetGeometryType()) {
            case eDetGeometry::kDetector_AMoRE_I: {
            } break;
            case eDetGeometry::kDetector_AMoRE200: {
                fEvtInfos->Branch("EdepOV", &fEvtInfo_EdepOV, "EdepOV[2]/D");
								fEvtInfos->Branch("InciAtOVC", &fEvtInfo_InciAtOVC, "InciAtOVC/I");
                fPrimAtOVC = new TNtupleD("PrimAtOVC", "Primaries at OVC", fgcVarListForPrim);
                fPrimAtOVC->SetDirectory(fOutputForPrim);
            } break;
            default:
                G4Exception(__PRETTY_FUNCTION__, "PRIM_NOTSUPPORT",
                            G4ExceptionSeverity::JustWarning,
                            "This geometry does not support recording primaries at borders.");
                break;
        }

        fROOTOutputFile->cd();
    }
}

void AmoreRootNtuple::CloseFile() {
    if (fRecordPrimary && fOutputForPrim != nullptr) {
        fOutputForPrim->Write();
        fOutputForPrim->Close();
        delete fOutputForPrim;
        fEvtInfos      = nullptr;
        fPrimAtCB      = nullptr;
				fPrimAtOVC     = nullptr;
        fOutputForPrim = nullptr;
    }
    CupRootNtuple::CloseFile();
}

void AmoreRootNtuple::ClearEvent() {
    CupRootNtuple::ClearEvent();
    fTIDListForPrimAtCB.clear();
    fEvtInfo_EdepOV[0]       = 0;
    fEvtInfo_EdepOV[1]       = 0;
    fEvtInfo_HittedCMONum = 0;
    fEvtInfo_InciAtCB     = 0;
		fEvtInfo_InciAtOVC    = 0;
}

int AmoreRootNtuple::CountHittedCMOs() {
    eDetGeometry DetectorType;
    DetectorType = AmoreDetectorConstruction::GetDetGeometryType();
    switch (DetectorType) {
        case eDetGeometry::kDetector_AMoRE200: {
            G4int NumberOfCMOs   = Ctgsd->GetNTotCell();
						cout << " DEBUG: The number of Total Cell = " << NumberOfCMOs << endl;
            G4int hittedCMOs     = 0;
            for (int i = 0; i < NumberOfCMOs; i++) {
                TCell *WaterCell = static_cast<TCell *>(Ctgsd->GetCell()->At(i));
                if (WaterCell->GetEdep() != 0) hittedCMOs++;
            }
						cout << "Hitted CMOs - " << hittedCMOs << endl;
            return hittedCMOs;
        } break;
        case eDetGeometry::kDetector_AMoREPilot: {
            TCell *nowCell;
            G4int hittedCMOs = 0;
            for (G4int i = 0; i < 6; i++) {
                nowCell = static_cast<TCell *>(Ctgsd->GetCell()->At(i));
                if (nowCell->GetEdep() > 0) hittedCMOs++;
            }
            return hittedCMOs;
        } break;
        case eDetGeometry::kDetector_AMoREPilotRUN5: {
            TCell *nowCell;
            G4int hittedCMOs = 0;
            for (G4int i = 0; i < 6; i++) {
                nowCell = static_cast<TCell *>(Ctgsd->GetCell()->At(i));
                if (nowCell->GetEdep() > 0) hittedCMOs++;
            }
            return hittedCMOs;
        } break;
        case eDetGeometry::kDetector_AMoRE_I: {
            G4int hittedCMOs = 0;
            if (fModuleArray == nullptr) return -1;
            for (G4int i = 0; i < fModuleArray->GetTotalNumOfDetectorModules(); i++) {
                DetectorModule_Amore &nowModule =
                    static_cast<DetectorModule_Amore &>((*fModuleArray)[i]);
                if (nowModule.GetCrystalEdep() > 0) hittedCMOs++;
            }
						cout << "Hitted CMOs - " << hittedCMOs << endl;
            return hittedCMOs;
        } break;
        default:
            return -1;
            break;
    }
}

bool AmoreRootNtuple::RecordCut() {
    using namespace std;
    eDetGeometry DetectorType;
    DetectorType = AmoreDetectorConstruction::GetDetGeometryType();
    switch (DetectorType) {
        case eDetGeometry::kDetector_AMoRE200: {
            G4int hittedCMOs = CountHittedCMOs();

            if (hittedCMOs >= 1)
                return true;
            else
                return false;
						/*
            if (hittedCMOs != 1)
                return false;
            else
                return true;
								*/
        } break;
        case eDetGeometry::kDetector_AMoREPilot:
        case eDetGeometry::kDetector_AMoREPilotRUN5:
        case eDetGeometry::kDetector_AMoRE_I: {
            G4int hittedCMOs = CountHittedCMOs();
            if (hittedCMOs >= 1)
                return true;
            else
                return false;
        } break;
        default:
            return true;
            break;
    }
}

void AmoreRootNtuple::SetMDSD(const G4Event *aEvent) {
    G4String colName                    = "MDSD/AmoreModuleSDColl";
    G4SDManager *sdMan                  = G4SDManager::GetSDMpointer();
    G4int mdsdID                        = sdMan->GetCollectionID(colName);
    AmoreModuleHitsCollection *nowHColl = nullptr;

    nowHColl = static_cast<AmoreModuleHitsCollection *>(aEvent->GetHCofThisEvent()->GetHC(mdsdID));

    if (nowHColl == nullptr)
        G4Exception(__PRETTY_FUNCTION__, "MDSD_NOTEXIST", G4ExceptionSeverity::FatalException,
                    "Module SD is not found.");

    size_t i;
    eDetGeometry DetectorType;
    DetectorType = AmoreDetectorConstruction::GetDetGeometryType();
    switch (DetectorType) {
        case eDetGeometry::kDetector_AMoRE_I: {
            for (i = 0; i < nowHColl->GetSize(); i++) {
                AmoreModuleHit *nowHit = (*nowHColl)[i];

                DetectorModule_Amore &nowModule = (*fModuleArray)[nowHit->GetModuleID()];

                nowModule.SetCrystalEdep(nowHit->GetCrystalEdep());
                nowModule.SetGeWaferEdep(nowHit->GetGeWaferEdep());
                nowModule.SetQuenchedCrystalEdep(nowHit->GetCrystalQEdep());
                nowModule.SetQuenchedGeWaferEdep(nowHit->GetGeWaferQEdep());

                nowModule.SetGoldFilmEdep(0, nowHit->GetCrystalGoldFilmEdep(0));
                nowModule.SetGoldFilmEdep(1, nowHit->GetGeWaferGoldFilmEdep(0));
                nowModule.SetGoldFilmEdep(2, nowHit->GetGeWaferGoldFilmEdep(1));
                nowModule.SetGoldFilmEdep(3, nowHit->GetGeWaferGoldFilmEdep(2));
            }
        }
        default:
            break;
    }
}

void AmoreRootNtuple::SetTGSD(const G4Event *a_event) {
    //////////////////////////
    // TG Sensitive Detector
    /////////////////////////
    G4String colName;
    G4SDManager *SDman           = G4SDManager::GetSDMpointer();
    TGSDHCID                     = SDman->GetCollectionID(colName = "TGSD/TGSDColl");
    G4HCofThisEvent *HCE         = a_event->GetHCofThisEvent();
    CupScintHitsCollection *ECHC = 0;

    if (TGSDHCID < 0) return;
    if (HCE) ECHC = (CupScintHitsCollection *)(HCE->GetHC(TGSDHCID));

    int iHit                 = 0;
    double totalE            = 0.;
    double totalEquenched    = 0.;
    int nTotCell             = 0;
    double tgTotEdep         = 0.;
    double tgTotEdepQuenched = 0.;
    G4String logVolName;
    // int		idxDetID = 0;
    // char		detName[300];
    int idxTag[1000] = {};

    int nCell = 0;
    TCell tcell;

    nTotCell = ECHC->entries();
    Ctgsd->SetNTotCell(nTotCell);
    double eDep;
    double eDepQuenched;

    eDetGeometry DetectorType;
    DetectorType = AmoreDetectorConstruction::GetDetGeometryType();
    G4cout << "Ntuple DetectorType= " << DetectorType << G4endl;
    switch (DetectorType) {
        case eDetGeometry::kDetector_AMoRE200: //// for AMoRE200
            for (int i1 = 0; i1 < nTotCell; i1++) {
                if (idxTag[i1] == -1) continue;
                tgcellEdepQuenched[i1] = 0.;
                CupScintHit *aHit      = (*ECHC)[i1];
                eDep                   = aHit->GetEdep();
                eDepQuenched           = aHit->GetEdepQuenched();
                tcell.SetEdep(eDep);
                tcell.SetEdepQuenched(eDepQuenched);
                if (eDep > 0.) {
                    iHit++;
                    totalE += eDep;
                    totalEquenched += eDepQuenched;

                    G4cout << "i1= " << i1 << ", volumeName= " << aHit->GetLogV()->GetName()
                           << ", tgcellEdep[idxDetID]= " << eDep
                           << ", tgcellEdepQuenched[idxDetID]= " << eDepQuenched << G4endl;
                }
                tgTotEdep         = totalE;
                tgTotEdepQuenched = totalEquenched;

                tcell.SetCellID(i1);
                new ((*tclcell)[nCell]) TCell(tcell);
                nCell++;
            }
            break;
        case eDetGeometry::kDetector_AMoRE10: //// for AMoRE10
            for (int i1 = 0; i1 < nTotCell; i1++) {
                if (idxTag[i1] == -1) continue;
                tgcellEdepQuenched[i1] = 0.;
                CupScintHit *aHit      = (*ECHC)[i1];
                eDep                   = aHit->GetEdep();
                eDepQuenched           = aHit->GetEdepQuenched();
                tcell.SetEdep(eDep);
                tcell.SetEdepQuenched(eDepQuenched);
                if (eDep > 0.) {
                    iHit++;
                    totalE += eDep;
                    totalEquenched += eDepQuenched;

                    G4cout << "i1= " << i1 << ", tgcellEdep[idxDetID]= " << eDep
                           << ", tgcellEdepQuenched[idxDetID]= " << eDepQuenched << G4endl;
                }
                tgTotEdep         = totalE;
                tgTotEdepQuenched = totalEquenched;

                tcell.SetCellID(i1);
                new ((*tclcell)[nCell]) TCell(tcell);
                nCell++;
            }
            break;
        case eDetGeometry::kDetector_AMoREPilot: // for AMoREPilot
            for (int i1 = 0; i1 < nTotCell; i1++) {
                if (idxTag[i1] == -1) continue;
                tgcellEdepQuenched[i1] = 0.;
                CupScintHit *aHit      = (*ECHC)[i1];
                eDep                   = aHit->GetEdep();
                eDepQuenched           = aHit->GetEdepQuenched();
                tcell.SetEdep(eDep);
                tcell.SetEdepQuenched(eDepQuenched);
                if (eDep > 0.) {
                    iHit++;
                    totalE += eDep;
                    totalEquenched += eDepQuenched;

                    G4cout << "i1= " << i1 << ", tgcellEdep[idxDetID]= " << eDep
                           << ", tgcellEdepQuenched[idxDetID]= " << eDepQuenched << G4endl;
                }
                tgTotEdep         = totalE;
                tgTotEdepQuenched = totalEquenched;

                tcell.SetCellID(i1);
                new ((*tclcell)[nCell]) TCell(tcell);
                nCell++;
            }
            break;
        case eDetGeometry::kDetector_AMoREPilotRUN5: // for AMoREPilot
            for (int i1 = 0; i1 < nTotCell; i1++) {
                if (idxTag[i1] == -1) continue;
                tgcellEdepQuenched[i1] = 0.;
                CupScintHit *aHit      = (*ECHC)[i1];
                eDep                   = aHit->GetEdep();
                eDepQuenched           = aHit->GetEdepQuenched();
                tcell.SetEdep(eDep);
                tcell.SetEdepQuenched(eDepQuenched);
                if (eDep > 0.) {
                    iHit++;
                    totalE += eDep;
                    totalEquenched += eDepQuenched;

                    G4cout << "i1= " << i1 << ", tgcellEdep[idxDetID]= " << eDep
                           << ", tgcellEdepQuenched[idxDetID]= " << eDepQuenched << G4endl;
                }
                tgTotEdep         = totalE;
                tgTotEdepQuenched = totalEquenched;

                tcell.SetCellID(i1);
                new ((*tclcell)[nCell]) TCell(tcell);
                nCell++;
            }
            break;
        case eDetGeometry::kDetector_MyDetector: 
            for (int i1 = 0; i1 < nTotCell; i1++) {
                if (idxTag[i1] == -1) continue;
                tgcellEdepQuenched[i1] = 0.;
                CupScintHit *aHit      = (*ECHC)[i1];
                eDep                   = aHit->GetEdep();
                eDepQuenched           = aHit->GetEdepQuenched();
                tcell.SetEdep(eDep);
                tcell.SetEdepQuenched(eDepQuenched);
                if (eDep > 0.) {
                    iHit++;
                    totalE += eDep;
                    totalEquenched += eDepQuenched;

                    G4cout << "i1= " << i1 << ", volumeName= " << aHit->GetLogV()->GetName()
                           << ", tgcellEdep[idxDetID]= " << eDep
                           << ", tgcellEdepQuenched[idxDetID]= " << eDepQuenched << G4endl;
                }
                tgTotEdep         = totalE;
                tgTotEdepQuenched = totalEquenched;

                tcell.SetCellID(i1);
                new ((*tclcell)[nCell]) TCell(tcell);
                nCell++;
            }
            break;
        default:
            G4cout << "### Detector type is not valid!!! Signal is not filled to Cell!!!\n";
            break;
    }
    (void)tgTotEdep;
    (void)tgTotEdepQuenched;
    Ctgsd->SetTotEdep(totalE);
    Ctgsd->SetTotEdepQuenched(totalEquenched);
    Ctgsd->SetNTotCell(nTotCell);
		G4cout << " nHit: " << iHit  << ", TotEdep: " << totalE << endl;
}

void AmoreRootNtuple::RecordET(const G4Track *a_track) {
    Int_t trid          = a_track->GetTrackID();
    Int_t prntid        = a_track->GetParentID();
    G4String volume     = a_track->GetVolume()->GetName();
    G4String pname      = a_track->GetDefinition()->GetParticleName();
    G4int atomicmass    = a_track->GetDefinition()->GetAtomicMass();
    G4int atomicnumber  = a_track->GetDefinition()->GetAtomicNumber();
    Double_t ke         = a_track->GetStep()->GetPreStepPoint()->GetKineticEnergy();
    Double_t globaltime = a_track->GetGlobalTime();
    Double_t localtime  = a_track->GetLocalTime();

    G4ThreeVector pos = a_track->GetPosition();
    Float_t xx        = (float)pos.x();
    Float_t yy        = (float)pos.y();
    Float_t zz        = (float)pos.z();

    G4String procname;
    const G4VProcess *EndProcess = a_track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep();
    if (EndProcess) procname = EndProcess->GetProcessName();

    TTrack ttr;
    Int_t cuppdgcode;

    char partname[100];
    char volname[100];
    char prcsname[100];

    sprintf(partname, "%s", (char *)pname.data());
    sprintf(volname, "%s", (char *)volume.data());
    sprintf(prcsname, "%s", (char *)procname.data());

    G4ParticleDefinition *pdef = a_track->GetDefinition();
    if (G4IonTable::IsIon(pdef)) {
        cuppdgcode = (1000 * atomicnumber + atomicmass) * 10;
    } else {
        cuppdgcode = pdef->GetPDGEncoding();
    }

    if (StatusTrack) {
        ttr.SetParticleName(partname);
        ttr.SetAtomicNumber(atomicnumber);
        ttr.SetAtomicMass(atomicmass);
        ttr.SetPDGcode(cuppdgcode);
        ttr.SetTrackID(trid);
        ttr.SetParentID(prntid);
        ttr.SetKineticEnergy(ke);
        ttr.SetX(xx);
        ttr.SetY(yy);
        ttr.SetZ(zz);
        ttr.SetGlobalTime(globaltime);
        ttr.SetLocalTime(localtime);
        ttr.SetProcessName(prcsname);
        ttr.SetVolumeName(volname);

        EndTrackList->push_back(new TTrack(ttr));
    }
}

void AmoreRootNtuple::RecordTrack(const G4Track *a_track) {
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // EJ: G4TrackStatus
    //     fAlive=0: Continue the tracking
    //     fStopButAlive=1: Invoke active rest physics processes and kill the current track
    //     afterward fStopAndKill=2: Kill the current track fKillTrackAndSecondaries=3: Kill the
    //     current track and also associated secondaries fSuspend=4: Suspend the current track
    //     fPostponeToNextEvent=5: Postpones the tracking of thecurrent track to the next event
    /////////////////////////////////////////////////////////////////////////////////////////////////

    Int_t motherCopyNo;
    Int_t trid          = a_track->GetTrackID();
    Int_t prntid        = a_track->GetParentID();
    G4String volume     = a_track->GetVolume()->GetName();
    G4String pname      = a_track->GetDefinition()->GetParticleName();
    G4int atomicmass    = a_track->GetDefinition()->GetAtomicMass();
    G4int atomicnumber  = a_track->GetDefinition()->GetAtomicNumber();
    Double_t ke         = a_track->GetKineticEnergy();
    Double_t globaltime = a_track->GetGlobalTime();
    Double_t localtime  = a_track->GetLocalTime();

    G4ThreeVector pos = a_track->GetPosition();
    Float_t xx        = (float)pos.x();
    Float_t yy        = (float)pos.y();
    Float_t zz        = (float)pos.z();

    G4String procname;
    const G4VProcess *creatorProcess = a_track->GetCreatorProcess();
    if (creatorProcess) procname = creatorProcess->GetProcessName();

    TTrack ttr;
    Int_t cuppdgcode;

    char partname[100];
    char volname[100];
    char prcsname[100];

    sprintf(partname, "%s", (char *)pname.data());
    sprintf(volname, "%s", (char *)volume.data());
    sprintf(prcsname, "%s", (char *)procname.data());

    if (!strcmp(volname, "physCMOCell")) {
        motherCopyNo = a_track->GetTouchableHandle()->GetCopyNumber(1);
        sprintf(volname, "%s%d", "physCMOCell", motherCopyNo);
    }

    G4ParticleDefinition *pdef = a_track->GetDefinition();
    if (G4IonTable::IsIon(pdef)) {
        cuppdgcode = (1000 * atomicnumber + atomicmass) * 10;
    } else {
        cuppdgcode = pdef->GetPDGEncoding();
    }

    if (StatusTrack) {
        ttr.SetParticleName(partname);
        ttr.SetAtomicNumber(atomicnumber);
        ttr.SetAtomicMass(atomicmass);
        ttr.SetPDGcode(cuppdgcode);
        ttr.SetTrackID(trid);
        ttr.SetParentID(prntid);
        ttr.SetKineticEnergy(ke);
        ttr.SetX(xx);
        ttr.SetY(yy);
        ttr.SetZ(zz);
        ttr.SetGlobalTime(globaltime);
        ttr.SetLocalTime(localtime);
        ttr.SetProcessName(prcsname);
        ttr.SetVolumeName(volname);

        new ((*tcltr)[nTrack]) TTrack(ttr);
        nTrack++;
    }
}

void AmoreRootNtuple::RecordStep(const G4Step *a_step) {
    Int_t trid         = a_step->GetTrack()->GetTrackID();
    Int_t prntid       = a_step->GetTrack()->GetParentID();
    G4String volume    = a_step->GetTrack()->GetVolume()->GetName();
    G4int istep        = a_step->GetTrack()->GetCurrentStepNumber();
    G4String pname     = a_step->GetTrack()->GetDefinition()->GetParticleName();
    G4int atomicmass   = a_step->GetTrack()->GetDefinition()->GetAtomicMass();
    G4int atomicnumber = a_step->GetTrack()->GetDefinition()->GetAtomicNumber();
    Double_t dep       = a_step->GetTotalEnergyDeposit();

    G4StepPoint *preStep  = a_step->GetPreStepPoint();
    G4StepPoint *postStep = a_step->GetPostStepPoint();

    G4String procname   = postStep->GetProcessDefinedStep()->GetProcessName();
    Double_t ke         = postStep->GetKineticEnergy();
    Double_t globaltime = postStep->GetGlobalTime();
    Double_t localtime  = postStep->GetLocalTime();
    G4ThreeVector pos   = postStep->GetPosition();
    Float_t xx          = (float)pos.x();
    Float_t yy          = (float)pos.y();
    Float_t zz          = (float)pos.z();

    G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
    G4int motherCopyNo             = -1;

    //	TTrack ttr;
    TStep tst;
    Int_t cuppdgcode;

    char volname[100];

    bool isWorldRegion;

    isWorldRegion =
        theTouchable->GetVolume(theTouchable->GetHistoryDepth())->GetLogicalVolume() ==
        theTouchable->GetVolume()->GetLogicalVolume()->GetRegion()->GetRootLogicalVolumeIterator().
        operator*();

    sprintf(volname, "%s", (char *)volume.data());
    if (!isWorldRegion) {
        switch (AmoreDetectorConstruction::GetDetGeometryType()) {
            case eDetGeometry::kDetector_AMoRE200:
                if (!strcmp(volname, "physCMOCell")) {
                    motherCopyNo = theTouchable->GetCopyNumber(1);
                    sprintf(volname, "%s%d", "physCMOCell", motherCopyNo);
                }
                break;
            case eDetGeometry::kDetector_AMoRE_I: {
                constexpr const char *pvCheckerString = "_Crystal_PV";
                //size_t checkerIndex                 = volume.index(pvCheckerString);
                size_t checkerIndex                 = volume.find(pvCheckerString);
                if (checkerIndex != G4String::npos) {
                    G4VPhysicalVolume *theMother;
                    size_t envelopeCopyNo = 1;
										size_t nowMotherLevel = 1;
                    theMother = theTouchable->GetVolume(nowMotherLevel++);
                    while (theMother != nullptr) {
                        if (theMother->GetLogicalVolume()->IsRootRegion()) {
                            if (theMother->GetMotherLogical() ==
                                nullptr) { // Reached the end of the world
                                G4Exception(__PRETTY_FUNCTION__, "MDSD_REGION_FAIL", FatalException,
                                            "Finding a root region for SD has been failed.");
                            }
                            envelopeCopyNo = theMother->GetCopyNo();
                            break;
                        }
                        theMother = theTouchable->GetVolume(nowMotherLevel++);
                    }

                    volume.insert(checkerIndex, "_");
                    volume.insert(checkerIndex + 1, std::to_string(envelopeCopyNo));
                    sprintf(volname, "%s", (char *)volume.data());
                }
            } break;
            default:
                break;
        }
    }

    G4ParticleDefinition *pdef = a_step->GetTrack()->GetDefinition();
    if (G4IonTable::IsIon(pdef)) {
        cuppdgcode = (1000 * atomicnumber + atomicmass) * 10;
    } else {
        cuppdgcode = pdef->GetPDGEncoding();
    }

    tst.SetParticleName(pname.data());
    tst.SetPDGcode(cuppdgcode);
    tst.SetTrackID(trid);
    tst.SetParentID(prntid);
    tst.SetKineticEnergy(ke);
    tst.SetEnergyDeposit(dep);
    tst.SetX(xx);
    tst.SetY(yy);
    tst.SetZ(zz);
    tst.SetGlobalTime(globaltime);
    tst.SetLocalTime(localtime);
    tst.SetProcessName(procname.data());
    tst.SetVolumeName(volname);
    tst.SetStepNo(istep);

    if (StatusStep) {
        new ((*tclst)[nStep]) TStep(tst);
        nStep++;
    }

    if (fRecordPrimary) RecordPrimaryAtBorder(a_step);
}

void AmoreRootNtuple::RecordPrimaryAtBorder(const G4Step *aStep) {
    G4StepPoint *postStep                   = aStep->GetPostStepPoint();

    const AmoreDetectorConstruction *tDetCons = static_cast<const AmoreDetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    eCavernType tNowCT    = AmoreDetectorConstruction::GetCavernType();
    //eVetoGeometry tNowVGT = AmoreDetectorConstruction::GetVetoGeometryType();

    G4bool fileFlushed = false;
    auto flushTuple    = [&](TNtupleD *aTuple, G4int &aFillCnt) {
        if (aTuple != nullptr && kEvtModForPrim != 0 && aFillCnt > kEvtModForPrim) {
            aTuple->FlushBaskets();
            aFillCnt = 0;
            if (!fileFlushed) {
                aTuple->GetCurrentFile()->Flush();
                fileFlushed = true;
            }
        }
    };

    auto judgeTID = [&](std::map<G4int, G4int> &aList, G4int aTID) {
        auto result = aList.find(aTID);
        if (result == aList.end()) {
            aList[aTID] = fValuesForPrim[9] = 0;
        } else {
            fValuesForPrim[9] = ++(result->second);
        }
    };

    fValuesForPrim[0] = fEvtInfo_EvtID;
    fValuesForPrim[1] = postStep->GetPosition().x();
    fValuesForPrim[2] = postStep->GetPosition().y();
    fValuesForPrim[3] = postStep->GetPosition().z();
    fValuesForPrim[4] = postStep->GetKineticEnergy();
    fValuesForPrim[5] = postStep->GetMomentum().x();
    fValuesForPrim[6] = postStep->GetMomentum().y();
    fValuesForPrim[7] = postStep->GetMomentum().z();
    fValuesForPrim[8] = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    fValuesForPrim[9] = -1;

    if (aStep->GetTrack()->GetParentID() != 0) {
        AmoreTrackInformation *aATI =
            static_cast<AmoreTrackInformation *>(aStep->GetTrack()->GetUserInformation());
        fValuesForPrim[10] = aATI->GetParentDefinition()->GetPDGEncoding();
        auto res_iter      = fEvtInfo_VolumeTbl->find(aATI->GetBirthPV()->GetName());
        if (res_iter == fEvtInfo_VolumeTbl->end()) {
            fValuesForPrim[11] = -100;
            G4Exception(__PRETTY_FUNCTION__, "PRIM_VOLTBL_NOEXIST",
                        G4ExceptionSeverity::JustWarning, "We couldn't find a PV in the PV table.");
        } else
            fValuesForPrim[11] = (*res_iter).second;
    } else {
        fValuesForPrim[10] = 0;
        fValuesForPrim[11] = -1;
    }

    switch (AmoreDetectorConstruction::GetDetGeometryType()) {
        case eDetGeometry::kDetector_AMoRE_I: {
            if (tDetCons->Judge_CavernBorder(aStep)) {
                judgeTID(fTIDListForPrimAtCB, aStep->GetTrack()->GetTrackID());
                fPrimAtCB->Fill(fValuesForPrim);
                fPrimFillCntAtCB++;
                fEvtInfo_InciAtCB++;
            }
            flushTuple(fPrimAtCB, fPrimFillCntAtCB);
        } break;
        case eDetGeometry::kDetector_AMoRE200: {
            switch (tNowCT) {
                case eCavernType::kCavern_Toy_HemiSphere:
                    if (tDetCons->Judge_CavernBorder(aStep)) {
                        judgeTID(fTIDListForPrimAtCB, aStep->GetTrack()->GetTrackID());
                        fPrimAtCB->Fill(fValuesForPrim);
                        fPrimFillCntAtCB++;
                        fEvtInfo_InciAtCB++;
                    } 
                    break;
                case eCavernType::kCavern_Toy_Cylinder:
                    if (AmoreDetectorConstruction::GetDetGeometryType() ==
                        eDetGeometry::kDetector_AMoRE200)
                        G4Exception(__PRETTY_FUNCTION__, "CAVERN", G4ExceptionSeverity::JustWarning,
                                    "Cavern type of cylinder toy model hasn't been implemented in "
                                    "AMoRE-II.");
                    return;
                    break;
                case eCavernType::kCavern_RealModel:
                    if (tDetCons->Judge_CavernBorder(aStep)) {
                        judgeTID(fTIDListForPrimAtCB, aStep->GetTrack()->GetTrackID());
                        fPrimAtCB->Fill(fValuesForPrim);
                        fPrimFillCntAtCB++;
                        fEvtInfo_InciAtCB++;
                    }
                    break;
                default:
                    G4Exception(__PRETTY_FUNCTION__, "CAVERN", G4ExceptionSeverity::JustWarning,
                                "Cavern type is wrong.");
                    return;
            }

						if(tDetCons->Judge_200_OVCBorder(aStep))
						{
                judgeTID(fTIDListForPrimAtOVC, aStep->GetTrack()->GetTrackID());
                fPrimAtOVC->Fill(fValuesForPrim);
                fPrimFillCntAtOVC++;
                fEvtInfo_InciAtOVC++;
						}

            flushTuple(fPrimAtCB, fPrimFillCntAtCB);
            flushTuple(fPrimAtOVC, fPrimFillCntAtOVC);
						
        } break;
        default:
            break;
    }
}

void AmoreRootNtuple::RecordBeginOfEvent(const G4Event *a_event) {
    CupRootNtuple::RecordBeginOfEvent(a_event);
    fEvtInfo_EvtID = a_event->GetEventID();

    // Genarate Volume Table
    if (fEvtInfo_VolumeTbl->size() == 0) {
        G4PhysicalVolumeStore *aPVStore = G4PhysicalVolumeStore::GetInstance();

        int idx = 0;
        for (auto &nowPV : *aPVStore) {
            auto &nowName = nowPV->GetName();
            if (fEvtInfo_VolumeTbl->find(nowName.c_str()) == fEvtInfo_VolumeTbl->end()) {
                (*fEvtInfo_VolumeTbl)[nowName.c_str()] = idx++;
            }
        }
    }
}

void AmoreRootNtuple::RecordPrimaryEvtInfos(const G4Event *) {
    switch (AmoreDetectorConstruction::GetDetGeometryType()) {
        case eDetGeometry::kDetector_AMoRE_I: {
            fEvtInfo_HittedCMONum = CountHittedCMOs();
            fEvtInfos->Fill();
        } break;
        case eDetGeometry::kDetector_AMoRE200: {
            Double_t edepOV[2]       = {-1,-1};
            const Int_t vetoIdx1   = Ctgsd->GetCell()->GetEntries() - 2;
            const TCell *vetoCell1 = static_cast<TCell *>((Ctgsd->GetCell())->At(vetoIdx1));
            const Int_t vetoIdx2   = Ctgsd->GetCell()->GetEntries() - 1;
            const TCell *vetoCell2 = static_cast<TCell *>((Ctgsd->GetCell())->At(vetoIdx2));
            edepOV[0]                = vetoCell1->GetEdep();
            edepOV[1]                = vetoCell2->GetEdep();
            fEvtInfo_EdepOV[0]       = edepOV[0];
            fEvtInfo_EdepOV[1]       = edepOV[1];
            fEvtInfo_HittedCMONum = CountHittedCMOs();
            fEvtInfos->Fill();
        } break;
        default:
            break;
    }
}

void AmoreRootNtuple::RecordEndOfEvent(const G4Event *a_event) {
    SetEventInfo(a_event);
    if (StatusPrimary) {
        SetPrimary(a_event);
    }
    if (StatusPhoton) {
        SetPhoton();
    }
    if (StatusScint) {
        SetScintillation();
    }

    if (fModuleArray != nullptr) {
        SetMDSD(a_event);
    } else {
        SetTGSD(a_event);
    }

    if (StatusTrack) {
        Ctrack->SetNtrack(nTrack);
    }
    if (StatusStep) {
        Cstep->SetNstep(nStep);
    }

    if (StatusMuon) {
        SetMuonSD(a_event);
    }

    if (fRecordPrimary) {
        RecordPrimaryEvtInfos(a_event);
    }

    // put to tree
    if (fRecordWithCut) {
        if (RecordCut()) {
            fROOTOutputTree->Fill();
            fRecordedEvt++;
        }
    } else {
        fROOTOutputTree->Fill();
        fRecordedEvt++;
    }

    if (kEvtMod != 0 && fRecordedEvt > kEvtMod) {
        fROOTOutputTree->FlushBaskets();
        fROOTOutputTree->GetCurrentFile()->Flush();
        fRecordedEvt = 0;
    }

    ClearEvent();

    G4cout << "///////////////////////////////// End of Event "
              "/////////////////////////////////////////"
           << G4endl;
}
