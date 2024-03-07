// AmoreRootNtuple.hh

#ifndef AmoreRootNtuple_h
#define AmoreRootNtuple_h 1

// The AmoreRecorderBase object.  We're implementing this abstract class
// with methods that fill histograms.
#include "CLHEP/Vector/ThreeVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "AmoreSim/AmoreDetectorConstruction.hh"
#include "AmoreSim/AmoreRootNtupleMessenger.hh"
#include "AmoreSim/AmoreTrajectoryPoint.hh"
#include "CupSim/CupRootNtuple.hh"

#include "MCObjs/DetectorArray_Amore.hh"
#include "MCObjs/DetectorModule_Amore.hh"
#include "MCObjs/EvtInfo.hh"
#include "MCObjs/EvtStep.hh"
#include "MCObjs/EvtTrack.hh"
#include "MCObjs/MLCS.hh"
#include "MCObjs/MuonSD.hh"
#include "MCObjs/PMTSD.hh"
#include "MCObjs/Photon.hh"
#include "MCObjs/Primary.hh"
#include "MCObjs/Scint.hh"
#include "MCObjs/TGSD.hh"
#include "MCObjs/Vertex.hh"

#include <map>
#include <vector>
// Forward declarations for ROOT.
// class TFile;
// class TNtuple;

// Forward declarations of G4 classes that are arguments to
// AmoreRecorderBase methods.
class G4Event;
class G4Step;
class AmoreRootNtupleMessenger;

// class CupRootNtuple: public CupRecorderBase
class AmoreRootNtuple : public CupRootNtuple {
  private:
    G4int fRecordedEvt;
    G4bool fRecordWithCut;
    G4bool fRecordPrimary;

    AmoreRootNtupleMessenger *myAmoreNtupleMessenger;

    TFile *fOutputForPrim;
    TTree *fEvtInfos;

    TNtupleD *fPrimAtCB;
		TNtupleD *fPrimAtOVC;
    G4int fPrimFillCntAtCB;
    G4int fPrimFillCntAtOVC;

    DetectorArray_Amore *fModuleArray;

  protected:
    using eDetGeometry  = AmoreDetectorConstruction::eDetGeometry;
    using eCavernType   = AmoreDetectorConstruction::eCavernType;
		using ePhaseAMoRE200 = AmoreDetectorConstruction::ePhaseAMoRE200;
    std::vector<TTrack *> *EndTrackList;

    G4int fEvtInfo_EvtID;
    G4double fEvtInfo_EdepOV[2];
    G4int fEvtInfo_HittedCMONum;
    G4int fEvtInfo_InciAtCB;
    G4int fEvtInfo_InciAtOVC;
    std::map<std::string, int> *fEvtInfo_VolumeTbl;

    static constexpr G4int fgcNVarForPrim = 12;
    static constexpr const char *fgcVarListForPrim =
        "EvtID:X:Y:Z:KE:Px:Py:Pz:PDG:InciNum:M_PDG:BirthPVIdx";
    G4double fValuesForPrim[fgcNVarForPrim];
    std::map<G4int, G4int> fTIDListForPrimAtCB;
    std::map<G4int, G4int> fTIDListForPrimAtOVC;

    void ClearEvent();

  public:
    AmoreRootNtuple();
    ~AmoreRootNtuple();

    virtual void RecordBeginOfEvent(const G4Event *);
    virtual void RecordEndOfEvent(const G4Event *);
    virtual void SetTGSD(const G4Event *a_event);
    // virtual void SetMuonSD(const G4Event *a_event);
    virtual void SetMDSD(const G4Event *a_event);
    virtual void OpenFile(const G4String filename, G4bool outputMode);
    virtual void CloseFile();

    virtual void CreateTree();

    virtual bool RecordCut();
    virtual int CountHittedCMOs();
    virtual void RecordStep(const G4Step *);
    virtual void RecordTrack(const G4Track *);
    virtual void RecordET(const G4Track *);
    virtual void RecordPrimaryAtBorder(const G4Step *aStep);
    virtual void RecordPrimaryEvtInfos(const G4Event *aEvent);
    void ClearET() {
        for (auto i : *EndTrackList)
            i->Delete();
        EndTrackList->clear();
    };

    inline void SetRecordCut(G4bool a) { fRecordWithCut = a; }
    inline G4bool GetRecordCut() { return fRecordWithCut; }

    inline void SetRecordPrim(G4bool a) {
        if (fROOTOutputFile != nullptr) {
            G4Exception(__PRETTY_FUNCTION__, "PRIM_OPEN_ERR", JustWarning,
                        " The file already has been opened!  Primary record setting will not be "
                        "changed...  ");
        } else
            fRecordPrimary = a;
    }
    inline G4bool GetRecordPrim() { return fRecordPrimary; }

    enum {
        max_primary_particles   = 16,
        max_hits_for_ROOT       = 200000,
        max_gammas              = 10000,
        max_opticalphotons      = 10000,
        max_protons             = 1000,
        max_nTrk                = 1000000,
        max_secondary_particles = 1000
    };
    enum { nLayer_vessel = 4, inout_muon = 2 };
    enum {
        max_muonSecondary           = 10000,
        max_neutronCapture          = 100,
        max_neutronCaptureSecondary = 100,
        max_correlBKG               = 1000000
    };
    enum { kEvtMod = 1000, kEvtModForPrim = 100000 };
};
#endif
