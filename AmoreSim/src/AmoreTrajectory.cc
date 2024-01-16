//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file runAndEvent/Amore/src/AmoreTrajectory.cc
/// \brief Implementation of the AmoreTrajectory class
//
// $Id: $
//
#include "AmoreSim/AmoreTrajectory.hh"
#include "AmoreSim/AmoreTrajectoryPoint.hh"

#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4ParticleTable.hh"
#include "G4SteppingManager.hh"
#include "G4TrackingManager.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4Version.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

#if G4VERSION_NUMBER <= 999
G4Allocator<AmoreTrajectory> *faTrajAllocator = nullptr;
#else
G4ThreadLocal G4Allocator<AmoreTrajectory> *faTrajAllocator = nullptr;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectory::AmoreTrajectory()
    : fPositionRecord(0), fTrackID(0), fParentID(0), fPDGEncoding(0), fPDGCharge(0.0),
      fParticleName(""), fInitialKineticEnergy(0.), fInitialMomentum(G4ThreeVector()),
      prevSecNum(0) {
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectory::AmoreTrajectory(const G4Track *aTrack, G4TrackingManager *aMan) {
    G4ParticleDefinition *fpParticleDefinition = aTrack->GetDefinition();
    fParticleName                              = fpParticleDefinition->GetParticleName();
    fPDGCharge                                 = fpParticleDefinition->GetPDGCharge();
    fPDGEncoding                               = fpParticleDefinition->GetPDGEncoding();
    fTrackID                                   = aTrack->GetTrackID();
    fParentID                                  = aTrack->GetParentID();
    fInitialKineticEnergy                      = aTrack->GetKineticEnergy();
    fInitialMomentum                           = aTrack->GetMomentum();
    fPositionRecord                            = new TrajectoryPointContainer();
    prevSecNum                                 = 0;
    fpMan                                      = aMan;
    // Following is for the first trajectory point
    fPositionRecord->push_back(new AmoreTrajectoryPoint(aTrack));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectory::AmoreTrajectory(AmoreTrajectory &right) : G4VTrajectory() {
    fParticleName         = right.fParticleName;
    fPDGCharge            = right.fPDGCharge;
    fPDGEncoding          = right.fPDGEncoding;
    fTrackID              = right.fTrackID;
    fParentID             = right.fParentID;
    fInitialKineticEnergy = right.fInitialKineticEnergy;
    fInitialMomentum      = right.fInitialMomentum;
    fPositionRecord       = new TrajectoryPointContainer();
    fpMan                 = right.fpMan;

    for (size_t i = 0; i < right.fPositionRecord->size(); i++) {
        AmoreTrajectoryPoint *rightPoint = (AmoreTrajectoryPoint *)((*(right.fPositionRecord))[i]);
        fPositionRecord->push_back(new AmoreTrajectoryPoint(*rightPoint));
    }

    prevSecNum = right.prevSecNum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectory::~AmoreTrajectory() {
    if (fPositionRecord) {
        //  fPositionRecord->clearAndDestroy();
        size_t i;
        for (i = 0; i < fPositionRecord->size(); i++) {
            delete (*fPositionRecord)[i];
        }
        fPositionRecord->clear();
        delete fPositionRecord;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreTrajectory::ShowTrajectory(std::ostream &os) const {
    // Invoke the default implementation in G4VTrajectory...
    G4VTrajectory::ShowTrajectory(os);
    // ... or override with your own code here.
}

/***
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreTrajectory::DrawTrajectory() const
{
// Invoke the default implementation in G4VTrajectory...
G4VTrajectory::DrawTrajectory();
// ... or override with your own code here.
}
 ***/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreTrajectory::DrawTrajectory(G4int i_mode) const {
    // Invoke the default implementation in G4VTrajectory...
#if G4VERSION_NUMBER <= 999
    G4VTrajectory::DrawTrajectory(i_mode);
#else
    G4VTrajectory::DrawTrajectory();
#endif
    // ... or override with your own code here.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::map<G4String, G4AttDef> *AmoreTrajectory::GetAttDefs() const {
    G4bool isNew;
    std::map<G4String, G4AttDef> *store = G4AttDefStore::GetInstance("AmoreTrajectory", isNew);
    if (isNew) {

        G4String id("ID");
        (*store)[id] = G4AttDef(id, "Track ID", "Physics", "", "G4int");

        G4String pid("PID");
        (*store)[pid] = G4AttDef(pid, "Parent ID", "Physics", "", "G4int");

        G4String pn("PN");
        (*store)[pn] = G4AttDef(pn, "Particle Name", "Physics", "", "G4String");

        G4String ch("Ch");
        (*store)[ch] = G4AttDef(ch, "Charge", "Physics", "e+", "G4double");

        G4String pdg("PDG");
        (*store)[pdg] = G4AttDef(pdg, "PDG Encoding", "Physics", "", "G4int");

        G4String ike("IKE");
        (*store)[ike] =
            G4AttDef(ike, "Initial kinetic energy", "Physics", "G4BestUnit", "G4double");

        G4String iMom("IMom");
        (*store)[iMom] =
            G4AttDef(iMom, "Initial momentum", "Physics", "G4BestUnit", "G4ThreeVector");

        G4String iMag("IMag");
        (*store)[iMag] =
            G4AttDef(iMag, "Initial momentum magnitude", "Physics", "G4BestUnit", "G4double");

        G4String ntp("NTP");
        (*store)[ntp] = G4AttDef(ntp, "No. of points", "Physics", "", "G4int");
    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4AttValue> *AmoreTrajectory::CreateAttValues() const {
    std::vector<G4AttValue> *values = new std::vector<G4AttValue>;

    values->push_back(G4AttValue("ID", G4UIcommand::ConvertToString(fTrackID), ""));

    values->push_back(G4AttValue("PID", G4UIcommand::ConvertToString(fParentID), ""));

    values->push_back(G4AttValue("PN", fParticleName, ""));

    values->push_back(G4AttValue("Ch", G4UIcommand::ConvertToString(fPDGCharge), ""));

    values->push_back(G4AttValue("PDG", G4UIcommand::ConvertToString(fPDGEncoding), ""));

    values->push_back(G4AttValue("IKE", G4BestUnit(fInitialKineticEnergy, "Energy"), ""));

    values->push_back(G4AttValue("IMom", G4BestUnit(fInitialMomentum, "Energy"), ""));

    values->push_back(G4AttValue("IMag", G4BestUnit(fInitialMomentum.mag(), "Energy"), ""));

    values->push_back(G4AttValue("NTP", G4UIcommand::ConvertToString(GetPointEntries()), ""));

#ifdef G4ATTDEBUG
    G4cout << G4AttCheck(values, GetAttDefs());
#endif

    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreTrajectory::AppendStep(const G4Step *aStep) {
    // fPositionRecord->push_back( new AmoreTrajectoryPoint(
    // aStep->GetPostStepPoint()->GetPosition(),
    // aStep->GetPreStepPoint()->GetMaterial() ));
    AmoreTrajectoryPoint *np = new AmoreTrajectoryPoint(aStep);
    fPositionRecord->push_back(np);
    const G4TrackVector *trkvec = aStep->GetSecondary();
    G4int size4this             = trkvec->size() - prevSecNum;
    for (G4int i = prevSecNum; i < trkvec->size(); i++)
        np->AddSecondary((*trkvec)[i]);
    prevSecNum += size4this;

    G4SteppingManager *aStepMan = fpMan->GetSteppingManager();
    np->SetNRest(aStepMan->GetfN2ndariesAtRestDoIt());
    np->SetNAlong(aStepMan->GetfN2ndariesAlongStepDoIt());
    np->SetNPost(aStepMan->GetfN2ndariesPostStepDoIt());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ParticleDefinition *AmoreTrajectory::GetParticleDefinition() {
    return (G4ParticleTable::GetParticleTable()->FindParticle(fParticleName));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreTrajectory::MergeTrajectory(G4VTrajectory *secondTrajectory) {
    if (!secondTrajectory) return;

    AmoreTrajectory *seco = (AmoreTrajectory *)secondTrajectory;
    G4int ent             = seco->GetPointEntries();
    for (G4int i = 1; i < ent; i++) // initial point of the second trajectory
                                    // should not be merged
    {
        fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
        //    fPositionRecord->push_back(seco->fPositionRecord->removeAt(1));
    }
    delete (*seco->fPositionRecord)[0];
    seco->fPositionRecord->clear();
}
