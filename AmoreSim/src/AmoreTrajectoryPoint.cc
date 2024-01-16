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
/// \file runAndEvent/Amore/src/AmoreTrajectoryPoint.cc
/// \brief Implementation of the AmoreTrajectoryPoint class
//
// $Id: $
//
#include "AmoreSim/AmoreTrajectoryPoint.hh"

#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "G4Version.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

#if G4VERSION_NUMBER <= 999
G4Allocator<AmoreTrajectoryPoint>* faTrajPointAllocator = nullptr;
#else
G4ThreadLocal G4Allocator<AmoreTrajectoryPoint>* faTrajPointAllocator = nullptr;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectoryPoint::AmoreTrajectoryPoint(const G4Track *aTrack) {
    fPosition     = aTrack->GetPosition();
    fMomentum     = aTrack->GetMomentum();
    fPolarization = aTrack->GetPolarization();
    fPDGCode      = aTrack->GetDynamicParticle()->GetPDGcode();
    fTrackID      = aTrack->GetTrackID();
    fKinE         = aTrack->GetKineticEnergy();
    fTime         = aTrack->GetGlobalTime();
    fMass         = aTrack->GetDynamicParticle()->GetMass();
    fProcName     = G4String("initStep");
    fNextVolName  = aTrack->GetNextVolume()->GetName();
    fpSec         = new G4TrackVector();
    fPD           = aTrack->GetParticleDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectoryPoint::AmoreTrajectoryPoint(const G4Step *aStep) {
    fPDGCode      = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
    fPosition     = aStep->GetPostStepPoint()->GetPosition();
    fMomentum     = aStep->GetPostStepPoint()->GetMomentum();
    fPolarization = aStep->GetTrack()->GetPolarization();
    fTrackID      = aStep->GetTrack()->GetTrackID();
    fKinE         = aStep->GetTrack()->GetKineticEnergy();
    fTime         = aStep->GetTrack()->GetGlobalTime();
    fMass         = aStep->GetPostStepPoint()->GetMass();
    fpSec         = new G4TrackVector();
    fPD           = aStep->GetTrack()->GetDefinition();

    if (aStep->GetTrack()->GetNextVolume() != 0)
        fNextVolName = aStep->GetTrack()->GetNextVolume()->GetName();
    else
        fNextVolName = G4String("OutOfWorld");

    if (aStep->GetPostStepPoint()->GetProcessDefinedStep() != 0)
        fProcName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    else {
        if (fPDGCode == 0)
            fProcName = G4String("initStep");
        else
            fProcName = G4String("User Limit");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectoryPoint::AmoreTrajectoryPoint(const AmoreTrajectoryPoint &right)
    : G4VTrajectoryPoint(), fPosition(right.fPosition), fMomentum(right.fMomentum),
      fPolarization(right.fPolarization), fPDGCode(right.fPDGCode), fTrackID(right.fTrackID),
      fKinE(right.fKinE), fTime(right.fTime), fMass(right.fMass), fProcName(right.fProcName),
      fNextVolName(right.fNextVolName), fPD(right.fPD) {
    fpSec = new G4TrackVector();
    for (int i = 0; right.fpSec->size(); i++)
        fpSec->push_back((*right.fpSec)[i]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AmoreTrajectoryPoint::~AmoreTrajectoryPoint() {
    for (auto i : *fpSec)
        delete i;
    delete fpSec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::map<G4String, G4AttDef> *AmoreTrajectoryPoint::GetAttDefs() const {
    G4bool isNew;
    std::map<G4String, G4AttDef> *store = G4AttDefStore::GetInstance("AmoreTrajectoryPoint", isNew);
    if (isNew) {
        G4String pos("Pos");
        (*store)[pos] = G4AttDef(pos, "Position", "Physics", "G4BestUnit", "G4ThreeVector");
        G4String mom("Mom");
        (*store)[mom] = G4AttDef(mom, "Momentum", "Physics", "G4BestUnit", "G4ThreeVector");
        G4String sn("Step");
        (*store)[sn] = G4AttDef(sn, "PDGCode", "Physics", "Count", "G4int");
        G4String kinE("KinE");
        (*store)[kinE] = G4AttDef(kinE, "Kinetic Energy", "Physics", "G4BestUnit", "G4double");
        G4String GT("Time");
        (*store)[GT] = G4AttDef(GT, "Global Time", "Physics", "G4BestUnit", "G4double");
        G4String sl("Mass");
        (*store)[sl] = G4AttDef(sl, "Mass", "Physics", "G4BestUnit", "G4double");
        G4String tl("Pol");
        (*store)[tl] = G4AttDef(tl, "Particle Polarization", "Physics", "G4BestUnit", "G4double");
        G4String pn("ProcName");
        (*store)[pn] = G4AttDef(pn, "Process Name", "Physics", "Name", "G4String");
        G4String nvn("NVName");
        (*store)[nvn] = G4AttDef(nvn, "Next Volume Name", "Physics", "Name", "G4String");
    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4AttValue> *AmoreTrajectoryPoint::CreateAttValues() const {
    std::vector<G4AttValue> *values = new std::vector<G4AttValue>;

    values->push_back(G4AttValue("Pos", G4BestUnit(fPosition, "Length"), ""));
    values->push_back(G4AttValue("Mom", G4BestUnit(fPosition, "Energy"), ""));
    values->push_back(G4AttValue("Pol", G4BestUnit(fPolarization, "Length"), ""));
    values->push_back(G4AttValue("PDGCode", fPDGCode, ""));
    values->push_back(G4AttValue("KinE", G4BestUnit(fKinE, "Energy"), ""));
    values->push_back(G4AttValue("Time", G4BestUnit(fTime, "Time"), ""));
    values->push_back(G4AttValue("Mass", G4BestUnit(fMass, "Energy"), ""));
    values->push_back(G4AttValue("ProcName", fProcName, ""));
    values->push_back(G4AttValue("NVName", fNextVolName, ""));

#ifdef G4ATTDEBUG
    G4cout << G4AttCheck(values, GetAttDefs());
#endif

    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AmoreTrajectoryPoint::AddSecondary(G4Track *aTrack) { fpSec->push_back(new G4Track(*aTrack)); }
