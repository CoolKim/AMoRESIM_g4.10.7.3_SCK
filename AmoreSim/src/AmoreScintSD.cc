#include "AmoreSim/AmoreScintSD.hh"
#include "CupSim/CupScintSD.hh"

#include "AmoreSim/AmoreScintillation.hh"
//#include "CupSim/CupScintillation.hh"
//#include "G4LossTableManager.hh"
//#include "G4EmSaturation.hh"

// Constructor /////////////////////////////////////////////////////////////
AmoreScintSD::AmoreScintSD(G4String name, int arg_max_tgs) : CupScintSD(name, arg_max_tgs) {}

// Destructor //////////////////////////////////////////////////////////////
AmoreScintSD::~AmoreScintSD() {}

G4bool AmoreScintSD::ProcessHits(G4Step *aStep, G4TouchableHistory * /*ROhist*/) {
    //  G4EmSaturation * emSaturation = G4LossTableManager::Instance()->EmSaturation();

    //  G4double edep_quenched = emSaturation->VisibleEnergyDeposition(aStep);
    //  G4double edep_quenched = CupScintillation::GetTotEdepQuenched();
    G4double edep_quenched             = AmoreScintillation::GetTotEdepQuenched();
    G4double edep                      = aStep->GetTotalEnergyDeposit();
    G4ParticleDefinition *particleType = aStep->GetTrack()->GetDefinition();
    G4String particleName              = particleType->GetParticleName();

    //  if(edep==0.) return true;
    if (edep == 0. || particleName == "opticalphoton") return true;
    G4cout << "EJ: edep= " << edep << ", visible= " << edep_quenched << G4endl;

    // EJ: for LSVetoFullDetector (20150910)
    G4StepPoint *preStepPoint            = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable       = preStepPoint->GetTouchableHandle();
    G4int copyNo                         = theTouchable->GetCopyNumber();
    G4int motherCopyNo                   = theTouchable->GetCopyNumber(1);
    G4VPhysicalVolume *thePhysical       = theTouchable->GetVolume();
    G4VPhysicalVolume *theMotherPhysical = theTouchable->GetVolume(1);
    G4String motherVolName               = theMotherPhysical->GetName();
    // G4cout << "EJ: physical volume name= " << thePhysical->GetName() << ", mother physical volume
    // name= " << motherVolName << G4endl; G4cout << "EJ: copyNo= " << copyNo << ", motherCopyNo= "
    // << motherCopyNo << G4endl;

    if ((strstr(motherVolName, "Envelope")) != NULL) copyNo = motherCopyNo;
    // G4cout << "EJ2: copyNo= " << copyNo << ", motherCopyNo= " << motherCopyNo << G4endl;

    CupScintHit *aHit = (*hitsCollection)[copyNo];
    // check if it is first touch
    if (!(aHit->GetLogV())) {
        // fill volume information
        aHit->SetLogV(thePhysical->GetLogicalVolume());
        G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
        aTrans.Invert();
        aHit->SetRot(aTrans.NetRotation());
        aHit->SetPos(aTrans.NetTranslation());
    }

    // add energy deposition
    aHit->AddEdep(edep);
    aHit->AddEdepQuenched(edep_quenched);

    return true;
}
