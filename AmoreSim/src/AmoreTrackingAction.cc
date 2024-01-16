#include "G4EventManager.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"

#include "AmoreSim/AmoreTrackInformation.hh"
#include "AmoreSim/AmoreTrackingAction.hh"
#include "AmoreSim/AmoreTrajectory.hh"
#include "CupSim/CupRecorderBase.hh"
#include "G4Trajectory.hh"

AmoreTrackingAction::AmoreTrackingAction(AmoreRootNtuple *r)
    : tracknum(0), CupTrackingAction(r), recorder(r) {}

void AmoreTrackingAction::PreUserTrackingAction(const G4Track *aTrack) {
    if (aTrack->GetTrackID() == 1)
        tracknum = 0;
    else if (aTrack->GetDefinition()->GetParticleName() != "opticalphoton")
        tracknum++;

    if (aTrack->GetParentID() != 0) {
        AmoreTrackInformation *aTrackInformation =
            static_cast<AmoreTrackInformation *>(aTrack->GetUserInformation());
        aTrackInformation->SetBirthPV(aTrack->GetVolume());
    }

    CupTrackingAction::PreUserTrackingAction(aTrack);
}

void AmoreTrackingAction::PostUserTrackingAction(const G4Track *aTrack) {
    CupTrackingAction::PostUserTrackingAction(aTrack);

    if (aTrack->GetDefinition()->GetParticleName() == "opticalphoton")
        return;
    else if (recorder)
        recorder->RecordET(aTrack);

    G4TrackVector *aSecondaries = fpTrackingManager->GimmeSecondaries();

    for (auto &now2nd : *aSecondaries)
        now2nd->SetUserInformation(new AmoreTrackInformation(aTrack));

    if (tracknum > 1000000) {
        G4EventManager::GetEventManager()->AbortCurrentEvent();
    }
}
