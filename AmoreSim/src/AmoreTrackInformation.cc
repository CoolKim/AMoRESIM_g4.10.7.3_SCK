#include "AmoreSim/AmoreTrackInformation.hh"

#ifdef G4MULTITHREADED
G4ThreadLocal G4Allocator<AmoreTrackInformation> *aTrackInformationAllocator = nullptr;
#else
G4Allocator<AmoreTrackInformation> *aTrackInformationAllocator = nullptr;
#endif

void AmoreTrackInformation::Print() const { G4cout << "There are no words to tell you!" << G4endl; }

AmoreTrackInformation::AmoreTrackInformation(const G4Track *aTrack)
    : fMotherDef(aTrack->GetDefinition()), fBirthPV(nullptr) {}

AmoreTrackInformation &AmoreTrackInformation::operator=(const AmoreTrackInformation &right) {
    fMotherDef = right.fMotherDef;
    fBirthPV   = right.fBirthPV;
    return *this;
}

AmoreTrackInformation::~AmoreTrackInformation() {}
