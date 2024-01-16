#ifndef AmoreTrackInformation_h
#define AmoreTrackInformation_h 1

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VUserTrackInformation.hh"
#include "globals.hh"

class G4VPhysicalVolume;

class AmoreTrackInformation : public G4VUserTrackInformation {
  public:
    AmoreTrackInformation() = delete;
    AmoreTrackInformation(const G4Track *aTrack);
    virtual ~AmoreTrackInformation();

    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);

    AmoreTrackInformation &operator=(const AmoreTrackInformation &right);

    virtual void Print() const;

    void SetBirthPV(const G4VPhysicalVolume *a) { fBirthPV = a; }
    const G4VPhysicalVolume *GetBirthPV() const { return fBirthPV; }

    void SetParentDefinition(G4ParticleDefinition *a) { fMotherDef = a; }
    const G4ParticleDefinition *GetParentDefinition() const { return fMotherDef; }

  private:
    const G4ParticleDefinition *fMotherDef;
    const G4VPhysicalVolume *fBirthPV;
};

#ifdef G4MULTITHREADED
extern G4ThreadLocal G4Allocator<AmoreTrackInformation> *aTrackInformationAllocator;
#else
extern G4Allocator<AmoreTrackInformation> *aTrackInformationAllocator;
#endif

inline void *AmoreTrackInformation::operator new(size_t) {
    if (!aTrackInformationAllocator)
        aTrackInformationAllocator = new G4Allocator<AmoreTrackInformation>;
    return (void *)aTrackInformationAllocator->MallocSingle();
}
inline void AmoreTrackInformation::operator delete(void *aTrackInfo) {
    aTrackInformationAllocator->FreeSingle((AmoreTrackInformation *)aTrackInfo);
}

#endif
