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
/// \file runAndEvent/Amore/include/AmoreTrajectoryPoint.hh
/// \brief Definition of the AmoreTrajectoryPoint class
//
// $Id: $
//
#ifndef AmoreTrajectoryPoint_h
#define AmoreTrajectoryPoint_h 1

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4TrackVector.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4Version.hh"
#include "globals.hh"
#include "G4Step.hh"

class G4Material;
class G4VProcess;
class G4VPhysicalVolume;

//
/// Trajectory point class
///
/// - new, delete and "==" operators are overwritten
///
/// - const G4ThreeVector GetPosition() const
///    gets the position of this trajectory
///
/// - const G4Material* GetMaterial() const
///    gets material that this trajectory has.
///
/// - const std::map<G4String,G4AttDef>* GetAttDefs() const
///    defines the position and the material as attiributes
///
/// - std::vector<G4AttValue>* CreateAttValues() const
///    sets and returns the attributes
//
////////////////////////
class AmoreTrajectoryPoint : public G4VTrajectoryPoint
////////////////////////
{

    //--------
  public:
    //--------

    // Constructor/Destructor
    AmoreTrajectoryPoint(const G4Track *aTrack);
    AmoreTrajectoryPoint(const G4Step *aStep);
    AmoreTrajectoryPoint(const AmoreTrajectoryPoint &right);
    virtual ~AmoreTrajectoryPoint();

    // Operators
    inline void *operator new(size_t);
    inline void operator delete(void *aTrajectoryPoint);
    inline int operator==(const AmoreTrajectoryPoint &right) const { return (this == &right); };

    // Get/Set functions
    inline virtual const G4ThreeVector GetPosition() const { return fPosition; };
    inline virtual const G4ThreeVector GetMomentum() const { return fMomentum; };
    inline G4int GetPDGcode() const { return fPDGCode; }
    inline const G4String GetProcessName() const { return fProcName; }
    inline const G4String GetNextVolName() const { return fNextVolName; }
    inline G4double GetKinE() const { return fKinE; }
    inline G4double GetGlobalTime() const { return fTime; }
    inline G4double GetMass() const { return fMass; }
    inline G4ThreeVector GetPolarization() const { return fPolarization; }
    inline const G4ParticleDefinition *GetPD() const { return fPD; }

    void AddSecondary(G4Track *);
    inline const G4Track *GetSecondary(G4int i) const { return (*fpSec)[i]; }
    inline G4int GetSecondarySize() const { return fpSec->size(); }

    inline void SetNRest(G4int a) { fNRest = a; }
    inline void SetNAlong(G4int a) { fNAlong = a; }
    inline void SetNPost(G4int a) { fNPost = a; }
    inline void SetTrackID(G4int a) { fTrackID = a; }

    inline G4int GetNRest() { return fNRest; }
    inline G4int GetNAlong() { return fNAlong; }
    inline G4int GetNPost() { return fNPost; }
    inline G4int GetTrackID() { return fTrackID; }

    // Get method for HEPRep style attributes
    virtual const std::map<G4String, G4AttDef> *GetAttDefs() const;
    virtual std::vector<G4AttValue> *CreateAttValues() const;

    //---------
  private:
    //---------

    // Member data
    G4ThreeVector fPosition;
    G4ThreeVector fMomentum;
    G4ThreeVector fPolarization;
    G4int fPDGCode;
    G4int fNRest;
    G4int fNAlong;
    G4int fNPost;
    G4int fTrackID;
    G4double fKinE;
    G4double fTime;
    G4double fMass;
    G4String fProcName;
    G4String fNextVolName;
    const G4ParticleDefinition *fPD;
    G4TrackVector *fpSec;
};

#if G4VERSION_NUMBER <= 999
extern G4Allocator<AmoreTrajectoryPoint> *faTrajPointAllocator;
#else
extern G4ThreadLocal G4Allocator<AmoreTrajectoryPoint> *faTrajPointAllocator;
#endif

inline void *AmoreTrajectoryPoint::operator new(size_t) {
    void *aTrajectoryPoint;
    if (faTrajPointAllocator == nullptr)
        faTrajPointAllocator = new G4Allocator<AmoreTrajectoryPoint>;
    aTrajectoryPoint = (void *)faTrajPointAllocator->MallocSingle();
    return aTrajectoryPoint;
}

inline void AmoreTrajectoryPoint::operator delete(void *aTrajectoryPoint) {
    faTrajPointAllocator->FreeSingle((AmoreTrajectoryPoint *)aTrajectoryPoint);
}

#endif
