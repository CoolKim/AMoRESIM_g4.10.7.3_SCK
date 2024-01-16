#include "AmoreSim/AmorePhysicsOp.hh"
#include "AmoreSim/AmoreScintillation.hh"
#include "CupSim/CupOpAttenuation.hh"
#include "CupSim/CupOpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4EmSaturation.hh"
//#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

#include "G4Version.hh"

AmorePhysicsOp::AmorePhysicsOp(const G4String &name) : G4VPhysicsConstructor(name) {
    OpVerbLevel = 0;
}

AmorePhysicsOp::~AmorePhysicsOp() {}

void AmorePhysicsOp::ConstructProcess() {
    // EJ: start
    // scintillation process
    AmoreScintillation *theScintProcessDef = new AmoreScintillation("Scintillation");
    // theScintProcessDef->DumpPhysicsTable();
    theScintProcessDef->SetTrackSecondariesFirst(true);
    theScintProcessDef->SetScintillationYieldFactor(1.0);     //
    theScintProcessDef->SetScintillationExcitationRatio(0.0); //
    theScintProcessDef->SetVerboseLevel(OpVerbLevel);

    G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
    theScintProcessDef->AddSaturation(emSaturation);

    // optical processes
    CupOpAttenuation *theAttenuationProcess = new CupOpAttenuation();
    theAttenuationProcess->UseTimeProfile("exponential");
    theAttenuationProcess->SetVerboseLevel(OpVerbLevel);

    //G4OpBoundaryProcess *theBoundaryProcess = new G4OpBoundaryProcess();
    CupOpBoundaryProcess *theBoundaryProcess = new CupOpBoundaryProcess();
    theBoundaryProcess->SetVerboseLevel(OpVerbLevel);

    // Cerenkov
    G4Cerenkov *theCerenkovProcess = new G4Cerenkov();
    theCerenkovProcess->SetTrackSecondariesFirst(true);

    G4ParticleTable::G4PTblDicIterator *theParticleIterator;
#if G4VERSION_NUMBER <= 999
    theParticleIterator = G4VPhysicsConstructor::theParticleIterator;
#else
    theParticleIterator = GetParticleIterator();
#endif

    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
        G4ParticleDefinition *particle = theParticleIterator->value();
        G4ProcessManager *pmanager     = particle->GetProcessManager();
        G4String particleName          = particle->GetParticleName();
        if (theScintProcessDef->IsApplicable(*particle)) {
            pmanager->AddProcess(theScintProcessDef);
            pmanager->SetProcessOrderingToLast(theScintProcessDef, idxAtRest);
            pmanager->SetProcessOrderingToLast(theScintProcessDef, idxPostStep);
        }
        if (theCerenkovProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theCerenkovProcess);
            pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
        }

        if (particleName == "opticalphoton") {
            pmanager->AddDiscreteProcess(theAttenuationProcess);
            pmanager->AddDiscreteProcess(theBoundaryProcess);
        }
    }
}
