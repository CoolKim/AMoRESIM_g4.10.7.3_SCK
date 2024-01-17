#include "G4Version.hh"
#if G4VERSION_NUMBER >= 1000

#include "AmoreSim/AmorePLManager.hh"
#include "AmoreSim/AmorePhysicsList.hh"

#include <cstdlib>
#include <iostream>
#include <vector>

#include "G4BuilderType.hh"
#include "G4OpticalPhysics.hh"
#include "G4ThermalNeutrons.hh"
#include "G4UIcommand.hh"

G4PhysListFactory *AmorePLManager::fgPLFactory = nullptr;
int AmorePLManager::fgPLManCnt                 = 0;

using namespace std;
AmorePLManager::AmorePLManager(const std::string &db_name)
    : fBuilt(false), fInitialized(false), fBuildOptical(false), fEnableScintillation(false),
      fEnableCerenkov(false), fOmitHadronPhys(false), fThermalNeutron(false), fUseCupPL(false),
      fPhysicsList(nullptr), fDB(CupStrParam::GetDB()) {
    cout << "AmorePLManager -- uses database file" << db_name << endl;
    OpenDBFile(db_name);
    Initialize();
    fgPLManCnt++;
}

AmorePLManager::AmorePLManager()
    : fBuilt(false), fInitialized(false), fBuildOptical(false), fEnableScintillation(false),
      fEnableCerenkov(false), fOmitHadronPhys(false), fThermalNeutron(false), fUseCupPL(false),
      fPhysicsList(nullptr), fDB(CupStrParam::GetDB()) {
    cout << "AmorePLManager -- uses default database file name PL_settings.dat" << endl;
    OpenDBFile("PL_settings.dat");
    Initialize();
    fgPLManCnt++;
}

AmorePLManager::~AmorePLManager() {
    cout << "Deleting AmorePLManager..." << endl;
    delete fPhysicsList;
    delete fDummyMessenger;
    if (--fgPLManCnt == 0) {
        cout << "There are no instantiated AmorePLManager. Deleting fgPLFactory" << endl;
        delete fgPLFactory;
        fgPLFactory = nullptr;
    } else {
        cout << "There are another AmorePLManager which instantiated. It will not delete "
                "fgPLFactory."
             << endl;
        cout << "fgPLManCnt: " << fgPLManCnt << endl;
    }
}

void AmorePLManager::OpenDBFile(const std::string &db_name) {
    if (db_name.find('/') != db_name.npos)
        fDB.ReadFile(db_name.c_str());
    else {
        char *lEnvStrAmoreData = getenv("AmoreDATA");
        char *lEnvStrCupData   = getenv("CupDATA");
        if (lEnvStrAmoreData != nullptr) {
            string lStringEnvAD = string(lEnvStrAmoreData);
            fDB.ReadFile(string(lStringEnvAD + "/" + db_name).c_str());
        } else if (lEnvStrCupData != nullptr) {
            string lStringEnvCD = string(lEnvStrCupData);
            fDB.ReadFile(string(lStringEnvCD + "/" + db_name).c_str());
        } else {
            string lStringData = string("data");
            fDB.ReadFile(string(lStringData + "/" + db_name).c_str());
        }
    }
}

void AmorePLManager::Initialize() {
    if (fBuilt) return;

    G4String useCupPLStr = fDB.GetWithDefault("UseCupPhysList", "false");
    fUseCupPL            = (useCupPLStr == "true");
    if (fUseCupPL) {
        cout << "CupPhysicsList will be used. All other settings will be ignored." << endl;
        fInitialized = true;
        return;
    }

    if (fgPLFactory == nullptr) {
        fgPLFactory = new G4PhysListFactory;
        cout << "PhysListFactory has been initiated." << endl;
    } else {
        cout << "PhysListFactory already has been initiated." << endl;
    }

    G4String enableOpticalStr = fDB.GetWithDefault("EnableOptical", "false");
    fBuildOptical             = (enableOpticalStr == "true");
    if (fBuildOptical) {
        cout << "Optical physics has been enabled." << endl;
        G4String enableScint = fDB.GetWithDefault("EnableScintillation", "true");
        fEnableScintillation = (enableScint == "true");
        if (fEnableScintillation) cout << "Optical physics -- Scintillation enabled." << endl;
        G4String enableCerenkov = fDB.GetWithDefault("EnableCerenkov", "true");
        fEnableCerenkov         = (enableCerenkov == "true");
        if (fEnableCerenkov) cout << "Optical physics -- Cerenkov enabled." << endl;
    } else
        cout << "Optical physics has been disabled." << endl;

    G4String omitHadPhysStr = fDB.GetWithDefault("OmitHadronPhysics", "false");
    fOmitHadronPhys         = (omitHadPhysStr == "true");
    if (fOmitHadronPhys)
        cout << "Hadron physics omitted." << endl;
    else {
        G4String enableThermalNeut = fDB.GetWithDefault("EnableThermalNeutron", "false");
        fThermalNeutron            = (enableThermalNeut == "true");
        if (fThermalNeutron)
            cout << "Thermal neutron scattering physics has been enabled." << endl;
        else
            cout << "Thermal neutron scattering physics has been disabled." << endl;
    }

    G4String refPhysName   = fDB.GetWithDefault("RefPhysListName", "QGSP_BERT_HP");
    G4String EMPhysicsName = fDB.GetWithDefault("EMPhysicsName", "");
    if (!SetRefPhysListName(refPhysName)) return;
    if (!SetEMPhysicsName(EMPhysicsName)) return;
    fInitialized = true;
}

void AmorePLManager::BuildPhysicsList() {
    if (!fInitialized) return;
    if (fBuilt) return;

    if (fUseCupPL) {
        fPhysicsList = new AmorePhysicsList();
        fBuilt       = true;
        return;
    }

    fPhysicsList    = fgPLFactory->GetReferencePhysList(fRefPLName + fEMPhysName);
    fDummyMessenger = new AmoreCPLDummyMessenger();

    if (fBuildOptical) {
        G4OpticalPhysics *nowOpticalPhysics = new G4OpticalPhysics();
/*
        if (fEnableScintillation)
            nowOpticalPhysics->Configure(kScintillation, true);
        else
            nowOpticalPhysics->Configure(kScintillation, false);
        if (fEnableCerenkov)
            nowOpticalPhysics->Configure(kCerenkov, true);
        else
            nowOpticalPhysics->Configure(kCerenkov, false);
*/

        fPhysicsList->RegisterPhysics(nowOpticalPhysics);
    }

    if (fOmitHadronPhys) {
        fPhysicsList->RemovePhysics(bHadronElastic);
        fPhysicsList->RemovePhysics(bHadronInelastic);
        fPhysicsList->RemovePhysics(bStopping);
        int size = fPhysicsList->GetSubInstanceManager().offset->physicsVector->size();
        vector<G4String> omitList;
        for (int i = 0; i < size; i++) {
            bool isOmittable        = false;
            const G4String &nowName = fPhysicsList->GetPhysics(i)->GetPhysicsName();
            isOmittable |= nowName.contains("hInelastic");
            isOmittable |= nowName.contains("hElastic");
            isOmittable |= nowName.contains("hadronic");
            isOmittable |= nowName.contains("Hadronic");
            isOmittable |= nowName.contains("ionInelastic");
            isOmittable |= nowName.contains("neutron");
            isOmittable |= nowName.contains("Stopping");
            isOmittable |= nowName.contains("stopping");
            cout << "Omitted hadronic process: " << nowName << endl;
            if (isOmittable) omitList.push_back(nowName);
        }
        for (auto nowName : omitList)
            fPhysicsList->RemovePhysics(nowName);
    } else {
        if (fThermalNeutron) fPhysicsList->RegisterPhysics(new G4ThermalNeutrons(0));
    }

    fBuilt = true;
}

#endif
