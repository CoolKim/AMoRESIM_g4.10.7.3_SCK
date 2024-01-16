#ifndef __h_AmorePLManager_
#define __h_AmorePLManager_

#include <string>
#include <vector>

#include "G4PhysListFactory.hh"
#include "G4UImessenger.hh"
#include "G4VModularPhysicsList.hh"
#include "G4Version.hh"

#include "CupSim/CupStrParam.hh"

class G4UIcommand;

#if G4VERSION_NUMBER >= 1000
class AmoreCPLDummyMessenger : public G4UImessenger {
  public:
    AmoreCPLDummyMessenger()
        : fGCC(nullptr), fECC(nullptr), fPCC(nullptr), fACC(nullptr), fMCC(nullptr), fPLC(nullptr),
          fLC(nullptr), fCRC(nullptr), fVC(nullptr), fYFC(nullptr), fCSCON(nullptr),
          fCSCOFF(nullptr), fCSCV(nullptr) {

        fGCC    = new G4UIcommand("/Cup/phys/CutGamma", this);
        fECC    = new G4UIcommand("/Cup/phys/CutEl", this);
        fPCC    = new G4UIcommand("/Cup/phys/CutPos", this);
        fACC    = new G4UIcommand("/Cup/phys/CutsAll", this);
        fMCC    = new G4UIcommand("/Cup/phys/DetectorCuts", this);
        fPLC    = new G4UIcommand("/Cup/phys/Physics", this);
        fLC     = new G4UIcommand("/Cup/phys/ListPhysics", this);
        fVC     = new G4UIcommand("/Cup/phys/verbose", this);
        fYFC    = new G4UIcommand("/Cup/phys/alphaYieldFactor", this);
        fCRC    = new G4UIcommand("/Cup/phys/EnableCrystalRegion", this);
        fCSCON  = new G4UIcommand("/cupscint/on", this);
        fCSCOFF = new G4UIcommand("/cupscint/off", this);
        fCSCV   = new G4UIcommand("/cupscint/verbose", this);
    };
    virtual ~AmoreCPLDummyMessenger() {
        delete fGCC;
        delete fECC;
        delete fPCC;
        delete fACC;
        delete fMCC;
        delete fPLC;
        delete fLC;
        delete fCRC;
        delete fVC;
        delete fYFC;
        delete fCSCON;
        delete fCSCOFF;
        delete fCSCV;
    };

  private:
    G4UIcommand *fGCC;
    G4UIcommand *fECC;
    G4UIcommand *fPCC;
    G4UIcommand *fACC;
    G4UIcommand *fMCC;
    G4UIcommand *fPLC;
    G4UIcommand *fLC;
    G4UIcommand *fCRC;
    G4UIcommand *fVC;
    G4UIcommand *fYFC;
    G4UIcommand *fCSCON;
    G4UIcommand *fCSCOFF;
    G4UIcommand *fCSCV;
};

class AmorePLManager {
  public:
    inline G4VUserPhysicsList *GetPhysicsList() {
        return (fBuilt == true) ? fPhysicsList : nullptr;
    }
    inline bool IsBuilt() const { return fBuilt; }
    void OpenDBFile(const std::string &db_name);
    AmorePLManager(const std::string &db_name);
    AmorePLManager();

    void BuildPhysicsList();
    inline bool SetRefPhysListName(const std::string &name) {
        const std::vector<G4String> &lRefPLList = fgPLFactory->AvailablePhysLists();
        for (auto nowPLName : lRefPLList) {
            if (name.length() == 0) {
                std::cout << "Empty reference physics list name" << std::endl;
                std::cout << "No reference physics list has been selected." << std::endl;
                return false;
            }
            if (name == nowPLName.c_str()) {
                std::cout << "Reference Physics List " << name << " has been selected."
                          << std::endl;
                fRefPLName = name;
                return true;
            }
        }
        std::cout << "Wrong reference physics list name: " << name << std::endl;
        std::cout << "No reference physics list has been selected." << std::endl;
        return false;
    }
    inline bool SetEMPhysicsName(const std::string &name) {
        const std::vector<G4String> &lRefEMPhysList = fgPLFactory->AvailablePhysListsEM();
        if (name == "default" || name.length() == 0) {
            std::cout << "Default EM physics has been selected." << std::endl;
            fEMPhysName = "";
            return true;
        }
        G4String EMPhysName;
        if (name[0] != '_')
            EMPhysName = "_" + name;
        else
            EMPhysName = name;
        for (auto nowEMPhysName : lRefEMPhysList) {
            if (EMPhysName == nowEMPhysName.c_str()) {
                std::cout << "EM physics " << name << " has been selected." << std::endl;
                fEMPhysName = EMPhysName;
                return true;
            }
        }
        std::cout << "Wrong EM physics name: " << name << std::endl;
        std::cout << "No EM physics has been selected." << std::endl;
        return false;
    }

    inline const std::string &GetRefPhysListName() const { return fRefPLName; }
    inline const std::string &GetEMPhysicsName() const { return fEMPhysName; }
    inline const std::vector<G4String> &GetAvailableRefPLList() const {
        return fgPLFactory->AvailablePhysLists();
    }
    inline const std::vector<G4String> &GetAvailableEMPhysList() const {
        return fgPLFactory->AvailablePhysListsEM();
    }

    ~AmorePLManager();

  private:
    bool fBuilt;
    bool fInitialized;

    bool fBuildOptical;
    bool fEnableScintillation;
    bool fEnableCerenkov;

    bool fOmitHadronPhys;
    bool fThermalNeutron;

    bool fUseCupPL;

    std::string fRefPLName;
    std::string fEMPhysName;

  protected:
    virtual void Initialize();
    static G4PhysListFactory *fgPLFactory;
    static int fgPLManCnt;
    G4VModularPhysicsList *fPhysicsList;
    CupStrParam &fDB;

  private:
    G4UImessenger *fDummyMessenger;
};
#endif
#endif
