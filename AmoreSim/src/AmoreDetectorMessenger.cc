////////////////////////////////////////////////////////////////
// AmoreDetectorMessenger
////////////////////////////////////////////////////////////////

#include "AmoreSim/AmoreDetectorMessenger.hh"
#include "AmoreSim/AmoreDetectorConstruction.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "fstream"  // for file streams
#include <stdlib.h> // for strtol

AmoreDetectorMessenger::AmoreDetectorMessenger(AmoreDetectorConstruction *Amoredetector)
    : AmoreDetector(Amoredetector) {
    // the AmoreDetector directory
    AmoreDetectorDir = new G4UIdirectory("/detGeometry/");
    AmoreDetectorDir->SetGuidance("Control the detector geometry options.");

    Amore200_DetectorDir = new G4UIdirectory("/detGeometry/200/");
    Amore200_DetectorDir->SetGuidance("Control the specific geometry options for AMoRE 200.");

    AmoreI_DetectorDir = new G4UIdirectory("/detGeometry/I/");
    AmoreI_DetectorDir->SetGuidance("Control the specific geometry options for AMoRE-I.");

    AmorePilot_DetectorDir = new G4UIdirectory("/detGeometry/Pilot/");
    AmorePilot_DetectorDir->SetGuidance("Control the specific geometry options for AMoRE Pilot.");

    AmorePilotRUN5_DetectorDir = new G4UIdirectory("/detGeometry/PilotRUN5/");
    AmorePilotRUN5_DetectorDir->SetGuidance("Control the specific geometry options for AMoRE PilotRUN5.");

    // the select command
    DetGeometrySelectCmd = new G4UIcommand("/detGeometry/select", this);
    DetGeometrySelectCmd->SetGuidance("Select which detector you want to build");
    DetGeometrySelectCmd->SetGuidance(
        "Use with no parameters to get list of available detector styles.");
    DetGeometrySelectCmd->AvailableForStates(G4State_PreInit);
    DetGeometrySelectCmd->SetParameter(new G4UIparameter("which", 's', true));

    NeutronModeCmd = new G4UIcommand("/detGeometry/NeutronMode", this);
    NeutronModeCmd->SetGuidance("Select enable neutron mode");
    NeutronModeCmd->AvailableForStates(G4State_PreInit);
    NeutronModeCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    RockgammaModeCmd = new G4UIcommand("/detGeometry/RockgammaMode", this);
    RockgammaModeCmd->SetGuidance("Select enable rockgamma simulation mode");
    RockgammaModeCmd->AvailableForStates(G4State_PreInit);
    RockgammaModeCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    AdditionalPECmd = new G4UIcommand("/detGeometry/AdditionalPE", this);
    AdditionalPECmd->SetGuidance("Select enable additional PE shield");
    AdditionalPECmd->AvailableForStates(G4State_PreInit);
    AdditionalPECmd->SetParameter(new G4UIparameter("enable", 'b', true));

    DebugModeCmd = new G4UIcommand("/detGeometry/DebugMode", this);
    DebugModeCmd->SetGuidance("Select enable debug mode");
    DebugModeCmd->AvailableForStates(G4State_PreInit);
    DebugModeCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    OverlapCheckCmd = new G4UIcommand("/detGeometry/OverlapCheck", this);
    OverlapCheckCmd->SetGuidance("Select enable overlap check");
    OverlapCheckCmd->AvailableForStates(G4State_PreInit);
    OverlapCheckCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    NeutShieldConfCmd = new G4UIcommand("/detGeometry/nShieldingToyConf", this);
    NeutShieldConfCmd->SetGuidance(
        "Select configurations for neutron shielding toy models in neutron mode");
    NeutShieldConfCmd->AvailableForStates(G4State_PreInit);
    NeutShieldConfCmd->SetParameter(new G4UIparameter("enable", 's', true));

    EnableOrigGeomCmd = new G4UIcommand("/detGeometry/EnableOrigGeom", this);
    EnableOrigGeomCmd->SetGuidance("Select enable original geometry");
    EnableOrigGeomCmd->AvailableForStates(G4State_PreInit);
    EnableOrigGeomCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableScintCmd = new G4UIcommand("/detGeometry/EnableScint", this);
    EnableScintCmd->SetGuidance("Select enable scintillator geometry");
    EnableScintCmd->AvailableForStates(G4State_PreInit);
    EnableScintCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableGantryCmd = new G4UIcommand("/detGeometry/EnableGantry", this);
    EnableGantryCmd->SetGuidance("Select enable gantry geometry");
    EnableGantryCmd->AvailableForStates(G4State_PreInit);
    EnableGantryCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableInnerDetCmd = new G4UIcommand("/detGeometry/EnableInnerDet", this);
    EnableInnerDetCmd->SetGuidance("Select enable inner detector geometry");
    EnableInnerDetCmd->AvailableForStates(G4State_PreInit);
    EnableInnerDetCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableInnermostCmd = new G4UIcommand("/detGeometry/EnableInnermost", this);
    EnableInnermostCmd->SetGuidance("Select enable innermost geometry in cryostat");
    EnableInnermostCmd->AvailableForStates(G4State_PreInit);
    EnableInnermostCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableRealConfCmd = new G4UIcommand("/detGeometry/EnableNeutronShield", this);
    EnableRealConfCmd->SetGuidance("Select enable real conf. for neutron shielding");
    EnableRealConfCmd->AvailableForStates(G4State_PreInit);
    EnableRealConfCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableRealConfCmd = new G4UIcommand("/detGeometry/EnableAdditionalPEShield", this);
    EnableRealConfCmd->SetGuidance("Select enable real conf. for neutron shielding");
    EnableRealConfCmd->AvailableForStates(G4State_PreInit);
    EnableRealConfCmd->SetParameter(new G4UIparameter("enable", 'b', true));

		AMoRE200PhaseSelectCmd = new G4UIcommand("/detGeometry/200/selectPhase", this);
		AMoRE200PhaseSelectCmd->SetGuidance("Select which phase you want to build");
		AMoRE200PhaseSelectCmd->AvailableForStates(G4State_PreInit);
		AMoRE200PhaseSelectCmd->SetParameter(new G4UIparameter("which", 's', true));

    CavernTypeSelectCmd = new G4UIcommand("/detGeometry/200/selectCavern", this);
    CavernTypeSelectCmd->SetGuidance("Select which cavern style you want to build");
    CavernTypeSelectCmd->AvailableForStates(G4State_PreInit);
    CavernTypeSelectCmd->SetParameter(new G4UIparameter("which", 's', true));

    EnableMumetalCmd = new G4UIcommand("/detGeometry/Pilot/EnableMumetal", this);
    EnableMumetalCmd->SetGuidance("Select enable mu-metal geometry");
    EnableMumetalCmd->AvailableForStates(G4State_PreInit);
    EnableMumetalCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableMumetalCmd = new G4UIcommand("/detGeometry/PilotRUN5/EnableMumetal", this);
    EnableMumetalCmd->SetGuidance("Select enable mu-metal geometry");
    EnableMumetalCmd->AvailableForStates(G4State_PreInit);
    EnableMumetalCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableTargRoomCmd = new G4UIcommand("/detGeometry/Pilot/EnableTargetRoom", this);
    EnableTargRoomCmd->SetGuidance("Select enable target room geometry");
    EnableTargRoomCmd->AvailableForStates(G4State_PreInit);
    EnableTargRoomCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableTargRoomCmd = new G4UIcommand("/detGeometry/PilotRUN5/EnableTargetRoom", this);
    EnableTargRoomCmd->SetGuidance("Select enable target room geometry");
    EnableTargRoomCmd->AvailableForStates(G4State_PreInit);
    EnableTargRoomCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableSuperMagneticShieldCmd =
        new G4UIcommand("/detGeometry/I/EnableSuperCMagneticShield", this);
    EnableSuperMagneticShieldCmd->SetGuidance(
        "Select enable superconducting magnet shielding of lead");
    EnableSuperMagneticShieldCmd->AvailableForStates(G4State_PreInit);
    EnableSuperMagneticShieldCmd->SetParameter(new G4UIparameter("enable", 'b', true));

    EnableCrystalArray = new G4UIcommand("/detGeometry/I/EnableCrystalArray", this);
    EnableCrystalArray->SetGuidance("Select enable crystal array for AMoRE-I");
    EnableCrystalArray->AvailableForStates(G4State_PreInit);
    EnableCrystalArray->SetParameter(new G4UIparameter("enable", 'b', true));
}

AmoreDetectorMessenger::~AmoreDetectorMessenger() {
    delete DetGeometrySelectCmd;
		delete AMoRE200PhaseSelectCmd;
    delete CavernTypeSelectCmd;

    delete EnableOrigGeomCmd;
    delete EnableScintCmd;
    delete EnableGantryCmd;
    delete EnableInnerDetCmd;
    delete EnableMumetalCmd;
    delete EnableInnermostCmd;
    delete EnableTargRoomCmd;
    delete EnableRealConfCmd;
    delete EnableSuperMagneticShieldCmd;
    delete EnableCrystalArray;
    delete NeutronModeCmd;
		delete RockgammaModeCmd;
    delete AdditionalPECmd;
    delete NeutShieldConfCmd;
		delete OverlapCheckCmd;
		delete DebugModeCmd;

    delete AmoreDetectorDir;
    delete Amore200_DetectorDir;
    delete AmorePilot_DetectorDir;
    delete AmorePilotRUN5_DetectorDir;
}

void AmoreDetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValues) {
    // GeometrySelectCmd
    if (command == DetGeometrySelectCmd) {
        if (newValues.length() == 0) {
            G4cout << "Available detector geometries: ";
            for (int i = 0; i < AmoreDetector->GetNumDetGeometryTypes(); i++) {
                AmoreDetectorConstruction::eDetGeometry nowEnum =
                    static_cast<AmoreDetectorConstruction::eDetGeometry>(i);
                G4cout << " " << AmoreDetector->GetDetGeometryTypeName(nowEnum);
            }
            G4cout << G4endl;
        } else {
            for (int i = 0; i < AmoreDetector->GetNumDetGeometryTypes(); i++) {
                AmoreDetectorConstruction::eDetGeometry nowEnum =
                    static_cast<AmoreDetectorConstruction::eDetGeometry>(i);
                if (newValues.compare(AmoreDetector->GetDetGeometryTypeName(nowEnum)) == 0) {
                    // EJ: start
                    G4cout << "detector/select " << AmoreDetector->GetDetGeometryTypeName(nowEnum)
                           << G4endl;
                    // EJ: end
                    AmoreDetector->SetWhichDetGeometry(nowEnum);
                    return;
                }
            }
            G4cerr << "Unknown detector geometry style " << newValues << G4endl;
        }
    } else if (command == AMoRE200PhaseSelectCmd) {
			if( newValues.length() == 0){
				G4cout << "Available AMoRE-II phase: " ;
				for(int i = 0; i < AmoreDetector->GetAMoRE200PhaseNumber(); i++)
					G4cout << " " << AmoreDetector->GetAMoRE200PhaseName( static_cast<AmoreDetectorConstruction::ePhaseAMoRE200>(i));
				G4cout << G4endl;
			} else {
				for(int i = 0; i < AmoreDetector->GetAMoRE200PhaseNumber(); i++){
					AmoreDetectorConstruction::ePhaseAMoRE200 nowEnum = static_cast<AmoreDetectorConstruction::ePhaseAMoRE200>(i);
					if(newValues.compare(AmoreDetector->GetAMoRE200PhaseName(nowEnum)) == 0){
						G4cout << "AMoRE200phase/select " << AmoreDetector->GetAMoRE200PhaseName(nowEnum) << G4endl;
						AmoreDetector->SetWhichAMoRE200Phase(nowEnum);
					}
				}
			}
		}
		else if (command == CavernTypeSelectCmd) {
        if (newValues.length() == 0) {
            G4cout << "Available cavern types(for AMoRE-II): ";
            for (int i = 0; i < AmoreDetector->GetNumCavernType(); i++)
                G4cout << " "
                       << AmoreDetector->GetCavernTypeName(
                              static_cast<AmoreDetectorConstruction::eCavernType>(i));
            G4cout << G4endl;
        } else {
            for (int i = 0; i < AmoreDetector->GetNumCavernType(); i++) {
                AmoreDetectorConstruction::eCavernType nowEnum =
                    static_cast<AmoreDetectorConstruction::eCavernType>(i);
                if (newValues.compare(AmoreDetector->GetCavernTypeName(nowEnum)) == 0) {
                    G4cout << "CavernType/select " << AmoreDetector->GetCavernTypeName(nowEnum)
                           << G4endl;
                    AmoreDetector->SetWhichCavernType(nowEnum);
                    return;
                }
            }
            G4cerr << "Unknown cavern type " << newValues << G4endl;
        }
    } else if (command == NeutShieldConfCmd) {
        if (newValues.length() == 0) {
            G4cout << "Available neutron shielding configurations(for AMoRE-Pilot): ";
            for (int i = 0; i < AmoreDetector->GetNumNShieldConfTypes(); i++)
                G4cout << " "
                       << AmoreDetector->GetNShieldConfTypeName(
                              static_cast<AmoreDetectorConstruction::eNShieldConf>(i));
            G4cout << G4endl;
        } else {
            for (int i = 0; i < AmoreDetector->GetNumNShieldConfTypes(); i++) {
                AmoreDetectorConstruction::eNShieldConf nowEnum =
                    static_cast<AmoreDetectorConstruction::eNShieldConf>(i);
                if (newValues.compare(AmoreDetector->GetNShieldConfTypeName(nowEnum)) == 0) {
                    G4cout << "NeutShieldConf/select "
                           << AmoreDetector->GetNShieldConfTypeName(nowEnum) << G4endl;
                    AmoreDetector->SetNShieldConfType(nowEnum);
                    return;
                }
            }
            G4cerr << "Unknown neutron shielding configuration " << newValues << G4endl;
        }
    } else if (command == EnableOrigGeomCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetEnableOriginalGeometry(inp);
    } else if (command == EnableScintCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetEnableScintillator(inp);
    } else if (command == EnableGantryCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetEnableGantry(inp);
    } else if (command == EnableInnerDetCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetEnableInnerDetector(inp);
    } else if (command == EnableMumetalCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->Set_Pilot_EnableMumetal(inp);
        AmoreDetector->Set_PilotRUN5_EnableMumetal(inp);
    } else if (command == EnableInnermostCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetEnableInnermost(inp);
    } else if (command == EnableTargRoomCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->Set_Pilot_EnableTargetRoom(inp);
        AmoreDetector->Set_PilotRUN5_EnableTargetRoom(inp);
    } else if (command == EnableRealConfCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetEnableRealConf(inp);
    } else if (command == NeutronModeCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetNeutronMode(inp);
    } else if (command == RockgammaModeCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetRockgammaMode(inp);
    } else if (command == AdditionalPECmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetAdditionalPE(inp);
    } else if (command == DebugModeCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetDebugMessage(inp);
    } else if (command == OverlapCheckCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->SetOverlapCheck(inp);
    } else if (command == EnableSuperMagneticShieldCmd) {
        G4bool inp = StoB(newValues);
        AmoreDetector->Set_I_EnableSuperConductingShield(inp);
    } else if (command == EnableCrystalArray) {
        G4bool inp = StoB(newValues);
        AmoreDetector->Set_I_EnableCrystalArray(inp);
    } else { // invalid command
        G4cerr << "invalid detector \"set\" command\n" << G4endl;
    }
}

G4String AmoreDetectorMessenger::GetCurrentValue(G4UIcommand *command) {
    // GeometrySelectCmd
    if (command == DetGeometrySelectCmd) {
        return AmoreDetector->GetDetGeometryTypeName(AmoreDetector->GetWhichDetGeometry());
    } else if (command == AMoRE200PhaseSelectCmd) {
			  return AmoreDetector->GetAMoRE200PhaseName(AmoreDetector->GetAMoRE200PhaseType());
		}else if (command == CavernTypeSelectCmd) {
        return AmoreDetector->GetCavernTypeName(AmoreDetector->GetCavernType());
    } else if (command == NeutShieldConfCmd) {
        return AmoreDetector->GetNShieldConfTypeName(AmoreDetector->GetNShieldConfType());
    } else if (command == EnableOrigGeomCmd) {
        return BtoS(AmoreDetector->GetEnableOriginalGeometry());
    } else if (command == EnableScintCmd) {
        return BtoS(AmoreDetector->GetEnableScintillator());
    } else if (command == EnableGantryCmd) {
        return BtoS(AmoreDetector->GetEnableGantry());
    } else if (command == EnableInnerDetCmd) {
        return BtoS(AmoreDetector->GetEnableInnerDetector());
    } else if (command == EnableMumetalCmd) {
        return BtoS(AmoreDetector->Get_Pilot_EnableMumetal());
        return BtoS(AmoreDetector->Get_PilotRUN5_EnableMumetal());
    } else if (command == EnableInnermostCmd) {
        return BtoS(AmoreDetector->GetEnableInnermost());
    } else if (command == EnableTargRoomCmd) {
        return BtoS(AmoreDetector->Get_Pilot_EnableTargetRoom());
        return BtoS(AmoreDetector->Get_PilotRUN5_EnableTargetRoom());
    } else if (command == EnableRealConfCmd) {
        return BtoS(AmoreDetector->GetEnableNeutronShield());
    } else if (command == NeutronModeCmd) {
        return BtoS(AmoreDetector->GetNeutronMode());
    } else if (command == RockgammaModeCmd) {
        return BtoS(AmoreDetector->GetRockgammaMode());
    } else if (command == AdditionalPECmd) {
        return BtoS(AmoreDetector->GetAdditionalPE());
    } else if (command == DebugModeCmd) {
        return BtoS(AmoreDetector->GetDebugMessage());
    } else if (command == OverlapCheckCmd) {
        return BtoS(AmoreDetector->GetOverlapCheck());
    } else if (command == EnableSuperMagneticShieldCmd) {
        return BtoS(AmoreDetector->Get_I_EnableSuperConductingShield());
    } else if (command == EnableCrystalArray) {
        return BtoS(AmoreDetector->Get_I_EnableCrystalArray());
    } else { // invalid command
        return G4String("invalid AmoreDetectorMessenger \"get\" command");
    }
}
