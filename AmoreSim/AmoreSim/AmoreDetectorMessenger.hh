// This file is part of the GenericLAND software library.
// $Id: AmoreDetectorMessenger.hh,v 1.3 2016/11/03 06:20:26 ysy Exp $
//
// AmoreDetectorMessenger.hh by Glenn Horton-Smith, Dec. 1999
#ifndef __AmoreDetectorMessenger_hh__
#define __AmoreDetectorMessenger_hh__ 1

#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class AmoreDetectorConstruction;

class AmoreDetectorMessenger : public G4UImessenger {
  public:
    AmoreDetectorMessenger(AmoreDetectorConstruction *AmoreDetector);
    ~AmoreDetectorMessenger();

    void SetNewValue(G4UIcommand *command, G4String newValues);
    G4String GetCurrentValue(G4UIcommand *command);

  private:
    AmoreDetectorConstruction *AmoreDetector;

    G4UIdirectory *AmoreDetectorDir;
    G4UIdirectory *AmoreI_DetectorDir;
    G4UIdirectory *Amore200_DetectorDir;
    G4UIdirectory *AmorePilot_DetectorDir;
    G4UIdirectory *AmorePilotRUN5_DetectorDir;

    G4UIcommand *DetGeometrySelectCmd;

    // For AMoRE 2
    G4UIcommand *VetoGeometrySelectCmd;
		G4UIcommand *AMoRE200PhaseSelectCmd;
		G4UIcommand *HatSizeSelectCmd;
    G4UIcommand *CavernTypeSelectCmd;
    G4UIcommand *NeutronModeCmd;
		G4UIcommand *RockgammaModeCmd;
		G4UIcommand *AdditionalPECmd;
    G4UIcommand *NeutShieldConfCmd;
		G4UIcommand *DebugModeCmd;
		G4UIcommand *OverlapCheckCmd;

    // For AMoRE I
    G4UIcommand *EnableSuperMagneticShieldCmd;
    G4UIcommand *EnableCrystalArray;

    // For AMoRE Pilot
    G4UIcommand *EnableOrigGeomCmd;
    G4UIcommand *EnableScintCmd;
    G4UIcommand *EnableGantryCmd;
    G4UIcommand *EnableInnerDetCmd;
    G4UIcommand *EnableMumetalCmd;
    G4UIcommand *EnableInnermostCmd;
    G4UIcommand *EnableTargRoomCmd;
    G4UIcommand *EnableRealConfCmd;

    class AmoreDetectorMessenger *myMessenger;
};

#endif
