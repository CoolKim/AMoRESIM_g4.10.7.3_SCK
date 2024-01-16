// AmoreCMOParameterisation source file
// 2017/9/11 by BaseHardware

#include "AmoreSim/AmoreCMOParameterisation.hh"

#include "G4ExtrudedSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AmoreCMOParameterisation::AmoreCMOParameterisation(G4double aCellsize_R, G4double aCellsize_H,
                                                   G4double aRSpacing, G4double aZSpacing,
                                                   eTranslationMethod aMethod)
    : G4VPVParameterisation(), fInitialized(false), fCSR(aCellsize_R), fCSH(aCellsize_H),
      fRS(aRSpacing), fHS(aZSpacing), fFN(1), fLN(0), fM(aMethod) {
    if (!CheckSettings()) {
        G4cout
            << "AmoreCMOParameterisation::AmoreCMOParameterisation says:" << G4endl
            << "This instance was constructed with optimizing option. The instance will not to be"
            << G4endl
            << "initialized in this time. Please call Optimize() then Initialize() to initilize "
               "this instance."
            << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AmoreCMOParameterisation::AmoreCMOParameterisation(G4double aCellsize_R, G4double aCellsize_H,
                                                   G4double aRSpacing, G4double aZSpacing,
                                                   G4int aFloorNum, G4int aLayerNum,
                                                   eTranslationMethod aMethod)
    : G4VPVParameterisation(), fInitialized(false), fCSR(aCellsize_R), fCSH(aCellsize_H),
      fRS(aRSpacing), fHS(aZSpacing), fFN(aFloorNum), fLN(aLayerNum), fM(aMethod) {
    Initialize();
}

AmoreCMOParameterisation::~AmoreCMOParameterisation() { delete[] fTranslationTable; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool AmoreCMOParameterisation::Optimize(G4int TargetNum, G4double TargetAR) {
    switch (fM) {
        case kChuhMethod: {
            std::vector<G4ThreeVector> candidates;
            G4int i, j;
            for (i = 0; (j = TargetNum / (2 * i * (i + 1) + 1)) > 0; i++) {
                G4double tEnvelL = (i * (2 * fCSR + fRS) + 2 * fCSR) * 2;
                G4double tTotH   = fCSH * j + fHS * (j - 1);
                G4int tTotNum    = (2 * i * (i + 1) + 1) * j;
                G4ThreeVector tCand(i, tTotNum, tEnvelL / tTotH);
                candidates.push_back(tCand);
            }

            G4double min   = DBL_MAX;
            G4int ResultLN = 0, ResultFN = 1;
            for (auto now : candidates) {
                if (std::abs(now.z() - TargetAR) + (TargetNum - now.y()) / TargetNum < min) {
                    min      = std::abs(now.z() - TargetAR) + (TargetNum - now.y()) / TargetNum;
                    ResultLN = now.x();
                    ResultFN = TargetNum / (2 * ResultLN * (ResultLN + 1) + 1);
                }
            }

            if (min == DBL_MAX) {
                G4cout << "Optimization failed." << G4endl;
                return true;
            }

            SetFloorNumber(ResultFN);
            SetLayerNumber(ResultLN);
            /*std::vector<G4ThreeVector> candidates;
            G4int i,j;
            for(i=0;(j=TargetNum/(2*i*(i+1)+1))>0;i++){
              G4double tEnvelL = i*(2*fCSR+fRS) + 2*fCSR;
              G4double tTotH = fCSH*j + fHS*(j-1);
              G4int tTotNum = (2*i*(i+1)+1)*j;
              G4ThreeVector tCand(i, tTotNum, tEnvelL/tTotH);
              candidates.push_back(tCand);
            }

            G4double min = DBL_MAX;
            G4int ResultLN=0, ResultFN=1;
            for(auto now : candidates){
              if(std::abs(now.z()-TargetAR) + (TargetNum-now.y())/TargetNum < min){
                min = std::abs(now.z()-TargetAR) + (TargetNum-now.y())/TargetNum;
                ResultLN = now.x();
                ResultFN = TargetNum/(2*ResultLN*(ResultLN+1)+1);
              }
            }

            if(min == DBL_MAX) {
              G4cout << "Optimization failed." << G4endl;
              return true;
            }

            SetFloorNumber(ResultFN);
            SetLayerNumber(ResultLN);*/
            break;
        }
        case kBasehwMethod: {
            std::vector<G4ThreeVector> candidates;
            G4int i, j;
            for (i = 0; (j = TargetNum / (3 * i * (i + 1) + 1)) > 0; i++) {
                G4double tEnvelL = (i * (2 * fCSR + fRS) + fCSR / std::cos(30 * deg)) * 2;
                G4double tTotH   = fCSH * j + fHS * (j - 1);
                G4int tTotNum    = (3 * i * (i + 1) + 1) * j;
                G4ThreeVector tCand(i, tTotNum, tEnvelL / tTotH);
                candidates.push_back(tCand);
            }

            G4double min   = DBL_MAX;
            G4int ResultLN = 0, ResultFN = 1;
            for (auto now : candidates) {
                if (std::abs(now.z() - TargetAR) + (TargetNum - now.y()) / TargetNum < min) {
                    min      = std::abs(now.z() - TargetAR) + (TargetNum - now.y()) / TargetNum;
                    ResultLN = now.x();
                    ResultFN = TargetNum / (3 * ResultLN * (ResultLN + 1) + 1);
                }
            }

            if (min == DBL_MAX) {
                G4cout << "Optimization failed." << G4endl;
                return true;
            }

            SetFloorNumber(ResultFN);
            SetLayerNumber(ResultLN);
            break;
        }
        default:
            return true;
    }
    return false;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool AmoreCMOParameterisation::Initialize() {
    G4bool retval = false;
    if ((retval = CheckSettings())) return retval;

    switch (fM) {
        case kChuhMethod:
            fEnvelL       = fLN * (2 * fCSR + fRS) + 2 * fCSR;
            fTotH         = fCSH * fFN + fHS * (fFN - 1);
            fFloorCellNum = 2 * (fLN + 2) * (fLN + 1) - (fLN + 1) * 4 + 1;
            fTotCellNum   = fFloorCellNum * fFN;
            fCellSolid    = new G4Box("Envelope-CMOCell", fCSR, fCSR, fCSH / 2.);

            fTranslationTable = new G4ThreeVector[fTotCellNum];
            if ((retval = GenerateTranslationTable())) return retval;
            if ((retval = BuildEnvelopeSolid())) return retval;
            break;

        case kBasehwMethod:
            fEnvelL       = fLN * (2 * fCSR + fRS) + fCSR / std::cos(30 * deg);
            fTotH         = fCSH * fFN + fHS * (fFN - 1);
            fFloorCellNum = 3 * (fLN + 2) * (fLN + 1) - (fLN + 1) * 6 + 1;
            fTotCellNum   = fFloorCellNum * fFN;
            fCellSolid    = new G4Tubs("Envelope-CMOCell", 0, fCSR, fCSH / 2., 0, 360 * deg);

            fTranslationTable = new G4ThreeVector[fTotCellNum];
            if ((retval = GenerateTranslationTable())) return retval;
            if ((retval = BuildEnvelopeSolid())) return retval;
            break;
    }
    fInitialized = true;
    return retval;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool AmoreCMOParameterisation::BuildEnvelopeSolid() {
    switch (fM) {
        case kChuhMethod: {
            std::vector<G4TwoVector> Rhombus_Templete;
            G4TwoVector Rhombus_Points[4] = {G4TwoVector(1, 0), G4TwoVector(0, 1),
                                             G4TwoVector(-1, 0), G4TwoVector(0, -1)};

            Rhombus_Templete.push_back(Rhombus_Points[0]);
            Rhombus_Templete.push_back(Rhombus_Points[1]);
            Rhombus_Templete.push_back(Rhombus_Points[2]);
            Rhombus_Templete.push_back(Rhombus_Points[3]);

            fEnvelopeSolid = new G4ExtrudedSolid("CMOParamEnvelope", Rhombus_Templete, fTotH / 2.,
                                                 0, fEnvelL, 0, fEnvelL);
            break;
        }
        case kBasehwMethod: {
            std::vector<G4TwoVector> Hexagon_Templete;
            G4TwoVector Hexagon_Points[6] = {G4TwoVector(1, 0), G4TwoVector(1, 0),
                                             G4TwoVector(1, 0), G4TwoVector(1, 0),
                                             G4TwoVector(1, 0), G4TwoVector(1, 0)};

            Hexagon_Points[0].rotate(0 * deg);
            Hexagon_Points[1].rotate(60 * deg);
            Hexagon_Points[2].rotate(120 * deg);
            Hexagon_Points[3].rotate(180 * deg);
            Hexagon_Points[4].rotate(240 * deg);
            Hexagon_Points[5].rotate(300 * deg); // make the hexagon points

            Hexagon_Templete.push_back(Hexagon_Points[0]);
            Hexagon_Templete.push_back(Hexagon_Points[1]);
            Hexagon_Templete.push_back(Hexagon_Points[2]);
            Hexagon_Templete.push_back(Hexagon_Points[3]);
            Hexagon_Templete.push_back(Hexagon_Points[4]);
            Hexagon_Templete.push_back(Hexagon_Points[5]); // make the hexagon

            fEnvelopeSolid = new G4ExtrudedSolid("CMOParamEnvelope", Hexagon_Templete, fTotH / 2.,
                                                 0, fEnvelL, 0, fEnvelL);
            break;
        }
    }
    return false;
}

G4bool AmoreCMOParameterisation::CheckSettings() {
    G4bool retval = false;

    if (fFN < 1) {
        G4Exception(
            "AmoreCMOParameterisation::AmoreCMOParameterisation()", "Argument is wrong",
            JustWarning,
            "aFloorNum must be above than 0, AmoreCMOParameterisation will not to be initialized");
        retval = true;
    }
    if (fLN < 0) {
        G4Exception("AmoreCMOParameterisation::AmoreCMOParameterisation()", "Argument is wrong",
                    JustWarning,
                    "aLayerNum must be a positive number or 0, AmoreCMOParameterisation will not "
                    "to be initialized");
        retval = true;
    }
    if (fCSH < 0 || fCSR < 0 || fRS < 0 || fHS < 0) {
        G4Exception("AmoreCMOParameterisation::AmoreCMOParameterisation()", "Argument is wrong",
                    JustWarning,
                    "Cellsize and Spacing must be a positive number, AmoreCMOParameterisation will "
                    "not to be initialized");
        retval = true;
    }
    switch (fM) {
        case kChuhMethod:
        case kBasehwMethod:
            break;
        default:
            G4Exception("AmoreCMOParameterisation::AmoreCMOParameterisation()", "Argument is wrong",
                        JustWarning,
                        "Placing Method must be the one of the value of eTranslationMethod, "
                        "AmoreCMOParameterisation will not to be initialized");
            retval = true;
            break;
    }

    return retval;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool AmoreCMOParameterisation::GenerateTranslationTable() {
    if (fTranslationTable == nullptr) {
        G4Exception(
            "AmoreCMOParameterisation::GenerateTranslationTable()", "Invalid settings", JustWarning,
            "fTranslationTable is nullptr!, AmoreCMOParameterisation will not to be initialized!!");
        return true;
    } else {
        switch (fM) {
            case kChuhMethod: {
                G4ThreeVector Displacer(std::sqrt(2) * (2. * fCSR + fRS));
                G4ThreeVector RaisingFloor(0, 0, fCSH + fHS);
                Displacer.rotateZ(45 * deg);
                for (int i = 0; i <= fLN; i++) {
                    for (int j = 0; j < i * 4 + (i == 0); j++) {
                        G4int nowindex = 2 * i * (i + 1) - 4 * i + j + 1 - (i == 0);

                        if (j == 0)
                            fTranslationTable[nowindex] =
                                G4ThreeVector((2 * fCSR + fRS) * i, 0, -fTotH / 2. + fCSH / 2.);
                        else
                            fTranslationTable[nowindex] =
                                fTranslationTable[nowindex - 1] + Displacer;

                        if (i != 0)
                            if (j % i == 0) Displacer.rotateZ(90 * deg);

                        for (int k = 1; k < fFN; k++)
                            fTranslationTable[nowindex + k * fFloorCellNum] =
                                fTranslationTable[nowindex + (k - 1) * fFloorCellNum] +
                                RaisingFloor;
                    }
                }
                break;
            }

            case kBasehwMethod: {
                G4ThreeVector Displacer(2 * fCSR + fRS);
                G4ThreeVector RaisingFloor(0, 0, fCSH + fHS);
                Displacer.rotateZ(60 * deg);

                for (int i = 0; i <= fLN; i++) {
                    for (int j = 0; j < i * 6 + (i == 0); j++) {
                        G4int nowindex = 3 * i * (i + 1) - 6 * i + j + 1 - (i == 0);

                        if (j == 0)
                            fTranslationTable[nowindex] =
                                G4ThreeVector((2 * fCSR + fRS) * i, 0, -fTotH / 2. + fCSH / 2.);
                        else
                            fTranslationTable[nowindex] =
                                fTranslationTable[nowindex - 1] + Displacer;

                        if (i != 0)
                            if (j % i == 0) Displacer.rotateZ(60 * deg);

                        for (int k = 1; k < fFN; k++)
                            fTranslationTable[nowindex + k * fFloorCellNum] =
                                fTranslationTable[nowindex + (k - 1) * fFloorCellNum] +
                                RaisingFloor;
                    }
                }
                break;
            }
        }
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AmoreCMOParameterisation::ComputeTransformation(const G4int copyNo,
                                                     G4VPhysicalVolume *physVol) const {
    if (!fInitialized) {
        G4Exception("AmoreCMOParameterisation::ComputeTransformation()", "Not initialized",
                    RunMustBeAborted, "AmoreCMOParameterisation is not initialized!!");
        return;
    }

    physVol->SetTranslation(fTranslationTable[copyNo]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
