#include "globals.hh"

#include "AmoreSim/AmoreDetectorConstruction.hh" // the DetectorConstruction class header
#include "AmoreSim/AmoreDetectorStaticInfo.hh"
#include "AmoreSim/AmoreModuleHit.hh"
#include "AmoreSim/AmoreModuleSD.hh"
#include "CupSim/CupPMTOpticalModel.hh"    // for same PMT optical model as main sim
#include "CupSim/CupPMTSD.hh"              // for making sensitive photocathodes
#include "CupSim/CupVetoSD.hh"             // for making sensitive photocathodes

#include "G4Element.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"

#include "G4Colour.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Version.hh"
#include "G4VisAttributes.hh"

#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <stdlib.h>

#include "G4Region.hh"
#include "G4Types.hh"

#include "G4NistManager.hh"

using namespace CLHEP;

using CrystalModuleInfo           = AmoreDetectorStaticInfo::CrystalModuleInfo;
//constexpr G4double toBeDetermined = AmoreDetectorStaticInfo::toBeDetermined;

// Crystal size from Hanbeom at Jan.2021.
const CrystalModuleInfo
    AmoreDetectorStaticInfo::AMoRE_I::crystalModuleInfoList[totalNumOfTower][maxModuleNumInTower] =
			{{{"SE09", "Copper", "CaMoO4", "Vm2000", "Gold",5.046 * cm, 4.8 * cm, 4.0 * cm, 354.2 * g},
				{"SE07", "Copper", "CaMoO4", "Vm2000", "Gold",5.05 * cm, 4.84 * cm, 4.1 * cm, 357.6 * g},
				{"SE03", "Copper", "CaMoO4", "Vm2000", "Silver",5.048 * cm, 5.2 * cm, 4.5 * cm, 426.0 * g},
				{"LMO5", "Copper", "Li2MoO4", "Vm2000", "Gold",6.365 * cm, 5.0 * cm, 5.0 * cm, 377.5 * g},
        {"EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY",-1, -1, -1, -1},
				{"EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY",-1, -1, -1,-1 }},
			 {{"SE06", "Copper", "CaMoO4", "Vm2000", "Gold",5.027 * cm, 4.8 * cm, 4.0 * cm, 355.5 * g},
        {"SE05", "Copper", "CaMoO4", "Vm2000", "Gold",5.041 * cm, 4.9 * cm, 4.2 * cm, 373.2 * g},
        {"SE04", "Copper", "CaMoO4", "Vm2000", "Gold",5.022 * cm, 5.4 * cm, 4.7 * cm, 473.2 * g},
        {"SB28", "Copper", "CaMoO4", "Vm2000", "Gold",2.6 * cm, 5.0 * cm, 4.30 * cm, 196.0 * g},
        {"LMO1", "Copper", "Li2MoO4", "Teflon", "Gold",4.917 * cm, 5.0 * cm, 5.0 * cm, 300.1 * g},
				{"EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY",-1, -1, -1, -1}},
			 {{"SS68", "Copper", "CaMoO4", "Vm2000", "Gold",4.0 * cm, 5.3 * cm, 4.12 * cm, 352.0 * g},
        {"S35", "Copper", "CaMoO4", "Teflon", "Gold",4.0 * cm, 4.4 * cm, 4.4 * cm, 256.0 * g},
        {"LMO3", "Copper", "Li2MoO4", "Vm2000", "Gold",5.053 * cm, 5.0 * cm, 5.0 * cm, 311.5 * g},
        {"LMO4", "Copper", "Li2MoO4", "Vm2000", "Gold",5.014 * cm, 5.0 * cm, 5.0 * cm, 308.1 * g},
        {"LMO2", "Copper", "Li2MoO4", "Vm2000", "Silver",5.017 * cm, 5.0 * cm, 5.0 * cm, 312.0 * g},
				{"EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY",-1, -1, -1, -1}},
			 {{"SE08", "Copper", "CaMoO4", "Vm2000", "Gold",5.08 * cm, 4.8 * cm, 4.1 * cm, 358.2 * g},
        {"SE02", "Copper", "CaMoO4", "Vm2000", "Gold",5.07 * cm, 4.76 * cm, 3.98 * cm, 340.0 * g},
        {"SE01", "Copper", "CaMoO4", "Vm2000", "Gold",5.10 * cm, 4.4 * cm, 4.4 * cm, 353.0 * g},
        {"SB29", "Copper", "CaMoO4", "Vm2000", "Gold",5.1 * cm, 5.0 * cm, 4.7 * cm, 390.0 * g},
        {"EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY",-1, -1, -1, -1},
        {"EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY",-1, -1, -1, -1}}};


/**
 * @author Mona
 * @date 2019/2/27
 * @brief SUPER CONDUCTING MAGNET SHIELDING
 * @param aWhere The location where the SMS will be located
 * @param aTlate A translation of SMS shield
 * @param aType Type of SMS shield (1: w/ lead shield, 2: w/o lead shield)
 * @param aFrameMaterial Material for copper frame
 * @param aSpaceMaterial Material of vacant space (usually set to be vacuum)
 * @param aShieldMaterial Material for shield at the surface of SMS geometry
 * */
G4LogicalVolume *AmoreDetectorConstruction::Build_I_SCMagnetShieldAt(
    G4LogicalVolume *aWhere, G4ThreeVector aTlate, G4int aType, G4Material *aFrameMaterial,
    G4Material *aSpaceMaterial, G4Material *aShieldMaterial) {
    G4LogicalVolume *resultLV;
    G4ThreeVector baseTlate = aTlate;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4double scShieldInnerSizeX;
    G4double scShieldInnerSizeY;
    G4double scShieldInnerSizeZ;

    G4double scsCuFrameBottomCutCoef;
    G4double scsCuFrameBottomCutSlope;

    G4double scsCuFrameSkeletonHori1Length;
    G4double scsCuFrameSkeletonHori2Length;
    G4double scsCuFrameSkeletonDiagLength;
    G4double scsCuFrameSkeletonVertLength;

    G4bool isSCSCuFrameCutted =
        !(scsCuFrameBottomCutAngle == 0. * deg || scsCuFrameBottomCutAngle == 90. * deg);

    G4VSolid *scShieldSolid;
    G4VSolid *scShieldInnerSolid;
    G4VSolid *scsCuFrameBottomSolid;
    G4Tubs *scsCuFrameBottomHole;
    G4Tubs *scsCuFrameTopDisk;

    G4Box *scsCuFrameTopBlock1;
    G4Box *scsCuFrameTopBlock2;
    G4Box *scsCuFrameTopLongBlock1;
    G4Box *scsCuFrameTopLongBlock2;
    G4Tubs *scsCuFrameHoleAtTopBlock;

    G4Box *scsCuFrameSkeletonHori1Box;
    G4Box *scsCuFrameSkeletonHori2Box;
    G4Box *scsCuFrameSkeletonDiagBox;
    G4Box *scsCuFrameSkeletonVert1Box;
    G4Box *scsCuFrameSkeletonVert2Box;

    G4LogicalVolume *scShieldLV = nullptr;
    G4LogicalVolume *scShieldInnerLV;
    G4LogicalVolume *scsCuFrameBottomLV;
    G4LogicalVolume *scsCuFrameBottomHoleLV;
    G4LogicalVolume *scsCuFrameTopDiskLV;

    G4LogicalVolume *scsCuFrameTopBlock1LV;
    G4LogicalVolume *scsCuFrameTopBlock2LV;
    G4LogicalVolume *scsCuFrameTopLongBlock1LV;
    G4LogicalVolume *scsCuFrameTopLongBlock2LV;
    G4LogicalVolume *scsCuFrameHoleAtTopBlockLV;

    G4LogicalVolume *scsCuFrameSkeletonHori1LV;
    G4LogicalVolume *scsCuFrameSkeletonHori2LV;
    G4LogicalVolume *scsCuFrameSkeletonDiagLV;
    G4LogicalVolume *scsCuFrameSkeletonVert1LV;
    G4LogicalVolume *scsCuFrameSkeletonVert2LV;

    G4ThreeVector scsCuFrameBottomHoleNowPos;
    G4ThreeVector scsCuFrameTopDiskPos;
    G4ThreeVector scsCuFrameTopBlocksNowPos;
    G4ThreeVector scsCuFrameSkeletonNowPos;

    G4RotationMatrix *scsCuFrameSkeletonRotMtx[2];

    std::vector<G4TwoVector> scShieldInnerVertices;
    std::vector<G4TwoVector> scsCuFrameBottomVertices;

    auto get2DDistanceOfTwoPosition = [](const G4TwoVector &a, const G4TwoVector &b) -> G4double {
        return (b - a).mag();
    };

    switch (aType) {
        case 1:
        case 2: {
            if (isSCSCuFrameCutted) {
                scsCuFrameBottomVertices.resize(8);
                scsCuFrameBottomCutSlope = tan(scsCuFrameBottomCutAngle + 90. * deg);
                scsCuFrameBottomCutCoef =
                    scsCuFrameBottomCutPosY - scsCuFrameBottomCutSlope * scsCuFrameBottomCutPosX;

                scsCuFrameBottomVertices[0]     = scsCuFrameBottomVertices[3] =
                    scsCuFrameBottomVertices[4] = scsCuFrameBottomVertices[7] = G4TwoVector(
                        scShieldSizeX / 2.,
                        (scShieldSizeX / 2.) * scsCuFrameBottomCutSlope + scsCuFrameBottomCutCoef);
                scsCuFrameBottomVertices[1]     = scsCuFrameBottomVertices[2] =
                    scsCuFrameBottomVertices[5] = scsCuFrameBottomVertices[6] = G4TwoVector(
                        (scShieldSizeY / 2. - scsCuFrameBottomCutCoef) / scsCuFrameBottomCutSlope,
                        scShieldSizeY / 2.);

                scsCuFrameBottomVertices[2][0] *= -1;
                scsCuFrameBottomVertices[3][0] *= -1;
                scsCuFrameBottomVertices[4][0] *= -1;
                scsCuFrameBottomVertices[4][1] *= -1;
                scsCuFrameBottomVertices[5][0] *= -1;
                scsCuFrameBottomVertices[5][1] *= -1;
                scsCuFrameBottomVertices[6][1] *= -1;
                scsCuFrameBottomVertices[7][1] *= -1;

                scShieldInnerVertices = scsCuFrameBottomVertices;

                scShieldInnerVertices[0][0] -= scShieldLeadThickness;
                scShieldInnerVertices[3][0] += scShieldLeadThickness;
                scShieldInnerVertices[4][0] += scShieldLeadThickness;
                scShieldInnerVertices[7][0] -= scShieldLeadThickness;

                scShieldInnerVertices[1][1] -= scShieldLeadThickness;
                scShieldInnerVertices[2][1] -= scShieldLeadThickness;
                scShieldInnerVertices[5][1] += scShieldLeadThickness;
                scShieldInnerVertices[6][1] += scShieldLeadThickness;

                scShieldInnerSizeZ =
                    scShieldSizeZ - scsCuFrameBottomThickness - scShieldLeadThickness;

                scShieldSolid = new G4ExtrudedSolid("SCMagnetOuterSolid", scsCuFrameBottomVertices,
                                                    scShieldSizeZ / 2., 0, 1., 0, 1.);
                scShieldInnerSolid =
                    new G4ExtrudedSolid("SCMagnetInnerSolid", scShieldInnerVertices,
                                        scShieldInnerSizeZ / 2., 0, 1., 0, 1.);

                scsCuFrameBottomSolid =
                    new G4ExtrudedSolid("SCMagnetCuFrameBottomSolid", scsCuFrameBottomVertices,
                                        scsCuFrameBottomThickness / 2., 0, 1., 0, 1.);
            } else {
                scShieldSolid = new G4Box("SCMagnetOuterSolid", scShieldSizeX / 2.,
                                          scShieldSizeY / 2., scShieldSizeZ / 2.);

                scShieldInnerSizeX = scShieldSizeX - scShieldLeadThickness * 2.;
                scShieldInnerSizeY = scShieldSizeY - scShieldLeadThickness * 2.;
                scShieldInnerSizeZ =
                    scShieldSizeZ - scsCuFrameBottomThickness - scShieldLeadThickness;

                scShieldInnerSolid = new G4Box("SCMagnetInnerSolid", scShieldInnerSizeX / 2.,
                                               scShieldInnerSizeY / 2., scShieldInnerSizeZ / 2.);

                scsCuFrameBottomSolid =
                    new G4Box("SCMagnetCuFrameBottomSolid", scShieldSizeX / 2., scShieldSizeY / 2.,
                              scsCuFrameBottomThickness / 2.);
            }

            if (aType == 1) {
                scShieldLV =
                    new G4LogicalVolume(scShieldSolid, aShieldMaterial, "SCMagnetOuterShield_LV");
                scShieldLV->SetVisAttributes(fI_LeadShield_VisAttr);
            } else {
                scShieldLV =
                    new G4LogicalVolume(scShieldSolid, aSpaceMaterial, "SCMagnetOuterShield_LV");
                scShieldLV->SetVisAttributes(fI_Invisible_VisAttr);
            }

            scShieldInnerLV =
                new G4LogicalVolume(scShieldInnerSolid, aSpaceMaterial, "SCMagnetInnerSpace_LV");
            scShieldInnerLV->SetVisAttributes(fI_Invisible_VisAttr);

            scsCuFrameBottomLV = new G4LogicalVolume(scsCuFrameBottomSolid, aFrameMaterial,
                                                     "SCMagnetBottomShield_LV");
            scsCuFrameBottomLV->SetVisAttributes(fI_CopperDefault_VisAttr);

            scsCuFrameBottomHole = new G4Tubs("SCMagnetBottomHole", 0, scsCuFrameBottomHoleRadius,
                                              scsCuFrameBottomThickness / 2., 0, 360. * deg);
            scsCuFrameBottomHoleLV =
                new G4LogicalVolume(scsCuFrameBottomHole, aSpaceMaterial, "SCShieldBottomHole_LV");

            new G4PVPlacement(nullptr, {0, 0, 0}, scsCuFrameBottomHoleLV, "SCMagnetCuBottomHole_PV",
                              scsCuFrameBottomLV, false, 0, OverlapCheck);

            scsCuFrameBottomHoleNowPos =
                G4ThreeVector(scsCuFrameBottomHoleDistX, scsCuFrameBottomHoleDistY, 0);
            new G4PVPlacement(nullptr, scsCuFrameBottomHoleNowPos, scsCuFrameBottomHoleLV,
                              "SCMagnetCuBottomHole_PV", scsCuFrameBottomLV, false, 1, OverlapCheck);
            scsCuFrameBottomHoleNowPos[1] *= -1.;
            new G4PVPlacement(nullptr, scsCuFrameBottomHoleNowPos, scsCuFrameBottomHoleLV,
                              "SCMagnetCuBottomHole_PV", scsCuFrameBottomLV, false, 2, OverlapCheck);
            scsCuFrameBottomHoleNowPos[0] *= -1.;
            new G4PVPlacement(nullptr, scsCuFrameBottomHoleNowPos, scsCuFrameBottomHoleLV,
                              "SCMagnetCuBottomHole_PV", scsCuFrameBottomLV, false, 3, OverlapCheck);
            scsCuFrameBottomHoleNowPos[1] *= -1.;
            new G4PVPlacement(nullptr, scsCuFrameBottomHoleNowPos, scsCuFrameBottomHoleLV,
                              "SCMagnetCuBottomHole_PV", scsCuFrameBottomLV, false, 4, OverlapCheck);

            new G4PVPlacement(nullptr, {0, 0, -scShieldSizeZ / 2. + scsCuFrameBottomThickness / 2.},
                              scsCuFrameBottomLV, "SCMagnetBottomFrame_PV", scShieldLV, false, 0, OverlapCheck);

            new G4PVPlacement(nullptr, {0, 0, scShieldSizeZ / 2. - scShieldInnerSizeZ / 2.},
                              scShieldInnerLV, "SCMagnetInnerSpace_PV", scShieldLV, false, 0, OverlapCheck);

            new G4PVPlacement(nullptr, baseTlate, scShieldLV, "SCMagnetOuterShield_PV", aWhere,
                              false, 0, OverlapCheck);
            resultLV = scShieldInnerLV;

            scsCuFrameTopDisk = new G4Tubs("SCMagnetTopFrameDisk", 0, scsCuFrameTopDiskRadius,
                                           scsCuFrameTopDiskThickness / 2., 0, 360 * deg);
            scsCuFrameTopDiskLV =
                new G4LogicalVolume(scsCuFrameTopDisk, aFrameMaterial, "SCMagnetTopFrameDisk_LV");
            scsCuFrameTopDiskLV->SetVisAttributes(fI_CopperDefault_VisAttr);

            scsCuFrameTopDiskPos = {0, 0, scShieldSizeZ / 2. + scsCuFrameTopDiskThickness / 2.};
            new G4PVPlacement(nullptr, scsCuFrameTopDiskPos + baseTlate, scsCuFrameTopDiskLV,
                              "SCMagnetTopFrame_LV", aWhere, false, 0, OverlapCheck);

            scsCuFrameTopBlock1 =
                new G4Box("SCMagnetFrameTopBlock1_Box", scsCuFrameTopBlock1SizeX / 2.,
                          scsCuFrameTopBlock1SizeY / 2., scsCuFrameTopBlocksThickness / 2.);
            scsCuFrameTopBlock2 =
                new G4Box("SCMagnetFrameTopBlock2_Box", scsCuFrameTopBlock2SizeX / 2.,
                          scsCuFrameTopBlock2SizeY / 2., scsCuFrameTopBlocksThickness / 2.);
            scsCuFrameTopBlock1LV = new G4LogicalVolume(scsCuFrameTopBlock1, aFrameMaterial,
                                                        "SCMagnetFrameTopBlock1_LV");
            scsCuFrameTopBlock1LV->SetVisAttributes(fI_CopperDefault_VisAttr);
            scsCuFrameTopBlock2LV = new G4LogicalVolume(scsCuFrameTopBlock2, aFrameMaterial,
                                                        "SCMagnetFrameTopBlock2_LV");
            scsCuFrameTopBlock2LV->SetVisAttributes(fI_CopperDefault_VisAttr);

            scsCuFrameTopLongBlock1 =
                new G4Box("SCMagnetFrameTopLongBlock1_Box", scsCuFrameTopLongBlock1SizeX / 2.,
                          scsCuFrameTopLongBlock1SizeY / 2., scsCuFrameTopBlocksThickness / 2.);
            scsCuFrameTopLongBlock2 =
                new G4Box("SCMagnetFrameTopLongBlock2_Box", scsCuFrameTopLongBlock2SizeX / 2.,
                          scsCuFrameTopLongBlock2SizeY / 2., scsCuFrameTopBlocksThickness / 2.);
            scsCuFrameTopLongBlock1LV = new G4LogicalVolume(scsCuFrameTopLongBlock1, aFrameMaterial,
                                                            "SCMagnetFrameTopLongBlock1_LV");
            scsCuFrameTopLongBlock1LV->SetVisAttributes(fI_CopperDefault_VisAttr);
            scsCuFrameTopLongBlock2LV = new G4LogicalVolume(scsCuFrameTopLongBlock2, aFrameMaterial,
                                                            "SCMagnetFrameTopLongBlock2_LV");
            scsCuFrameTopLongBlock2LV->SetVisAttributes(fI_CopperDefault_VisAttr);

            scsCuFrameHoleAtTopBlock =
                new G4Tubs("SCMagnetFrameHoleAtTopBlock_Tub", 0, scsCuFrameHoleAtTopBlockRadius,
                           scsCuFrameHoleAtTopBlockLength / 2., 0, 360. * deg);
            scsCuFrameHoleAtTopBlockLV = new G4LogicalVolume(
                scsCuFrameHoleAtTopBlock, aFrameMaterial, "SCMagnetFrameHoleAtTopBlock_LV");

            new G4PVPlacement(nullptr, {-2. * scsCuFrameHoleAtTopBlockDist, 0, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock2_PV",
                              scsCuFrameTopBlock2LV, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, {-1. * scsCuFrameHoleAtTopBlockDist, 0, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock2_PV",
                              scsCuFrameTopBlock2LV, false, 1, OverlapCheck);
            new G4PVPlacement(nullptr, {0, 0, 0}, scsCuFrameHoleAtTopBlockLV,
                              "SCMagnetFrameHoleAtTopBlock2_PV", scsCuFrameTopBlock2LV, false, 2, OverlapCheck);
            new G4PVPlacement(nullptr, {1. * scsCuFrameHoleAtTopBlockDist, 0, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock2_PV",
                              scsCuFrameTopBlock2LV, false, 3, OverlapCheck);
            new G4PVPlacement(nullptr, {2. * scsCuFrameHoleAtTopBlockDist, 0, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock2_PV",
                              scsCuFrameTopBlock2LV, false, 4, OverlapCheck);

            new G4PVPlacement(nullptr, {0, -2. * scsCuFrameHoleAtTopBlockDist, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock1_PV",
                              scsCuFrameTopBlock1LV, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, {0, -1. * scsCuFrameHoleAtTopBlockDist, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock1_PV",
                              scsCuFrameTopBlock1LV, false, 1, OverlapCheck);
            new G4PVPlacement(nullptr, {0, 0, 0}, scsCuFrameHoleAtTopBlockLV,
                              "SCMagnetFrameHoleAtTopBlock1_PV", scsCuFrameTopBlock1LV, false, 2, OverlapCheck);
            new G4PVPlacement(nullptr, {0, 1. * scsCuFrameHoleAtTopBlockDist, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock1_PV",
                              scsCuFrameTopBlock1LV, false, 3, OverlapCheck);
            new G4PVPlacement(nullptr, {0, 2. * scsCuFrameHoleAtTopBlockDist, 0},
                              scsCuFrameHoleAtTopBlockLV, "SCMagnetFrameHoleAtTopBlock1_PV",
                              scsCuFrameTopBlock1LV, false, 4, OverlapCheck);

            scsCuFrameTopBlocksNowPos = {-scShieldSizeX / 2. - scsCuFrameTopLongBlock1SizeX / 2., 0,
                                         scShieldSizeZ / 2. - scsCuFrameTopBlocksThickness / 2.};
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos,
                              scsCuFrameTopLongBlock1LV, "SCMagnetFrameTopLongBlock1_PV", aWhere,
                              false, 0, OverlapCheck);
            scsCuFrameTopBlocksNowPos[0] *= -1;
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos,
                              scsCuFrameTopLongBlock1LV, "SCMagnetFrameTopLongBlock1_PV", aWhere,
                              false, 1, OverlapCheck);

            scsCuFrameTopBlocksNowPos = {0, -scShieldSizeY / 2. - scsCuFrameTopLongBlock2SizeY / 2.,
                                         scShieldSizeZ / 2. - scsCuFrameTopBlocksThickness / 2.};
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos,
                              scsCuFrameTopLongBlock2LV, "SCMagnetFrameTopLongBlock2_PV", aWhere,
                              false, 0, OverlapCheck);
            scsCuFrameTopBlocksNowPos[1] *= -1;
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos,
                              scsCuFrameTopLongBlock2LV, "SCMagnetFrameTopLongBlock2_PV", aWhere,
                              false, 1, OverlapCheck);

            scsCuFrameTopBlocksNowPos = {-scShieldSizeX / 2. - scsCuFrameTopLongBlock1SizeX -
                                             scsCuFrameTopBlock1SizeX / 2.,
                                         0, scShieldSizeZ / 2. - scsCuFrameTopBlocksThickness / 2.};
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos, scsCuFrameTopBlock1LV,
                              "SCMagnetTopBlock1_PV", aWhere, false, 0, OverlapCheck);
            scsCuFrameTopBlocksNowPos[0] *= -1.;
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos, scsCuFrameTopBlock1LV,
                              "SCMagnetTopBlock1_PV", aWhere, false, 1, OverlapCheck);

            scsCuFrameTopBlocksNowPos = {0,
                                         -scShieldSizeY / 2. - scsCuFrameTopLongBlock2SizeY -
                                             scsCuFrameTopBlock2SizeY / 2.,
                                         scShieldSizeZ / 2. - scsCuFrameTopBlocksThickness / 2.};
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos, scsCuFrameTopBlock2LV,
                              "SCMagnetTopBlock2_PV", aWhere, false, 0, OverlapCheck);
            scsCuFrameTopBlocksNowPos[1] *= -1.;
            new G4PVPlacement(nullptr, baseTlate + scsCuFrameTopBlocksNowPos, scsCuFrameTopBlock2LV,
                              "SCMagnetTopBlock2_PV", aWhere, false, 1, OverlapCheck);

            if (isSCSCuFrameCutted) {
                scsCuFrameSkeletonDiagLength = get2DDistanceOfTwoPosition(
                    scsCuFrameBottomVertices[0], scsCuFrameBottomVertices[1]);
                scsCuFrameSkeletonHori1Length = get2DDistanceOfTwoPosition(
                    scsCuFrameBottomVertices[3], scsCuFrameBottomVertices[4]);
                scsCuFrameSkeletonHori2Length = get2DDistanceOfTwoPosition(
                    scsCuFrameBottomVertices[1], scsCuFrameBottomVertices[2]);

                scsCuFrameSkeletonDiagBox =
                    new G4Box("SCMagnetFrameSkeletonDiag_Box", scsCuFrameSkeletonDiagLength / 2.,
                              scsCuFrameSkeletonThickness / 2., scsCuFrameSkeletonWidth / 2.);

                scsCuFrameSkeletonDiagLV = new G4LogicalVolume(
                    scsCuFrameSkeletonDiagBox, aFrameMaterial, "SCMagnetFrameSkeletonDiag_LV");
                scsCuFrameSkeletonDiagLV->SetVisAttributes(fI_CopperDefault_VisAttr);
                G4TwoVector scsCuFrameSkeletonDiagDisplacer(1, 0);
                scsCuFrameSkeletonDiagDisplacer.rotate(scsCuFrameBottomCutAngle);
                scsCuFrameSkeletonDiagDisplacer *= scsCuFrameSkeletonThickness / 2.;
                scsCuFrameSkeletonDiagDisplacer +=
                    ((scsCuFrameBottomVertices[0] + scsCuFrameBottomVertices[1]) / 2.);

                scsCuFrameSkeletonNowPos = {scsCuFrameSkeletonDiagDisplacer[0],
                                            scsCuFrameSkeletonDiagDisplacer[1], 0};

                scsCuFrameSkeletonRotMtx[0] = new G4RotationMatrix;
                scsCuFrameSkeletonRotMtx[0]->rotateZ(-scsCuFrameBottomCutAngle + 90 * deg);
                scsCuFrameSkeletonRotMtx[1] = new G4RotationMatrix;
                scsCuFrameSkeletonRotMtx[1]->rotateZ(scsCuFrameBottomCutAngle - 90 * deg);

                for (G4int i = 0; i < 2; i++) {
                    scsCuFrameSkeletonNowPos[2] = scShieldSizeZ / 2. - scsCuFrameSkeletonWidth / 2.;
                    for (G4int j = 0; j < 4; j++) {
                        new G4PVPlacement(scsCuFrameSkeletonRotMtx[i],
                                          baseTlate + scsCuFrameSkeletonNowPos,
                                          scsCuFrameSkeletonDiagLV, "SCMagnetFrameSkeletonDiag_PV",
                                          aWhere, false, i * 4 + j, OverlapCheck);
                        scsCuFrameSkeletonNowPos[2] -= scsCuFrameSkeletonDistance;
                    }
                    scsCuFrameSkeletonNowPos[0] *= -1;
                    scsCuFrameSkeletonNowPos[1] *= -1;
                    scsCuFrameSkeletonNowPos[2] = scShieldSizeZ / 2. - scsCuFrameSkeletonWidth / 2.;
                    for (G4int j = 0; j < 4; j++) {
                        new G4PVPlacement(scsCuFrameSkeletonRotMtx[i],
                                          baseTlate + scsCuFrameSkeletonNowPos,
                                          scsCuFrameSkeletonDiagLV, "SCMagnetFrameSkeletonDiag_PV",
                                          aWhere, false, i * 4 + j + 8, OverlapCheck);
                        scsCuFrameSkeletonNowPos[2] -= scsCuFrameSkeletonDistance;
                    }
                    scsCuFrameSkeletonNowPos[1] *= -1;
                }
            } else {
                scsCuFrameSkeletonHori1Length = scShieldSizeY;
                scsCuFrameSkeletonHori2Length = scShieldSizeX;
            }
            scsCuFrameSkeletonVertLength = scShieldSizeZ - scsCuFrameTopBlocksThickness;

            scsCuFrameSkeletonHori1Box =
                new G4Box("SCMagnetFrameSkeletonHori1_Box", scsCuFrameSkeletonThickness / 2.,
                          scsCuFrameSkeletonHori1Length / 2., scsCuFrameSkeletonWidth / 2.);
            scsCuFrameSkeletonHori2Box =
                new G4Box("SCMagnetFrameSkeletonHori2_Box", scsCuFrameSkeletonHori2Length / 2.,
                          scsCuFrameSkeletonThickness / 2., scsCuFrameSkeletonWidth / 2.);
            scsCuFrameSkeletonVert1Box =
                new G4Box("SCMagnetFrameSkeletonVert1_Box", scsCuFrameSkeletonThickness / 2.,
                          scsCuFrameSkeletonWidth / 2., scsCuFrameSkeletonVertLength / 2.);
            scsCuFrameSkeletonVert2Box =
                new G4Box("SCMagnetFrameSkeletonVert2_Box", scsCuFrameSkeletonWidth / 2.,
                          scsCuFrameSkeletonThickness / 2., scsCuFrameSkeletonVertLength / 2.);

            scsCuFrameSkeletonHori1LV = new G4LogicalVolume(
                scsCuFrameSkeletonHori1Box, aFrameMaterial, "SCMagnetFrameSkeletonHori1_LV");
            scsCuFrameSkeletonHori2LV = new G4LogicalVolume(
                scsCuFrameSkeletonHori2Box, aFrameMaterial, "SCMagnetFrameSkeletonHori2_LV");
            scsCuFrameSkeletonVert1LV = new G4LogicalVolume(
                scsCuFrameSkeletonVert1Box, aFrameMaterial, "SCMagnetFrameSkeletonVert1_LV");
            scsCuFrameSkeletonVert2LV = new G4LogicalVolume(
                scsCuFrameSkeletonVert2Box, aFrameMaterial, "SCMagnetFrameSkeletonVert2_LV");
            scsCuFrameSkeletonHori1LV->SetVisAttributes(fI_CopperDefault_VisAttr);
            scsCuFrameSkeletonHori2LV->SetVisAttributes(fI_CopperDefault_VisAttr);
            scsCuFrameSkeletonVert1LV->SetVisAttributes(fI_CopperDefault_VisAttr);
            scsCuFrameSkeletonVert2LV->SetVisAttributes(fI_CopperDefault_VisAttr);

            scsCuFrameSkeletonNowPos = {-scShieldSizeX / 2. - scsCuFrameSkeletonThickness / 2., 0,
                                        scShieldSizeZ / 2. - scsCuFrameSkeletonWidth / 2.};
            for (G4int i = 0; i < 3; i++) {
                scsCuFrameSkeletonNowPos[2] -= scsCuFrameSkeletonDistance;
                for (G4int j = 0; j < 2; j++) {
                    new G4PVPlacement(nullptr, baseTlate + scsCuFrameSkeletonNowPos,
                                      scsCuFrameSkeletonHori1LV, "SCMagnetFrameSkeletonHori1_PV",
                                      aWhere, false, i * 2 + j, OverlapCheck);
                    scsCuFrameSkeletonNowPos[0] *= -1.;
                }
            }

            scsCuFrameSkeletonNowPos = {
                -scShieldSizeX / 2. - scsCuFrameSkeletonThickness * 3 / 2.,
                scsCuFrameSkeletonHori1Length / 2. -
                    ((isSCSCuFrameCutted) ? 0 : (scsCuFrameSkeletonWidth / 2.)),
                scShieldSizeZ / 2. - scsCuFrameTopBlocksThickness -
                    scsCuFrameSkeletonVertLength / 2.};
            for (G4int i = 0; i < 2; i++) {
                for (G4int j = 0; j < 2; j++) {
                    new G4PVPlacement(nullptr, baseTlate + scsCuFrameSkeletonNowPos,
                                      scsCuFrameSkeletonVert1LV, "SCMagnetFrameSkeletonVert1_PV",
                                      aWhere, false, i * 2 + j, OverlapCheck);
                    scsCuFrameSkeletonNowPos[1] *= -1;
                }
                scsCuFrameSkeletonNowPos[0] *= -1;
            }

            scsCuFrameSkeletonNowPos = {0, -scShieldSizeY / 2. - scsCuFrameSkeletonThickness / 2.,
                                        scShieldSizeZ / 2. - scsCuFrameSkeletonWidth / 2.};
            for (G4int i = 0; i < 3; i++) {
                scsCuFrameSkeletonNowPos[2] -= scsCuFrameSkeletonDistance;
                for (G4int j = 0; j < 2; j++) {
                    new G4PVPlacement(nullptr, baseTlate + scsCuFrameSkeletonNowPos,
                                      scsCuFrameSkeletonHori2LV, "SCMagnetFrameSkeletonHori2_PV",
                                      aWhere, false, i * 2 + j, OverlapCheck);
                    scsCuFrameSkeletonNowPos[1] *= -1.;
                }
            }

            scsCuFrameSkeletonNowPos = {
                scsCuFrameSkeletonHori2Length / 2. -
                    ((isSCSCuFrameCutted) ? 0 : (scsCuFrameSkeletonWidth / 2.)),
                -scShieldSizeY / 2. - scsCuFrameSkeletonThickness * 3 / 2.,
                scShieldSizeZ / 2. - scsCuFrameTopBlocksThickness -
                    scsCuFrameSkeletonVertLength / 2.};
            for (G4int i = 0; i < 2; i++) {
                for (G4int j = 0; j < 2; j++) {
                    new G4PVPlacement(nullptr, baseTlate + scsCuFrameSkeletonNowPos,
                                      scsCuFrameSkeletonVert2LV, "SCMagnetFrameSkeletonVert2_PV",
                                      aWhere, false, i * 2 + j, OverlapCheck);
                    scsCuFrameSkeletonNowPos[1] *= -1;
                }
                scsCuFrameSkeletonNowPos[0] *= -1;
            }
        } break;
        default: {
            resultLV = nullptr;
            G4Exception(__PRETTY_FUNCTION__, "AMORE_I_TYPEERROR",
                        G4ExceptionSeverity::FatalErrorInArgument,
                        "Wrong type for Superconducting Magnet Shield");
            break;
        }
    }

    return resultLV;
}

// modified by Jeewon Seo at 2021. Mar.
// Two option (if aRealistic is turned off, simple box will be installed.)
G4LogicalVolume *AmoreDetectorConstruction::Build_I_DetectorArrayUpperPanel(
    G4bool aRealistic, G4Material *aFrameMaterial, G4Material *aSpaceMaterial) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4bool isCornerCutted    = !(cornerCuttedAngle == 0. * deg || cornerCuttedAngle == 90. * deg);
    G4double cornerCuttedSlope, cornerCuttedCoef;
		G4double subPanelSizeX = type3_innerPartSizeX+type3_panelFrameThick*2;
		G4double subPanelSizeY = type3_innerPartSizeY+type3_panelFrameThick*2;

    G4VSolid *currentSolid;
    G4Box *simpleUpperPanel;
    G4ExtrudedSolid *realisticUpperPanel;

    G4Box *innerCarvingBox;
    G4Box *innerMajorHoleBox;
    G4Box *innerMinorHoleBox1;
    G4Box *innerMinorHoleBox2;
    G4Tubs *edgeHole;

    G4LogicalVolume *edgeHoleLV;
    G4LogicalVolume *innerMajorHoleLV;
    G4LogicalVolume *innerMinorHoleLV1;
    G4LogicalVolume *innerMinorHoleLV2;

    G4ThreeVector holePosition;

		G4Box *subUpperPanel;

    std::vector<G4TwoVector> realisticUpperPanelVertices;

    if (aRealistic) {
        resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
            "Realistic_DetectorArrayUpperPanel_LV", false);
        if (resultLV == nullptr) {
            G4cout << __PRETTY_FUNCTION__
                   << ": Realistic CMO detector array upper panel selected (It could be very slow)"
                   << G4endl;
            if (isCornerCutted) {
                G4cout << __PRETTY_FUNCTION__ << ": Cutting the edges..." << G4endl;
                cornerCuttedSlope = tan(cornerCuttedAngle + 90. * deg);

                realisticUpperPanelVertices.resize(8);
                cornerCuttedCoef = innerPartSizeY / 2. - cornerCuttedSlope * innerPartSizeX / 2.;

                realisticUpperPanelVertices[0]     = realisticUpperPanelVertices[3] =
                    realisticUpperPanelVertices[4] = realisticUpperPanelVertices[7] =
                        G4TwoVector(upperPanelSizeX / 2.,
                                    (upperPanelSizeX / 2.) * cornerCuttedSlope + cornerCuttedCoef);
                realisticUpperPanelVertices[1]     = realisticUpperPanelVertices[2] =
                    realisticUpperPanelVertices[5] = realisticUpperPanelVertices[6] =
                        G4TwoVector((upperPanelSizeY / 2. - cornerCuttedCoef) / cornerCuttedSlope,
                                    upperPanelSizeY / 2.);

                realisticUpperPanelVertices[2][0] *= -1;
                realisticUpperPanelVertices[3][0] *= -1;
                realisticUpperPanelVertices[4][0] *= -1;
                realisticUpperPanelVertices[4][1] *= -1;
                realisticUpperPanelVertices[5][0] *= -1;
                realisticUpperPanelVertices[5][1] *= -1;
                realisticUpperPanelVertices[6][1] *= -1;
                realisticUpperPanelVertices[7][1] *= -1;

                realisticUpperPanel =
                    new G4ExtrudedSolid("CMOArrayUpperPanelSolid", realisticUpperPanelVertices,
                                        type3_panelFrameThick / 2., 0, 1., 0, 1.);

								subUpperPanel = new G4Box("subUpperPanelBox",subPanelSizeX/2., subPanelSizeY/2., type3_subPanelThick/2.);

								currentSolid = new G4UnionSolid("CMOArrayUPSolid",realisticUpperPanel,subUpperPanel,nullptr,
										{0,0,-type3_panelFrameThick/2.-type3_subPanelThick/2.});
            } else {
                simpleUpperPanel = new G4Box("CMOArrayUpperPanelBox", upperPanelSizeX / 2.,
                                             upperPanelSizeY / 2., upperPanelThick / 2.);
                currentSolid     = simpleUpperPanel;
            }

            G4cout << __PRETTY_FUNCTION__ << ": Carving square well at the center..." << G4endl;
            innerCarvingBox = new G4Box("InnerCarvingBox", type3_innerPartSizeX / 2., type3_innerPartSizeY / 2.,
                                        type3_innerPartDepth / 2. + innerCarvingTol / 2.);
            currentSolid    = new G4SubtractionSolid(
                "CMOAUP_Carved_Panel", currentSolid, innerCarvingBox, nullptr,
                {0, 0, type3_panelFrameThick/2.-type3_subPanelThick+type3_innerPartDepth/2. + innerCarvingTol/2.});

            resultLV = new G4LogicalVolume(currentSolid, aFrameMaterial,
                                           "Realistic_DetectorArrayUpperPanel_LV");
            resultLV->SetVisAttributes(fI_CopperDefault_VisAttr);

            G4cout << __PRETTY_FUNCTION__ << ": Punching inner square holes..." << G4endl;
            innerMajorHoleBox  = new G4Box("InnerMajorHoleBox", type3_innerMajorHoleSizeX / 2.,
                                          type3_innerMajorHoleSizeY / 2.,type3_panelFrameThick / 2.);
            innerMinorHoleBox1 = new G4Box("InnerMinorHole1Box", type3_innerMinorHole1SizeX / 2.,
                                           type3_innerMinorHole1SizeY / 2., type3_panelFrameThick / 2.);
            innerMinorHoleBox2 = new G4Box("InnerMinorHole2Box", type3_innerMinorHole2SizeX / 2.,
                                           type3_innerMinorHole2SizeY / 2., type3_panelFrameThick / 2.);

            innerMajorHoleLV =
                new G4LogicalVolume(innerMajorHoleBox, aSpaceMaterial, "InnerMajorHole_LV");
            innerMinorHoleLV1 =
                new G4LogicalVolume(innerMinorHoleBox1, aSpaceMaterial, "InnerMinorHole1_LV");
            innerMinorHoleLV2 =
                new G4LogicalVolume(innerMinorHoleBox2, aSpaceMaterial, "InnerMinorHole2_LV");

            holePosition =
                G4ThreeVector(
                    -type3_innerPartSizeX / 2. + type3_innerMajorHoleSizeX / 2. + type3_innerMajorHoleDistanceX,
                    -type3_innerPartSizeY / 2. + type3_innerMajorHoleSizeY / 2. + type3_innerMajorHoleDistanceY) +
                G4ThreeVector(0, 0, -type3_subPanelThick);
            new G4PVPlacement(nullptr, holePosition, innerMajorHoleLV, "InnerMajorHole_PV",
                              resultLV, false, 0, OverlapCheck);
            holePosition[0] *= -1.;
            new G4PVPlacement(nullptr, holePosition, innerMajorHoleLV, "InnerMajorHole_PV",
                              resultLV, false, 1, OverlapCheck);
            holePosition[0] *= -1.;
            holePosition[1] *= -1.;
            new G4PVPlacement(nullptr, holePosition, innerMajorHoleLV, "InnerMajorHole_PV",
                              resultLV, false, 2, OverlapCheck);
            holePosition[0] *= -1.;
            new G4PVPlacement(nullptr, holePosition, innerMajorHoleLV, "InnerMajorHole_PV",
                              resultLV, false, 3, OverlapCheck);

            new G4PVPlacement(nullptr, {0, 0, -type3_subPanelThick},
                              innerMinorHoleLV1, "InnerMinorHole1_PV", resultLV, false, 0, OverlapCheck);

            holePosition = G4ThreeVector(0, type3_innerMinorHole2SizeY / 2. + type3_innerMinorHoleDistance +
                                                type3_innerMinorHole1SizeY / 2.) +
                           G4ThreeVector(0, 0, -type3_subPanelThick);
            new G4PVPlacement(nullptr, holePosition, innerMinorHoleLV2, "InnerMinorHole2_PV",
                              resultLV, false, 0, OverlapCheck);
            holePosition[1] *= -1.;
            new G4PVPlacement(nullptr, holePosition, innerMinorHoleLV2, "InnerMinorHole2_PV",
                              resultLV, false, 1, OverlapCheck);

            G4cout << __PRETTY_FUNCTION__ << ": Punching the holes at the four corners..."
                   << G4endl;
            edgeHole =
                new G4Tubs("EdgeHole", 0, edgeHoleRadius, type3_panelFrameThick / 2., 0, 360. * deg);
            edgeHoleLV = new G4LogicalVolume(edgeHole, aSpaceMaterial, "EdgeHole_LV");

            holePosition =
                G4ThreeVector(upperPanelSizeX / 2. - type3_edgeHoleRadius / 2. - type3_edgeHoleDistanceEX);
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 4, OverlapCheck);
            holePosition[1] += type3_edgeHoleDistanceY;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 5, OverlapCheck);
            holePosition[1] -= 2. * type3_edgeHoleDistanceY;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 3, OverlapCheck);

            holePosition[0] *= -1.;
            holePosition[1] += type3_edgeHoleDistanceY;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 1, OverlapCheck);
            holePosition[1] += type3_edgeHoleDistanceY;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 2, OverlapCheck);
            holePosition[1] -= 2. * type3_edgeHoleDistanceY;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 0, OverlapCheck);

            holePosition =
                G4ThreeVector(0, upperPanelSizeY / 2. - type3_edgeHoleRadius / 2. - type3_edgeHoleDistanceEY);
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 10, OverlapCheck);
            holePosition[0] += type3_edgeHoleDistanceX;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 11, OverlapCheck);
            holePosition[0] -= 2. * type3_edgeHoleDistanceX;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 9, OverlapCheck);

            holePosition[1] *= -1.;
            holePosition[0] += type3_edgeHoleDistanceX;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 7, OverlapCheck);
            holePosition[0] += type3_edgeHoleDistanceX;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 8, OverlapCheck);
            holePosition[0] -= 2. * type3_edgeHoleDistanceX;
            new G4PVPlacement(nullptr, holePosition, edgeHoleLV, "EdgeHole_PV", resultLV, true, 6, OverlapCheck);
        }

    } else {
        G4cout << __PRETTY_FUNCTION__ << ": Simplified CMO detector array upper panel selected"
               << G4endl;
        resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
            "Simple_DetectorArrayUpperPanel_LV", false);
        if (resultLV == nullptr) {
            simpleUpperPanel = new G4Box("CMOArrayUpperPanelBox", upperPanelSizeX / 2.,
                                         upperPanelSizeY / 2., upperPanelThick / 2.);
            resultLV         = new G4LogicalVolume(simpleUpperPanel, aFrameMaterial,
                                           "Simple_DetectorArrayUpperPanel_LV");
        }
    }
    return resultLV;
}

G4LogicalVolume *AmoreDetectorConstruction::Build_I_DetectorArrayBottomPanel(
    G4bool aRealistic, G4Material *aFrameMaterial, G4Material *aSpaceMaterial) {
    G4LogicalVolume *resultLV;
    using namespace std;
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4Box *simpleBottomPanel;

		//for realistic geometry
		G4Box *topSmallBox;
		G4Box *topSmallBoxHole;
		G4Box *topBigBox;
		G4Box *topBigBoxHole;
		G4Box *baseSmallBox;
		G4Box *baseSmallBoxHole;
		G4Box *baseBigBox;
		G4Box *baseBigBoxHole;
		G4Box *longFrameBox;
		G4Box *longFrameHole;

		G4ExtrudedSolid *longFrameCutter;

		G4VSolid *topSmallSolid;
		G4VSolid *topBigSolid;
		G4VSolid *baseSmallSolid;
		G4VSolid *baseBigSolid;
		G4VSolid *longFrameSolid;
		G4VSolid *cuttedlongFrameSolid;

		G4LogicalVolume *topSmallFrameLV;
		G4LogicalVolume *topBigFrameLV;
		G4LogicalVolume *baseSmallFrameLV;
		G4LogicalVolume *baseBigFrameLV;
		G4LogicalVolume *longFrameLV;
		G4LogicalVolume *cuttedlongFrameLV;

		G4ThreeVector topSmallFramePos;
		G4ThreeVector topBigFramePos;
		G4ThreeVector baseSmallFramePos;
		G4ThreeVector baseBigFramePos;
		G4ThreeVector longFramePos;

		G4double topSmallFrameDis;
		G4double baseBoxSizeY = topSmallDis + 2 * bottomPanelFrameThick;


    std::vector<G4TwoVector> cutterSolidVertices;

    simpleBottomPanel = new G4Box("Simple_DetectorArrayBottomPanelBox", bottomPanelSizeX / 2.,
                                  bottomPanelSizeY / 2., bottomPanelThick / 2.);
    if (aRealistic) {
				if(fDbgMsgOn) {
					G4cout << __PRETTY_FUNCTION__
               << ": Realistic CMO detector array bottom panel selected (It could be very slow)"
               << G4endl;
				}
        resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
            "Realistic_DetectorArrayBottomPanel_LV", false);
        if (resultLV == nullptr) {
            resultLV = new G4LogicalVolume(simpleBottomPanel, aSpaceMaterial,
                                           "Realistic_DetectorArrayBottomPanel_LV");
            resultLV->SetVisAttributes(fI_Invisible_VisAttr);

						topSmallBox     = new G4Box("topSmall_Box",bottomPanelFrameWidth/2.,
								                topBoxSizeY/2.,topBoxThick/2.);
						topSmallBoxHole = new G4Box("topSmall_Hole",bottomPanelFrameWidth/2. + solidBooleanTol, 
								                topHoleSizeY/2., topHoleThick/2. + solidBooleanTol);
						topSmallSolid   = new G4SubtractionSolid("topSmallSolid",topSmallBox,topSmallBoxHole,
								                nullptr,{0,0, -bottomPanelFrameThick/2. - solidBooleanTol});

						topSmallFrameLV = new G4LogicalVolume(topSmallSolid, aFrameMaterial, "bottomPanelSmallFrame1_LV");
						topSmallFrameLV->SetVisAttributes(fI_CopperDefault_VisAttr);

						topBigBox     = new G4Box("topBig_Box", bottomPanelFrameBigWidth / 2.,
								              topBoxSizeY / 2., topBoxThick / 2.);
						topBigBoxHole = new G4Box("topBig_Hole", bottomPanelFrameBigWidth / 2.,
								              topHoleSizeY / 2., topHoleThick / 2. + solidBooleanTol);
						topBigSolid   = new G4SubtractionSolid("topBigSolid", topBigBox, topBigBoxHole,
								              nullptr, {0, 0, -bottomPanelFrameThick / 2. - solidBooleanTol});

						topBigFrameLV = new G4LogicalVolume(topBigSolid, aFrameMaterial, "bottomPanelBigFrame1_LV");
						topBigFrameLV->SetVisAttributes(fI_CopperDefault_VisAttr);

						baseSmallBox     = new G4Box("baseSmall_Box",bottomPanelFrameWidth/2., 
															   baseBoxSizeY/2., bottomPanelFrameThick/2.);
						baseSmallBoxHole = new G4Box("baseSmall_Hole", baseSmallHoleSizeX / 2.,
								                 baseHoleSizeY / 2., bottomPanelFrameThick / 2. + solidBooleanTol);
						baseSmallSolid   = new G4SubtractionSolid("baseSmallSolid", baseSmallBox, baseSmallBoxHole,
								                 nullptr, {0, 0, 0});

						baseSmallFrameLV = new G4LogicalVolume(baseSmallSolid, aFrameMaterial, "bottomPanelSmallFrame2_LV");
						baseSmallFrameLV->SetVisAttributes(fI_CopperDefault_VisAttr);

						baseBigBox     = new G4Box("baseBic_Box", bottomPanelFrameBigWidth / 2.,
								               baseBoxSizeY / 2., bottomPanelFrameThick / 2. );
						baseBigBoxHole = new G4Box("baseBig_Hole", baseBigHoleSizeX / 2.,
								               baseHoleSizeY / 2., bottomPanelFrameThick / 2. + solidBooleanTol);
						baseBigSolid   = new G4SubtractionSolid("baseBigSolid", baseBigBox, baseBigBoxHole,
								               nullptr, {0, 0, 0});

						baseBigFrameLV = new G4LogicalVolume(baseBigSolid, aFrameMaterial, "bottomPanelBigFrame2_LV");
						baseBigFrameLV->SetVisAttributes(fI_CopperDefault_VisAttr);

						longFrameBox    = new G4Box("longFrame_Box",bottomPanelSizeX / 2.,
								             bottomPanelFrameWidth / 2., bottomPanelFrameThick / 2.);
						longFrameHole   = new G4Box("longFrame_Hole", longFrameHoleSizeX / 2.,
								             longFrameHoleSizeY / 2., bottomPanelFrameThick / 2. + solidBooleanTol);

						cutterSolidVertices.resize(4);
						cutterSolidVertices[0] = {cutterSizeMinorX/2., -bottomPanelFrameThick};
						cutterSolidVertices[1] = {-cutterSizeMinorX/2., -bottomPanelFrameThick};
						cutterSolidVertices[2] = {-cutterSizeMajorX/2., bottomPanelFrameThick+solidBooleanTol};
						cutterSolidVertices[3] = {cutterSizeMajorX/2.,bottomPanelFrameThick+solidBooleanTol};

						longFrameCutter = new G4ExtrudedSolid("longFrameCutterSolid", cutterSolidVertices,
								                bottomPanelFrameThick / 2. + solidBooleanTol, 0, 1., 0, 1.);

						longFrameSolid = new G4SubtractionSolid("longFrameSolid",longFrameBox, longFrameHole, 
								               nullptr, {-bottomPanelSizeX / 4., 0, 0});
						longFrameSolid = new G4SubtractionSolid("longFrameSolid",longFrameSolid, longFrameHole,
								               nullptr, {bottomPanelSizeX / 4., 0, 0});
						longFrameSolid = new G4SubtractionSolid("longFrameSolid", longFrameSolid, baseBigBoxHole,
								               nullptr, {0, 0, 0});

						longFrameLV = new G4LogicalVolume(longFrameSolid, aFrameMaterial, "bottomPanelLongFrame1_LV");
						longFrameLV->SetVisAttributes(fI_CopperDefault_VisAttr);

						cuttedlongFrameSolid = new G4SubtractionSolid("cuttedlongFrameSolid",longFrameSolid, longFrameCutter,
								                     nullptr, {0,bottomPanelFrameThick/2.,0});
						cuttedlongFrameSolid = new G4SubtractionSolid("cuttedlongFrameSolid",
										                 cuttedlongFrameSolid, longFrameCutter, nullptr, 
																		 {bottomPanelSizeX/2. - bottomPanelFrameWidth + cutterSizeMinorX/2.,
																		 bottomPanelFrameThick/2.,0});
						cuttedlongFrameSolid = new G4SubtractionSolid("cuttedlongFrameSolid",
										                 cuttedlongFrameSolid, longFrameCutter, nullptr, 
																		 {-bottomPanelSizeX/2. + bottomPanelFrameWidth - cutterSizeMinorX/2.,
																		 bottomPanelFrameThick/2.,0});

						cuttedlongFrameLV = new G4LogicalVolume(cuttedlongFrameSolid, aFrameMaterial, "bottomPanelLongFrame2_LV");
						cuttedlongFrameLV->SetVisAttributes(fI_CopperDefault_VisAttr);

						// positioning
						new G4PVPlacement(nullptr, {0,0,-bottomPanelThick/2. + bottomPanelFrameThick/2.}, longFrameLV,
								"bottomPanelLongFrame1_PV", resultLV, false, 0, OverlapCheck);

						G4RotationMatrix *FrameRotationMtx = new G4RotationMatrix();
						FrameRotationMtx->rotateZ(180. * deg);
						new G4PVPlacement(nullptr, {0,-bottomPanelSizeY/2. + bottomPanelFrameWidth/2.,
								-bottomPanelThick/2. + bottomPanelFrameThick/2.},
								cuttedlongFrameLV, "bottomPaneLongFrame2_PV", resultLV, false, 0, OverlapCheck);
						new G4PVPlacement(FrameRotationMtx, {0,bottomPanelSizeY/2. - bottomPanelFrameWidth/2., 
								-bottomPanelThick/2. + bottomPanelFrameThick/2.},
								cuttedlongFrameLV, "bottomPaneLongFrame2_PV", resultLV, false, 0, OverlapCheck);

						topSmallFramePos = {bottomPanelSizeX/2. - bottomPanelFrameWidth/2., 
							                 bottomPanelSizeY/2. - topBoxSizeY / 2., 
															 bottomPanelThick/2.- topBoxThick/2.};
						topSmallFrameDis = topBoxSizeY + topSmallDis;

						topBigFramePos   = {0, bottomPanelSizeY/2. - topBoxSizeY / 2.,
							                  bottomPanelThick/2. - topBoxThick / 2. };

						for(int i = 0; i < 2; i++)
						{
							baseSmallFramePos = topSmallFramePos;
							baseSmallFramePos[1] -= topSmallFrameDis / 2.;
							baseSmallFramePos[2] -= topBoxThick/2. + bottomPanelFrameThick /2.; 

							new G4PVPlacement(nullptr, baseSmallFramePos, baseSmallFrameLV,
									"bottomPanelSmallFrame2_PV", resultLV, false, 0, OverlapCheck);
							baseSmallFramePos[0] *= -1;
							new G4PVPlacement(nullptr, baseSmallFramePos, baseSmallFrameLV,
									"bottomPanelSmallFrame2_PV", resultLV, false, 0, OverlapCheck);

							baseBigFramePos = topBigFramePos;
							baseBigFramePos[1] -= topSmallFrameDis / 2.;
							baseBigFramePos[2] -= topBoxThick / 2. + bottomPanelFrameThick / 2.;

							new G4PVPlacement(nullptr, baseBigFramePos, baseBigFrameLV,
									"bottomPanelBigFrame2_PV", resultLV, false, 0, OverlapCheck);

							for(int j = 0; j < 2; j++)
							{
								new G4PVPlacement(nullptr, topSmallFramePos, topSmallFrameLV, 
										"bottomPanelSmallFrame1_PV", resultLV, false, 0, OverlapCheck);	
								topSmallFramePos[0] *= -1;
								new G4PVPlacement(nullptr, topSmallFramePos, topSmallFrameLV,
										"bottomPanelSmallFrame1_PV", resultLV, false, 0, OverlapCheck);

								new G4PVPlacement(nullptr, topBigFramePos, topBigFrameLV,
										"bottomPanelBigFrame1_PV", resultLV, false, 0, OverlapCheck);

								if(j!=1) 
								{
									topSmallFramePos[1] -= topSmallFrameDis;
									topBigFramePos[1] -= topSmallFrameDis;
								}
							}
							topSmallFramePos[1] *= -1;
							topBigFramePos[1] *= -1;
						}
        }
    } else {
        G4cout << __PRETTY_FUNCTION__ << ": Simplified CMO detector array bottom panel selected"
               << G4endl;
        resultLV = new G4LogicalVolume(simpleBottomPanel, aFrameMaterial,
                                       "Simple_DetectorArrayBottomPanel_LV");
        resultLV->SetVisAttributes(fI_CopperDefault_VisAttr);
    }
    return resultLV;
}

// There are 3 type ( type1: run5, type2: run6, type3: real AmoreI )
G4LogicalVolume *AmoreDetectorConstruction::Build_I_DetectorArray(G4int aType,
                                                                  G4Material *aArrayFrameMaterial,
                                                                  G4Material *aSpaceMaterial) {
    using namespace std;
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4LogicalVolume *resultLV;

    // For Type 1
    G4double moduleMarginFromXEndXSide;
    G4double moduleMarginFromXEndYSide;

    G4double towerRodLength;
    G4double towerRodBottomLength;

    G4Box *arrayBox;
    G4Tubs *towerRodTopFullTub;
    G4Tubs *nowTowerRodTopTub;
    G4Tubs *towerRodBottomTub;
    G4VSolid *towerRodSolid;
    G4LogicalVolume *towerRodLV;

    G4ThreeVector nowModulePos;
    G4ThreeVector nowModulePosInXYPlane;
    G4ThreeVector arrayUpperEndPos;
    G4ThreeVector arrayBottomEndPos;
    G4ThreeVector nowTowerRodPos;
    G4ThreeVector nowTowerRodDisplace;

		// For Type 3
		G4Box *weightForUpperPanelBox;
		G4Box *AdditionalHeightBox;

		G4LogicalVolume *weightForUpperPanelLV;
		G4LogicalVolume *AdditionalHeightLV;

		G4ThreeVector nowWeightPos;
		G4ThreeVector WeightPosInXYPlane;

    switch (aType) {
        case 1: // Run5 configuration
            resultLV =
                G4LogicalVolumeStore::GetInstance()->GetVolume("DetectorArray_Type1_LV", false);
            if (resultLV == nullptr) {
                arrayBox =
                    new G4Box("DetectorArray_Type1_Box", arraySizeX / 2. + type1_detModuleMarginX,
                              arraySizeY / 2. + type1_detModuleMarginY, arrayHeight / 2.);
                resultLV = new G4LogicalVolume(arrayBox, aSpaceMaterial, "DetectorArray_Type1_LV");
                resultLV->SetVisAttributes(fI_Invisible_VisAttr);

                new G4PVPlacement(
                    nullptr, {0, 0, arrayHeight / 2. - upperPanelThick / 2.},
                    Build_I_DetectorArrayUpperPanel(true, aArrayFrameMaterial, aSpaceMaterial),
                    "DetectorArray_UpperPanel_Type1_PV", resultLV, false, 0, OverlapCheck);
                new G4PVPlacement(
                    nullptr, {0, 0, -arrayHeight / 2. + bottomPanelThick / 2.},
                    Build_I_DetectorArrayBottomPanel(true, aArrayFrameMaterial, aSpaceMaterial),
                    "DetectorArray_UpperPanel_Type1_PV", resultLV, false, 0, OverlapCheck);

                towerRodTopFullTub = new G4Tubs("TowerRodTopFull_Tub", 0, towerRodTopRadius,
                                                towerRodTopLength / 2., 0, 360. * deg);

                arrayUpperEndPos  = G4ThreeVector(0, 0, arrayHeight / 2. - upperPanelThick);
                arrayBottomEndPos = G4ThreeVector(0, 0, -arrayHeight / 2. + bottomPanelThick);

                moduleMarginFromXEndXSide =
                    (bottomPanelSizeX - 2. * copperFrameSizeX - type1_distanceBetweenTower) / 2.;
                moduleMarginFromXEndYSide =
                    (bottomPanelSizeY - 2. * copperFrameSizeY - type1_distanceBetweenTower) / 2.;

                nowModulePosInXYPlane = G4ThreeVector(
                    -bottomPanelSizeX / 2. + copperFrameSizeX / 2. + moduleMarginFromXEndXSide,
                    bottomPanelSizeY / 2. - copperFrameSizeY / 2. - moduleMarginFromXEndYSide, 0);
                G4int nowCopyNo = 0;

                for (G4int i = 0; i < totalNumOfTower; i++) {

                    nowModulePos = arrayBottomEndPos;
                    nowModulePos += nowModulePosInXYPlane;
                    for (G4int j = 0; j < maxModuleNumInTower; j++) {
                        const CrystalModuleInfo &nowInfo = crystalModuleInfoList[i][j];
                        if (nowInfo.fCrystalHeight < 0.) continue;
                        if (nowModulePos.getZ() + nowInfo.fCrystalHeight >
                            arrayUpperEndPos.getZ()) {
                            G4Exception(
                                __PRETTY_FUNCTION__, "AMORE_I_DETARRAY",
                                G4ExceptionSeverity::JustWarning,
                                "Total height of a tower has reached the upper end of array. "
                                "Remaining modules of this tower will be truncated.");
                            break;
                        }

                        nowModulePos[nowModulePos.Z] +=
                            nowInfo.fCrystalHeight / 2. + smallBlock1SizeZ;

                        G4int nowPosID[2] = {i, j};
                        AmoreModuleSDInfo nowSDInfo;
                        nowSDInfo.fModuleName =
                            (nowInfo.fName + "[" + to_string(i) + "][" + to_string(j) + "]");
                        nowSDInfo.fModuleID         = nowCopyNo;
                        nowSDInfo.fCrystalPosIdx[0] = nowPosID[0];
                        nowSDInfo.fCrystalPosIdx[1] = nowPosID[1];

                        G4LogicalVolume *nowDetectorModuleLV = Build_I_SingleDetectorModule(
                            false, 1, nowInfo, aSpaceMaterial, &nowSDInfo, (j == 0));
                        nowDetectorModuleLV->SetVisAttributes(fI_Invisible_VisAttr);
                        fI_DetectorModuleRegion->AddRootLogicalVolume(nowDetectorModuleLV);

                        G4VPhysicalVolume *nowDetPV =
                            new G4PVPlacement(nullptr, nowModulePos, nowDetectorModuleLV,
                                              ("DetEnvl_" + nowInfo.fName + "[" + to_string(i) +
                                               "][" + to_string(j) + "]_PV")
                                                  .c_str(),
                                              resultLV, false, nowCopyNo++, OverlapCheck);
                        nowSDInfo.fModulePV = nowDetPV;
                        fModuleSDInfos.insert(nowSDInfo);

                        nowModulePos[nowModulePos.Z] +=
                            nowInfo.fCrystalHeight / 2. + smallBlock1SizeZ;
                    }
                    Place_I_GeometriesInRooftopForTowerOf(resultLV, 1, nowModulePos);

                    towerRodLength = fabs(nowModulePos.getZ() - arrayUpperEndPos.getZ());
                    if (towerRodLength > towerRodTopLength) {
                        towerRodBottomLength = towerRodLength - towerRodTopLength;
                        towerRodBottomTub =
                            new G4Tubs("Tower[" + to_string(i) + "]_BottomTub", 0, towerRodRadius,
                                       towerRodBottomLength / 2. + solidBooleanTol, 0, 360 * deg);
                        towerRodSolid = new G4UnionSolid(
                            "Tower[" + to_string(i) + "]_RodSolid", towerRodBottomTub,
                            towerRodTopFullTub, nullptr,
                            {0, 0,
                             towerRodBottomLength / 2. + towerRodTopLength / 2. - solidBooleanTol});
                        towerRodLV = new G4LogicalVolume(towerRodSolid, aArrayFrameMaterial,
                                                         "Tower[" + to_string(i) + "]_Rod_LV");
                        towerRodLV->SetVisAttributes(fI_CopperDefault_VisAttr);

                        nowTowerRodPos =
                            G4ThreeVector(0, 0, towerRodBottomLength / 2. + solidBooleanTol);

                    } else {
                        nowTowerRodTopTub =
                            new G4Tubs("Tower[" + to_string(i) + "]Tub", 0, towerRodRadius,
                                       towerRodLength / 2., 0, 360 * deg);
                        towerRodLV = new G4LogicalVolume(nowTowerRodTopTub, aArrayFrameMaterial,
                                                         "Tower[" + to_string(i) + "]_Rod_LV");
                        towerRodLV->SetVisAttributes(fI_CopperDefault_VisAttr);

                        nowTowerRodPos = G4ThreeVector(0, 0, towerRodLength / 2.);
                    }
                    nowTowerRodDisplace =
                        G4ThreeVector(copperFrameSizeX / 2. + smallBlock1DistX,
                                      copperFrameSizeY / 2. + smallBlock1DistY, 0);

                    nowTowerRodPos += nowModulePos;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 0, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 1, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodDisplace[0] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 2, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 3, OverlapCheck);

                    nowModulePosInXYPlane[0] *= -1.;
                    if (i == 1) nowModulePosInXYPlane[1] *= -1.;
                }
            }
            break;
        case 2: // Run6 configuration
            resultLV =
                G4LogicalVolumeStore::GetInstance()->GetVolume("DetectorArray_Type2_LV", false);
            if (resultLV == nullptr) {
                arrayBox =
                    new G4Box("DetectorArray_Type2_Box", arraySizeX / 2. + type2_detModuleMarginX,
                              arraySizeY / 2. + type2_detModuleMarginY, arrayHeight / 2.);
                resultLV = new G4LogicalVolume(arrayBox, aSpaceMaterial, "DetectorArray_Type2_LV");
                resultLV->SetVisAttributes(fI_Invisible_VisAttr);

                new G4PVPlacement(
                    nullptr, {0, 0, arrayHeight / 2. - upperPanelThick / 2.},
                    Build_I_DetectorArrayUpperPanel(true, aArrayFrameMaterial, aSpaceMaterial),
                    "DetectorArray_UpperPanel_Type2_PV", resultLV, false, 0, OverlapCheck);
                new G4PVPlacement(
                    nullptr, {0, 0, -arrayHeight / 2. + bottomPanelThick / 2.},
                    Build_I_DetectorArrayBottomPanel(true, aArrayFrameMaterial, aSpaceMaterial),
                    "DetectorArray_UpperPanel_Type2_PV", resultLV, false, 0, OverlapCheck);

                towerRodTopFullTub = new G4Tubs("TowerRodTopFullTub", 0, towerRodTopRadius,
                                                towerRodTopLength / 2., 0, 360. * deg);

                arrayUpperEndPos  = G4ThreeVector(0, 0, arrayHeight / 2. - upperPanelThick);
                arrayBottomEndPos = G4ThreeVector(0, 0, -arrayHeight / 2. + bottomPanelThick);

                moduleMarginFromXEndXSide =
                    (bottomPanelSizeX - 2. * copperFrameSizeX - type2_distanceBetweenTower) / 2.;
                moduleMarginFromXEndYSide =
                    (bottomPanelSizeY - 2. * copperFrameSizeY - type2_distanceBetweenTower) / 2.;

                nowModulePosInXYPlane = G4ThreeVector(
                    -bottomPanelSizeX / 2. + copperFrameSizeX / 2. + moduleMarginFromXEndXSide,
                    bottomPanelSizeY / 2. - copperFrameSizeY / 2. - moduleMarginFromXEndYSide, 0);
                G4int nowCopyNo = 0;

                for (G4int i = 0; i < totalNumOfTower; i++) {

                    nowModulePos = arrayBottomEndPos;
                    nowModulePos += nowModulePosInXYPlane;
                    for (G4int j = 0; j < maxModuleNumInTower; j++) {
                        const CrystalModuleInfo &nowInfo = crystalModuleInfoList[i][j];
                        if (nowInfo.fCrystalHeight < 0.) continue;
                        if (nowModulePos.getZ() + nowInfo.fCrystalHeight >
                            arrayUpperEndPos.getZ()) {
                            G4Exception(
                                __PRETTY_FUNCTION__, "AMORE_I_DETARRAY",
                                G4ExceptionSeverity::JustWarning,
                                "Total height of a tower has reached the upper end of array. "
                                "Remaining modules of this tower will be truncated.");
                            break;
                        }

                        nowModulePos[nowModulePos.Z] +=
                            nowInfo.fCrystalHeight / 2. + smallBlock1SizeZ;

                        G4int nowPosID[2] = {i, j};
                        AmoreModuleSDInfo nowSDInfo;
                        nowSDInfo.fModuleName =
                            (nowInfo.fName + "[" + to_string(i) + "][" + to_string(j) + "]");
                        nowSDInfo.fModuleID         = nowCopyNo;
                        nowSDInfo.fCrystalPosIdx[0] = nowPosID[0];
                        nowSDInfo.fCrystalPosIdx[1] = nowPosID[1];

                        G4LogicalVolume *nowDetectorModuleLV = Build_I_SingleDetectorModule(
                            false, 2, nowInfo, aSpaceMaterial, &nowSDInfo, (j == 0));
                        nowDetectorModuleLV->SetVisAttributes(fI_Invisible_VisAttr);
                        fI_DetectorModuleRegion->AddRootLogicalVolume(nowDetectorModuleLV);

                        G4VPhysicalVolume *nowDetPV =
                            new G4PVPlacement(nullptr, nowModulePos, nowDetectorModuleLV,
                                              ("DetEnvl_" + nowInfo.fName + "[" + to_string(i) +
                                               "][" + to_string(j) + "]_PV")
                                                  .c_str(),
                                              resultLV, false, nowCopyNo++, OverlapCheck);
                        nowSDInfo.fModulePV = nowDetPV;
                        fModuleSDInfos.insert(nowSDInfo);

                        nowModulePos[nowModulePos.Z] +=
                            nowInfo.fCrystalHeight / 2. + smallBlock1SizeZ;
                    }
                    Place_I_GeometriesInRooftopForTowerOf(resultLV, 2, nowModulePos);

                    towerRodLength = fabs(nowModulePos.getZ() - arrayUpperEndPos.getZ());
                    if (towerRodLength > towerRodTopLength) {
                        towerRodBottomLength = towerRodLength - towerRodTopLength;
                        towerRodBottomTub =
                            new G4Tubs("Tower[" + to_string(i) + "]_BottomTub", 0, towerRodRadius,
                                       towerRodBottomLength / 2. + solidBooleanTol, 0, 360 * deg);
                        towerRodSolid = new G4UnionSolid(
                            "Tower[" + to_string(i) + "]_RodSolid", towerRodBottomTub,
                            towerRodTopFullTub, nullptr,
                            {0, 0,
                             towerRodBottomLength / 2. + towerRodTopLength / 2. - solidBooleanTol});
                        towerRodLV = new G4LogicalVolume(towerRodSolid, aArrayFrameMaterial,
                                                         "Tower[" + to_string(i) + "]_Rod_LV");
                        towerRodLV->SetVisAttributes(fI_CopperDefault_VisAttr);

                        nowTowerRodPos =
                            G4ThreeVector(0, 0, towerRodBottomLength / 2. + solidBooleanTol);

                    } else {
                        nowTowerRodTopTub =
                            new G4Tubs("Tower[" + to_string(i) + "]Tub", 0, towerRodRadius,
                                       towerRodLength / 2., 0, 360 * deg);
                        towerRodLV = new G4LogicalVolume(nowTowerRodTopTub, aArrayFrameMaterial,
                                                         "Tower[" + to_string(i) + "]_Rod_LV");
                        towerRodLV->SetVisAttributes(fI_CopperDefault_VisAttr);

                        nowTowerRodPos = G4ThreeVector(0, 0, towerRodLength / 2.);
                    }
                    nowTowerRodDisplace =
                        G4ThreeVector(copperFrameSizeX / 2. + smallBlock1DistX,
                                      copperFrameSizeY / 2. + smallBlock1DistY, 0);

                    nowTowerRodPos += nowModulePos;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 0, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 1, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodDisplace[0] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 2, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 3, OverlapCheck);

                    nowModulePosInXYPlane[0] *= -1.;
                    if (i == 1) nowModulePosInXYPlane[1] *= -1.;
                }
            }
            break;
        case 3: // REAL configuration
            resultLV =
                G4LogicalVolumeStore::GetInstance()->GetVolume("DetectorArray_Type3_LV", false);
            if (resultLV == nullptr) {
                arrayBox =
                    new G4Box("DetectorArray_Type3_Box", arraySizeX / 2. + type2_detModuleMarginX,
                              arraySizeY / 2. + type2_detModuleMarginY, arrayHeight / 2.);
                resultLV = new G4LogicalVolume(arrayBox, aSpaceMaterial, "DetectorArray_Type3_LV");
                resultLV->SetVisAttributes(fI_Invisible_VisAttr);

                new G4PVPlacement(
                    nullptr, {0, 0, arrayHeight / 2. - type3_panelFrameThick / 2.},
                    Build_I_DetectorArrayUpperPanel(true, aArrayFrameMaterial, aSpaceMaterial),
                    "DetectorArray_UpperPanel_Type3_PV", resultLV, false, 0, OverlapCheck);

                new G4PVPlacement(
                    nullptr, {0, 0, -arrayHeight / 2. + bottomPanelThick / 2. + type3_additionalHeightThick},
                    Build_I_DetectorArrayBottomPanel(true, aArrayFrameMaterial, aSpaceMaterial),
                    "DetectorArray_BottomPanel_Type3_PV", resultLV, false, 0, OverlapCheck);

								AdditionalHeightBox = new G4Box("AdditionalHeight_Box",
										additionalHeightSizeX / 2., additionalHeightSizeY / 2., type3_additionalHeightThick / 2.);

								AdditionalHeightLV = new G4LogicalVolume(AdditionalHeightBox, aArrayFrameMaterial, "AdditionalHeight_LV");
								AdditionalHeightLV->SetVisAttributes(fI_CopperDefault_VisAttr);

							  new G4PVPlacement(
										nullptr, {0, 0, -arrayHeight / 2. + type3_additionalHeightThick / 2.},
										AdditionalHeightLV, "AdditionalHeight_PV",resultLV,false,0,OverlapCheck);

                arrayUpperEndPos  = G4ThreeVector(0, 0, arrayHeight / 2. - type3_upperPanelThick);
                arrayBottomEndPos = G4ThreeVector(0, 0, -arrayHeight / 2. + bottomPanelThick + type3_additionalHeightThick);
                moduleMarginFromXEndXSide =
                    (bottomPanelSizeX - 2. * type3_copperFrameSizeX - type2_distanceBetweenTower) / 2.;
                moduleMarginFromXEndYSide =
                    (bottomPanelSizeY - 2. * type3_copperFrameSizeY - type2_distanceBetweenTower) / 2.;

                nowModulePosInXYPlane = G4ThreeVector(
                    -bottomPanelSizeX / 2. + type3_copperFrameSizeX / 2. + moduleMarginFromXEndXSide,
                    bottomPanelSizeY / 2. - type3_copperFrameSizeY / 2. - moduleMarginFromXEndYSide, 0);
                G4int nowCopyNo = 0;

								WeightPosInXYPlane = G4ThreeVector(
											-type3_innerPartSizeX / 2. + type3_innerMajorHoleSizeX / 2. + type3_innerMajorHoleDistanceX,
											-type3_innerPartSizeY / 2. + type3_innerMajorHoleSizeY / 2. + type3_innerMajorHoleDistanceY,0	);

                for (G4int i = 0; i < totalNumOfTower; i++) {

                    nowModulePos = arrayBottomEndPos;
                    nowModulePos += nowModulePosInXYPlane;
                    for (G4int j = 0; j < maxModuleNumInTower; j++) {
                        const CrystalModuleInfo &nowInfo = crystalModuleInfoList[i][j];
                        if (nowInfo.fCrystalHeight < 0.) continue;
                        if (nowModulePos.getZ() + nowInfo.fCrystalHeight >
                            arrayUpperEndPos.getZ()) {
                            G4Exception(
                                __PRETTY_FUNCTION__, "AMORE_I_DETARRAY",
                                G4ExceptionSeverity::JustWarning,
                                "Total height of a tower has reached the upper end of array. "
                                "Remaining modules of this tower will be truncated.");
                            break;
                        }

                        nowModulePos[nowModulePos.Z] +=
                            nowInfo.fCrystalHeight / 2. + type3_smallBlock3SizeZ / 2. + type3_smallBlock1SizeZ / 2.;

                        G4int nowPosID[2] = {i, j};
                        AmoreModuleSDInfo nowSDInfo;
                        nowSDInfo.fModuleName =
                            (nowInfo.fName + "[" + to_string(i) + "][" + to_string(j) + "]");
                        nowSDInfo.fModuleID         = nowCopyNo;
                        nowSDInfo.fCrystalPosIdx[0] = nowPosID[0];
                        nowSDInfo.fCrystalPosIdx[1] = nowPosID[1];

                        G4LogicalVolume *nowDetectorModuleLV = Build_I_SingleDetectorModule(
                            true, 3, nowInfo, aSpaceMaterial, &nowSDInfo, (j == 0));
                        nowDetectorModuleLV->SetVisAttributes(fI_Invisible_VisAttr);
                        fI_DetectorModuleRegion->AddRootLogicalVolume(nowDetectorModuleLV);

                        G4VPhysicalVolume *nowDetPV =
                            new G4PVPlacement(nullptr, nowModulePos, nowDetectorModuleLV,
                                              ("DetEnvl_" + nowInfo.fName + "[" + to_string(i) +
                                               "][" + to_string(j) + "]_PV")
                                                  .c_str(),
                                              resultLV, false, nowCopyNo++, OverlapCheck);
                        nowSDInfo.fModulePV = nowDetPV;
                        fModuleSDInfos.insert(nowSDInfo);

                        nowModulePos[nowModulePos.Z] +=
                            nowInfo.fCrystalHeight / 2. + type3_smallBlock1SizeZ /2. + type3_smallBlock3SizeZ/2.;

												if(fDbgMsgOn) {
													std::cout << "crystal module position: " << nowDetPV->GetTranslation() << std::endl; 
												}
                    }
										// bolts for upper crystal
										Place_I_GeometriesInRooftopForTowerOf(resultLV, 3, nowModulePos);

										// weight for upper panel
										nowWeightPos = G4ThreeVector(
											0, 0, arrayUpperEndPos.getZ() - weightForUpperPanelThick[i] / 2.);
                    nowWeightPos += WeightPosInXYPlane;

										if(i!=3){
											weightForUpperPanelBox = new G4Box("weightForUpperPanel_Box",
													weightForUpperPanelSizeX / 2.,
													weightForUpperPanelSizeY / 2.,
													weightForUpperPanelThick[i] / 2.);

											weightForUpperPanelLV = new G4LogicalVolume(weightForUpperPanelBox, aArrayFrameMaterial,
																								"Weight[" + to_string(i) + "]_LV");
											weightForUpperPanelLV->SetVisAttributes(fI_CopperDefault_VisAttr);

											new G4PVPlacement(nullptr, nowWeightPos, weightForUpperPanelLV, "weightForUpperPanel_PV",
                              resultLV, false, 0, OverlapCheck);
										}

										// tower Rod	
                    towerRodLength = fabs(nowModulePos.getZ() - arrayUpperEndPos.getZ());

										towerRodSolid = 
                        new G4Tubs("Tower[" + to_string(i) + "]_RodSolid", 0, towerRodRadius,
                                    towerRodLength / 2., 0, 360 * deg);
                    towerRodLV = new G4LogicalVolume(towerRodSolid, aArrayFrameMaterial,
                                                      "Tower[" + to_string(i) + "]_Rod_LV");
                    towerRodLV->SetVisAttributes(fI_CopperDefault_VisAttr);

                    nowTowerRodPos = G4ThreeVector(0, 0, towerRodLength / 2.);
                    nowTowerRodDisplace = G4ThreeVector(type3_copperFrameSizeX / 2. + smallBlock1DistX,
                                      type3_copperFrameSizeY / 2. + smallBlock1DistY, 0);

                    nowTowerRodPos += nowModulePos;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 0, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 1, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodDisplace[0] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 2, OverlapCheck);
                    nowTowerRodPos -= nowTowerRodDisplace;
                    nowTowerRodDisplace[1] *= -1.;
                    nowTowerRodPos += nowTowerRodDisplace;
                    new G4PVPlacement(nullptr, nowTowerRodPos, towerRodLV,
                                      "Tower[" + to_string(i) + "]_Rod_PV", resultLV, false, 3, OverlapCheck);

                    nowModulePosInXYPlane[0] *= -1.;
										WeightPosInXYPlane[0] *= -1.;
                    if (i == 1) {
											nowModulePosInXYPlane[1] *= -1.;
											WeightPosInXYPlane[1] *= -1.;
										}
                }
            }
            break;
        default:
            resultLV = nullptr;
            G4Exception(__PRETTY_FUNCTION__, "AMORE_I_TYPEERROR",
                        G4ExceptionSeverity::FatalErrorInArgument,
                        "Wrong type for an array of detector modules");
            break;
    }
    return resultLV;
}

// Modified by Jeewon at 2021. Mar.
// Type3 added for real AMoRE-I geometry
G4LogicalVolume *AmoreDetectorConstruction::Build_I_SingleDetectorModule(
    G4bool aRealistic, G4int aType, const CrystalModuleInfo &aCrystal, G4Material *aSpaceMaterial,
    AmoreModuleSDInfo *aMSDForThisModule, G4bool aIsFloor) {

    G4LogicalVolume *resultLV = nullptr;
    const G4String lvSuffix   = aIsFloor ? "_Floor_LV" : "_LV";
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4double reflectorLengthZ;

    G4Material *nowCrystalMaterial;
    G4Material *nowFrameMaterial;
    G4Material *nowReflectorMaterial;

    G4VSolid *moduleEnvelopeSolid;
    G4VSolid *reflectorSolid;
    G4VSolid *reflectorCarverSolid;
    G4Tubs *moduleSupportRodSolid;
    G4Tubs *crystalBottomGoldFilmTub;
    G4EllipticalTube *crystalSolid;

    G4LogicalVolume *crystalLV;
    G4LogicalVolume *reflectorLV;
    G4LogicalVolume *moduleSupportRodLV;
    G4LogicalVolume *crystalBottomGoldFilmLV;
    G4ThreeVector nowSupportRodDisplace;

		// For Type3
		G4double type3_moduleSizeZ = aCrystal.fCrystalHeight + type3_smallBlock1SizeZ + type3_smallBlock3SizeZ ;

		// For Type1 and Type2
    G4double reflectorCerverSizeX = reflectorLengthX - reflectorThick * 2.;
    G4double reflectorCerverSizeY = reflectorLengthY - reflectorThick * 2.;
    G4double reflectorCarverSizeZ;

    G4String nowEnvelopeLVName = "Type" + std::to_string(aType) + "_" + aCrystal.fName + lvSuffix;
    switch (aType) {
        case 1:
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(nowEnvelopeLVName, false);
            if (resultLV == nullptr) {
                G4ThreeVector nowCopperFramePos;
                nowCrystalMaterial   = G4Material::GetMaterial(aCrystal.fCrystalMaterialName);
                nowFrameMaterial     = G4Material::GetMaterial(aCrystal.fFrameMaterialName);
                nowReflectorMaterial = G4Material::GetMaterial(aCrystal.fReflectorMaterialName);
                if (nowCrystalMaterial == nullptr || nowFrameMaterial == nullptr) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_001",
                                G4ExceptionSeverity::FatalException,
                                "Material pointers for detector module are not found.");
                    nowFrameMaterial = nowCrystalMaterial = _vacuum;
                }
                moduleEnvelopeSolid = new G4Box(aCrystal.fName + "_DetModuleBox",
                                                copperFrameSizeX / 2. + type1_detModuleMarginX,
                                                copperFrameSizeY / 2. + type1_detModuleMarginY,
                                                aCrystal.fCrystalHeight / 2. + smallBlock1SizeZ);

                resultLV =
                    new G4LogicalVolume(moduleEnvelopeSolid, aSpaceMaterial, nowEnvelopeLVName);

                crystalSolid = new G4EllipticalTube(
                    aCrystal.fName + "_CrystalSolid", aCrystal.fCrystalMajorDiameter / 2.,
                    aCrystal.fCrystalMinorDiameter / 2., aCrystal.fCrystalHeight / 2.);
                crystalLV = new G4LogicalVolume(crystalSolid, nowCrystalMaterial,
                                                aCrystal.fName + "_Crystal_LV");
                crystalLV->SetVisAttributes(fI_Crystal_VisAttr);
                fI_CrystalLVs.insert(crystalLV);
                aMSDForThisModule->fCrystalLV = crystalLV;

                crystalBottomGoldFilmLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "GoldFilmBelowCrystal_Type1_LV", false);
                if (crystalBottomGoldFilmLV == nullptr) {
                    crystalBottomGoldFilmTub =
                        new G4Tubs("GoldFilmBelowCrystal_Tub", 0, crystalBottomGoldRadius,
                                   crystalBottomGoldThickness / 2., 0, 360. * deg);

                    crystalBottomGoldFilmLV = new G4LogicalVolume(crystalBottomGoldFilmTub, _gold,
                                                                  "GoldFilmBelowCrystal_Type1_LV");

                    crystalBottomGoldFilmLV->SetVisAttributes(fI_GoldDefault_VisAttr);

                    fI_CrystalGoldFilmLVs.insert(crystalBottomGoldFilmLV);
                }
                aMSDForThisModule->fCrystalGoldFilmLV = crystalBottomGoldFilmLV;

                new G4PVPlacement(
                    nullptr,
                    {0, 0, -aCrystal.fCrystalHeight / 2. + crystalBottomGoldThickness / 2.},
                    crystalBottomGoldFilmLV, aCrystal.fName + "_BottomGoldFilm_PV", crystalLV,
                    false, 0, OverlapCheck);

                moduleSupportRodSolid =
                    new G4Tubs(aCrystal.fName + "_SupportRodTub", 0, moduleSupportRodRadius,
                               aCrystal.fCrystalHeight / 2. - copperFrameSizeZ - smallBlock2SizeZ,
                               0, 360. * deg);

                moduleSupportRodLV = new G4LogicalVolume(moduleSupportRodSolid, nowFrameMaterial,
                                                         aCrystal.fName + "_SupportRod_LV");
                moduleSupportRodLV->SetVisAttributes(fI_CopperFrame_VisAttr);

                nowSupportRodDisplace = G4ThreeVector(copperFrameSizeX / 2. - smallBlock2DistX,
                                                      copperFrameSizeY / 2. - smallBlock2DistY, 0);
                for (G4int i = 0; i < 2; i++) {
                    for (G4int j = 0; j < 2; j++) {
                        new G4PVPlacement(nullptr, nowSupportRodDisplace, moduleSupportRodLV,
                                          aCrystal.fName + "_SupportRod_PV", resultLV, false,
                                          i * 2 + j, OverlapCheck);
                        nowSupportRodDisplace[0] *= -1.;
                    }
                    nowSupportRodDisplace[1] *= -1.;
                }

                nowCopperFramePos =
                    G4ThreeVector(0, 0, (aCrystal.fCrystalHeight - copperFrameSizeZ) / 2.);

                new G4PVPlacement(nullptr, nowCopperFramePos, Build_I_CopperFrameUpper(1),
                                  aCrystal.fName + "_UpperCuFrame_PV", resultLV, false, 0, OverlapCheck);
                new G4PVPlacement(nullptr, -nowCopperFramePos, Build_I_CopperFrameBottom(1),
                                  aCrystal.fName + "_BottomCuFrame_PV", resultLV, false, 0, OverlapCheck);
                Place_I_GeometriesInDetModuleOf(resultLV, 1, aCrystal, !aIsFloor);

                new G4PVPlacement(
                    nullptr, {0, 0, (aCrystal.fCrystalHeight + photonFrameSizeZ) / 2.},
                    Build_I_PhotonDetector(false,1, nowFrameMaterial, aSpaceMaterial, aMSDForThisModule),
                    aCrystal.fName + "_PhotonDetector_PV", resultLV, false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {0, 0, 0}, crystalLV, aCrystal.fName + "_Crystal_PV",
                                  resultLV, false, 0, OverlapCheck);

                if (nowReflectorMaterial != nullptr) {
                    reflectorLengthZ     = aCrystal.fCrystalHeight + reflectorThick;
                    reflectorCarverSizeZ = aCrystal.fCrystalHeight + solidBooleanTol * 2.;
                    reflectorSolid =
                        new G4Box(aCrystal.fName + "_ReflectorOutterBox", reflectorLengthX / 2.,
                                  reflectorLengthY / 2., reflectorLengthZ / 2.);
                    reflectorCarverSolid =
                        new G4Box(aCrystal.fName + "_ReflectorCarverBox", reflectorCerverSizeX / 2.,
                                  reflectorCerverSizeY / 2., reflectorCarverSizeZ / 2.);
                    reflectorSolid = new G4SubtractionSolid(
                        aCrystal.fName + "_ReflectorSolid", reflectorSolid, reflectorCarverSolid,
                        nullptr, {0, 0, reflectorThick / 2. + solidBooleanTol});

                    reflectorLV = new G4LogicalVolume(reflectorSolid, nowReflectorMaterial,
                                                      aCrystal.fName + "_Reflector_LV");
                    reflectorLV->SetVisAttributes(fI_Reflector_VisAttr);
                    new G4PVPlacement(nullptr, {0, 0, -reflectorThick / 2.}, reflectorLV,
                                      aCrystal.fName + "_Reflector_PV", resultLV, false, 0, OverlapCheck);
                } else {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_NOREFLECTOR",
                                G4ExceptionSeverity::JustWarning,
                                ("There is no matching material for a crystal (" + aCrystal.fName +
                                 "). The reflector will not be constructed.")
                                    .c_str());
                }
            } else {
                G4bool found = false;
                for (auto &nowSDInfo : fModuleSDInfos) {
                    if (nowSDInfo.fModulePV->GetLogicalVolume()->GetName() == nowEnvelopeLVName) {
                        aMSDForThisModule->fCrystalLV         = nowSDInfo.fCrystalLV;
                        aMSDForThisModule->fGeWaferLV         = nowSDInfo.fGeWaferLV;
                        aMSDForThisModule->fCrystalGoldFilmLV = nowSDInfo.fCrystalGoldFilmLV;
                        aMSDForThisModule->fGeWaferGoldFilmLV = nowSDInfo.fGeWaferGoldFilmLV;

                        found = true;
                        break;
                    }
                }
                if (found == false) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_HITINFO_ERROR",
                                G4ExceptionSeverity::RunMustBeAborted,
                                "LV is found but SD Hit info for that has not found.");
                }
            }
            break;
        case 2:
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(nowEnvelopeLVName, false);
            if (resultLV == nullptr) {
                G4ThreeVector nowCopperFramePos;
                nowCrystalMaterial   = G4Material::GetMaterial(aCrystal.fCrystalMaterialName);
                nowFrameMaterial     = G4Material::GetMaterial(aCrystal.fFrameMaterialName);
                nowReflectorMaterial = G4Material::GetMaterial(aCrystal.fReflectorMaterialName);
                if (nowCrystalMaterial == nullptr || nowFrameMaterial == nullptr) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_001",
                                G4ExceptionSeverity::FatalException,
                                "Material pointers for detector module are not found.");
                    nowFrameMaterial = nowCrystalMaterial = _vacuum;
                }
                moduleEnvelopeSolid = new G4Box(aCrystal.fName + "_DetModuleBox",
                                                copperFrameSizeX / 2. + type2_detModuleMarginX,
                                                copperFrameSizeY / 2. + type2_detModuleMarginY,
                                                aCrystal.fCrystalHeight / 2. + smallBlock1SizeZ);

                resultLV =
                    new G4LogicalVolume(moduleEnvelopeSolid, aSpaceMaterial, nowEnvelopeLVName);

                crystalSolid = new G4EllipticalTube(
                    aCrystal.fName + "_CrystalSolid", aCrystal.fCrystalMajorDiameter / 2.,
                    aCrystal.fCrystalMinorDiameter / 2., aCrystal.fCrystalHeight / 2.);
                crystalLV = new G4LogicalVolume(crystalSolid, nowCrystalMaterial,
                                                aCrystal.fName + "_Crystal_LV");
                crystalLV->SetVisAttributes(fI_Crystal_VisAttr);
                fI_CrystalLVs.insert(crystalLV);
                aMSDForThisModule->fCrystalLV = crystalLV;

                crystalBottomGoldFilmLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "GoldFilmBelowCrystal_Type2_LV", false);
                if (crystalBottomGoldFilmLV == nullptr) {
                    crystalBottomGoldFilmTub =
                        new G4Tubs("GoldFilmBelowCrystal_Tub", 0, crystalBottomGoldRadius,
                                   crystalBottomGoldThickness / 2., 0, 360. * deg);

                    crystalBottomGoldFilmLV = new G4LogicalVolume(crystalBottomGoldFilmTub, _gold,
                                                                  "GoldFilmBelowCrystal_Type2_LV");

                    crystalBottomGoldFilmLV->SetVisAttributes(fI_GoldDefault_VisAttr);

                    fI_CrystalGoldFilmLVs.insert(crystalBottomGoldFilmLV);
                }
                aMSDForThisModule->fCrystalGoldFilmLV = crystalBottomGoldFilmLV;

                new G4PVPlacement(
                    nullptr,
                    {0, 0, -aCrystal.fCrystalHeight / 2. + crystalBottomGoldThickness / 2.},
                    crystalBottomGoldFilmLV, aCrystal.fName + "_BottomGoldFilm_PV", crystalLV,
                    false, 0, OverlapCheck);

                moduleSupportRodSolid =
                    new G4Tubs(aCrystal.fName + "_SupportRodTub", 0, moduleSupportRodRadius,
                               aCrystal.fCrystalHeight / 2. - copperFrameSizeZ - smallBlock2SizeZ,
                               0, 360. * deg);

                moduleSupportRodLV = new G4LogicalVolume(moduleSupportRodSolid, nowFrameMaterial,
                                                         aCrystal.fName + "_SupportRod_LV");
                moduleSupportRodLV->SetVisAttributes(fI_CopperFrame_VisAttr);

                nowSupportRodDisplace = G4ThreeVector(copperFrameSizeX / 2. - smallBlock2DistX,
                                                      copperFrameSizeY / 2. - smallBlock2DistY, 0);
                for (G4int i = 0; i < 2; i++) {
                    for (G4int j = 0; j < 2; j++) {
                        new G4PVPlacement(nullptr, nowSupportRodDisplace, moduleSupportRodLV,
                                          aCrystal.fName + "_SupportRod_PV", resultLV, false,
                                          i * 2 + j, OverlapCheck);
                        nowSupportRodDisplace[0] *= -1.;
                    }
                    nowSupportRodDisplace[1] *= -1.;
                }

                nowCopperFramePos =
                    G4ThreeVector(0, 0, (aCrystal.fCrystalHeight - copperFrameSizeZ) / 2.);

                new G4PVPlacement(nullptr, nowCopperFramePos, Build_I_CopperFrameUpper(1),
                                  aCrystal.fName + "_UpperCuFrame_PV", resultLV, false, 0, OverlapCheck);
                new G4PVPlacement(nullptr, -nowCopperFramePos, Build_I_CopperFrameBottom(1),
                                  aCrystal.fName + "_BottomCuFrame_PV", resultLV, false, 0, OverlapCheck);
                Place_I_GeometriesInDetModuleOf(resultLV, 2, aCrystal, !aIsFloor);

                new G4PVPlacement(
                    nullptr, {0, 0, (aCrystal.fCrystalHeight + photonFrameSizeZ) / 2.},
                    Build_I_PhotonDetector(false,2, nowFrameMaterial, aSpaceMaterial, aMSDForThisModule),
                    aCrystal.fName + "_PhotonDetector_PV", resultLV, false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {0, 0, 0}, crystalLV, aCrystal.fName + "_Crystal_PV",
                                  resultLV, false, 0, OverlapCheck);

                if (nowReflectorMaterial != nullptr) {
                    reflectorLengthZ =
                        aCrystal.fCrystalHeight + crystalHeaterSpacing + reflectorThick;
                    reflectorCarverSizeZ =
                        aCrystal.fCrystalHeight + crystalHeaterSpacing + solidBooleanTol * 2.;
                    reflectorSolid =
                        new G4Box(aCrystal.fName + "_ReflectorOutterBox", reflectorLengthX / 2.,
                                  reflectorLengthY / 2., reflectorLengthZ / 2.);
                    reflectorCarverSolid =
                        new G4Box(aCrystal.fName + "_ReflectorCarverBox", reflectorCerverSizeX / 2.,
                                  reflectorCerverSizeY / 2., reflectorCarverSizeZ / 2.);
                    reflectorSolid = new G4SubtractionSolid(
                        aCrystal.fName + "_ReflectorSolid", reflectorSolid, reflectorCarverSolid,
                        nullptr, {0, 0, reflectorThick / 2. + solidBooleanTol});

                    reflectorLV = new G4LogicalVolume(reflectorSolid, nowReflectorMaterial,
                                                      aCrystal.fName + "_Reflector_LV");
                    reflectorLV->SetVisAttributes(fI_Reflector_VisAttr);
                    new G4PVPlacement(
                        nullptr, {0, 0, -reflectorThick / 2. - crystalHeaterSpacing / 2.},
                        reflectorLV, aCrystal.fName + "_Reflector_PV", resultLV, false, 0, OverlapCheck);
                } else {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_NOREFLECTOR",
                                G4ExceptionSeverity::JustWarning,
                                ("There is no matching material for a crystal (" + aCrystal.fName +
                                 "). The reflector will not be constructed.")
                                    .c_str());
                }
            } else {
                G4bool found = false;
                for (auto &nowSDInfo : fModuleSDInfos) {
                    if (nowSDInfo.fModulePV->GetLogicalVolume()->GetName() == nowEnvelopeLVName) {
                        aMSDForThisModule->fCrystalLV         = nowSDInfo.fCrystalLV;
                        aMSDForThisModule->fGeWaferLV         = nowSDInfo.fGeWaferLV;
                        aMSDForThisModule->fCrystalGoldFilmLV = nowSDInfo.fCrystalGoldFilmLV;
                        aMSDForThisModule->fGeWaferGoldFilmLV = nowSDInfo.fGeWaferGoldFilmLV;

                        found = true;
                        break;
                    }
                }
                if (found == false) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_HITINFO_ERROR",
                                G4ExceptionSeverity::RunMustBeAborted,
                                "LV is found but SD Hit info for that has not found.");
                }
            }
            break;
				case 3: // realistic geometry
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(nowEnvelopeLVName, false);
            if (resultLV == nullptr) 
						{
							nowCrystalMaterial   = G4Material::GetMaterial(aCrystal.fCrystalMaterialName);
              nowFrameMaterial     = G4Material::GetMaterial(aCrystal.fFrameMaterialName);
              nowReflectorMaterial = G4Material::GetMaterial(aCrystal.fReflectorMaterialName);
							G4Material *nowPhononCollectorMaterial 
									                 = G4Material::GetMaterial(aCrystal.fPhononCollectorMaterialName);
              if (nowCrystalMaterial == nullptr || nowFrameMaterial == nullptr) 
							{
								G4Exception(__PRETTY_FUNCTION__, "AMORE_I_001",
                                G4ExceptionSeverity::FatalException,
                                "Material pointers for detector module are not found.");
                nowFrameMaterial = nowCrystalMaterial = _vacuum;
              }

							// module envelope (MOTHER VOLUME) ----------------------------
              moduleEnvelopeSolid = new G4Box(aCrystal.fName + "_DetModuleBox",
                                                type3_copperFrameSizeX / 2. + type2_detModuleMarginX,
                                                type3_copperFrameSizeY / 2. + type2_detModuleMarginY,
																								type3_moduleSizeZ /2. );

              resultLV = new G4LogicalVolume(moduleEnvelopeSolid, aSpaceMaterial, nowEnvelopeLVName);

							// crystal -----------------------------------------------------
							double temp_diameter_M = aCrystal.fCrystalMajorDiameter/2.;
							double temp_diameter_m = aCrystal.fCrystalMinorDiameter/2.;
							double temp_density = nowCrystalMaterial->GetDensity();
							double crystal_mass = aCrystal.fCrystalMass-0.15*g;
							double temp_crystal_mass = 3.14 * temp_diameter_M * temp_diameter_m * aCrystal.fCrystalHeight * temp_density;
							double temp_factor = sqrt(crystal_mass/temp_crystal_mass);
							
              crystalSolid = new G4EllipticalTube(
                    aCrystal.fName + "_CrystalSolid", (aCrystal.fCrystalMajorDiameter / 2.)*temp_factor,
                    (aCrystal.fCrystalMinorDiameter / 2.)*temp_factor, aCrystal.fCrystalHeight / 2.);

              crystalLV = new G4LogicalVolume(crystalSolid, nowCrystalMaterial, aCrystal.fName + "_Crystal_LV");
              crystalLV->SetVisAttributes(fI_Crystal_VisAttr);
              fI_CrystalLVs.insert(crystalLV);
              aMSDForThisModule->fCrystalLV = crystalLV;

							fI_crystalsRegion->AddRootLogicalVolume(crystalLV);

							// Crystal Positioning -----------------------------------------
              new G4PVPlacement(nullptr, {0, 0, type3_moduleSizeZ/2. - type3_smallBlock1SizeZ - aCrystal.fCrystalHeight/2.}, 
										crystalLV, aCrystal.fName + "_Crystal_PV", resultLV, false, 0, OverlapCheck);
				
							if(fDbgMsgOn) {
								std::cout << aCrystal.fName << " mass: " << crystalLV->GetMass()/g << " g" << std::endl;
							}

							// Bottom Gold Film --------------------------------------------
              crystalBottomGoldFilmLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "GoldFilmBelowCrystal_Type3_LV", false);
              if (crystalBottomGoldFilmLV == nullptr) 
							{
								crystalBottomGoldFilmTub = new G4Tubs("GoldFilmBelowCrystal_Tub", 0, 
										crystalBottomGoldRadius, crystalBottomGoldThickness / 2., 0, 360. * deg);

                crystalBottomGoldFilmLV = new G4LogicalVolume( crystalBottomGoldFilmTub, 
										nowPhononCollectorMaterial, "GoldFilmBelowCrystal_Type3_LV");
                crystalBottomGoldFilmLV->SetVisAttributes(fI_GoldDefault_VisAttr);
                fI_CrystalGoldFilmLVs.insert(crystalBottomGoldFilmLV);
							}
              aMSDForThisModule->fCrystalGoldFilmLV = crystalBottomGoldFilmLV;

              new G4PVPlacement( nullptr,
                    {0, 0, -aCrystal.fCrystalHeight / 2. + crystalBottomGoldThickness / 2.},
                    crystalBottomGoldFilmLV, aCrystal.fName + "_BottomGoldFilm_PV", crystalLV,
                    false, 0, OverlapCheck);

							// Mudule Support Rod ---------------------------------------------
              moduleSupportRodSolid = new G4Tubs(aCrystal.fName + "_SupportRodTub", 0, 
									type3_moduleSupportRodRadius, 
									aCrystal.fCrystalHeight / 2. - type3_copperFrameSizeZ, 0, 360. * deg);

              moduleSupportRodLV = new G4LogicalVolume( moduleSupportRodSolid, 
									nowFrameMaterial, aCrystal.fName + "_SupportRod_LV");
              moduleSupportRodLV->SetVisAttributes(fI_CopperFrame_VisAttr);

              nowSupportRodDisplace = G4ThreeVector(type3_copperFrameSizeX / 2. - smallBlock2DistX,
                              type3_copperFrameSizeY / 2. - smallBlock2DistY, 
															type3_moduleSizeZ/2 - type3_smallBlock1SizeZ - aCrystal.fCrystalHeight /2. );
              for (G4int i = 0; i < 2; i++) 
							{
								for (G4int j = 0; j < 2; j++) 
									{
                    new G4PVPlacement(nullptr, nowSupportRodDisplace, moduleSupportRodLV,
												aCrystal.fName + "_SupportRod_PV", resultLV, false, i * 2 + j, OverlapCheck);
                    nowSupportRodDisplace[0] *= -1.;
                  }
                nowSupportRodDisplace[1] *= -1.;
              }

							// Upper Copper Frame ------------------------------------------
              new G4PVPlacement(nullptr, 
									{0,0,type3_moduleSizeZ/2.-type3_smallBlock1SizeZ - type3_copperFrameSizeZ/2.},
									Build_I_CopperFrameUpper(3), aCrystal.fName + "_UpperCuFrame_PV", 
									resultLV, false, 0, OverlapCheck);

							// Lower Copper Frame ------------------------------------------
              new G4PVPlacement(nullptr, 
									{0,0,-type3_moduleSizeZ/2 + type3_smallBlock3SizeZ + type3_copperFrameSizeZ/2.},
									Build_I_CopperFrameBottom(3), aCrystal.fName + "_BottomCuFrame_PV", 
									resultLV, false, 0, OverlapCheck);

							// PhotonDetector and Photon PCBs ------------------------------
              new G4PVPlacement( nullptr, 
									{0, 0, type3_moduleSizeZ/2. - type3_smallBlock1SizeZ + 
									  (type3_photonFrameSizeZ + type3_photonFrameTopSpace)/ 2.},
                  Build_I_PhotonDetector(aRealistic, 3, nowFrameMaterial, aSpaceMaterial, aMSDForThisModule),
                  aCrystal.fName + "_PhotonDetector_PV", resultLV, false, 0, OverlapCheck);

							// Phonon detector and PCBs ------------------------------------ 
							Place_I_GeometriesInDetModuleOf(resultLV, 3, aCrystal, !aIsFloor);

							if(aRealistic)
							{ // Bolts and Clamps
								Place_I_GeometriesInDetModuleOf(resultLV, 4, aCrystal, !aIsFloor);

								// Bolts for Post
								Place_I_GeometriesInDetModuleOf(resultLV, 5, aCrystal, !aIsFloor);
							}

							// Reflector ---------------------------------------------------
              if (nowReflectorMaterial != nullptr) 
							{
								reflectorLengthZ = aCrystal.fCrystalHeight + crystalHeaterSpacing;
								double type3_reflectorLengthX = aCrystal.fCrystalMajorDiameter * temp_factor + reflectorThick;
								double type3_reflectorLengthY = aCrystal.fCrystalMinorDiameter * temp_factor + reflectorThick;

								G4VSolid *reflectorBottom = new G4Box(aCrystal.fName + "_ReflectorBottomBox", 
										type3_reflectorLengthX / 2., type3_reflectorLengthY / 2., reflectorThick / 2.);

								G4Tubs *reflectorHole = new G4Tubs(aCrystal.fName + "_ReflectorHoleTub",
										0, 10., reflectorThick + solidBooleanTol*2., 0,360. * deg);

								G4VSolid *reflectorBottomSolid = new G4SubtractionSolid(
												aCrystal.fName + "_ReflectorBottomSolid", 
												reflectorBottom, reflectorHole, nullptr, {0, 0, 0});

								reflectorSolid = new G4EllipticalTube(aCrystal.fName+"_ReflectorOutterTub",
										type3_reflectorLengthX / 2., type3_reflectorLengthY / 2., reflectorLengthZ / 2.);
								reflectorCarverSolid = new G4EllipticalTube(aCrystal.fName+"_ReflectorCarverTub",
										(type3_reflectorLengthX - reflectorThick) / 2.,
										(type3_reflectorLengthY - reflectorThick) / 2.,
										reflectorLengthZ);

                reflectorSolid = new G4SubtractionSolid( aCrystal.fName + "_ReflectorBarrelSolid", 
										reflectorSolid, reflectorCarverSolid, nullptr, {0, 0, 0});
								reflectorSolid = new G4UnionSolid( aCrystal.fName + "_ReflectorSolid",
										reflectorSolid, reflectorBottomSolid,	nullptr, 
										{0, 0, -reflectorLengthZ / 2. - reflectorThick / 2.});

                reflectorLV = new G4LogicalVolume(reflectorSolid, 
										nowReflectorMaterial, aCrystal.fName + "_Reflector_LV");
                reflectorLV->SetVisAttributes(fI_Reflector_VisAttr);

                new G4PVPlacement( nullptr, 
										{0, 0, type3_moduleSizeZ/2 -type3_smallBlock1SizeZ - aCrystal.fCrystalHeight/2.
												-reflectorThick / 2. - crystalHeaterSpacing / 2.},
                    reflectorLV, aCrystal.fName + "_Reflector_PV", resultLV, false, 0, OverlapCheck);
              } else {
                G4Exception(__PRETTY_FUNCTION__, "AMORE_I_NOREFLECTOR",
                              G4ExceptionSeverity::JustWarning,
                    ("There is no matching material for a crystal (" + aCrystal.fName +
                          "). The reflector will not be constructed.").c_str());
              }
            } else {
              G4bool found = false;
							for (auto &nowSDInfo : fModuleSDInfos) 
							{
                if (nowSDInfo.fModulePV->GetLogicalVolume()->GetName() == nowEnvelopeLVName) 
								{
                  aMSDForThisModule->fCrystalLV         = nowSDInfo.fCrystalLV;
                  aMSDForThisModule->fGeWaferLV         = nowSDInfo.fGeWaferLV;
                  aMSDForThisModule->fCrystalGoldFilmLV = nowSDInfo.fCrystalGoldFilmLV;
                  aMSDForThisModule->fGeWaferGoldFilmLV = nowSDInfo.fGeWaferGoldFilmLV;
                  found = true;
                  break;
                }
							}

              if (found == false) 
							{
                G4Exception(__PRETTY_FUNCTION__, "AMORE_I_HITINFO_ERROR",
                                G4ExceptionSeverity::RunMustBeAborted,
                                "LV is found but SD Hit info for that has not found.");
              }
            }
          break;
      default:
        resultLV = nullptr;
        G4Exception(__PRETTY_FUNCTION__, "AMORE_I_TYPEERROR",
                        G4ExceptionSeverity::FatalErrorInArgument,
                        "Wrong type for a single CMO module");
          break;
    }
    return resultLV;
}

// Written by Jeewon Seo @ 2021.03.01.
// This function for constructing the PCB for photon detector
// There are realistic option. 
// If realistic option is turned off, simple shape of PCB will be installed.
// There are two cases for bolts ( Head of bolts is in ....)
// case1: w/  bolts
// case2: w/o bolts
void AmoreDetectorConstruction::Place_I_PCBs(G4bool aRealistic, G4LogicalVolume *aWhere, 
		                                          G4int aType, G4ThreeVector aTlate){

		using namespace AmoreDetectorStaticInfo;
		using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4VSolid *photonDetectorPCBTopSolid = nullptr;
		G4VSolid *photonDetectorPCBBaseSolid = nullptr;
		G4VSolid *photonDetectorPCBInnerSolid = nullptr;
		G4VSolid *photonDetectorPCBSolderSolid = nullptr;
		G4VSolid *photonDetectorPCBSolid = nullptr;
		G4Box *photonDetectorPCBHole = nullptr;

		G4LogicalVolume *photonDetectorPCBLV;
		G4LogicalVolume *photonDetectorPCBSolderLV = nullptr;

		G4ThreeVector nowPCBPosition;
		G4ThreeVector nowSolderPosition;
		G4ThreeVector nowHolePosition = {5,-5,0};

		G4TwoVector nowVertexDisplacer;
		std::vector<G4TwoVector> photonDetectorPCBVertices;

		G4double type3_PCBthick = type3_photDetPCBBaseBoxThick-type3_photDetPCBBaseInnerThick;

		// For case 1 only
		G4ThreeVector nowBoltsPosition = {type3_photDetPCBBaseBoxSizeX / 2. - photDetPCBBoltsHoleDistX, 
										                  type3_photDetPCBBaseBoxSizeY / 2. - photDetPCBBoltsHoleDistY,  
																			0};

		G4LogicalVolume *logicalVolumeChecker = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV",false);
		switch (aType) { 
			case 1: { // photon Detector PCB with bolts
				if(logicalVolumeChecker == nullptr){
					// PCB Brass bolts
					G4Tubs *photonDetectorPCBBoltsHole = new G4Tubs("PCBBoltsHoleForPhotonDetectorTubs",0,
							photDetPCBBoltsHoleRadius, photDetPCBBoltsHoleHeight / 2 + solidBooleanTol, 0, 360. * deg);

					// PCB For photon detector ----------------
          photonDetectorPCBVertices.resize(6);
					photonDetectorPCBVertices[0] = {type3_photDetPCBBaseBoxSizeX / 2.,
																					type3_photDetPCBBaseBoxSizeY / 2.};

          nowVertexDisplacer = {type3_photDetPCBBaseBoxSizeX1,0};
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[0] - nowVertexDisplacer;

          nowVertexDisplacer = {0,type3_photDetPCBBaseBoxSizeY1};
					photonDetectorPCBVertices[2] = photonDetectorPCBVertices[1] - nowVertexDisplacer;

					photonDetectorPCBVertices[4] = {-type3_photDetPCBBaseBoxSizeX / 2.,
																					-type3_photDetPCBBaseBoxSizeY / 2.};

          nowVertexDisplacer = {0,type3_photDetPCBBaseBoxSizeY2};
					photonDetectorPCBVertices[3] = photonDetectorPCBVertices[4] + nowVertexDisplacer;

					photonDetectorPCBVertices[5] = {type3_photDetPCBBaseBoxSizeX / 2.,
																					-type3_photDetPCBBaseBoxSizeY / 2.};

          photonDetectorPCBTopSolid =
              new G4ExtrudedSolid("PCBForPhotonDetectorTopSolid", photonDetectorPCBVertices,
                                  type3_PCBthick / 2., 0, 1., 0, 1.);

					photonDetectorPCBTopSolid = new G4SubtractionSolid("Type3_PhotonFrameUpper_BoltsHole1",
											photonDetectorPCBTopSolid, photonDetectorPCBBoltsHole, nullptr, nowBoltsPosition);

					nowBoltsPosition[1] *= -1.;
					photonDetectorPCBTopSolid = new G4SubtractionSolid("Type3_PhotonFrameUpper_BoltsHole2",
											photonDetectorPCBTopSolid, photonDetectorPCBBoltsHole, nullptr, nowBoltsPosition);

					nowBoltsPosition[0] *= -1.;
					photonDetectorPCBTopSolid = new G4SubtractionSolid("Type3_PhotonFrameUpper_BoltsHole2",
											photonDetectorPCBTopSolid, photonDetectorPCBBoltsHole, nullptr, nowBoltsPosition);

					// PCB inner space for photon detector -----------
					photonDetectorPCBVertices.resize(7);
					photonDetectorPCBVertices[0] = {type3_photDetPCBBaseBoxSizeX / 2. -
																						type3_photDetPCBBaseBoxSizeX1,
																					type3_photDetPCBBaseBoxSizeY / 2. -
																						type3_photDetPCBBaseBoxSizeY1 + solidBooleanTol};
					photonDetectorPCBVertices[2] = {-type3_photDetPCBBaseBoxSizeX / 2. +
									                            type3_photDetPCBBaseBoxSizeX1,
																					-type3_photDetPCBBaseBoxSizeY / 2.-solidBooleanTol};
					photonDetectorPCBVertices[3] = {type3_photDetPCBBaseBoxSizeX / 2. -
									                            type3_photDetPCBBaseBoxSizeX1,
									                        -type3_photDetPCBBaseBoxSizeY / 2.-solidBooleanTol};
					photonDetectorPCBVertices[5] = {type3_photDetPCBBaseBoxSizeX / 2. -
									                            type3_photDetPCBBaseInnerSizeX1,
																					-type3_photDetPCBBaseBoxSizeY / 2. +
																							type3_photDetPCBBaseInnerSizeY3};
					photonDetectorPCBVertices[6] = { type3_photDetPCBBaseBoxSizeX / 2. -
									                            type3_photDetPCBBaseInnerSizeX1,
																						type3_photDetPCBBaseBoxSizeY / 2. -
																							type3_photDetPCBBaseInnerSizeY4}; 

					nowVertexDisplacer = {0, type3_photDetPCBBaseInnerSizeY1+solidBooleanTol};
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[2] + nowVertexDisplacer;

					nowVertexDisplacer = {0, type3_photDetPCBBaseInnerSizeY2};
					photonDetectorPCBVertices[4] = photonDetectorPCBVertices[3] + nowVertexDisplacer;

					photonDetectorPCBInnerSolid = new G4ExtrudedSolid("PCBForPhotonDetectorInner", 
							                               photonDetectorPCBVertices,
												                      type3_photDetPCBBaseInnerThick / 2., 0, 1., 0, 1.);

					// PCB solder ---------------------------
					nowVertexDisplacer = {0, solidBooleanTol};
					photonDetectorPCBVertices[3] = photonDetectorPCBVertices[3] + nowVertexDisplacer;
					photonDetectorPCBVertices[2] = photonDetectorPCBVertices[2] + nowVertexDisplacer;
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[1] - nowVertexDisplacer;
					photonDetectorPCBVertices[0] = photonDetectorPCBVertices[0] - nowVertexDisplacer;
          photonDetectorPCBSolderSolid = new G4ExtrudedSolid("SolderForPhotonDetector", 
							                                photonDetectorPCBVertices,
							                                 photDetPCBThick / 2., 0, 1., 0, 1.);

					// PCB Hole -----------------------------
					photonDetectorPCBHole = new G4Box ("PCBForPhotonDetectorHoleBox", photDetPCBHoleSizeX / 2.,
							                                photDetPCBHoleSizeY / 2., photDetPCBHoleSizeZ / 2.);
					photonDetectorPCBSolderSolid = new G4SubtractionSolid("SolderForPhotonDetectorSolid",
							                                photonDetectorPCBSolderSolid, photonDetectorPCBHole, 
																							nullptr, nowHolePosition);

					// PCB base for photon detector ----------
					nowVertexDisplacer = {type3_photDetPCBBaseInnerSizeX1-solidBooleanTol,0};
					photonDetectorPCBVertices[6] = photonDetectorPCBVertices[6] + nowVertexDisplacer;
					photonDetectorPCBVertices[5] = photonDetectorPCBVertices[5] + nowVertexDisplacer;
					photonDetectorPCBVertices[4] = photonDetectorPCBVertices[4] + nowVertexDisplacer;
					photonDetectorPCBVertices[3] = photonDetectorPCBVertices[3] + nowVertexDisplacer;
					photonDetectorPCBVertices[2] = photonDetectorPCBVertices[2] - nowVertexDisplacer;

					nowVertexDisplacer = {0, type3_photDetPCBBaseInnerSizeY5};
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[2] + nowVertexDisplacer;

					photonDetectorPCBBaseSolid = 
							new G4ExtrudedSolid("PCBForPhotonDetectorBaseSolid", photonDetectorPCBVertices,
												            type3_photDetPCBBaseInnerThick / 2., 0, 1., 0, 1.);

					// Making one volume of PCB for photon detector--------
					photonDetectorPCBSolid = new G4UnionSolid("PCBForPhotonDetectorSolid1",
							photonDetectorPCBTopSolid, photonDetectorPCBBaseSolid, nullptr,
							{0,0, -type3_PCBthick/2. - type3_photDetPCBBaseInnerThick/2.});

					photonDetectorPCBSolid = new G4SubtractionSolid("PCBForPhotonDetectorSolid2",
							photonDetectorPCBSolid, photonDetectorPCBInnerSolid, nullptr, 
							{0,0, type3_PCBthick/2. - type3_photDetPCBBaseInnerThick/2.});

					photonDetectorPCBSolid = new G4SubtractionSolid("PCBForPhotonDetectorSolid",
							                                photonDetectorPCBSolid, photonDetectorPCBHole,
																							nullptr, nowHolePosition);

					if(aRealistic){
						photonDetectorPCBLV = new G4LogicalVolume(photonDetectorPCBSolid, 
								                                       _copper, "PCBForPhotonDetector_LV");
						photonDetectorPCBSolderLV = new G4LogicalVolume(photonDetectorPCBSolderSolid,
						                                                 _solder, "PhotonDetectorPCBSolder_LV");
						photonDetectorPCBSolderLV->SetVisAttributes(fI_LeadDefault_VisAttr);

					} else {
						photonDetectorPCBLV = new G4LogicalVolume(photonDetectorPCBTopSolid, 
								                                       _copper, "PCBForPhotonDetector_LV");
					}
          photonDetectorPCBLV->SetVisAttributes(fI_PCB_VisAttr);

				} else {
					photonDetectorPCBLV = 
						G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV");
				}

				nowPCBPosition = { 0, 0, - type3_PCBthick / 2.};
				nowPCBPosition += aTlate;
				new G4PVPlacement(nullptr, nowPCBPosition, photonDetectorPCBLV,
						              "PCBForPhotonDetector_PV", aWhere, false, 0, OverlapCheck);

				nowSolderPosition = { 0, 0, - type3_photDetPCBBaseBoxThick + type3_PCBthick + photDetPCBThick/2.};
				new G4PVPlacement(nullptr, aTlate + nowSolderPosition, photonDetectorPCBSolderLV,
						              "PhotonDetectorPCBSolder_PV", aWhere, false, 0, OverlapCheck);

			break;
			}
			case 2: { // photon Detector PCB without bolts
				if(logicalVolumeChecker == nullptr){
					// PCB For photon detector ----------------
          photonDetectorPCBVertices.resize(6);
					photonDetectorPCBVertices[0] = {type3_photDetPCBBaseBoxSizeX / 2.,
																					type3_photDetPCBBaseBoxSizeY / 2.};

          nowVertexDisplacer = {type3_photDetPCBBaseBoxSizeX1,0};
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[0] - nowVertexDisplacer;

          nowVertexDisplacer = {0,type3_photDetPCBBaseBoxSizeY1};
					photonDetectorPCBVertices[2] = photonDetectorPCBVertices[1] - nowVertexDisplacer;

					photonDetectorPCBVertices[4] = {-type3_photDetPCBBaseBoxSizeX / 2.,
																					-type3_photDetPCBBaseBoxSizeY / 2.};

          nowVertexDisplacer = {0,type3_photDetPCBBaseBoxSizeY2};
					photonDetectorPCBVertices[3] = photonDetectorPCBVertices[4] + nowVertexDisplacer;

					photonDetectorPCBVertices[5] = {type3_photDetPCBBaseBoxSizeX / 2.,
																					-type3_photDetPCBBaseBoxSizeY / 2.};

          photonDetectorPCBTopSolid =
              new G4ExtrudedSolid("PCBForPhotonDetectorTopSolid", photonDetectorPCBVertices,
                                  type3_PCBthick / 2., 0, 1., 0, 1.);

					// PCB inner space for photon detector -----------
					photonDetectorPCBVertices.resize(7);
					photonDetectorPCBVertices[0] = {type3_photDetPCBBaseBoxSizeX / 2. -
																						type3_photDetPCBBaseBoxSizeX1,
																					type3_photDetPCBBaseBoxSizeY / 2. -
																						type3_photDetPCBBaseBoxSizeY1 + solidBooleanTol};
					photonDetectorPCBVertices[2] = {-type3_photDetPCBBaseBoxSizeX / 2. +
									                            type3_photDetPCBBaseBoxSizeX1,
																					-type3_photDetPCBBaseBoxSizeY / 2.-solidBooleanTol};
					photonDetectorPCBVertices[3] = {type3_photDetPCBBaseBoxSizeX / 2. -
									                            type3_photDetPCBBaseBoxSizeX1,
									                        -type3_photDetPCBBaseBoxSizeY / 2.-solidBooleanTol};
					photonDetectorPCBVertices[5] = {type3_photDetPCBBaseBoxSizeX / 2. -
									                            type3_photDetPCBBaseInnerSizeX1,
																					-type3_photDetPCBBaseBoxSizeY / 2. +
																							type3_photDetPCBBaseInnerSizeY3};
					photonDetectorPCBVertices[6] = { type3_photDetPCBBaseBoxSizeX / 2. -
									                            type3_photDetPCBBaseInnerSizeX1,
																						type3_photDetPCBBaseBoxSizeY / 2. -
																							type3_photDetPCBBaseInnerSizeY4}; 

					nowVertexDisplacer = {0, type3_photDetPCBBaseInnerSizeY1+solidBooleanTol};
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[2] + nowVertexDisplacer;

					nowVertexDisplacer = {0, type3_photDetPCBBaseInnerSizeY2};
					photonDetectorPCBVertices[4] = photonDetectorPCBVertices[3] + nowVertexDisplacer;

					photonDetectorPCBInnerSolid = new G4ExtrudedSolid("PCBForPhotonDetectorInner", 
							                               photonDetectorPCBVertices,
												                      type3_photDetPCBBaseInnerThick / 2., 0, 1., 0, 1.);

					// PCB solder ---------------------------
					nowVertexDisplacer = {0, solidBooleanTol};
					photonDetectorPCBVertices[3] = photonDetectorPCBVertices[3] + nowVertexDisplacer;
					photonDetectorPCBVertices[2] = photonDetectorPCBVertices[2] + nowVertexDisplacer;
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[1] - nowVertexDisplacer;
					photonDetectorPCBVertices[0] = photonDetectorPCBVertices[0] - nowVertexDisplacer;
          photonDetectorPCBSolderSolid = new G4ExtrudedSolid("SolderForPhotonDetector", 
							                                photonDetectorPCBVertices,
							                                 photDetPCBThick / 2., 0, 1., 0, 1.);

					// PCB Hole -----------------------------
					photonDetectorPCBHole = new G4Box ("PCBForPhotonDetectorHoleBox", photDetPCBHoleSizeX / 2.,
							                                photDetPCBHoleSizeY / 2., photDetPCBHoleSizeZ / 2.);
					photonDetectorPCBSolderSolid = new G4SubtractionSolid("SolderForPhotonDetectorSolid",
							                                photonDetectorPCBSolderSolid, photonDetectorPCBHole, 
																							nullptr, nowHolePosition);

					// PCB base for photon detector ----------
					nowVertexDisplacer = {type3_photDetPCBBaseInnerSizeX1-solidBooleanTol,0};
					photonDetectorPCBVertices[6] = photonDetectorPCBVertices[6] + nowVertexDisplacer;
					photonDetectorPCBVertices[5] = photonDetectorPCBVertices[5] + nowVertexDisplacer;
					photonDetectorPCBVertices[4] = photonDetectorPCBVertices[4] + nowVertexDisplacer;
					photonDetectorPCBVertices[3] = photonDetectorPCBVertices[3] + nowVertexDisplacer;
					photonDetectorPCBVertices[2] = photonDetectorPCBVertices[2] - nowVertexDisplacer;

					nowVertexDisplacer = {0, type3_photDetPCBBaseInnerSizeY5};
					photonDetectorPCBVertices[1] = photonDetectorPCBVertices[2] + nowVertexDisplacer;

					photonDetectorPCBBaseSolid = 
							new G4ExtrudedSolid("PCBForPhotonDetectorBaseSolid", photonDetectorPCBVertices,
												            type3_photDetPCBBaseInnerThick / 2., 0, 1., 0, 1.);

					// Making one volume of PCB for photon detector--------
					photonDetectorPCBSolid = new G4UnionSolid("PCBForPhotonDetectorSolid1",
							photonDetectorPCBTopSolid, photonDetectorPCBBaseSolid, nullptr,
							{0,0, -type3_PCBthick/2. - type3_photDetPCBBaseInnerThick/2.});

					photonDetectorPCBSolid = new G4SubtractionSolid("PCBForPhotonDetectorSolid2",
							photonDetectorPCBSolid, photonDetectorPCBInnerSolid, nullptr, 
							{0,0, type3_PCBthick/2. - type3_photDetPCBBaseInnerThick/2.});

					photonDetectorPCBSolid = new G4SubtractionSolid("PCBForPhotonDetectorSolid",
							                                photonDetectorPCBSolid, photonDetectorPCBHole,
																							nullptr, nowHolePosition);

					if(aRealistic){
						photonDetectorPCBLV = new G4LogicalVolume(photonDetectorPCBSolid, 
								                                       _copper, "PCBForPhotonDetector_LV");
						photonDetectorPCBSolderLV = new G4LogicalVolume(photonDetectorPCBSolderSolid,
						                                                 _solder, "PhotonDetectorPCBSolder_LV");
						photonDetectorPCBSolderLV->SetVisAttributes(fI_LeadDefault_VisAttr);
					} else {
						photonDetectorPCBLV = new G4LogicalVolume(photonDetectorPCBTopSolid, 
								                                       _copper, "PCBForPhotonDetector_LV");
					}
          photonDetectorPCBLV->SetVisAttributes(fI_PCB_VisAttr);

				} else {
					photonDetectorPCBLV = 
						G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV");
				}

				nowPCBPosition = { 0, 0, - type3_PCBthick / 2.};
				nowPCBPosition += aTlate;
				new G4PVPlacement(nullptr, nowPCBPosition, photonDetectorPCBLV,
						              "PCBForPhotonDetector_PV", aWhere, false, 0, OverlapCheck);

				nowSolderPosition = { 0, 0, - type3_photDetPCBBaseBoxThick + type3_PCBthick + photDetPCBThick/2.};
				new G4PVPlacement(nullptr, aTlate + nowSolderPosition, photonDetectorPCBSolderLV,
						              "PhotonDetectorPCBSolder_PV", aWhere, false, 0, OverlapCheck);
			break;
			}
		}
}

// CopperFrame by Daehoon, HA!!! and Mona /////////////////////////
// Modified by Jeewon at 2021.Mar.--------------------------------
// case 3 added for real AMORE-I geometry
// with realistic option,
//      holes for bolts and clamps will be included.
// ---------------------------------------------------------------
G4LogicalVolume *
    AmoreDetectorConstruction::Build_I_PhotonDetector(G4bool aRealistic, G4int aType, G4Material *aFrameMaterial,
                                                      G4Material *aSpaceMaterial,
                                                      AmoreModuleSDInfo *aMSDForThisModule) {

    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    G4LogicalVolume *resultLV;
    G4String lvName;

    G4Box *photonFrameBase;
    G4Tubs *photonFrameHole;
    G4VSolid *photonFrame;
    G4VSolid *photonFrameHoleAndSlot;

    G4Tubs *germaniumWaferTub;
    G4Tubs *geWaferUpperGoldFilmTub;

    G4LogicalVolume *holeAndSlotLV;
    G4LogicalVolume *germaniumWaferLV;
    G4LogicalVolume *geWaferUpperGoldFilmLV;

    G4ThreeVector nowUpperGoldFilmPos;

		// For Type1,
    G4double type1_photonFrameHoleLength = photonFrameSizeZ - type1_photonFrameBaseThick;
    std::vector<G4TwoVector> type1_holderVertices;

		// For Type2,
    G4double type2_photonFrameHoleLength =
        photonFrameSizeZ - type2_photonFrameBaseThick - type2_photonFrameRoofThick;
    G4double type2_photonFrameRoofHoleRadius =
        type2_photonFrameHoleRadius - type2_photonFrameRoofWidth;
    G4double type2_photonFrameHeartLeapMajorRadius = type2_photonFrameHoleRadius / 2.;
    G4double type2_PEEKClampRadialDistance =
        type2_photonFrameHoleRadius - type2_photonFramePeekRdaius / 2. +
        type2_photonFrameClampBox2SizeY + type2_photonFrameClampBox1SizeXY / 2.;

    G4ThreeVector type2_photonFrameHeartLeapPosition;
    G4ThreeVector type2_photonFrameBaseSupportNowPosition;
    G4ThreeVector type2_PEEKNowPosition;
    G4ThreeVector type2_PEEKClampNowPosition;
    std::vector<G4TwoVector> type2_PEEKClampVertices;

    G4RotationMatrix *type2_photonFrameHeartLeapRotMtx;
    G4RotationMatrix *type2_PEEKClampRotMtx[4];

		// For Type3,
		G4double type3_photonFrameHoleBoxSizeY = type3_photonFrameSizeY /2.;
		G4double type3_photonFrameHoleBoxPositionY = -type3_photonFrameHoleBoxSizeY / 2.;
		G4RotationMatrix *type3_photonFrameSupportRotMtx[4];

    switch (aType) {
        case 1:
						if(fDbgMsgOn) {
							G4cout << __PRETTY_FUNCTION__
                   << ": Realistic Photon frame type 1 selected (It could be very slow)" << G4endl;
						}
            lvName   = "PhotonDetector_Type1_LV";
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(lvName, false);
            if (resultLV == nullptr) {
                photonFrameBase = new G4Box("PhotonFrame_BaseBox", photonFrameSizeX / 2.,
                                            photonFrameSizeY / 2., photonFrameSizeZ / 2.);

                G4cout << __PRETTY_FUNCTION__ << ": Carving the well at the base ..." << G4endl;
                G4Box *photonFrameBaseCarver = new G4Box(
                    "Type1_PhotonFrame_BaseCarverXBox", copperFrameSizeX / 2 + solidBooleanTol,
                    copperFrameSizeY / 2 - type1_photonFrameBaseWidth, type1_photonFrameBaseThick);
                photonFrame = new G4SubtractionSolid("Type1_PhotonFrame_CarvingStage1Solid",
                                                     photonFrameBase, photonFrameBaseCarver,
                                                     nullptr, {0, 0, -photonFrameSizeZ / 2.});

                photonFrameBaseCarver = new G4Box(
                    "Type1_PhotonFrame_BaseCarverYBox", type1_photonFrameBaseSpaceWidth / 2,
                    photonFrameSizeY / 2. + solidBooleanTol, type1_photonFrameBaseThick);
                photonFrame = new G4SubtractionSolid("Type1_PhotonFrame_CarvingStage2Solid",
                                                     photonFrame, photonFrameBaseCarver, nullptr,
                                                     {0, 0, -photonFrameSizeZ / 2.});
                resultLV    = new G4LogicalVolume(photonFrame, aFrameMaterial, lvName);
                resultLV->SetVisAttributes(fI_OpticalFrame_VisAttr);

                G4cout << __PRETTY_FUNCTION__ << ": Assembling holder for photon detector ..."
                       << G4endl;

                type1_holderVertices.resize(9);
                type1_holderVertices[0] = {waferHolderWidth / 2., waferHolderLength};
                type1_holderVertices[3] = type1_holderVertices[6] = type1_holderVertices[0];
                type1_holderVertices[1] = {-waferHolderWidth / 2., waferHolderLength};
                type1_holderVertices[4] = type1_holderVertices[7] = type1_holderVertices[1];
                type1_holderVertices[3].rotate(120. * deg);
                type1_holderVertices[4].rotate(120. * deg);
                type1_holderVertices[6].rotate(240. * deg);
                type1_holderVertices[7].rotate(240. * deg);

                type1_holderVertices[2] = {-waferHolderWidth / 2.,
                                           (waferHolderWidth / 2.) * sin(120 * deg) +
                                               (waferHolderWidth / 2.) * (1 + cos(120 * deg)) *
                                                   tan((90 - 120) * deg)};
                type1_holderVertices[5] = type1_holderVertices[8] = type1_holderVertices[2];
                type1_holderVertices[5].rotate(120. * deg);
                type1_holderVertices[8].rotate(240. * deg);

                G4ThreeVector type1_waferHolderTlate =
                    G4ThreeVector(0, 0, photonFrameSizeZ / 2. - waferHolderThick / 2.);
                G4VSolid *type1_waferHolder =
                    new G4ExtrudedSolid("PhotonFrame_WaferHolderBaseSolid", type1_holderVertices,
                                        waferHolderThick / 2., 0, 1., 0, 1.);

                G4cout << __PRETTY_FUNCTION__ << ": Punching hole at the center..." << G4endl;
                photonFrameHole =
                    new G4Tubs("PhotonFrame_MainHoleTub", 0, type1_photonFrameHoleRadius,
                               type1_photonFrameHoleLength / 2., 0, 360 * deg);
                G4Tubs *type1_photonDetectorSlot = new G4Tubs(
                    "PhotonFrame_SlotTub", 0, germaniumWaferRadius + type1_germaniumWaferSlotDepth,
                    type1_germaniumWaferSlotThick / 2., 0, 360 * deg);
                photonFrameHoleAndSlot = new G4UnionSolid(
                    "PhotonFrame_HoleAndSlot", photonFrameHole, type1_photonDetectorSlot);

                type1_waferHolder = new G4IntersectionSolid("PhotonFrame_WaferHolderSolid",
                                                            type1_waferHolder, photonFrameHole);

                G4LogicalVolume *waferHolderLV =
                    new G4LogicalVolume(type1_waferHolder, aFrameMaterial, "WaferHolder_LV");
                waferHolderLV->SetVisAttributes(fI_OpticalFrame_VisAttr);
                holeAndSlotLV = new G4LogicalVolume(photonFrameHoleAndSlot, aSpaceMaterial,
                                                    "FrameHoleAndSlot_LV");

                new G4PVPlacement(nullptr,
                                  {0, 0, photonFrameSizeZ / 2. - type1_photonFrameHoleLength / 2.},
                                  holeAndSlotLV, "PhotonFrame_HoleAndSlot_PV", resultLV, false, 0, OverlapCheck);

                new G4PVPlacement(
                    nullptr, {0, 0, type1_photonFrameHoleLength / 2. - waferHolderThick / 2.},
                    waferHolderLV, "PhotonDetector_Holder_PV", holeAndSlotLV, false, 0, OverlapCheck);

                G4cout << __PRETTY_FUNCTION__ << ": Adding a germanium detector..." << G4endl;
                germaniumWaferTub = new G4Tubs("PhotonFrame_GeWafer", 0, germaniumWaferRadius,
                                               germaniumWaferThick / 2., 0, 360. * deg);
                //germaniumWaferLV  = new G4LogicalVolume(germaniumWaferTub, _gewafer,
                germaniumWaferLV  = new G4LogicalVolume(germaniumWaferTub, _SiWafer,
                                                       "PhotonDetector_GeWafer_Type1_LV");
                germaniumWaferLV->SetVisAttributes(fI_GeWafer_VisAttr);

                aMSDForThisModule->fGeWaferLV = germaniumWaferLV;
                fI_GeWaferLVs.insert(germaniumWaferLV);

                new G4PVPlacement(nullptr, {0, 0, germaniumWaferVacThick / 2.}, germaniumWaferLV,
                                  "Type1_PhotonDetector_GeWafer_PV", holeAndSlotLV, false, 0, OverlapCheck);

                geWaferUpperGoldFilmTub =
                    new G4Tubs("GoldFilmAboveGeWafer_Tub", 0, geWaferUpperGoldRadius,
                               geWaferUpperGoldThickness / 2., 0, 360. * deg);
                geWaferUpperGoldFilmLV = new G4LogicalVolume(geWaferUpperGoldFilmTub, _gold,
                                                             "GoldFilmAboveGeWafer_Type1_LV");
                fI_GeWaferGoldFilmLVs.insert(geWaferUpperGoldFilmLV);

                geWaferUpperGoldFilmLV->SetVisAttributes(fI_GoldDefault_VisAttr);
                aMSDForThisModule->fGeWaferGoldFilmLV = geWaferUpperGoldFilmLV;

                nowUpperGoldFilmPos =
                    G4ThreeVector(0, -geWaferUpperGoldRadialDistance,
                                  germaniumWaferThick / 2. + geWaferUpperGoldThickness / 2. +
                                      germaniumWaferVacThick / 2.);

                new G4PVPlacement(nullptr, nowUpperGoldFilmPos, geWaferUpperGoldFilmLV,
                                  "GeWaferUpperGoldFilm_PV", holeAndSlotLV, false, 0, OverlapCheck);
                nowUpperGoldFilmPos.rotateZ(120. * deg);
                new G4PVPlacement(nullptr, nowUpperGoldFilmPos, geWaferUpperGoldFilmLV,
                                  "GeWaferUpperGoldFilm_PV", holeAndSlotLV, false, 1, OverlapCheck);
                nowUpperGoldFilmPos.rotateZ(120. * deg);
                new G4PVPlacement(nullptr, nowUpperGoldFilmPos, geWaferUpperGoldFilmLV,
                                  "GeWaferUpperGoldFilm_PV", holeAndSlotLV, false, 2, OverlapCheck);
            } else {
                germaniumWaferLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "PhotonDetector_GeWafer_Type1_LV", false);
                geWaferUpperGoldFilmLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "GoldFilmAboveGeWafer_Type1_LV", false);
                aMSDForThisModule->fGeWaferLV         = germaniumWaferLV;
                aMSDForThisModule->fGeWaferGoldFilmLV = geWaferUpperGoldFilmLV;
                if (germaniumWaferLV == nullptr || geWaferUpperGoldFilmLV == nullptr) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_HITINFO_ERROR",
                                G4ExceptionSeverity::RunMustBeAborted,
                                "LV is found but SD Hit info for that has not found.");
                }
            }
            break;
        case 2:
						if(fDbgMsgOn) {
							G4cout << __PRETTY_FUNCTION__
                   << ": Realistic Photon frame type 2 selected (It could be very slow)" << G4endl;
						}
            lvName   = "PhotonDetector_Type2_LV";
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(lvName, false);
            if (resultLV == nullptr) {
                photonFrame = new G4Box("Type2_PhotonFrame_Box", photonFrameSizeX / 2.,
                                        photonFrameSizeY / 2., photonFrameSizeZ / 2.);
                resultLV    = new G4LogicalVolume(photonFrame, aSpaceMaterial, lvName);
                resultLV->SetVisAttributes(fI_Invisible_VisAttr);

                photonFrameBase =
                    new G4Box("Type2_PhotonFrame_BaseSuppportBox", type2_photonFrameBaseWidthX / 2.,
                              type2_photonFrameBaseWidthY / 2., type2_photonFrameBaseThick / 2.);

                G4LogicalVolume *type2_photonFrameBaseLV = new G4LogicalVolume(photonFrameBase, aFrameMaterial,
                                                              "PhotonFrameBaseSuppport_Type2_LV");
                type2_photonFrameBaseLV->SetVisAttributes(fI_OpticalFrame_VisAttr);

                type2_photonFrameBaseSupportNowPosition =
                    G4ThreeVector(-photonFrameSizeX / 2. + type2_photonFrameBaseWidthX / 2.,
                                  photonFrameSizeY / 2. - type2_photonFrameBaseWidthY / 2.,
                                  -photonFrameSizeZ / 2. + type2_photonFrameBaseThick / 2.);

                new G4PVPlacement(nullptr, type2_photonFrameBaseSupportNowPosition,
                                  type2_photonFrameBaseLV, "PhotonFrameBaseSupport_Type2_PV",
                                  resultLV, false, 0, OverlapCheck);
                type2_photonFrameBaseSupportNowPosition[0] *= -1.;
                new G4PVPlacement(nullptr, type2_photonFrameBaseSupportNowPosition,
                                  type2_photonFrameBaseLV, "PhotonFrameBaseSupport_Type2_PV",
                                  resultLV, false, 1, OverlapCheck);
                type2_photonFrameBaseSupportNowPosition[1] *= -1.;
                type2_photonFrameBaseSupportNowPosition[0] *= -1.;
                new G4PVPlacement(nullptr, type2_photonFrameBaseSupportNowPosition,
                                  type2_photonFrameBaseLV, "PhotonFrameBaseSupport_Type2_PV",
                                  resultLV, false, 2, OverlapCheck);
                type2_photonFrameBaseSupportNowPosition[0] *= -1.;
                new G4PVPlacement(nullptr, type2_photonFrameBaseSupportNowPosition,
                                  type2_photonFrameBaseLV, "PhotonFrameBaseSupport_Type2_PV",
                                  resultLV, false, 3, OverlapCheck);

                G4Box *type2_photonFrameUpperBox = new G4Box(
                    "Type2_PhotonFrameUpperBox", photonFrameSizeX / 2., photonFrameSizeY / 2.,
                    (photonFrameSizeZ - type2_photonFrameBaseThick) / 2.);
                G4LogicalVolume *type2_photonFrameUpperLV = new G4LogicalVolume(
                    type2_photonFrameUpperBox, aFrameMaterial, "PhotonFrameUpper_Type2_LV");
                type2_photonFrameUpperLV->SetVisAttributes(fI_OpticalFrame_VisAttr);

                new G4PVPlacement(nullptr, {0, 0, type2_photonFrameBaseThick / 2.},
                                  type2_photonFrameUpperLV, "PhotonFrameUpper_Type2_PV", resultLV,
                                  false, 0, OverlapCheck);

                photonFrameHole =
                    new G4Tubs("Type2_PhotonFrame_MainHoleTub", 0, type2_photonFrameHoleRadius,
                               type2_photonFrameHoleLength / 2., 0, 360. * deg);

                G4VSolid *type2_photonFrameRoofHole = new G4Tubs(
                    "Type2_PhotonFrame_RoofHoleTub", 0, type2_photonFrameRoofHoleRadius,
                    type2_photonFrameRoofThick / 2. + solidBooleanTolSmall,
                    -type2_photonFrameHeartAngle, type2_photonFrameHeartAngle * 2. + 180. * deg);

                G4EllipticalTube *type2_photonFrameHoleEllipsoid = new G4EllipticalTube(
                    "Type2_PhotonFrame_RoofHoleHeartLeafTube",
                    type2_photonFrameHeartLeapMinorRadius, type2_photonFrameHeartLeapMajorRadius,
                    type2_photonFrameRoofThick / 2. + solidBooleanTolSmall);

                type2_photonFrameHeartLeapPosition = G4ThreeVector(
                    type2_photonFrameRoofHoleRadius - type2_photonFrameHeartLeapMajorRadius, 0, 0);
                type2_photonFrameHeartLeapPosition.rotateZ(-type2_photonFrameHeartAngle +
                                                           solidBooleanTolAngle);

                type2_photonFrameHeartLeapRotMtx = new G4RotationMatrix();
                type2_photonFrameHeartLeapRotMtx->rotateZ(-90 * deg + type2_photonFrameHeartAngle -
                                                          solidBooleanTolAngle);

                type2_photonFrameRoofHole = new G4UnionSolid(
                    "Type2_PhotonFrame_RoofHoleHeartStage1Solid", type2_photonFrameRoofHole,
                    type2_photonFrameHoleEllipsoid, type2_photonFrameHeartLeapRotMtx,
                    type2_photonFrameHeartLeapPosition);

                type2_photonFrameHeartLeapPosition = G4ThreeVector(
                    -type2_photonFrameRoofHoleRadius + type2_photonFrameHeartLeapMajorRadius, 0, 0);
                type2_photonFrameHeartLeapPosition.rotateZ(type2_photonFrameHeartAngle -
                                                           solidBooleanTolAngle);

                type2_photonFrameHeartLeapRotMtx->rotateZ(-2. * type2_photonFrameHeartAngle +
                                                          2. * solidBooleanTolAngle);

                type2_photonFrameRoofHole = new G4UnionSolid(
                    "Type2_PhotonFrame_RoofHoleHeartStage2Solid", type2_photonFrameRoofHole,
                    type2_photonFrameHoleEllipsoid, type2_photonFrameHeartLeapRotMtx,
                    type2_photonFrameHeartLeapPosition);

                G4Tubs *type2_PEEKHoleTub =
                    new G4Tubs("Type2_PhotonFrame_PEEKHoleTube", 0, type2_photonFramePeekHoleRadius,
                               type2_photonFramePeekHoleThick / 2., 0., 360. * deg);

                type2_PEEKNowPosition = G4ThreeVector(type2_photonFrameHoleRadius, 0.,
                                                      -type2_photonFrameHoleLength / 2. +
                                                          type2_photonFramePeekHoleThick / 2.);
                type2_PEEKNowPosition.rotateZ(45. * deg);

                photonFrameHoleAndSlot =
                    new G4UnionSolid("Type2_PhotonFrameHole_PEEKHoleStage1Solid", photonFrameHole,
                                     type2_PEEKHoleTub, nullptr, type2_PEEKNowPosition);
                type2_PEEKNowPosition.rotateZ(90. * deg);
                photonFrameHoleAndSlot = new G4UnionSolid(
                    "Type2_PhotonFrameHole_PEEKHoleStage2Solid", photonFrameHoleAndSlot,
                    type2_PEEKHoleTub, nullptr, type2_PEEKNowPosition);

                type2_PEEKNowPosition.rotateZ(90. * deg);
                photonFrameHoleAndSlot = new G4UnionSolid(
                    "Type2_PhotonFrameHole_PEEKHoleStage3Solid", photonFrameHoleAndSlot,
                    type2_PEEKHoleTub, nullptr, type2_PEEKNowPosition);

                type2_PEEKNowPosition.rotateZ(90. * deg);
                photonFrameHoleAndSlot = new G4UnionSolid(
                    "Type2_PhotonFrameHole_PEEKHoleStage4Solid", photonFrameHoleAndSlot,
                    type2_PEEKHoleTub, nullptr, type2_PEEKNowPosition);

                photonFrameHoleAndSlot =
                    new G4UnionSolid("Type2_PhotonFrameHole_HeartAttachedSolid",
                                     photonFrameHoleAndSlot, type2_photonFrameRoofHole, nullptr,
                                     {0, 0,
                                      type2_photonFrameHoleLength / 2 +
                                          type2_photonFrameRoofThick / 2. - solidBooleanTolSmall});

                holeAndSlotLV = new G4LogicalVolume(photonFrameHoleAndSlot, aSpaceMaterial,
                                                    "PhotonFrameHoleAndSlot_Type2_LV");

                G4Tubs *type2_PEEKTub =
                    new G4Tubs("Type2_PhotonFramePEEKTubs", 0, type2_photonFramePeekRdaius,
                               type2_photonFramePeekThick / 2., 0, 360 * deg);
                G4LogicalVolume *type2_PEEKLV =
                    new G4LogicalVolume(type2_PEEKTub, _peek1, "PhotonFramePEEK_Type2_LV");
                type2_PEEKLV->SetVisAttributes(fI_PEEK_VisAttr);

                type2_PEEKNowPosition = G4ThreeVector(type2_photonFrameHoleRadius, 0.,
                                                      -type2_photonFrameHoleLength / 2. +
                                                          type2_photonFramePeekThick / 2.);
                type2_PEEKNowPosition.rotateZ(45. * deg);

                for (G4int i = 0; i < 4; i++) {
                    new G4PVPlacement(nullptr, type2_PEEKNowPosition, type2_PEEKLV,
                                      "PhotonFramePEEK_Type2_PV", holeAndSlotLV, false, i, OverlapCheck);
                    type2_PEEKNowPosition.rotateZ(90. * deg);
                }

                new G4PVPlacement(nullptr, {0, 0, -type2_photonFrameRoofThick / 2.}, holeAndSlotLV,
                                  "PhotonFrameHoleAndSlot_Type2_PV", type2_photonFrameUpperLV,
                                  false, 0, OverlapCheck);

                G4Tubs *type2_germaniumWaferSlotTub =
                    new G4Tubs("Type2_PhotonFrame_GeWaferSlotTube", 0, germaniumWaferRadius,
                               type2_germaniumWaferSlotThick / 2., 0, 360. * deg);
                G4LogicalVolume *type2_germaniumWaferSlotLV =
                    new G4LogicalVolume(type2_germaniumWaferSlotTub, aSpaceMaterial,
                                        "PhotonDetector_GeWaferSlot_Type2_LV");

                germaniumWaferTub = new G4Tubs("PhotonFrame_GeWafer", 0, germaniumWaferRadius,
                                               germaniumWaferThick / 2., 0, 360. * deg);
                //germaniumWaferLV  = new G4LogicalVolume(germaniumWaferTub, _gewafer,
                germaniumWaferLV  = new G4LogicalVolume(germaniumWaferTub, _SiWafer,
                                                       "PhotonDetector_GeWafer_Type2_LV");
                germaniumWaferLV->SetVisAttributes(fI_GeWafer_VisAttr);
                aMSDForThisModule->fGeWaferLV = germaniumWaferLV;
                fI_GeWaferLVs.insert(germaniumWaferLV);

                new G4PVPlacement(nullptr, {}, germaniumWaferLV, "PhotonDetector_GeWafer_Type2_PV",
                                  type2_germaniumWaferSlotLV, false, 0, OverlapCheck);

                geWaferUpperGoldFilmTub =
                    new G4Tubs("Type2_GoldFilmAboveGeWafer_Tub", 0, geWaferUpperGoldRadius,
                               geWaferUpperGoldThickness / 2., 0, 360. * deg);
                geWaferUpperGoldFilmLV = new G4LogicalVolume(geWaferUpperGoldFilmTub, _gold,
                                                             "GoldFilmAboveGeWafer_Type2_LV");
                geWaferUpperGoldFilmLV->SetVisAttributes(fI_GoldDefault_VisAttr);
                aMSDForThisModule->fGeWaferGoldFilmLV = geWaferUpperGoldFilmLV;
                fI_GeWaferGoldFilmLVs.insert(geWaferUpperGoldFilmLV);

                nowUpperGoldFilmPos =
                    G4ThreeVector(0, -geWaferUpperGoldRadialDistance,
                                  germaniumWaferThick / 2. + geWaferUpperGoldThickness / 2.);

                for (G4int i = 0; i < 3; i++) {
                    new G4PVPlacement(nullptr, nowUpperGoldFilmPos, geWaferUpperGoldFilmLV,
                                      "GeWaferUpperGoldFilm_Type2_PV", type2_germaniumWaferSlotLV,
                                      false, i, OverlapCheck);
                    nowUpperGoldFilmPos.rotateZ(120. * deg);
                }

                new G4PVPlacement(nullptr,
                                  {0, 0,
                                   -photonFrameSizeZ / 2. + type2_photonFrameBaseThick -
                                       type2_germaniumWaferSlotThick / 2.},
                                  type2_germaniumWaferSlotLV, "PhotonDetector_GeWaferSlot_Type2_PV",
                                  resultLV, false, 0, OverlapCheck);

                type2_PEEKClampVertices.resize(8);
                type2_PEEKClampVertices[0] = {type2_photonFrameClampBox1SizeXY / 2.,
                                              -type2_photonFrameClampBox1SizeXY / 2.};
                type2_PEEKClampVertices[1] = {type2_photonFrameClampBox1SizeXY / 2.,
                                              type2_photonFrameClampBox1SizeXY / 2.};
                type2_PEEKClampVertices[2] = {type2_photonFrameClampBox2SizeX / 2.,
                                              type2_photonFrameClampBox1SizeXY / 2.};
                type2_PEEKClampVertices[3] = {type2_photonFrameClampBox2SizeX / 2.,
                                              type2_photonFrameClampBox1SizeXY / 2. +
                                                  type2_photonFrameClampBox2SizeY};
                type2_PEEKClampVertices[4] = {-type2_photonFrameClampBox2SizeX / 2.,
                                              type2_photonFrameClampBox1SizeXY / 2. +
                                                  type2_photonFrameClampBox2SizeY};
                type2_PEEKClampVertices[5] = {-type2_photonFrameClampBox2SizeX / 2.,
                                              type2_photonFrameClampBox1SizeXY / 2.};
                type2_PEEKClampVertices[6] = {-type2_photonFrameClampBox1SizeXY / 2.,
                                              type2_photonFrameClampBox1SizeXY / 2.};
                type2_PEEKClampVertices[7] = {-type2_photonFrameClampBox1SizeXY / 2.,
                                              -type2_photonFrameClampBox1SizeXY / 2.};

                G4ExtrudedSolid *type2_PEEKClampSolid =
                    new G4ExtrudedSolid("Type2_PhotonFramePEEKClampSolid", type2_PEEKClampVertices,
                                        type2_photonFrameClampThick / 2., 0, 1., 0., 1.);

                G4LogicalVolume *type2_PEEKClampLV = new G4LogicalVolume(type2_PEEKClampSolid, _teflon,
                                                        "PhotonFramePEEKClamp_Type2_LV");
                type2_PEEKClampLV->SetVisAttributes(fI_Clamp_VisAttr);

                G4Tubs *type2_PEEKClampBoltTub =
                    new G4Tubs("Type2_PEEKClampBoltTub", 0, type2_photonFrameClampBoltRadius,
                               type2_photonFrameClampThick / 2., 0, 360. * deg);
                G4LogicalVolume *type2_PEEKClampBoltLV = new G4LogicalVolume(type2_PEEKClampBoltTub, _stainless,
                                                            "PEEKClampBolt_Type2_LV");
                type2_PEEKClampBoltLV->SetVisAttributes(fI_ClampBolt_VisAttr);

                new G4PVPlacement(nullptr, {}, type2_PEEKClampBoltLV,
                                  "PhotonFramePEEKClampBolt_Type2_PV", type2_PEEKClampLV, false, 0, OverlapCheck);

                type2_PEEKClampNowPosition = {type2_PEEKClampRadialDistance, 0,
                                              -photonFrameSizeZ / 2. + type2_photonFrameBaseThick -
                                                  type2_germaniumWaferSlotThick -
                                                  type2_photonFrameClampThick / 2.};

                type2_PEEKClampNowPosition.rotateZ(45. * deg);

                type2_PEEKClampRotMtx[0] = new G4RotationMatrix();
                type2_PEEKClampRotMtx[0]->rotateZ(-45. * deg - 90. * deg);

                new G4PVPlacement(type2_PEEKClampRotMtx[0], type2_PEEKClampNowPosition,
                                  type2_PEEKClampLV, "PhotonFramePEEKClamp_Type2_PV", resultLV,
                                  false, 0, OverlapCheck);

                for (G4int i = 1; i < 4; i++) {
                    type2_PEEKClampRotMtx[i] = new G4RotationMatrix(*type2_PEEKClampRotMtx[i - 1]);
                    type2_PEEKClampRotMtx[i]->rotateZ(-90. * deg);

                    type2_PEEKClampNowPosition.rotateZ(90. * deg);

                    new G4PVPlacement(type2_PEEKClampRotMtx[i], type2_PEEKClampNowPosition,
                                      type2_PEEKClampLV, "PhotonFramePEEKClamp_Type2_PV", resultLV,
                                      false, i, OverlapCheck);
                }
            } else {
                germaniumWaferLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "PhotonDetector_GeWafer_Type2_LV", false);
                geWaferUpperGoldFilmLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "GoldFilmAboveGeWafer_Type2_LV", false);
                aMSDForThisModule->fGeWaferLV         = germaniumWaferLV;
                aMSDForThisModule->fGeWaferGoldFilmLV = geWaferUpperGoldFilmLV;
                if (germaniumWaferLV == nullptr || geWaferUpperGoldFilmLV == nullptr) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_HITINFO_ERROR",
                                G4ExceptionSeverity::RunMustBeAborted,
                                "LV is found but SD Hit info for that has not found.");
                }
            }
            break;
        case 3: // Added for Real AMoRE-I configuration 
						if(fDbgMsgOn) {
							G4cout << __PRETTY_FUNCTION__
                   << ": Realistic Photon frame type 3 selected (It could be very slow)" << G4endl;
						}
            lvName   = "PhotonDetector_Type3_LV";
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(lvName, false);
            if (resultLV == nullptr) {
                photonFrame = new G4Box("Type3_PhotonFrame_Box", 
										type3_photonFrameSizeX / 2.,
                    type3_photonFrameSizeY / 2., 
										(type3_photonFrameSizeZ + type3_photonFrameTopSpace)/ 2.);

								// PhotonFrame Upper ----------------------------------
                G4Box *type3_photonFrameUpperBox = new G4Box("Type3_PhotonFrameUpperBox", 
										type3_photonFrameSizeX / 2., 
										type3_photonFrameSizeY / 2.,
                    (type3_photonFrameSizeZ - type3_photonFrameBaseThick) / 2.);

                G4Tubs *type3_photonFrameHoleTubs = new G4Tubs("Type3_PhotonFrame_MainHoleTub", 0, 
										type3_photonFrameHoleRadius, 
										type3_photonFrameHoleThick / 2. + solidBooleanTolBig, 
										0, 360. *deg);

								G4Box *type3_photonFrameHoleBox = new G4Box("Type3_PhotonFrame_HoleBox", 
										type3_photonFrameHoleBoxSizeX / 2., 
										type3_photonFrameHoleBoxSizeY /2. , 
										type3_photonFrameHoleThick / 2. + solidBooleanTolBig);

                G4UnionSolid *type3_photonFrameHole = new G4UnionSolid("Type3_PhotonFrameHole", 
										type3_photonFrameHoleTubs, type3_photonFrameHoleBox, nullptr, 
										{0, type3_photonFrameHoleBoxPositionY-solidBooleanTolBig,0});

								G4SubtractionSolid *type3_photonFrameUpper = new G4SubtractionSolid("Type3_PhotonFrameUpper",
										type3_photonFrameUpperBox, type3_photonFrameHole, nullptr, {0,0,0});

								// PhotonFrame Base Support 4ea ---------------------
								std::vector<G4TwoVector> Pentagon_Template;
								G4TwoVector Pentagon_Points[5] = {G4TwoVector(0.5, 0.5), G4TwoVector(-0.5, 0.5), G4TwoVector(-0.5,-0.5),
																										G4TwoVector(0,-0.5), G4TwoVector(0.5,0)};

								Pentagon_Template.push_back(Pentagon_Points[0]);
								Pentagon_Template.push_back(Pentagon_Points[1]);
								Pentagon_Template.push_back(Pentagon_Points[2]);
								Pentagon_Template.push_back(Pentagon_Points[3]);
								Pentagon_Template.push_back(Pentagon_Points[4]);

								G4ExtrudedSolid *type3_photonFrameBase = new G4ExtrudedSolid(
										"Type3_PhotonFrame_BaseSupportBox", Pentagon_Template, type3_photonFrameBaseThick / 2., 
										0, type3_photonFrameBaseWidthX, 0, type3_photonFrameBaseWidthY);

                G4ThreeVector type3_photonFrameBaseSupportNowPosition = G4ThreeVector(
										-type3_photonFrameSizeX / 2. + type3_photonFrameBaseWidthX / 2.,
                    type3_photonFrameSizeY / 2. - type3_photonFrameBaseWidthY / 2.,
										- (type3_photonFrameSizeZ - type3_photonFrameBaseThick + type3_photonFrameBaseThick) / 2.);

								type3_photonFrameSupportRotMtx[0] = new G4RotationMatrix();
								for(int i = 0; i < 4; i++)
								{
									if(i!=0) 
									{
										type3_photonFrameSupportRotMtx[i] = new G4RotationMatrix(*type3_photonFrameSupportRotMtx[i - 1]);
										type3_photonFrameSupportRotMtx[i] -> rotateZ(90.* deg);
										type3_photonFrameBaseSupportNowPosition[0] *= -1.;
										if(i==2)
										{
											type3_photonFrameBaseSupportNowPosition[0] *= -1.;
											type3_photonFrameBaseSupportNowPosition[1] *= -1.;
										}

										photonFrameHoleAndSlot = new G4UnionSolid("Type3_photonFrameHoleAndSlot",
												photonFrameHoleAndSlot, type3_photonFrameBase, 
												type3_photonFrameSupportRotMtx[i], type3_photonFrameBaseSupportNowPosition);
									} else{
										photonFrameHoleAndSlot = new G4UnionSolid("Type3_photonFrameHoleAndSlot",
											type3_photonFrameUpper, type3_photonFrameBase,
											type3_photonFrameSupportRotMtx[i], type3_photonFrameBaseSupportNowPosition);
									}
								}

								// PhotonFrame PEEK 3ea -------------------------------------------
                G4Tubs *type3_PEEKTub = new G4Tubs("Type3_PhotonFramePEEKTubs", 0, 
										type3_photonFramePeekRdaius, 
										type3_photonFramePeekThick / 2., 0, 360 * deg);
								G4Tubs *type3_PEEKHoleTub = new G4Tubs("Type3_PhotonFramePEEKthinHole", 0,
										type3_photonFramePeekHoleRadius, 
										type3_photonFramePeekHoleThick / 2., 0, 360 * deg);
								G4SubtractionSolid *type3_PEEKs = new G4SubtractionSolid("Type3_PhotonFramePEEK",
										type3_PEEKTub, type3_PEEKHoleTub, nullptr, 
										{0,0,type3_photonFramePeekThick / 2. - type3_photonFramePeekHoleThick / 2.});

								G4Tubs *type3_PEEKBoltHead = new G4Tubs("Type3_PEEKBoltHead", 0, 
										type3_photonFramePeekBoltHeadRadius,
										(type3_photonFrameTopSpace - solidBooleanTol) / 2., 0, 360. * deg);
								G4Tubs *type3_PEEKBoltBody = new G4Tubs("Type3_PEEKBoltBody", 0,
										type3_photonFramePeekBoltRadius,
										(type3_photonFrameSizeZ - type3_photonFrameBaseThick) / 2., 0, 360 * deg);
								G4UnionSolid *type3_PEEKBolts = new G4UnionSolid("Type3_PhotonFramePEEKBolts",
										type3_PEEKBoltBody, type3_PEEKBoltHead, nullptr, {0, 0, 
										(type3_photonFrameSizeZ - type3_photonFrameBaseThick)/2. + 
											(type3_photonFrameTopSpace - solidBooleanTol) / 2.});
								type3_PEEKBolts = new G4UnionSolid("Type3_photonFramePEEKBolts",
										type3_PEEKBolts, type3_PEEKHoleTub, nullptr, {0, 0, 
										-(type3_photonFrameSizeZ - type3_photonFrameBaseThick) / 2. -
											type3_photonFramePeekHoleThick / 2.});
								
								G4Tubs *type3_PEEKBoltHole = new G4Tubs("Type3_PEEKBoltHole", 0,
										type3_photonFramePeekBoltRadius, type3_photonFrameSizeZ, 0, 360 * deg);

                G4ThreeVector type3_PEEKNowPosition = G4ThreeVector(0, 
										type3_photonFrameHoleRadius + type3_photonFrameCuBoltsRadius, 0);

								if(aRealistic)
								{
									// Hole for PEEK Bolts (3ea) --------------------------
									for (G4int i = 0; i < 3; i++) {
										photonFrameHoleAndSlot = new G4SubtractionSolid("Type3_photonFrameHoleAndSlot_PEEKHole",
												photonFrameHoleAndSlot, type3_PEEKBoltHole, nullptr,	type3_PEEKNowPosition);
                    type3_PEEKNowPosition.rotateZ(120. * deg);
									}

									// Hole for PCB Bolts (3ea) ----------------------------
									G4Tubs *type3_PCBBoltsHole = new G4Tubs("type3_PhotonFrameBoltsHole", 0,
											photDetPCBBoltsHoleRadius, 
											photDetPCBBoltsHoleHeight / 2 + solidBooleanTol, 
											0, 360. * deg);

									G4ThreeVector type3_PCBBoltsPosition = {
										type3_photDetPCBBaseBoxSizeX / 2. - photDetPCBBoltsHoleDistX, 
										type3_photDetPCBBaseBoxSizeY / 2. - photDetPCBBoltsHoleDistY,  
										(type3_photonFrameSizeZ - type3_photonFrameBaseThick - photDetPCBBoltsHoleHeight) / 2.
									    + type3_photonFrameTopSpace};

									// Make PCB bolts hole on photonframe
									photonFrameHoleAndSlot = new G4SubtractionSolid("Type3_PhotonFrameHoleAndSlot_BoltsHole1",
											photonFrameHoleAndSlot, type3_PCBBoltsHole, nullptr, type3_PCBBoltsPosition);

									type3_PCBBoltsPosition[1] *= -1.;
									photonFrameHoleAndSlot = new G4SubtractionSolid("Type3_PhotonFrameHoleAndSlot_BoltsHole2",
											photonFrameHoleAndSlot, type3_PCBBoltsHole, nullptr, type3_PCBBoltsPosition);

									type3_PCBBoltsPosition[0] *= -1.;
									photonFrameHoleAndSlot = new G4SubtractionSolid("Type3_PhotonFrameHoleAndSlot_BoltsHole3",
											photonFrameHoleAndSlot, type3_PCBBoltsHole, nullptr, type3_PCBBoltsPosition);

									// Make PCB bolts hole on photonframe mothervolume
									type3_PCBBoltsPosition[2] = (type3_photonFrameSizeZ + type3_photonFrameTopSpace 
											                          - photDetPCBBoltsHoleHeight)/ 2. ;
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForPCBBolts",
											photonFrame, type3_PCBBoltsHole, nullptr, type3_PCBBoltsPosition);
									type3_PCBBoltsPosition[0] *=-1.;
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForPCBBolts",
											photonFrame, type3_PCBBoltsHole, nullptr, type3_PCBBoltsPosition);
									type3_PCBBoltsPosition[1] *= -1.;
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForPCBBolts",
											photonFrame, type3_PCBBoltsHole, nullptr, type3_PCBBoltsPosition);

									// Hole for photon-phonon connection bolts (4ea) -----------------------
									G4Tubs *type3_Bolts_hole = new G4Tubs("type3_PhotonPhononBolts_hole", 0,
											type3_photonPhononBoltsRadius + solidBooleanTol,
											type3_photonPhononBoltsThick / 2. + solidBooleanTol/2., 0, 360. * deg);

									G4Tubs *type3_BoltsHead_hole = new G4Tubs("type3_photonPhononBoltsHead_hole",0,
											type3_photonPhononBoltsHeadRadius + solidBooleanTol,
											type3_photonPhononBoltsHeadThick / 2. + solidBooleanTol/2., 0, 360. * deg);

									G4UnionSolid *type3_BoltsHole = new G4UnionSolid("type3_photonPhononBoltsHole",
											type3_Bolts_hole,type3_BoltsHead_hole, nullptr, {0, 0, 
											(type3_photonPhononBoltsThick + type3_photonPhononBoltsHeadThick)/2. + solidBooleanTol});

									for(int i = 0; i < 4; i++)
									{
										if(i==0 || i==2)
										{
											type3_photonFrameBaseSupportNowPosition[0] *= -1.;
											type3_photonFrameBaseSupportNowPosition[1] *= -1.;
										} else type3_photonFrameBaseSupportNowPosition[0] *= -1.;

										// make bolts hole on photonframe 
										type3_photonFrameBaseSupportNowPosition[2] = (type3_photonFrameSizeZ - type3_photonFrameBaseThick) / 2.
										                                               - (type3_photonPhononBoltsThick + solidBooleanTol) / 2.;
										photonFrameHoleAndSlot = new G4SubtractionSolid("Type3_PhotonFrameHoleAndSlot_BoltsHole",
												photonFrameHoleAndSlot, type3_BoltsHole, nullptr, type3_photonFrameBaseSupportNowPosition);

										// make bolts hole on photonframe mothervolume
										type3_photonFrameBaseSupportNowPosition[2] = (type3_photonFrameSizeZ - type3_photonFrameTopSpace) / 2.
											                                            - (type3_photonPhononBoltsThick + solidBooleanTol)/ 2.;
										photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForBolts",
												photonFrame, type3_BoltsHole, nullptr, type3_photonFrameBaseSupportNowPosition);
									}

									// Hole for clamp (4ea) ---------------------------------------------
									G4double clampSizeX = type3_copperFrameSizeX / 2. - type3_smallBlockSizeX - type3_photonFrameSizeX / 2.;
									G4double clampPosX = type3_copperFrameSizeX / 2. - type3_smallBlockSizeX - clampSizeX * 2;
									G4double clampPosY = type3_copperFrameSizeY / 2. - type3_smallBlockSizeY - clampSizeX * 2;
									G4Tubs *type3_ClampTubs = new G4Tubs("type3_Clamp_Tub", 0,
											type3_photonFrameClampRadius, 
											type3_photonFrameClampThick / 2., 90 * deg, 180. * deg);
									G4Box *type3_ClampBox = new G4Box("type3_ClampBox",
										clampSizeX, type3_photonFrameClampRadius, type3_photonFrameClampSizeZ / 2.);
									G4UnionSolid *type3_ClampSolid = new G4UnionSolid("type3_ClampSolid", 
										type3_ClampTubs, type3_ClampBox, nullptr, {clampSizeX/2., 0, 
										type3_photonFrameClampThick / 2. - type3_photonFrameClampSizeZ / 2.});

									G4RotationMatrix *type3_ClampRotMtx = new G4RotationMatrix();
									G4ThreeVector nowClampPosition = G4ThreeVector(clampPosX, 0, 
											(type3_photonFrameSizeZ - type3_photonFrameTopSpace)/2. 
											- type3_photonFrameHoleThick - type3_photonFramePeekThick 
											- type3_germaniumWaferThick - type3_photonFrameClampThick / 2. );

									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForClamp",
											photonFrame, type3_ClampSolid, nullptr, nowClampPosition);

									nowClampPosition[0] -= clampPosX ;
									nowClampPosition[1] += clampPosY ;
									type3_ClampRotMtx->rotateZ(-90.* deg);
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForClamp",
											photonFrame, type3_ClampSolid, type3_ClampRotMtx, nowClampPosition);

									nowClampPosition[0] -= clampPosX ;
									nowClampPosition[1] -= clampPosY ;
									type3_ClampRotMtx->rotateZ(-90.*deg);
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForClamp",
											photonFrame, type3_ClampSolid, type3_ClampRotMtx, nowClampPosition);

									nowClampPosition[0] += clampPosX ;
									nowClampPosition[1] -= clampPosY ;
									type3_ClampRotMtx->rotateZ(-90.*deg);
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForClamp",
											photonFrame, type3_ClampSolid, type3_ClampRotMtx, nowClampPosition);

									// Hole for clamp bolts (2ea/4ea) ---------------------------------
									G4Tubs *clampBoltsBodyHole = new G4Tubs("clampBoltsBodyHole",0,
											type3_photonPhononBoltsRadius, type3_copperFrameSizeZ / 2., 0, 360. * deg);
									G4Tubs *clampBoltsHeadHole = new G4Tubs("clampBoltsHeadHole", 0,
											type3_photonPhononBoltsHeadRadius, type3_photonPhononBoltsHeadThick / 2., 0, 360. * deg);
									G4UnionSolid *clampBoltsHole = new G4UnionSolid("clampBoltsHole",
											clampBoltsBodyHole, clampBoltsHeadHole, nullptr, {0,0,
											type3_copperFrameSizeZ/2+type3_photonPhononBoltsHeadThick / 2.});

									nowClampPosition = {0, (type3_copperFrameSizeY - type3_copperFrameYWidth)/2.,
										-(type3_copperFrameSizeZ + type3_photonFrameSizeZ + type3_photonFrameTopSpace)/2};
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForClampBolts",
											photonFrame, clampBoltsHole, nullptr, nowClampPosition);
									nowClampPosition[1] *= -1;
									photonFrame = new G4SubtractionSolid("Type3_photonFrame_HoleForClampBolts",
											photonFrame, clampBoltsHole, nullptr, nowClampPosition);

								}
								
								// mother physical volume ( PHOTON FRAME ) -----
                resultLV    = new G4LogicalVolume(photonFrame, aSpaceMaterial, lvName);
                resultLV->SetVisAttributes(fI_Invisible_VisAttr);

								// photon frame physical volume -----
                G4LogicalVolume *photonFrameHoleAndSlotLV = new G4LogicalVolume(
                    photonFrameHoleAndSlot, aFrameMaterial, "PhotonFrameHoleAndSlot_Type3_LV");
                photonFrameHoleAndSlotLV->SetVisAttributes(fI_OpticalFrame_VisAttr);

                new G4PVPlacement(nullptr, 
										{0, 0, (type3_photonFrameBaseThick - type3_photonFrameTopSpace) / 2.}, 
										photonFrameHoleAndSlotLV, "PhotonFrameHoleAndSlot_Type3_PV", 
										resultLV, false, 0, OverlapCheck);

								// PEEK physical volume -----
                G4LogicalVolume *type3_PEEKLV = new G4LogicalVolume(
										type3_PEEKs, _peek1, "PhotonFramePEEK_Type3_LV");
                type3_PEEKLV->SetVisAttributes(fI_PEEK_VisAttr);

								G4LogicalVolume *type3_PEEKBoltLV = new G4LogicalVolume(
										type3_PEEKBolts, _copper, "PhotonFramePEEKBolt_Type3_LV");
								type3_PEEKBoltLV->SetVisAttributes(fI_Bolt_VisAttr);

                type3_PEEKNowPosition = G4ThreeVector(0, 
										type3_photonFrameHoleRadius + type3_photonFrameCuBoltsRadius,
										(type3_photonFrameSizeZ-type3_photonFrameTopSpace)/2. 
										- type3_photonFrameHoleThick - type3_photonFramePeekThick/2.);

                for (G4int i = 0; i < 3; i++) {
										type3_PEEKNowPosition[2] = (type3_photonFrameSizeZ-type3_photonFrameTopSpace) / 2.
											                          - type3_photonFrameHoleThick / 2.;
										new G4PVPlacement(nullptr, type3_PEEKNowPosition, type3_PEEKBoltLV,
												              "PhotonFramePEEKBolt_Type3_PV", resultLV, false, i, OverlapCheck);

										type3_PEEKNowPosition[2] -= (type3_photonFrameHoleThick + type3_photonFramePeekThick) / 2.;
                    new G4PVPlacement(nullptr, type3_PEEKNowPosition, type3_PEEKLV,
                                      "PhotonFramePEEK_Type3_PV", resultLV, false, i, OverlapCheck);
                    type3_PEEKNowPosition.rotateZ(120. * deg);
                }

								// PhotonFrame Silicon wafer : wafer changed for AMORE-I (AMORE-pilot: ge, AMORE-I: si)
                germaniumWaferTub = new G4Tubs("PhotonFrame_GeWafer", 0, 
										type3_germaniumWaferRadius, type3_germaniumWaferThick / 2., 0, 360. * deg);
                germaniumWaferLV  = new G4LogicalVolume(
										//germaniumWaferTub, _gewafer, "PhotonDetector_GeWafer_Type3_LV");
										germaniumWaferTub, _SiWafer, "PhotonDetector_GeWafer_Type3_LV");
                germaniumWaferLV->SetVisAttributes(fI_GeWafer_VisAttr);
                aMSDForThisModule->fGeWaferLV = germaniumWaferLV;
                fI_GeWaferLVs.insert(germaniumWaferLV);

                new G4PVPlacement(nullptr, 
										{0,0, type3_photonFrameSizeZ/2. - type3_photonFrameTopSpace/2.
										- type3_photonFrameHoleThick - type3_photonFramePeekThick - type3_germaniumWaferThick/2.},
										germaniumWaferLV, "PhotonDetector_GeWafer_Type3_PV", resultLV, false, 0, OverlapCheck);
								
								// PhotonFrame Gold Film 3ea --------------------------------------------
                geWaferUpperGoldFilmTub = new G4Tubs("Type3_GoldFilmAboveGeWafer_Tub", 0, 
										type3_geWaferUpperGoldRadius, 
										type3_geWaferUpperGoldThickness / 2., 
										0, 360. * deg);
                geWaferUpperGoldFilmLV = new G4LogicalVolume(
										geWaferUpperGoldFilmTub, _gold, "GoldFilmAboveGeWafer_Type3_LV");
                geWaferUpperGoldFilmLV->SetVisAttributes(fI_GoldDefault_VisAttr);
                aMSDForThisModule->fGeWaferGoldFilmLV = geWaferUpperGoldFilmLV;
                fI_GeWaferGoldFilmLVs.insert(geWaferUpperGoldFilmLV);

                nowUpperGoldFilmPos = G4ThreeVector(0, -geWaferUpperGoldRadialDistance,
										type3_photonFrameSizeZ/2. - type3_photonFrameTopSpace/2.
										- type3_photonFrameHoleThick - type3_photonFramePeekThick + type3_geWaferUpperGoldThickness/2.);

                for (G4int i = 0; i < 3; i++) {
                    new G4PVPlacement(nullptr, nowUpperGoldFilmPos, 
												geWaferUpperGoldFilmLV, "GeWaferUpperGoldFilm_Type3_PV", 
												resultLV, false, i, OverlapCheck);
                    nowUpperGoldFilmPos.rotateZ(120. * deg);
                }

								// Photon PCB -----------------------------
								G4ThreeVector nowPhotonPCBPos = G4ThreeVector(0, 0, (type3_photonFrameSizeZ + type3_photonFrameTopSpace)/2.);
								if(aRealistic)
								{	// Photon PCBs with hole for bolts 
									Place_I_PCBs(true, resultLV, 1, nowPhotonPCBPos);
								} else { // photon PCB without bolts
									Place_I_PCBs(true, resultLV, 2, nowPhotonPCBPos);
								}
            } else {
                germaniumWaferLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "PhotonDetector_GeWafer_Type3_LV", false);
                geWaferUpperGoldFilmLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
                    "GoldFilmAboveGeWafer_Type3_LV", false);
                aMSDForThisModule->fGeWaferLV         = germaniumWaferLV;
                aMSDForThisModule->fGeWaferGoldFilmLV = geWaferUpperGoldFilmLV;
                if (germaniumWaferLV == nullptr || geWaferUpperGoldFilmLV == nullptr) {
                    G4Exception(__PRETTY_FUNCTION__, "AMORE_I_HITINFO_ERROR",
                                G4ExceptionSeverity::RunMustBeAborted,
                                "LV is found but SD Hit info for that has not found.");
                }
            }
            break;
    }
    return resultLV;
}

///// CopperFrame by Daehoon, HA!!! /////////////////////////
///// Modified by Jeewon at 2021. Mar.
///// type=3, real AMoRE-I geometry added

G4LogicalVolume *AmoreDetectorConstruction::Build_I_CopperFrameUpper(G4int aType) {
    G4LogicalVolume *resultLV;

    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4Box *sol_CopPlatehole;
    G4Box *sol_CopPlate;
    G4VSolid *sol_Copframe;
    G4Box *sol_CopBoxA;
    G4ThreeVector smallBlockPosition;

		// for type1
    G4double copperFramePunchSizeZ = copperFrameSizeZ * 1.1;
    G4double smallBlock1DistZ      = smallBlock1SizeZ / 2.;
    G4double smallBlock2DistZ      = smallBlock2SizeZ / 2.;
    G4Tubs *sol_Cophole = nullptr;
    G4Box *sol_CopBoxB = nullptr;

		// for type3, type2
		G4double type3_copperFramePunchSizeZ = type3_copperFrameSizeZ * 1.1;
		G4double type3_smallBlock1DistZ      = type3_smallBlock1SizeZ / 2.;
		G4Box *sol_CopPlatehole1 = nullptr;
		G4Box *sol_CopPlatehole2 = nullptr;
		G4Box *sol_CopPlatehole3 = nullptr;

    switch (aType) {
        case 1:
            resultLV =
                G4LogicalVolumeStore::GetInstance()->GetVolume("CopperFrame_Upper_Type1_LV", false);
            if (resultLV == nullptr) {
                sol_CopPlate = new G4Box("CopperFrame_Upper_Base", copperFrameSizeX / 2,
                                         copperFrameSizeY / 2, copperFrameSizeZ / 2);

                sol_CopPlatehole = new G4Box(
                    "CopperFrame_Upper_Punch", copperFrameSizeX / 2 - copperFrameXWidth,
                    copperFrameSizeY / 2. - copperFrameYWidth, copperFramePunchSizeZ / 2.);

                sol_Cophole = new G4Tubs("CopperFrame_Upper_hole", 0, copperFrameHoleDiameter / 2,
                                         copperFrameSizeZ * 2.2, 0, 360 * deg);

                sol_CopBoxA =
                    new G4Box("CopperFrame_Upper_SmallBoxA", smallBlock1SizeX / 2.,
                              smallBlock1SizeY / 2., smallBlock1SizeZ / 2. + solidBooleanTol);
                sol_CopBoxB =
                    new G4Box("CopperFrame_Upper_SmallBoxB", smallBlock2SizeX / 2.,
                              smallBlock2SizeY / 2., smallBlock2SizeZ / 2. + solidBooleanTol);

                smallBlockPosition =
                    G4ThreeVector(copperFrameSizeX / 2. + smallBlock1DistX,
                                  copperFrameSizeY / 2. + smallBlock1DistY,
                                  copperFrameSizeZ / 2. + smallBlock1DistZ - solidBooleanTol);

                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_TopBlockAtached1Solid", sol_CopPlate,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_TopBlockAtached2Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_TopBlockAtached3Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_TopBlockAtached4Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);

                smallBlockPosition =
                    G4ThreeVector(copperFrameSizeX / 2. - smallBlock2DistX,
                                  copperFrameSizeY / 2. - smallBlock2DistY,
                                  -copperFrameSizeZ / 2. - smallBlock2DistZ + solidBooleanTol);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_UnderBlockAtached1Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_UnderBlockAtached2Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_UnderBlockAtached3Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Upper_UnderBlockAtached4Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);

                smallBlockPosition = G4ThreeVector(copperFrameSizeX / 2. - copperFrameHoleDistance,
                                                   copperFrameSizeY / 2. - copperFrameHoleDistance);
                smallBlockPosition[0] *= -1;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Upper_HolePunched1Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Upper_HolePunched2Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Upper_HolePunched3Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Upper_HolePunched4Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);

                sol_Copframe = new G4SubtractionSolid("CopperFrame_Upper_Solid", sol_Copframe,
                                                      sol_CopPlatehole);

                resultLV = new G4LogicalVolume(sol_Copframe, _copper, "CopperFrame_Upper_Type1_LV");
                resultLV->SetVisAttributes(fI_CopperFrame_VisAttr);
            }
            break;

				case 2: //Real AMoRE-I geometry w/o holes
						resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
								"CopperFrame_Upper_Type3_LV", false);
						if(resultLV == nullptr)
						{
							// Main Copper Frame -------------------------------
                sol_CopPlate = new G4Box("CopperFrame_Upper_Base", 
										type3_copperFrameSizeX / 2, type3_copperFrameSizeY / 2, type3_copperFrameSizeZ / 2);

							// Small Box (4ea) ----------------------------------
                sol_CopBoxA = new G4Box("CopperFrame_Upper_SmallBoxA", 
										type3_smallBlockSizeX / 2., type3_smallBlockSizeY / 2., 
										type3_smallBlock1SizeZ / 2. + solidBooleanTol);

                smallBlockPosition = G4ThreeVector(
										type3_copperFrameSizeX / 2. + smallBlock1DistX,
                    type3_copperFrameSizeY / 2. + smallBlock1DistY,
                    type3_copperFrameSizeZ / 2. + type3_smallBlock1DistZ - solidBooleanTol);

                smallBlockPosition[0] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached1Solid", 
										sol_CopPlate, sol_CopBoxA, nullptr, smallBlockPosition);
								
                smallBlockPosition[0] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached2Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);

                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached3Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);

                smallBlockPosition[0] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached4Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);

                resultLV = new G4LogicalVolume(sol_Copframe, _copper, "CopperFrame_Upper_Type1_LV");
                resultLV->SetVisAttributes(fI_CopperFrame_VisAttr);
						}
						break;
        case 3: // REAL CopperFrame Upper with holes
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(
								"CopperFrame_Upper_Type3_LV", false);
            if (resultLV == nullptr) 
						{
							// Main Copper Frame -------------------------------
                sol_CopPlate = new G4Box("CopperFrame_Upper_Base", 
										type3_copperFrameSizeX / 2, type3_copperFrameSizeY / 2, type3_copperFrameSizeZ / 2);

						  // Frame hole --------------------------------------
                sol_CopPlatehole = new G4Box("CopperFrame_Upper_Punch", 
										type3_copperFrameSizeX / 2 - type3_copperFrameXWidth,
                    type3_copperFrameSizeY / 2. - type3_copperFrameYWidth, 
										type3_copperFramePunchSizeZ / 2.);
								sol_CopPlatehole1 = new G4Box("copperFrame_hole_small2", 
										type3_copperFrameHoleWidth / 2., 
										type3_copperFrameHole1SizeY / 2.,
										type3_copperFrameSizeZ / 2. + solidBooleanTol);
								sol_CopPlatehole2 = new G4Box("CopperFrame_hole_small1",
										type3_copperFrameHole2SizeX / 2., 
										type3_copperFrameHoleWidth /2.,
										type3_copperFrameSizeZ / 2. + solidBooleanTol);
								sol_CopPlatehole3 = new G4Box("copperFrame_hole_small3", 
										type3_copperFrameHole3SizeX / 2.,
										type3_copperFrameYWidth / 2.,
										type3_copperFrameSizeZ / 2. + solidBooleanTol);

								G4ThreeVector frameHole1Position = 
									G4ThreeVector(type3_copperFrameSizeX / 2. + smallBlock1DistX,
											          type3_copperFrameSizeY / 2. + smallBlock1DistY + type3_frameHole1Dist,0);
								sol_Copframe =  new G4SubtractionSolid("CopperFrame_Upper_withHole1",
											sol_CopPlate, sol_CopPlatehole1, nullptr, frameHole1Position);
								for(int i = 0; i < 3; i++)
								{
									frameHole1Position[0] *= -1.;
									sol_Copframe =  new G4SubtractionSolid("CopperFrame_Upper_withHole1",
											sol_Copframe, sol_CopPlatehole1, nullptr, frameHole1Position);
									if(i == 0) frameHole1Position[1] *= -1.;
								}

								G4ThreeVector frameHole2Position = G4ThreeVector(
										type3_copperFrameSizeX / 2. - smallBlock2DistX - type3_frameHole2Dist,
										type3_copperFrameSizeY / 2. - smallBlock2DistY,	0);
								for(int i = 0; i<4; i++)
								{
									sol_Copframe =  new G4SubtractionSolid("CopperFrame_Upper_withHole2",
											sol_Copframe, sol_CopPlatehole2, nullptr, frameHole2Position);
									frameHole2Position[0] *= -1.;
									if(i == 1)	frameHole2Position[1] *= -1.;
								}

                sol_Copframe = new G4SubtractionSolid("CopperFrame_Upper_Solid", 
										sol_Copframe, sol_CopPlatehole);
								
							// Small Box (4ea) ----------------------------------
                sol_CopBoxA = new G4Box("CopperFrame_Upper_SmallBoxA", 
										type3_smallBlockSizeX / 2., type3_smallBlockSizeY / 2., 
										type3_smallBlock1SizeZ / 2. + solidBooleanTol);

                smallBlockPosition = G4ThreeVector(
										type3_copperFrameSizeX / 2. + smallBlock1DistX,
                    type3_copperFrameSizeY / 2. + smallBlock1DistY,
                    type3_copperFrameSizeZ / 2. + type3_smallBlock1DistZ - solidBooleanTol);

                smallBlockPosition[0] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached1Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);
								
                smallBlockPosition[0] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached2Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);

                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached3Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);

                smallBlockPosition[0] *= -1.;
                sol_Copframe = new G4UnionSolid("CopperFrame_Upper_TopBlockAtached4Solid", 
										sol_Copframe, sol_CopBoxA, nullptr, smallBlockPosition);

                resultLV = new G4LogicalVolume(sol_Copframe, _copper, "CopperFrame_Upper_Type1_LV");
                resultLV->SetVisAttributes(fI_CopperFrame_VisAttr);

								// Clamp bolts(4ea) -----------------------

								G4Tubs *clampBoltsBodyTubs = new G4Tubs("clampBoltsBodyTubs",0,
										type3_photonPhononBoltsRadius - solidBooleanTol, 
										(type3_copperFrameSizeZ - solidBooleanTol) / 2., 0, 360. * deg);
								G4LogicalVolume *clampBoltsBodyLV = new G4LogicalVolume(clampBoltsBodyTubs, _brass, "clampBoltsBody_LV");
								clampBoltsBodyLV->SetVisAttributes(fI_BrassBolt_VisAttr);

								G4ThreeVector nowSolidPosition = {type3_copperFrameSizeX / 2. - type3_copperFrameXWidth / 2, 0,0};
									//-type3_copperFrameSizeZ / 2. - type3_photonFrameSizeZ - type3_photonFrameTopSpace};
								new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsBodyLV,
										"ClampBoltsForCrystal_PV", resultLV, false, 0, OverlapCheck);

								nowSolidPosition[0] *=-1;
								new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsBodyLV,
										"ClampBoltsForCrystal_PV", resultLV, false, 1, OverlapCheck);

								nowSolidPosition[0] += (type3_copperFrameSizeX - type3_copperFrameXWidth) / 2.;
								nowSolidPosition[1] += (type3_copperFrameSizeY - type3_copperFrameYWidth) / 2.;
								new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsBodyLV,
										"ClampBoltsForCrystal_PV", resultLV, false, 2, OverlapCheck);

								nowSolidPosition[1] *= -1;
								new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsBodyLV,
										"ClampBoltsForCrystal_PV", resultLV, false, 3, OverlapCheck);

                // photon-phonon bolts body(4ea) ----------------
								G4Tubs *photonPhononBoltsBody = new G4Tubs("photonPhononBoltsBodyTubs", 0,
										type3_photonPhononBoltsRadius, 
										type3_copperFrameSizeZ/ 2., 0, 360. * deg);
								G4LogicalVolume *photonPhononBoltsBodyLV = new G4LogicalVolume(photonPhononBoltsBody,
										_brass, "photonPhononBoltsBody_LV");
								photonPhononBoltsBodyLV->SetVisAttributes(fI_BrassBolt_VisAttr);

								nowSolidPosition = {(type3_photonFrameBaseWidthX - type3_photonFrameSizeX)/2.,
									(type3_photonFrameSizeY - type3_photonFrameBaseWidthY) / 2.,0};
									//- type3_photonFrameSizeZ - type3_photonFrameTopSpace - type3_copperFrameSizeZ /2 };
								new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsBodyLV,
										"photonPhononBoltsBody_PV", resultLV, false, 0, OverlapCheck);
								
								nowSolidPosition[0] *= -1;
								new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsBodyLV,
										"photonPhononBoltsBody_PV", resultLV, false, 1, OverlapCheck);
								nowSolidPosition[1] *= -1;
								new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsBodyLV,
										"photonPhononBoltsBody_PV", resultLV, false, 2, OverlapCheck);
								nowSolidPosition[0] *= -1;
								new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsBodyLV,
										"photonPhononBoltsBody_PV", resultLV, false, 3, OverlapCheck);
								
            }
            break;
    }
    return resultLV;
}

/////CopperFrame by Daehoon, HA!!! /////////////////////////
///// Modified by Jeewon at 2021. Mar. (real AMoRE_I geometry added)
///// type=2, w/o holes
///// type=3, with holes
G4LogicalVolume *AmoreDetectorConstruction::Build_I_CopperFrameBottom(G4int aType) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4double copperFramePunchSizeZ = (copperFrameSizeZ + bottomSupportDepth) * 1.1;
    G4double smallBlock1DistZ      = smallBlock1SizeZ / 2.;
    G4double smallBlock2DistZ      = smallBlock2SizeZ / 2.;

    G4double type3_copperFramePunchSizeZ = (type3_copperFrameSizeZ + bottomSupportDepth) * 1.1;
    G4double type3_smallBlock3DistZ      = type3_smallBlock3SizeZ / 2.;

    G4Tubs *sol_Cophole;
    G4Box *sol_CopPlatehole;
    G4Box *sol_CopPlate;
    G4VSolid *sol_Copframe;
    G4Box *sol_CopBoxA;
    G4Box *sol_CopBoxB ;
    G4Box *sol_CopDownDiskA;
    G4Box *sol_CopDownDiskB;

    G4ThreeVector smallBlockPosition;

    switch (aType) {
        case 1:
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume("CopperFrame_Bottom_Type1_LV",
                                                                      false);
            if (resultLV == nullptr) {
                sol_CopPlate = new G4Box("CopperFrame_Bottom_Base", copperFrameSizeX / 2,
                                         copperFrameSizeY / 2, copperFrameSizeZ / 2);

                sol_CopPlatehole = new G4Box(
                    "CopperFrame_Bottom_Punch", copperFrameSizeX / 2 - copperFrameXWidth,
                    copperFrameSizeY / 2. - copperFrameYWidth, copperFramePunchSizeZ / 2.);

                sol_Cophole = new G4Tubs("CopperFrame_Bottom_hole", 0, copperFrameHoleDiameter / 2,
                                         copperFrameSizeZ * 2.2, 0, 360 * deg);

                sol_CopBoxA =
                    new G4Box("CopperFrame_Bottom_SmallBoxA", smallBlock1SizeX / 2.,
                              smallBlock1SizeY / 2., smallBlock1SizeZ / 2. + solidBooleanTol);
                sol_CopBoxB =
                    new G4Box("CopperFrame_Bottom_SmallBoxB", smallBlock2SizeX / 2.,
                              smallBlock2SizeY / 2., smallBlock2SizeZ / 2. + solidBooleanTol);

                smallBlockPosition =
                    G4ThreeVector(copperFrameSizeX / 2. + smallBlock1DistX,
                                  copperFrameSizeY / 2. + smallBlock1DistY,
                                  -copperFrameSizeZ / 2. - smallBlock1DistZ + solidBooleanTol);

                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached1Solid", sol_CopPlate,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached2Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached3Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached4Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);

                smallBlockPosition =
                    G4ThreeVector(copperFrameSizeX / 2. - smallBlock2DistX,
                                  copperFrameSizeY / 2. - smallBlock2DistY,
                                  copperFrameSizeZ / 2. + smallBlock2DistZ - solidBooleanTol);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_UnderBlockAtached1Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_UnderBlockAtached2Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_UnderBlockAtached3Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_UnderBlockAtached4Solid", sol_Copframe,
                                     sol_CopBoxB, nullptr, smallBlockPosition);

                smallBlockPosition = G4ThreeVector(copperFrameSizeX / 2. - copperFrameHoleDistance,
                                                   copperFrameSizeY / 2. - copperFrameHoleDistance);
                smallBlockPosition[0] *= -1;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Bottom_HolePunched1Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Bottom_HolePunched2Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Bottom_HolePunched3Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4SubtractionSolid("CopperFrame_Bottom_HolePunched4Solid", sol_Copframe,
                                           sol_Cophole, nullptr, smallBlockPosition);

                sol_CopDownDiskA =
                    new G4Box("CopperFrame_Bottom_SupportPlateA_Solid", bottomSupportASizeX / 2.,
                              bottomSupportASizeY / 2., bottomSupportASizeZ / 2. + solidBooleanTol);
                sol_CopDownDiskB =
                    new G4Box("CopperFrame_Bottom_SupportPlateB_Solid", bottomSupportBSizeX / 2.,
                              bottomSupportBSizeY / 2., bottomSupportBSizeZ / 2. + solidBooleanTol);

                sol_Copframe = new G4UnionSolid(
                    "CopperFrame_Bottom_SupportA_Added_Solid", sol_Copframe, sol_CopDownDiskA,
                    nullptr,
                    {bottomSupportASizeX / 2 + bottomSupportADistX, 0,
                     -bottomSupportASizeZ / 2. - copperFrameSizeZ / 2. + solidBooleanTol});
                sol_Copframe = new G4UnionSolid(
                    "CopperFrame_Bottom_SupportB_Added_Solid", sol_Copframe, sol_CopDownDiskB,
                    nullptr,
                    {-bottomSupportBSizeX / 2 - bottomSupportBDistX, 0,
                     -bottomSupportBSizeZ / 2. - copperFrameSizeZ / 2. + solidBooleanTol});
                sol_Copframe = new G4SubtractionSolid(
                    "CopperFrame_Bottom_Solid", sol_Copframe, sol_CopPlatehole, nullptr,
                    {0, 0,
                     copperFramePunchSizeZ / 2. - copperFrameSizeZ / 2. - bottomSupportDepth});

                resultLV =
                    new G4LogicalVolume(sol_Copframe, _copper, "CopperFrame_Bottom_Type1_LV");
                resultLV->SetVisAttributes(fI_CopperFrame_VisAttr);
            }
            break;
        case 3:
            resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume("CopperFrame_Bottom_Type3_LV",
                                                                      false);
            if (resultLV == nullptr) {
                sol_CopPlate = new G4Box("CopperFrame_Bottom_Base", type3_copperFrameSizeX / 2,
                                         type3_copperFrameSizeY / 2, type3_copperFrameSizeZ / 2);

                sol_CopPlatehole = new G4Box(
                    "CopperFrame_Bottom_Punch", type3_copperFrameSizeX / 2 - type3_copperFrameXWidth,
                    type3_copperFrameSizeY / 2. - type3_copperFrameYWidth, type3_copperFramePunchSizeZ / 2.);

                sol_Cophole = new G4Tubs("CopperFrame_Bottom_hole", 0, type3_copperFrameHoleDiameter / 2,
                                         type3_copperFrameSizeZ * 2.2, 0, 360 * deg);

                sol_CopBoxA =
                    new G4Box("CopperFrame_Bottom_SmallBoxA", type3_smallBlockSizeX / 2.,
                              type3_smallBlockSizeY / 2., type3_smallBlock3SizeZ / 2. + solidBooleanTol);
                //sol_CopBoxB =
                    //new G4Box("CopperFrame_Bottom_SmallBoxB", smallBlock2SizeX / 2.,
                              //smallBlock2SizeY / 2., smallBlock2SizeZ / 2. + solidBooleanTol);

                smallBlockPosition =
                    G4ThreeVector(type3_copperFrameSizeX / 2. + smallBlock1DistX,
                                  type3_copperFrameSizeY / 2. + smallBlock1DistY,
                                  -type3_copperFrameSizeZ / 2. - type3_smallBlock3DistZ + solidBooleanTol);

                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached1Solid", sol_CopPlate,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached2Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                smallBlockPosition[1] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached3Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);
                smallBlockPosition[0] *= -1.;
                sol_Copframe =
                    new G4UnionSolid("CopperFrame_Bottom_TopBlockAtached4Solid", sol_Copframe,
                                     sol_CopBoxA, nullptr, smallBlockPosition);

                sol_CopDownDiskA =
                    new G4Box("CopperFrame_Bottom_SupportPlateA_Solid", type3_bottomSupportASizeX / 2.,
                              type3_bottomSupportASizeY / 2., type3_bottomSupportASizeZ / 2. + solidBooleanTol);
                sol_CopDownDiskB =
                    new G4Box("CopperFrame_Bottom_SupportPlateB_Solid", type3_bottomSupportBSizeX / 2.,
                              type3_bottomSupportBSizeY / 2., type3_bottomSupportBSizeZ / 2. + solidBooleanTol);

                sol_Copframe = new G4UnionSolid(
                    "CopperFrame_Bottom_SupportA_Added_Solid", sol_Copframe, sol_CopDownDiskA,
                    nullptr,
                    {type3_copperFrameSizeX / 2. - type3_smallBlockSizeX - type3_bottomSupportADistX, 0,
                     -type3_bottomSupportASizeZ / 2. - type3_copperFrameSizeZ / 2. + solidBooleanTol});

                sol_Copframe = new G4UnionSolid(
                    "CopperFrame_Bottom_SupportB_Added_Solid", sol_Copframe, sol_CopDownDiskB,
                    nullptr,
                    {-type3_copperFrameSizeX / 2. + type3_smallBlockSizeX + type3_bottomSupportBDistX, 0,
                     -type3_bottomSupportBSizeZ / 2. - type3_copperFrameSizeZ / 2. + solidBooleanTol});

                sol_Copframe = new G4SubtractionSolid(
                    "CopperFrame_Bottom_Solid", sol_Copframe, sol_CopPlatehole, nullptr,
                    {0, 0,
                     type3_copperFramePunchSizeZ / 2. - type3_copperFrameSizeZ / 2. - bottomSupportDepth});

                resultLV =
                    new G4LogicalVolume(sol_Copframe, _copper, "CopperFrame_Bottom_Type1_LV");
                resultLV->SetVisAttributes(fI_CopperFrame_VisAttr);
            }
            break;
    }

    return resultLV;
}

G4LogicalVolume *
    AmoreDetectorConstruction::Build_I_TopScintillator(G4Material *aScintMaterial,
                                                       G4Material *aReflectorMaterial,
                                                       G4SurfaceProperty *aSurfaceProperty) {

    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    using namespace AmoreDetectorStaticInfo;

    constexpr const char *resultLVName = "MuonTopScintillator_LV";

    G4double topScint_trapFlatL_half  = topScint_trapZ_half * scint_FlatL_ratio;
    G4double topScint_trapThick1_half = scint_Thick_half;
    G4double topScint_trapThick2_half = topScint_trapX1_half;
    G4double topScint_boxThick_half   = scint_Thick_half;
    G4double topScint_TrapSlope_double =
        (topScint_trapX2_half - topScint_trapX1_half) / topScint_trapZ_half;

    G4VSolid *topScint_MotherSolid;
    G4Box *topScint_Mother_Box;
    G4Box *topScint_Box;
    G4VSolid *topScint_MotherTrap;
    G4Trd *topScint_MotherPMTTrap;
    G4Trd *topScint_FlatTrap;
    G4Trd *topScint_PMTTrap;

    G4RotationMatrix *scint_alignMtx_Uside = new G4RotationMatrix;
    G4RotationMatrix *scint_alignMtx_Dside = new G4RotationMatrix;

    resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(resultLVName, false);

    if (resultLV == nullptr) {
        scint_alignMtx_Uside->rotateX(-90 * deg);
        scint_alignMtx_Dside->rotateX(90 * deg);

        topScint_Mother_Box =
            new G4Box("Muon_Top_Scint_Mother_Box", topScint_boxWidth_half + scint_reflector_thick,
                      topScint_boxHeight_half + scint_reflector_thick,
                      topScint_boxThick_half + scint_reflector_thick);

        topScint_MotherTrap = new G4Trd(
            "Muon_Top_Scintillator_MS_PreTrap1",
            topScint_trapX2_half -
                topScint_TrapSlope_double * (topScint_trapFlatL_half + solidBooleanTolSmall / 2.),
            topScint_trapX2_half + topScint_TrapSlope_double * solidBooleanTolSmall / 2.,
            topScint_trapThick1_half + scint_reflector_thick,
            topScint_trapThick1_half + scint_reflector_thick,
            topScint_trapFlatL_half + solidBooleanTolSmall);
        topScint_MotherPMTTrap =
            new G4Trd("Muon_Top_Scintillator_MS_PreTrap2",
                      topScint_trapX1_half + topScint_TrapSlope_double * scint_reflector_thick / 2.,
                      topScint_trapX2_half - topScint_TrapSlope_double * topScint_trapFlatL_half,
                      topScint_trapThick2_half + scint_reflector_thick,
                      topScint_trapThick1_half + scint_reflector_thick,
                      topScint_trapZ_half - topScint_trapFlatL_half - scint_reflector_thick / 2.);
        topScint_MotherTrap = new G4UnionSolid(
            "Muon_Top_Scintillator_MS_Trap", topScint_MotherTrap, topScint_MotherPMTTrap, 0,
            G4ThreeVector(0, 0, -topScint_trapZ_half + scint_reflector_thick / 2.));

        topScint_MotherSolid = new G4UnionSolid(
            "Muon_Top_Mother_Solid_Temp1", topScint_Mother_Box, topScint_MotherTrap,
            scint_alignMtx_Uside,
            G4ThreeVector(-topScint_trapX2_half,
                          topScint_boxHeight_half + scint_reflector_thick + topScint_trapFlatL_half,
                          0));
        topScint_MotherSolid = new G4UnionSolid(
            "Muon_Top_Mother_Solid_Temp2", topScint_MotherSolid, topScint_MotherTrap,
            scint_alignMtx_Uside,
            G4ThreeVector(+topScint_trapX2_half,
                          topScint_boxHeight_half + scint_reflector_thick + topScint_trapFlatL_half,
                          0));
        topScint_MotherSolid = new G4UnionSolid(
            "Muon_Top_Mother_Mother_Solid_Temp3", topScint_MotherSolid, topScint_MotherTrap,
            scint_alignMtx_Dside,
            G4ThreeVector(
                -topScint_trapX2_half,
                -topScint_boxHeight_half - scint_reflector_thick - topScint_trapFlatL_half, 0));
        topScint_MotherSolid = new G4UnionSolid(
            "Muon_Top_Mother_Solid", topScint_MotherSolid, topScint_MotherTrap,
            scint_alignMtx_Dside,
            G4ThreeVector(
                +topScint_trapX2_half,
                -topScint_boxHeight_half - scint_reflector_thick - topScint_trapFlatL_half, 0));

        topScint_Box = new G4Box("Muon_Top_Scintillator_Box", topScint_boxWidth_half,
                                 topScint_boxHeight_half, topScint_boxThick_half);
        topScint_FlatTrap =
            new G4Trd("Muon_Top_Scintillator_Flat_Trd",
                      topScint_trapX2_half - topScint_TrapSlope_double * topScint_trapFlatL_half,
                      topScint_trapX2_half, topScint_trapThick1_half, topScint_trapThick1_half,
                      topScint_trapFlatL_half);
        topScint_PMTTrap =
            new G4Trd("Muon_Top_Scintillator_PMT_Trd", topScint_trapX1_half,
                      topScint_trapX2_half - (topScint_trapX2_half - topScint_trapX1_half) /
                                                 topScint_trapZ_half * topScint_trapFlatL_half,
                      topScint_trapThick2_half, topScint_trapThick1_half,
                      topScint_trapZ_half - topScint_trapFlatL_half);

        resultLV = new G4LogicalVolume(topScint_MotherSolid, aReflectorMaterial, resultLVName);
        new G4LogicalSkinSurface("Muon_Top_Scintillator_reflector_opsurf", resultLV,
                                 aSurfaceProperty);

        fI_TopScint_BoxLogical =
            new G4LogicalVolume(topScint_Box, aScintMaterial, "Muon_Top_Scintillator_Box_LV");
        fI_TopScint_PMTTrapLogical  = new G4LogicalVolume(topScint_PMTTrap, aScintMaterial,
                                                         "Muon_Top_Scintillator_PMT_Trd_LV");
        fI_TopScint_FlatTrapLogical = new G4LogicalVolume(topScint_FlatTrap, aScintMaterial,
                                                          "Muon_Top_Scintillator_Flat_Trd_LV");

        new G4PVPlacement(scint_alignMtx_Uside,
                          G4ThreeVector(-topScint_trapX2_half,
                                        topScint_boxHeight_half + topScint_trapFlatL_half, 0),
                          fI_TopScint_FlatTrapLogical, "Muon_Top_Scintillator_FlatTrd_PV", resultLV,
                          false, 0, OverlapCheck);
        new G4PVPlacement(
            scint_alignMtx_Uside,
            G4ThreeVector(-topScint_trapX2_half,
                          topScint_boxHeight_half + topScint_trapZ_half + topScint_trapFlatL_half,
                          0),
            fI_TopScint_PMTTrapLogical, "Muon_Top_Scintillator_PMT_Trd_PV", resultLV, false, 0,
            OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Uside,
                          G4ThreeVector(topScint_trapX2_half,
                                        topScint_boxHeight_half + topScint_trapFlatL_half, 0),
                          fI_TopScint_FlatTrapLogical, "Muon_Top_Scintillator_FlatTrd_PV", resultLV,
                          false, 1, OverlapCheck);
        new G4PVPlacement(
            scint_alignMtx_Uside,
            G4ThreeVector(topScint_trapX2_half,
                          topScint_boxHeight_half + topScint_trapZ_half + topScint_trapFlatL_half,
                          0),
            fI_TopScint_PMTTrapLogical, "Muon_Top_Scintillator_PMT_Trd_PV", resultLV, false, 1,
            OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Dside,
                          G4ThreeVector(-topScint_trapX2_half,
                                        -topScint_boxHeight_half - topScint_trapFlatL_half, 0),
                          fI_TopScint_FlatTrapLogical, "Muon_Top_Scintillator_FlatTrd_PV", resultLV,
                          false, 2, OverlapCheck);
        new G4PVPlacement(
            scint_alignMtx_Dside,
            G4ThreeVector(-topScint_trapX2_half,
                          -topScint_boxHeight_half - topScint_trapZ_half - topScint_trapFlatL_half,
                          0),
            fI_TopScint_PMTTrapLogical, "Muon_Top_Scintillator_PMT_Trd_PV", resultLV, false, 2,
            OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Dside,
                          G4ThreeVector(topScint_trapX2_half,
                                        -topScint_boxHeight_half - topScint_trapFlatL_half, 0),
                          fI_TopScint_FlatTrapLogical, "Muon_Top_Scintillator_FlatTrd_PV", resultLV,
                          false, 3, OverlapCheck);
        new G4PVPlacement(
            scint_alignMtx_Dside,
            G4ThreeVector(topScint_trapX2_half,
                          -topScint_boxHeight_half - topScint_trapZ_half - topScint_trapFlatL_half,
                          0),
            fI_TopScint_PMTTrapLogical, "Muon_Top_Scintillator_PMT_Trd_PV", resultLV, false, 3,
            OverlapCheck);
    }

    return resultLV;
}

G4LogicalVolume *
    AmoreDetectorConstruction::Build_I_Side_FBScintillator(G4Material *aScintMaterial,
                                                           G4Material *aReflectorMaterial,
                                                           G4SurfaceProperty *aSurfaceProperty) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    using namespace AmoreDetectorStaticInfo;
    constexpr const char *resultLVName = "MuonSideScintillatorFB_LV";

    G4double sideFBScint_trapFlatL_half  = sideFBScint_trapZ_half * scint_FlatL_ratio;
    G4double sideFBScint_trapThick1_half = scint_Thick_half;
    G4double sideFBScint_trapThick2_half = sideFBScint_trapX1_half;
    G4double sideFBScint_boxThick_half   = scint_Thick_half;
    G4double sideFBScint_TrapSlope_double =
        (sideFBScint_trapX2_half - sideFBScint_trapX1_half) / sideFBScint_trapZ_half;

    G4VSolid *sideFBScint_MotherSolid;
    G4Box *sideFBScint_Mother_Box;
    G4VSolid *sideFBScint_MotherTrap;
    G4Trd *sideFBScint_MotherPMTTrap;
    G4Box *sideFBScint_Box;
    G4Trd *sideFBScint_FlatTrap;
    G4Trd *sideFBScint_PMTTrap;

    G4RotationMatrix *scintMother_alignMtx_FBSide = new G4RotationMatrix;
    G4RotationMatrix *scint_alignMtx_Uside        = new G4RotationMatrix;
    resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(resultLVName, false);
    if (resultLV == nullptr) {
        scint_alignMtx_Uside->rotateX(-90 * deg);
        scintMother_alignMtx_FBSide->rotateY(-90 * deg);
        scintMother_alignMtx_FBSide->rotateZ(-90 * deg);

        sideFBScint_Mother_Box = new G4Box("Muon_FBSide_Scintillator_Mother_Box",
                                           sideFBScint_boxWidth_half + scint_reflector_thick,
                                           sideFBScint_boxHeight_half + scint_reflector_thick,
                                           sideFBScint_boxThick_half + scint_reflector_thick);
        sideFBScint_MotherTrap = new G4Trd(
            "Muon_FBSide_Scintillator_MS_PreTrap1",
            sideFBScint_trapX2_half - sideFBScint_TrapSlope_double *
                                          (sideFBScint_trapFlatL_half + solidBooleanTolSmall / 2.),
            sideFBScint_trapX2_half + sideFBScint_TrapSlope_double * solidBooleanTolSmall / 2.,
            sideFBScint_trapThick1_half + scint_reflector_thick,
            sideFBScint_trapThick1_half + scint_reflector_thick,
            sideFBScint_trapFlatL_half + solidBooleanTolSmall);
        sideFBScint_MotherPMTTrap = new G4Trd(
            "Muon_FBSide_Scintillator_MS_PreTrap2",
            sideFBScint_trapX1_half + sideFBScint_TrapSlope_double * scint_reflector_thick / 2.,
            sideFBScint_trapX2_half - sideFBScint_TrapSlope_double * sideFBScint_trapFlatL_half,
            sideFBScint_trapThick2_half + scint_reflector_thick,
            sideFBScint_trapThick1_half + scint_reflector_thick,
            sideFBScint_trapZ_half - sideFBScint_trapFlatL_half - scint_reflector_thick / 2.);
        sideFBScint_MotherTrap = new G4UnionSolid(
            "Muon_FBSide_Scintillator_MS_Trap", sideFBScint_MotherTrap, sideFBScint_MotherPMTTrap,
            0, G4ThreeVector(0, 0, -sideFBScint_trapZ_half + scint_reflector_thick / 2.));

        sideFBScint_MotherSolid =
            new G4UnionSolid("Muon_FBSide_Mother_Solid_Temp", sideFBScint_Mother_Box,
                             sideFBScint_MotherTrap, scint_alignMtx_Uside,
                             G4ThreeVector(sideFBScint_trapX2_half,
                                           sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half +
                                               scint_reflector_thick,
                                           0));
        sideFBScint_MotherSolid =
            new G4UnionSolid("Muon_FBSide_Mother_Solid", sideFBScint_MotherSolid,
                             sideFBScint_MotherTrap, scint_alignMtx_Uside,
                             G4ThreeVector(-sideFBScint_trapX2_half,
                                           sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half +
                                               scint_reflector_thick,
                                           0));

        sideFBScint_Box      = new G4Box("Muon_Side_Scintillator_FB_Box", sideFBScint_boxWidth_half,
                                    sideFBScint_boxHeight_half, sideFBScint_boxThick_half);
        sideFBScint_FlatTrap = new G4Trd("Muon_Side_Scintillator_FB_Flat_Trd",
                                         sideFBScint_trapX2_half - sideFBScint_TrapSlope_double *
                                                                       sideFBScint_trapFlatL_half,
                                         sideFBScint_trapX2_half, sideFBScint_trapThick1_half,
                                         sideFBScint_trapThick1_half, sideFBScint_trapFlatL_half);
        sideFBScint_PMTTrap  = new G4Trd(
            "Muon_Side_Scintillator_FB_PMT_Trd", sideFBScint_trapX1_half,
            sideFBScint_trapX2_half - sideFBScint_TrapSlope_double * sideFBScint_trapFlatL_half,
            sideFBScint_trapThick2_half, sideFBScint_trapThick1_half,
            sideFBScint_trapZ_half - sideFBScint_trapFlatL_half);

        resultLV = new G4LogicalVolume(sideFBScint_MotherSolid, aReflectorMaterial, resultLVName);
        new G4LogicalSkinSurface("Muon_FBSide_Scintillator_reflector_opsurf", resultLV,
                                 aSurfaceProperty);

        fI_SideFBScint_BoxLogical      = new G4LogicalVolume(sideFBScint_Box, aScintMaterial,
                                                        "Muon_Side_Scintillator_FB_Box_LV");
        fI_SideFBScint_PMTTrapLogical  = new G4LogicalVolume(sideFBScint_PMTTrap, aScintMaterial,
                                                            "Muon_Side_Scintillator_FB_PMT_Trd_LV");
        fI_SideFBScint_FlatTrapLogical = new G4LogicalVolume(
            sideFBScint_FlatTrap, aScintMaterial, "Muon_Side_Scintillator_FB_Flat_Trd_LV");

        new G4PVPlacement(nullptr, {}, fI_SideFBScint_BoxLogical, "Muon_Scintillator_FB_Box_PV",
                          resultLV, false, 0, OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Uside,
                          G4ThreeVector(-sideFBScint_trapX2_half,
                                        sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half, 0),
                          fI_SideFBScint_FlatTrapLogical, "Muon_Scintillator_FB_FlatTrd_PV",
                          resultLV, false, 0, OverlapCheck);
        new G4PVPlacement(scint_alignMtx_Uside,
                          G4ThreeVector(-sideFBScint_trapX2_half,
                                        sideFBScint_boxHeight_half + sideFBScint_trapZ_half +
                                            sideFBScint_trapFlatL_half,
                                        0),
                          fI_SideFBScint_PMTTrapLogical, "Muon_Scintillator_FB_PMT_Trd_PV",
                          resultLV, false, 0, OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Uside,
                          G4ThreeVector(sideFBScint_trapX2_half,
                                        sideFBScint_boxHeight_half + sideFBScint_trapFlatL_half, 0),
                          fI_SideFBScint_FlatTrapLogical, "Muon_Scintillator_FB_FlatTrd_PV",
                          resultLV, false, 1, OverlapCheck);
        new G4PVPlacement(scint_alignMtx_Uside,
                          G4ThreeVector(sideFBScint_trapX2_half,
                                        sideFBScint_boxHeight_half + sideFBScint_trapZ_half +
                                            sideFBScint_trapFlatL_half,
                                        0),
                          fI_SideFBScint_PMTTrapLogical, "Muon_Scintillator_FB_PMT_Trd_PV",
                          resultLV, false, 1, OverlapCheck);
    }
    return resultLV;
}

G4LogicalVolume *
    AmoreDetectorConstruction::Build_I_Side_LRScintillator(G4Material *aScintMaterial,
                                                           G4Material *aReflectorMaterial,
                                                           G4SurfaceProperty *aSurfaceProperty) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    constexpr const char *resultLVName = "MuonSideScintillatorLR_LV";

    G4Box *sideLRScint_Mother_Box;
    G4VSolid *sideLRScint_MotherTrap;
    G4VSolid *sideLRScint_MotherSolid;
    G4Trd *sideLRScint_MotherPMTTrap;
    G4Box *sideLRScint_Box;
    G4Trd *sideLRScint_FlatTrap;
    G4Trd *sideLRScint_PMTTrap;

    G4double sideLRScint_trapFlatL_half  = sideLRScint_trapZ_half * scint_FlatL_ratio;
    G4double sideLRScint_trapThick1_half = scint_Thick_half;
    G4double sideLRScint_trapThick2_half = sideLRScint_trapX1_half;
    G4double sideLRScint_boxThick_half   = scint_Thick_half;
    G4double sideLRScint_TrapSlope_double =
        (sideLRScint_trapX2_half - sideLRScint_trapX1_half) / sideLRScint_trapZ_half;

    G4RotationMatrix *scint_alignMtx_Rside;

    resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(resultLVName, false);
    if (resultLV == nullptr) {
        scint_alignMtx_Rside = new G4RotationMatrix;

        scint_alignMtx_Rside->rotateX(90 * deg);
        scint_alignMtx_Rside->rotateY(-90 * deg);

        sideLRScint_Mother_Box = new G4Box("Muon_LRSide_Scintillator_Mother_Box",
                                           sideLRScint_boxWidth_half + scint_reflector_thick,
                                           sideLRScint_boxHeight_half + scint_reflector_thick,
                                           sideLRScint_boxThick_half + scint_reflector_thick);
        sideLRScint_MotherTrap = new G4Trd(
            "Muon_LRSide_Scintillator_MS_PreTrap1",
            sideLRScint_trapX2_half - sideLRScint_TrapSlope_double *
                                          (sideLRScint_trapFlatL_half + solidBooleanTolSmall / 2.),
            sideLRScint_trapX2_half + sideLRScint_TrapSlope_double * solidBooleanTolSmall / 2.,
            sideLRScint_trapThick1_half + scint_reflector_thick,
            sideLRScint_trapThick1_half + scint_reflector_thick,
            sideLRScint_trapFlatL_half + solidBooleanTolSmall);
        sideLRScint_MotherPMTTrap = new G4Trd(
            "Muon_LRSide_Scintillator_MS_PreTrap2",
            sideLRScint_trapX1_half + sideLRScint_TrapSlope_double * scint_reflector_thick / 2.,
            sideLRScint_trapX2_half - sideLRScint_TrapSlope_double * sideLRScint_trapFlatL_half,
            sideLRScint_trapThick2_half + scint_reflector_thick,
            sideLRScint_trapThick1_half + scint_reflector_thick,
            sideLRScint_trapZ_half - sideLRScint_trapFlatL_half - scint_reflector_thick / 2.);
        sideLRScint_MotherTrap = new G4UnionSolid(
            "Muon_LRSide_Scintillator_MS_Trap", sideLRScint_MotherTrap, sideLRScint_MotherPMTTrap,
            0, G4ThreeVector(0, 0, -sideLRScint_trapZ_half + scint_reflector_thick / 2.));

        sideLRScint_MotherSolid = new G4UnionSolid(
            "Muon_LRSide_Scintillator_Solid_Temp1", sideLRScint_Mother_Box, sideLRScint_MotherTrap,
            scint_alignMtx_Rside,
            G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                              scint_reflector_thick,
                          -sideLRScint_boxHeight_half + 1. * sideLRScint_trapX2_half, 0));
        sideLRScint_MotherSolid = new G4UnionSolid(
            "Muon_LRSide_Scintillator_Solid_Temp2", sideLRScint_MotherSolid, sideLRScint_MotherTrap,
            scint_alignMtx_Rside,
            G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                              scint_reflector_thick,
                          -sideLRScint_boxHeight_half + 3. * sideLRScint_trapX2_half, 0));
        sideLRScint_MotherSolid = new G4UnionSolid(
            "Muon_LRSide_Scintillator_Solid", sideLRScint_MotherSolid, sideLRScint_MotherTrap,
            scint_alignMtx_Rside,
            G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                              scint_reflector_thick,
                          -sideLRScint_boxHeight_half + 5. * sideLRScint_trapX2_half, 0));

        sideLRScint_Box      = new G4Box("Muon_Side_Scintillator_LR_Box", sideLRScint_boxWidth_half,
                                    sideLRScint_boxHeight_half, sideLRScint_boxThick_half);
        sideLRScint_FlatTrap = new G4Trd("Muon_Side_Scintillator_LR_Flat_Trd",
                                         sideLRScint_trapX2_half - sideLRScint_TrapSlope_double *
                                                                       sideLRScint_trapFlatL_half,
                                         sideLRScint_trapX2_half, sideLRScint_trapThick1_half,
                                         sideLRScint_trapThick1_half, sideLRScint_trapFlatL_half);
        sideLRScint_PMTTrap  = new G4Trd(
            "Muon_Side_Scintillator_LR_PMT_Trd", sideLRScint_trapX1_half,
            sideLRScint_trapX2_half - sideLRScint_TrapSlope_double * sideLRScint_trapFlatL_half,
            sideLRScint_trapThick2_half, sideLRScint_trapThick1_half,
            sideLRScint_trapZ_half - sideLRScint_trapFlatL_half);

        resultLV = new G4LogicalVolume(sideLRScint_MotherSolid, aReflectorMaterial, resultLVName);
        new G4LogicalSkinSurface("Muon_LRSide_Scintillator_reflector_opsurf", resultLV,
                                 aSurfaceProperty);

        fI_SideLRScint_BoxLogical      = new G4LogicalVolume(sideLRScint_Box, aScintMaterial,
                                                        "Muon_Side_Scintillator_LR_Box_LV");
        fI_SideLRScint_PMTTrapLogical  = new G4LogicalVolume(sideLRScint_PMTTrap, aScintMaterial,
                                                            "Muon_Side_Scintillator_LR_PMT_Trd_LV");
        fI_SideLRScint_FlatTrapLogical = new G4LogicalVolume(
            sideLRScint_FlatTrap, aScintMaterial, "Muon_Side_Scintillator_LR_Flat_Trd_LV");

        new G4PVPlacement(nullptr, {}, fI_SideLRScint_BoxLogical, "Muon_Scintillator_LR_Box_PV",
                          resultLV, false, 0, OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Rside,
                          G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half,
                                        -sideLRScint_boxHeight_half + sideLRScint_trapX2_half, 0),
                          fI_SideLRScint_FlatTrapLogical, "Muon_Scintillator_LR_FlatTrd_PV",
                          resultLV, false, 0, OverlapCheck);
        new G4PVPlacement(scint_alignMtx_Rside,
                          G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                                            sideLRScint_trapZ_half,
                                        -sideLRScint_boxHeight_half + sideLRScint_trapX2_half, 0),
                          fI_SideLRScint_PMTTrapLogical, "Muon_Scintillator_LR_PMT_Trd_PV",
                          resultLV, false, 0, OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Rside,
                          G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half,
                                        -sideLRScint_boxHeight_half + 3. * sideLRScint_trapX2_half,
                                        0),
                          fI_SideLRScint_FlatTrapLogical, "Muon_Scintillator_LR_FlatTrd_PV",
                          resultLV, false, 1, OverlapCheck);
        new G4PVPlacement(scint_alignMtx_Rside,
                          G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                                            sideLRScint_trapZ_half,
                                        -sideLRScint_boxHeight_half + 3. * sideLRScint_trapX2_half,
                                        0),
                          fI_SideLRScint_PMTTrapLogical, "Muon_Scintillator_LR_PMT_Trd_PV",
                          resultLV, false, 1, OverlapCheck);

        new G4PVPlacement(scint_alignMtx_Rside,
                          G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half,
                                        -sideLRScint_boxHeight_half + 5. * sideLRScint_trapX2_half,
                                        0),
                          fI_SideLRScint_FlatTrapLogical, "Muon_Scintillator_LR_FlatTrd_PV",
                          resultLV, false, 2, OverlapCheck);
        new G4PVPlacement(scint_alignMtx_Rside,
                          G4ThreeVector(-sideLRScint_boxWidth_half - sideLRScint_trapFlatL_half -
                                            sideLRScint_trapZ_half,
                                        -sideLRScint_boxHeight_half + 5. * sideLRScint_trapX2_half,
                                        0),
                          fI_SideLRScint_PMTTrapLogical, "Muon_Scintillator_LR_PMT_Trd_PV",
                          resultLV, false, 2, OverlapCheck);
    }
    return resultLV;
}

G4LogicalVolume *AmoreDetectorConstruction::Build_I_Muffler_FBScintillator(
    G4Material *aScintMaterial, G4Material *aReflectorMaterial, G4Material *aSideReflectorMaterial,
    G4SurfaceProperty *aSurfaceProperty, G4SurfaceProperty *aSideReflectorSurfProp) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    using namespace AmoreDetectorStaticInfo;
    constexpr const char *resultLVName = "MuonMufflerScintillatorFB_LV";

    G4double mufflerFBScint_trapFlatL_half  = mufflerFBScint_trapZ_half * scint_muffler_FlatL_ratio;
    G4double mufflerFBScint_trapThick2_half = scint_muffler_Thick_half;
    G4double mufflerFBScint_boxThick_half   = scint_muffler_Thick_half;
    G4double mufflerFBScint_TrapSlope_double =
        (mufflerFBScint_trapX2_half - mufflerFBScint_trapX1_half) / mufflerFBScint_trapZ_half;
    G4bool isNotAllFlat = (scint_muffler_FlatL_ratio != 1);
    G4double mufflerFBScint_trapThick1_half =
        (isNotAllFlat) ? mufflerFBScint_trapX1_half : scint_muffler_Thick_half;
    G4bool isMufflerFBScint_SideReflector_attached;

    G4VSolid *mufflerFBScint_MotherSolid;
    G4Box *mufflerFBScint_Mother_Box;
    G4VSolid *mufflerFBScint_MotherTrap;
    G4Trd *mufflerFBScint_MotherPMTTrap;
    G4Box *mufflerFBScint_Box;
    G4Trd *mufflerFBScint_FlatTrap;
    G4Trd *mufflerFBScint_FlatTrapSideRefl;
    G4Trd *mufflerFBScint_PMTTrap;
    G4Trd *mufflerFBScint_PMTTrapSideRefl;
    G4Box *mufflerFBScint_boxSideReflector;

    G4LogicalVolume *mufflerFBScint_boxSideReflLV;
    G4LogicalVolume *mufflerFBScint_flatTrapSideReflLV;
    G4LogicalVolume *mufflerFBScint_pmtTrapSideReflLV;

    if (scint_muffler_sideRefl_thick > scint_reflector_thick) {
        G4Exception(__PRETTY_FUNCTION__, "SCINT_SIDEREFL_01", JustWarning,
                    "The thickness of a side reflector in muon detector at muffler FB side is too "
                    "thick. The reflector will not be attached.");
        isMufflerFBScint_SideReflector_attached = false;
    } else
        isMufflerFBScint_SideReflector_attached = true;

    G4RotationMatrix *scintMother_alignMtx_FBSide = new G4RotationMatrix;
    G4RotationMatrix *scint_alignMtx_Uside        = new G4RotationMatrix;
    resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(resultLVName, false);
    if (resultLV == nullptr) {
        scint_alignMtx_Uside->rotateX(-90 * deg);
        scintMother_alignMtx_FBSide->rotateY(-90 * deg);
        scintMother_alignMtx_FBSide->rotateZ(-90 * deg);

        mufflerFBScint_Mother_Box = new G4Box("Muon_FBMuffler_Scintillator_Mother_Box",
                                              mufflerFBScint_boxWidth_half + scint_reflector_thick,
                                              mufflerFBScint_boxHeight_half + scint_reflector_thick,
                                              mufflerFBScint_boxThick_half + scint_reflector_thick);
        if (isNotAllFlat) {
            mufflerFBScint_MotherTrap = new G4Trd(
                "Muon_FBMuffler_Scintillator_MS_PreTrap1",
                mufflerFBScint_trapX2_half -
                    mufflerFBScint_TrapSlope_double * // Slope(기울기) 값이 원래 값의 2배
                                                      // (Double)임에 주의할 것
                        (mufflerFBScint_trapFlatL_half +
                         solidBooleanTolSmall / 2.), // 해당 기울기 값 (2배임) 에 맞게 계산된 것
                mufflerFBScint_trapX2_half + mufflerFBScint_TrapSlope_double *
                                                 solidBooleanTolSmall /
                                                 2., // 해당 기울기 값 (2배임) 에 맞게 계산된 것
                mufflerFBScint_trapThick2_half + scint_reflector_thick,
                mufflerFBScint_trapThick2_half + scint_reflector_thick,
                mufflerFBScint_trapFlatL_half + solidBooleanTolSmall);
            mufflerFBScint_MotherPMTTrap =
                new G4Trd("Muon_FBMuffler_Scintillator_MS_PreTrap2",
                          mufflerFBScint_trapX1_half +
                              mufflerFBScint_TrapSlope_double * scint_reflector_thick / 2.,
                          mufflerFBScint_trapX2_half -
                              mufflerFBScint_TrapSlope_double * mufflerFBScint_trapFlatL_half,
                          mufflerFBScint_trapThick1_half + scint_reflector_thick,
                          mufflerFBScint_trapThick2_half + scint_reflector_thick,
                          mufflerFBScint_trapZ_half - mufflerFBScint_trapFlatL_half -
                              scint_reflector_thick / 2.);
            mufflerFBScint_MotherTrap = new G4UnionSolid(
                "Muon_FBMuffler_Scintillator_MS_Trap", mufflerFBScint_MotherTrap,
                mufflerFBScint_MotherPMTTrap, nullptr,
                G4ThreeVector(0, 0, -mufflerFBScint_trapZ_half + scint_reflector_thick / 2.));

            mufflerFBScint_MotherSolid =
                new G4UnionSolid("Muon_FBMuffler_Mother_Solid", mufflerFBScint_Mother_Box,
                                 mufflerFBScint_MotherTrap, scint_alignMtx_Uside,
                                 {0,
                                  mufflerFBScint_boxHeight_half + mufflerFBScint_trapFlatL_half +
                                      scint_reflector_thick,
                                  0});
        } else {
            mufflerFBScint_MotherTrap =
                new G4Trd("Muon_FBMuffler_Scintillator_MS_PreTrap1",
                          mufflerFBScint_trapX1_half +
                              mufflerFBScint_TrapSlope_double * // Slope(기울기) 값이 원래 값의 2배
                                                                // (Double)임에 주의할 것
                                  (scint_reflector_thick + solidBooleanTolSmall) / 2.,
                          mufflerFBScint_trapX2_half +
                              mufflerFBScint_TrapSlope_double * solidBooleanTolSmall / 2.,
                          mufflerFBScint_trapThick1_half +
                              scint_reflector_thick * mufflerFBScint_TrapSlope_double / 2. *
                                  (solidBooleanTolSmall + scint_reflector_thick),
                          mufflerFBScint_trapThick2_half + scint_reflector_thick,
                          mufflerFBScint_trapFlatL_half + solidBooleanTolSmall / 2. -
                              scint_reflector_thick / 2.);

            mufflerFBScint_MotherSolid =
                new G4UnionSolid("Muon_FBMuffler_Mother_Solid_Temp", mufflerFBScint_Mother_Box,
                                 mufflerFBScint_MotherTrap, scint_alignMtx_Uside,
                                 {0,
                                  mufflerFBScint_boxHeight_half + mufflerFBScint_trapFlatL_half +
                                      scint_reflector_thick / 2. - solidBooleanTolSmall / 2.,
                                  0});
        }

        resultLV =
            new G4LogicalVolume(mufflerFBScint_MotherSolid, aReflectorMaterial, resultLVName);
        new G4LogicalSkinSurface("Muon_FBMuffler_Scintillator_reflector_opsurf", resultLV,
                                 aSurfaceProperty);

        mufflerFBScint_Box =
            new G4Box("Muon_Muffler_Scintillator_FB_Box", mufflerFBScint_boxWidth_half,
                      mufflerFBScint_boxHeight_half, mufflerFBScint_boxThick_half);
        mufflerFBScint_FlatTrap =
            new G4Trd("Muon_Muffler_Scintillator_FB_Flat_Trd",
                      mufflerFBScint_trapX2_half -
                          mufflerFBScint_TrapSlope_double * mufflerFBScint_trapFlatL_half,
                      mufflerFBScint_trapX2_half, mufflerFBScint_trapThick2_half,
                      mufflerFBScint_trapThick2_half, mufflerFBScint_trapFlatL_half);

        fI_MufflerFBScint_BoxLogical      = new G4LogicalVolume(mufflerFBScint_Box, aScintMaterial,
                                                           "Muon_Muffler_Scintillator_FB_Box_LV");
        fI_MufflerFBScint_FlatTrapLogical = new G4LogicalVolume(
            mufflerFBScint_FlatTrap, aScintMaterial, "Muon_Muffler_Scintillator_FB_Flat_Trd_LV");

        if (isMufflerFBScint_SideReflector_attached) {
            mufflerFBScint_boxSideReflector =
                new G4Box("Muon_Muffler_Scintillator_FB_Box Side_Reflector_Box",
                          mufflerFBScint_boxWidth_half + scint_muffler_sideRefl_thick,
                          mufflerFBScint_boxHeight_half + scint_muffler_sideRefl_thick / 2.,
                          mufflerFBScint_boxThick_half);
            mufflerFBScint_FlatTrapSideRefl =
                new G4Trd("Muon_Muffler_Scintillator_FB_Flat Side_Reflector_Trd",
                          mufflerFBScint_trapX2_half -
                              mufflerFBScint_TrapSlope_double * mufflerFBScint_trapFlatL_half +
                              scint_muffler_sideRefl_thick,
                          mufflerFBScint_trapX2_half + scint_muffler_sideRefl_thick,
                          mufflerFBScint_trapThick2_half, mufflerFBScint_trapThick2_half,
                          mufflerFBScint_trapFlatL_half);

            mufflerFBScint_boxSideReflLV =
                new G4LogicalVolume(mufflerFBScint_boxSideReflector, aSideReflectorMaterial,
                                    "Muon_Muffler_Scintillator_FB_Box_Side_Reflector_LV");
            mufflerFBScint_flatTrapSideReflLV =
                new G4LogicalVolume(mufflerFBScint_FlatTrapSideRefl, aSideReflectorMaterial,
                                    "Muon_Muffler_Scintillator_FB_Flat_Trap_Side_Reflector_LV");

            new G4LogicalSkinSurface("Muon_FBMuffler_Scintillator_BoxSideReflector_opsurf",
                                     mufflerFBScint_boxSideReflLV, aSideReflectorSurfProp);
            new G4LogicalSkinSurface("Muon_FBMuffler_Scintillator_FlatTrapSideReflector_opsurf",
                                     mufflerFBScint_flatTrapSideReflLV, aSideReflectorSurfProp);

            mufflerFBScint_boxSideReflLV->SetVisAttributes(fI_Vikuiti_VisAttr);
            mufflerFBScint_flatTrapSideReflLV->SetVisAttributes(fI_Vikuiti_VisAttr);

            new G4PVPlacement(nullptr, {0, scint_muffler_sideRefl_thick / 2., 0},
                              fI_MufflerFBScint_BoxLogical, "Muon_Muffler_Scintillator_FB_Box_PV",
                              mufflerFBScint_boxSideReflLV, false, 0, OverlapCheck);

            new G4PVPlacement(nullptr, {}, fI_MufflerFBScint_FlatTrapLogical,
                              "Muon_Muffler_Scintillator_FB_FlatTrd_PV",
                              mufflerFBScint_flatTrapSideReflLV, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr, {0, -scint_muffler_sideRefl_thick / 2., 0}, mufflerFBScint_boxSideReflLV,
                "Muon_Muffler_Scintillator_FB_Box_Side_Reflector_PV", resultLV, false, 0, OverlapCheck);
            new G4PVPlacement(scint_alignMtx_Uside,
                              {0, mufflerFBScint_boxHeight_half + mufflerFBScint_trapFlatL_half, 0},
                              mufflerFBScint_flatTrapSideReflLV,
                              "Muon_Muffler_Scintillator_FB_FlatTrd_Side_Reflector_PV", resultLV,
                              false, 0, OverlapCheck);
        } else {
            new G4PVPlacement(nullptr, {}, fI_MufflerFBScint_BoxLogical,
                              "Muon_Muffler_Scintillator_FB_Box_PV", resultLV, false, 0, OverlapCheck);

            new G4PVPlacement(
                scint_alignMtx_Uside,
                G4ThreeVector(0, mufflerFBScint_boxHeight_half + mufflerFBScint_trapFlatL_half, 0),
                fI_MufflerFBScint_FlatTrapLogical, "Muon_Muffler_Scintillator_FB_FlatTrd_PV",
                resultLV, false, 0, OverlapCheck);
        }

        if (isNotAllFlat) {
            mufflerFBScint_PMTTrap =
                new G4Trd("Muon_Muffler_Scintillator_FB_PMT_Trd", mufflerFBScint_trapX1_half,
                          mufflerFBScint_trapX2_half -
                              mufflerFBScint_TrapSlope_double * mufflerFBScint_trapFlatL_half,
                          mufflerFBScint_trapThick1_half, mufflerFBScint_trapThick2_half,
                          mufflerFBScint_trapZ_half - mufflerFBScint_trapFlatL_half);
            fI_MufflerFBScint_PMTTrapLogical = new G4LogicalVolume(
                mufflerFBScint_PMTTrap, aScintMaterial, "Muon_Muffler_Scintillator_FB_PMT_Trd_LV");

            if (isMufflerFBScint_SideReflector_attached) {
                mufflerFBScint_PMTTrapSideRefl =
                    new G4Trd("Muon_Muffler_Scintillator_FB_PMT_Side_Reflector_Trd",
                              mufflerFBScint_trapX1_half + scint_muffler_sideRefl_thick,
                              mufflerFBScint_trapX2_half -
                                  mufflerFBScint_TrapSlope_double * mufflerFBScint_trapFlatL_half +
                                  scint_muffler_sideRefl_thick,
                              mufflerFBScint_trapThick1_half, mufflerFBScint_trapThick2_half,
                              mufflerFBScint_trapZ_half - mufflerFBScint_trapFlatL_half);

                mufflerFBScint_pmtTrapSideReflLV =
                    new G4LogicalVolume(mufflerFBScint_PMTTrapSideRefl, aSideReflectorMaterial,
                                        "Muon_Muffler_Scintillator_FB_PMT_Side_Reflector_LV");
                new G4LogicalSkinSurface("Muon_FBMuffler_Scintillator_PMTTrapSideReflector_opsurf",
                                         mufflerFBScint_pmtTrapSideReflLV, aSideReflectorSurfProp);
                mufflerFBScint_pmtTrapSideReflLV->SetVisAttributes(fI_Vikuiti_VisAttr);

                new G4PVPlacement(nullptr, {}, fI_MufflerFBScint_PMTTrapLogical,
                                  "Muon_Muffler_Scintillator_FB_PMT_Trd_PV",
                                  mufflerFBScint_pmtTrapSideReflLV, false, 0, OverlapCheck);

                new G4PVPlacement(scint_alignMtx_Uside,
                                  {0,
                                   mufflerFBScint_boxHeight_half + mufflerFBScint_trapZ_half +
                                       mufflerFBScint_trapFlatL_half,
                                   0},
                                  mufflerFBScint_pmtTrapSideReflLV,
                                  "Muon_Muffler_Scintillator_FB_PMT_Side_Reflector_PV", resultLV,
                                  false, 0, OverlapCheck);
            } else {
                new G4PVPlacement(scint_alignMtx_Uside,
                                  {0,
                                   mufflerFBScint_boxHeight_half + mufflerFBScint_trapZ_half +
                                       mufflerFBScint_trapFlatL_half,
                                   0},
                                  fI_MufflerFBScint_PMTTrapLogical,
                                  "Muon_Muffler_Scintillator_FB_PMT_Trd_PV", resultLV, false, 0,
                                  OverlapCheck);
            }
        }
    }
    return resultLV;
}

G4LogicalVolume *AmoreDetectorConstruction::Build_I_Muffler_LRScintillator(
    G4Material *aScintMaterial, G4Material *aReflectorMaterial, G4Material *aSideReflectorMaterial,
    G4SurfaceProperty *aSurfaceProperty, G4SurfaceProperty *aSideReflectorSurfProp) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    using namespace AmoreDetectorStaticInfo;
    constexpr const char *resultLVName = "MuonMufflerScintillatorLR_LV";

    G4double mufflerLRScint_trapFlatL_half  = mufflerLRScint_trapZ_half * scint_muffler_FlatL_ratio;
    G4double mufflerLRScint_trapThick2_half = scint_muffler_Thick_half;
    G4double mufflerLRScint_boxThick_half   = scint_muffler_Thick_half;
    G4double mufflerLRScint_TrapSlope_double =
        (mufflerLRScint_trapX2_half - mufflerLRScint_trapX1_half) / mufflerLRScint_trapZ_half;
    G4bool isNotAllFlat = (scint_muffler_FlatL_ratio != 1);
    G4double mufflerLRScint_trapThick1_half =
        (isNotAllFlat) ? mufflerLRScint_trapX1_half : scint_muffler_Thick_half;
    G4bool isMufflerLRScint_SideReflector_attached;

    G4VSolid *mufflerLRScint_MotherSolid;
    G4Box *mufflerLRScint_Mother_Box;
    G4VSolid *mufflerLRScint_MotherTrap;
    G4Trd *mufflerLRScint_MotherPMTTrap;
    G4Box *mufflerLRScint_Box;
    G4Trd *mufflerLRScint_FlatTrap;
    G4Trd *mufflerLRScint_FlatTrapSideRefl;
    G4Trd *mufflerLRScint_PMTTrap;
    G4Trd *mufflerLRScint_PMTTrapSideRefl;
    G4Box *mufflerLRScint_boxSideReflector;

    G4LogicalVolume *mufflerLRScint_boxSideReflLV;
    G4LogicalVolume *mufflerLRScint_flatTrapSideReflLV;
    G4LogicalVolume *mufflerLRScint_pmtTrapSideReflLV;

    if (scint_muffler_sideRefl_thick > scint_reflector_thick) {
        G4Exception(__PRETTY_FUNCTION__, "SCINT_SIDEREFL_01", JustWarning,
                    "The thickness of a side reflector in muon detector at muffler LR side is too "
                    "thick. The reflector will not be attached.");
        isMufflerLRScint_SideReflector_attached = false;
    } else
        isMufflerLRScint_SideReflector_attached = true;

    G4RotationMatrix *scintMother_alignMtx_LRSide = new G4RotationMatrix;
    G4RotationMatrix *scint_alignMtx_Uside        = new G4RotationMatrix;
    resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(resultLVName, false);
    if (resultLV == nullptr) {
        scint_alignMtx_Uside->rotateX(-90 * deg);
        scintMother_alignMtx_LRSide->rotateY(-90 * deg);
        scintMother_alignMtx_LRSide->rotateZ(-90 * deg);

        mufflerLRScint_Mother_Box = new G4Box("Muon_LRMuffler_Scintillator_Mother_Box",
                                              mufflerLRScint_boxWidth_half + scint_reflector_thick,
                                              mufflerLRScint_boxHeight_half + scint_reflector_thick,
                                              //mufflerLRScint_boxHeight_half,
                                              mufflerLRScint_boxThick_half + scint_reflector_thick);
        if (isNotAllFlat) {
            mufflerLRScint_MotherTrap = new G4Trd(
                "Muon_LRMuffler_Scintillator_MS_PreTrap1",
                mufflerLRScint_trapX2_half -
                    mufflerLRScint_TrapSlope_double * // Slope(기울기) 값이 원래 값의 2배
                                                      // (Double)임에 주의할 것
                        (mufflerLRScint_trapFlatL_half +
                         solidBooleanTolSmall / 2.), // 해당 기울기 값 (2배임) 에 맞게 계산된 것
                mufflerLRScint_trapX2_half + mufflerLRScint_TrapSlope_double *
                                                 solidBooleanTolSmall /
                                                 2., // 해당 기울기 값 (2배임) 에 맞게 계산된 것
                mufflerLRScint_trapThick2_half + scint_reflector_thick,
                mufflerLRScint_trapThick2_half + scint_reflector_thick,
                mufflerLRScint_trapFlatL_half + solidBooleanTolSmall);
            mufflerLRScint_MotherPMTTrap =
                new G4Trd("Muon_LRMuffler_Scintillator_MS_PreTrap2",
                          mufflerLRScint_trapX1_half +
                              mufflerLRScint_TrapSlope_double * scint_reflector_thick / 2.,
                          mufflerLRScint_trapX2_half -
                              mufflerLRScint_TrapSlope_double * mufflerLRScint_trapFlatL_half,
                          mufflerLRScint_trapThick1_half + scint_reflector_thick,
                          mufflerLRScint_trapThick2_half + scint_reflector_thick,
                          mufflerLRScint_trapZ_half - mufflerLRScint_trapFlatL_half -
                              scint_reflector_thick / 2.);
            mufflerLRScint_MotherTrap = new G4UnionSolid(
                "Muon_LRMuffler_Scintillator_MS_Trap", mufflerLRScint_MotherTrap,
                mufflerLRScint_MotherPMTTrap, nullptr,
                G4ThreeVector(0, 0, -mufflerLRScint_trapZ_half + scint_reflector_thick / 2.));

            mufflerLRScint_MotherSolid =
                new G4UnionSolid("Muon_LRMuffler_Mother_Solid", mufflerLRScint_Mother_Box,
                                 mufflerLRScint_MotherTrap, scint_alignMtx_Uside,
                                 {0,
                                  mufflerLRScint_boxHeight_half + mufflerLRScint_trapFlatL_half +
                                      scint_reflector_thick,
                                  0});

        } else {
            mufflerLRScint_MotherTrap =
                new G4Trd("Muon_LRMuffler_Scintillator_MS_PreTrap1",
                          mufflerLRScint_trapX1_half +
                              mufflerLRScint_TrapSlope_double * // Slope(기울기) 값이 원래 값의 2배
                                                                // (Double)임에 주의할 것
                                  (scint_reflector_thick + solidBooleanTolSmall) / 2.,
                          mufflerLRScint_trapX2_half +
                              mufflerLRScint_TrapSlope_double * solidBooleanTolSmall / 2.,
                          mufflerLRScint_trapThick1_half +
                              scint_reflector_thick * mufflerLRScint_TrapSlope_double / 2. *
                                  (solidBooleanTolSmall + scint_reflector_thick),
                          mufflerLRScint_trapThick2_half + scint_reflector_thick,
                          mufflerLRScint_trapFlatL_half + solidBooleanTolSmall / 2. -
                              scint_reflector_thick / 2.);

            mufflerLRScint_MotherSolid =
                new G4UnionSolid("Muon_LRMuffler_Mother_Solid_Temp", mufflerLRScint_Mother_Box,
                                 mufflerLRScint_MotherTrap, scint_alignMtx_Uside,
                                 {0,
                                  mufflerLRScint_boxHeight_half + mufflerLRScint_trapFlatL_half +
                                      scint_reflector_thick / 2. - solidBooleanTolSmall / 2.,
                                      //- solidBooleanTolSmall / 2.,
                                  0});

						scint_alignMtx_Uside->rotateX(180 * deg);
            mufflerLRScint_MotherSolid =
                new G4UnionSolid("Muon_LRMuffler_Mother_Solid", mufflerLRScint_MotherSolid,
                                 mufflerLRScint_MotherTrap, scint_alignMtx_Uside,
                                 {0,
                                  -mufflerLRScint_boxHeight_half - mufflerLRScint_trapFlatL_half +
                                      -scint_reflector_thick /2.  +solidBooleanTolSmall /2.,
                                  0});
        }

        resultLV =
            new G4LogicalVolume(mufflerLRScint_MotherSolid, aReflectorMaterial, resultLVName);
        new G4LogicalSkinSurface("Muon_LRMuffler_Scintillator_reflector_opsurf", resultLV,
                                 aSurfaceProperty);

        mufflerLRScint_Box =
            new G4Box("Muon_Muffler_Scintillator_LR_Box", mufflerLRScint_boxWidth_half,
                      mufflerLRScint_boxHeight_half, mufflerLRScint_boxThick_half);
        mufflerLRScint_FlatTrap =
            new G4Trd("Muon_Muffler_Scintillator_LR_Flat_Trd",
                      mufflerLRScint_trapX2_half -
                          mufflerLRScint_TrapSlope_double * mufflerLRScint_trapFlatL_half,
                      mufflerLRScint_trapX2_half, mufflerLRScint_trapThick2_half,
                      mufflerLRScint_trapThick2_half, mufflerLRScint_trapFlatL_half);

        fI_MufflerLRScint_BoxLogical      = new G4LogicalVolume(mufflerLRScint_Box, aScintMaterial,
                                                           "Muon_Muffler_Scintillator_LR_Box_LV");
        fI_MufflerLRScint_FlatTrapLogical = new G4LogicalVolume(
            mufflerLRScint_FlatTrap, aScintMaterial, "Muon_Muffler_Scintillator_LR_Flat_Trd_LV");

        if (isMufflerLRScint_SideReflector_attached) {
            mufflerLRScint_boxSideReflector =
                new G4Box("Muon_Muffler_Scintillator_LR_Box Side_Reflector_Box",
                          mufflerLRScint_boxWidth_half + scint_muffler_sideRefl_thick,
                          mufflerLRScint_boxHeight_half,
                          //mufflerLRScint_boxHeight_half + scint_muffler_sideRefl_thick / 2.,
                          mufflerLRScint_boxThick_half);
            mufflerLRScint_FlatTrapSideRefl =
                new G4Trd("Muon_Muffler_Scintillator_LR_Flat Side_Reflector_Trd",
                          mufflerLRScint_trapX2_half -
                              mufflerLRScint_TrapSlope_double * mufflerLRScint_trapFlatL_half +
                              scint_muffler_sideRefl_thick,
                          mufflerLRScint_trapX2_half + scint_muffler_sideRefl_thick,
                          mufflerLRScint_trapThick2_half, mufflerLRScint_trapThick2_half,
                          mufflerLRScint_trapFlatL_half);

            mufflerLRScint_boxSideReflLV =
                new G4LogicalVolume(mufflerLRScint_boxSideReflector, aSideReflectorMaterial,
                                    "Muon_Muffler_Scintillator_LR_Box_Side_Reflector_LV");
            mufflerLRScint_flatTrapSideReflLV =
                new G4LogicalVolume(mufflerLRScint_FlatTrapSideRefl, aSideReflectorMaterial,
                                    "Muon_Muffler_Scintillator_LR_Flat_Trap_Side_Reflector_LV");

            new G4LogicalSkinSurface("Muon_LRMuffler_Scintillator_BoxSideReflector_opsurf",
                                     mufflerLRScint_boxSideReflLV, aSideReflectorSurfProp);
            new G4LogicalSkinSurface("Muon_LRMuffler_Scintillator_FlatTrapSideReflector_opsurf",
                                     mufflerLRScint_flatTrapSideReflLV, aSideReflectorSurfProp);

            mufflerLRScint_boxSideReflLV->SetVisAttributes(fI_Vikuiti_VisAttr);
            mufflerLRScint_flatTrapSideReflLV->SetVisAttributes(fI_Vikuiti_VisAttr);

            //new G4PVPlacement(nullptr, {0, scint_muffler_sideRefl_thick / 2., 0},
            new G4PVPlacement(nullptr, {0, 0, 0},
                              fI_MufflerLRScint_BoxLogical, "Muon_Muffler_Scintillator_LR_Box_PV",
                              mufflerLRScint_boxSideReflLV, false, 0, OverlapCheck);

            new G4PVPlacement(nullptr, {}, fI_MufflerLRScint_FlatTrapLogical,
                              "Muon_Muffler_Scintillator_LR_FlatTrd_PV",
                              mufflerLRScint_flatTrapSideReflLV, false, 0, OverlapCheck);

            new G4PVPlacement(
                //nullptr, {0, -scint_muffler_sideRefl_thick / 2., 0}, mufflerLRScint_boxSideReflLV,
                nullptr, {0,0, 0}, mufflerLRScint_boxSideReflLV,
                "Muon_Muffler_Scintillator_LR_Box_Side_Reflector_PV", resultLV, false, 0, OverlapCheck);

            new G4PVPlacement(G4Transform3D(*scint_alignMtx_Uside,
                              {0, mufflerLRScint_boxHeight_half + mufflerLRScint_trapFlatL_half, 0}),
                              mufflerLRScint_flatTrapSideReflLV,
                              "Muon_Muffler_Scintillator_LR_FlatTrd_Side_Reflector1_PV", resultLV,
                              false, 0, OverlapCheck);

						scint_alignMtx_Uside->rotateX(180 * deg);
            new G4PVPlacement(G4Transform3D(*scint_alignMtx_Uside,
                              {0, - mufflerLRScint_boxHeight_half - mufflerLRScint_trapFlatL_half, 0}),
                              mufflerLRScint_flatTrapSideReflLV,
                              "Muon_Muffler_Scintillator_LR_FlatTrd_Side_Reflector2_PV", resultLV,
                              false, 0, OverlapCheck);
        } else {
            new G4PVPlacement(nullptr, {}, fI_MufflerLRScint_BoxLogical,
                              "Muon_Muffler_Scintillator_LR_Box_PV", resultLV, false, 0, OverlapCheck);

            new G4PVPlacement(
                scint_alignMtx_Uside,
                G4ThreeVector(0, mufflerLRScint_boxHeight_half + mufflerLRScint_trapFlatL_half, 0),
                fI_MufflerLRScint_FlatTrapLogical, "Muon_Muffler_Scintillator_LR_FlatTrd_PV",
                resultLV, false, 0, OverlapCheck);
        }

        if (isNotAllFlat) {
            mufflerLRScint_PMTTrap =
                new G4Trd("Muon_Muffler_Scintillator_LR_PMT_Trd", mufflerLRScint_trapX1_half,
                          mufflerLRScint_trapX2_half -
                              mufflerLRScint_TrapSlope_double * mufflerLRScint_trapFlatL_half,
                          mufflerLRScint_trapThick1_half, mufflerLRScint_trapThick2_half,
                          mufflerLRScint_trapZ_half - mufflerLRScint_trapFlatL_half);
            fI_MufflerLRScint_PMTTrapLogical = new G4LogicalVolume(
                mufflerLRScint_PMTTrap, aScintMaterial, "Muon_Muffler_Scintillator_LR_PMT_Trd_LV");

            if (isMufflerLRScint_SideReflector_attached) {
                mufflerLRScint_PMTTrapSideRefl =
                    new G4Trd("Muon_Muffler_Scintillator_LR_PMT_Side_Reflector_Trd",
                              mufflerLRScint_trapX1_half + scint_muffler_sideRefl_thick,
                              mufflerLRScint_trapX2_half -
                                  mufflerLRScint_TrapSlope_double * mufflerLRScint_trapFlatL_half +
                                  scint_muffler_sideRefl_thick,
                              mufflerLRScint_trapThick1_half, mufflerLRScint_trapThick2_half,
                              mufflerLRScint_trapZ_half - mufflerLRScint_trapFlatL_half);

                mufflerLRScint_pmtTrapSideReflLV =
                    new G4LogicalVolume(mufflerLRScint_PMTTrapSideRefl, aSideReflectorMaterial,
                                        "Muon_Muffler_Scintillator_LR_PMT_Side_Reflector_LV");
                new G4LogicalSkinSurface("Muon_LRMuffler_Scintillator_PMTTrapSideReflector_opsurf",
                                         mufflerLRScint_pmtTrapSideReflLV, aSideReflectorSurfProp);
                mufflerLRScint_pmtTrapSideReflLV->SetVisAttributes(fI_Vikuiti_VisAttr);

                new G4PVPlacement(nullptr, {}, fI_MufflerLRScint_PMTTrapLogical,
                                  "Muon_Muffler_Scintillator_LR_PMT_Trd_PV",
                                  mufflerLRScint_pmtTrapSideReflLV, false, 0, OverlapCheck);

                new G4PVPlacement(scint_alignMtx_Uside,
                                  {0,
                                   mufflerLRScint_boxHeight_half + mufflerLRScint_trapZ_half +
                                       mufflerLRScint_trapFlatL_half,
                                   0},
                                  mufflerLRScint_pmtTrapSideReflLV,
                                  "Muon_Muffler_Scintillator_LR_PMT_Side_Reflector_PV", resultLV,
                                  false, 0, OverlapCheck);
            } else {
                new G4PVPlacement(G4Transform3D(*scint_alignMtx_Uside,
                                  {0,
                                   -mufflerLRScint_boxHeight_half - mufflerLRScint_trapZ_half -
                                       mufflerLRScint_trapFlatL_half,
                                   0}),
                                  fI_MufflerLRScint_PMTTrapLogical,
                                  "Muon_Muffler_Scintillator_LR_PMT_Trd_PV", resultLV, false, 0,
                                  OverlapCheck);
								scint_alignMtx_Uside->rotateX(180 * deg);
                new G4PVPlacement(G4Transform3D(*scint_alignMtx_Uside,
                                  {0,
                                   mufflerLRScint_boxHeight_half + mufflerLRScint_trapZ_half +
                                       mufflerLRScint_trapFlatL_half,
                                   0}),
                                  fI_MufflerLRScint_PMTTrapLogical,
                                  "Muon_Muffler_Scintillator_LR_PMT_Trd_PV", resultLV, false, 0,
                                  OverlapCheck);

            }
        }
    }
    return resultLV;
}

void AmoreDetectorConstruction::Place_I_GeometriesInRooftopForTowerOf(
    G4LogicalVolume *aWhere, G4int aType, G4ThreeVector aPosOfRooftopFloor) {

    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4Tubs *towerBoltM3S;

    G4Box *detModuleTopPinConnectorBox;

    G4Box *adhesiveFT2850ABox;
    G4Box *adhesiveFT2850BBox;
    G4Box *adhesiveFT2850CBox;

    G4Box *photonDetectorPCBEpoxyBox;
    G4Box *photonDetectorPCBSolder1Box;
    G4Box *photonDetectorPCBSolder2Box;
    G4ExtrudedSolid *photonDetectorPCBSolid;

    G4LogicalVolume *logicalVolumeChecker;
    G4LogicalVolume *towerBoltM3SLV;

    G4LogicalVolume *detModuleTopPinConnectorLV;

    G4LogicalVolume *adhesiveFT2850ALV;
    G4LogicalVolume *adhesiveFT2850BLV;
    G4LogicalVolume *adhesiveFT2850CLV;

    G4LogicalVolume *photonDetectorPCBLV;
    G4LogicalVolume *photonDetectorPCBEpoxyLV;
    G4LogicalVolume *photonDetectorPCBSolder1LV;
    G4LogicalVolume *photonDetectorPCBSolder2LV;

    G4ThreeVector nowTowerBoltPos;
    G4ThreeVector nowTowerBoltDisplace;
    G4ThreeVector nowSolidPosition;
    G4ThreeVector nowPCBPosition;

    G4TwoVector nowVertexDisplacer;
    std::vector<G4TwoVector> photonDetectorPCBVertices;
    std::vector<G4TwoVector> photonDetectorPCBInnerVertices;

    G4RotationMatrix *nowAdhesiveFT2850RotMtx;
    G4RotationMatrix *photonDetectorEpoxyRotMtx;

    switch (aType) {
        case 1: {
            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("BoltM3S_LV", false);
            if (logicalVolumeChecker == nullptr) {
                towerBoltM3S   = new G4Tubs("BoltM3S_Tub", 0, detModuleBoltM3SRadius,
                                          detModuleBoltM3SLength / 2., 0, 360 * deg);
                towerBoltM3SLV = new G4LogicalVolume(towerBoltM3S, _stainless, "BoltM3S_LV");
                towerBoltM3SLV->SetVisAttributes(fI_Bolt_VisAttr);
            } else {
                towerBoltM3SLV = G4LogicalVolumeStore::GetInstance()->GetVolume("BoltM3S_LV");
            }

            nowTowerBoltPos      = G4ThreeVector(0, 0, detModuleBoltM3SLength / 2.);
            nowTowerBoltDisplace = G4ThreeVector(
                -copperFrameSizeX / 2. + detModuleBoltM3SRadius + dMBM3SDistanceXFromXEnd,
                copperFrameSizeY / 2. - copperFrameHoleDistance, 0);

            nowTowerBoltPos += aPosOfRooftopFloor;
            nowTowerBoltPos += nowTowerBoltDisplace;
            new G4PVPlacement(nullptr, nowTowerBoltPos, towerBoltM3SLV, "BoltM3S_AtTower_PV",
                              aWhere, false, 0, OverlapCheck);
            nowTowerBoltPos -= nowTowerBoltDisplace;
            nowTowerBoltDisplace[0] *= -1.;
            nowTowerBoltPos += nowTowerBoltDisplace;
            new G4PVPlacement(nullptr, nowTowerBoltPos, towerBoltM3SLV, "BoltM3S_AtTower_PV",
                              aWhere, false, 1, OverlapCheck);
            nowTowerBoltPos -= nowTowerBoltDisplace;
            nowTowerBoltDisplace[1] *= -1.;
            nowTowerBoltDisplace[0] *= -1.;
            nowTowerBoltPos += nowTowerBoltDisplace;
            new G4PVPlacement(nullptr, nowTowerBoltPos, towerBoltM3SLV, "BoltM3S_AtTower_PV",
                              aWhere, false, 2, OverlapCheck);
            nowTowerBoltPos -= nowTowerBoltDisplace;
            nowTowerBoltDisplace[0] *= -1.;
            nowTowerBoltPos += nowTowerBoltDisplace;
            new G4PVPlacement(nullptr, nowTowerBoltPos, towerBoltM3SLV, "BoltM3S_AtTower_PV",
                              aWhere, false, 3, OverlapCheck);

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("AdhesiveFT2850_A_Box", false);
            if (logicalVolumeChecker == nullptr) {
                adhesiveFT2850ABox =
                    new G4Box("AdhesiveFT2850_A_Box", detModuleFT2850ASizeX / 2.,
                              detModuleFT2850ASizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850BBox =
                    new G4Box("AdhesiveFT2850_B_Box", detModuleFT2850BSizeX / 2.,
                              detModuleFT2850BSizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850CBox =
                    new G4Box("AdhesiveFT2850_C_Box", detModuleFT2850CSizeX / 2.,
                              detModuleFT2850CSizeY / 2., detModuleFT2850Thickness / 2.);

                adhesiveFT2850ALV =
                    new G4LogicalVolume(adhesiveFT2850ABox, _stycast, "AdhesiveFT2850_A_LV");
                adhesiveFT2850BLV =
                    new G4LogicalVolume(adhesiveFT2850BBox, _stycast, "AdhesiveFT2850_B_LV");
                adhesiveFT2850CLV =
                    new G4LogicalVolume(adhesiveFT2850CBox, _stycast, "AdhesiveFT2850_C_LV");

                adhesiveFT2850ALV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850BLV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850CLV->SetVisAttributes(fI_FT2850_VisAttr);

                detModuleTopPinConnectorBox =
                    new G4Box("PinConnectorTop_Box", detModuleTopPinConnectorSizeX / 2.,
                              detModuleTopPinConnectorSizeY / 2., detModulePinConnectorSizeZ / 2.);

                detModuleTopPinConnectorLV =
                    new G4LogicalVolume(detModuleTopPinConnectorBox, _kevlar, "PinConnectorTop_LV");
            } else {
                G4LogicalVolumeStore *nowInstance = G4LogicalVolumeStore::GetInstance();

                adhesiveFT2850ALV = nowInstance->GetVolume("AdhesiveFT2850_A_LV", false);
                adhesiveFT2850BLV = nowInstance->GetVolume("AdhesiveFT2850_B_LV", false);
                adhesiveFT2850CLV = nowInstance->GetVolume("AdhesiveFT2850_C_LV", false);

                detModuleTopPinConnectorLV = nowInstance->GetVolume("PinConnectorTop_LV");
            }

            // Adhesive parts for optical detector cabling
            nowSolidPosition =
                G4ThreeVector(0, dMFT2850AYDistanceFromCenter, detModuleFT2850Thickness / 2.);
            nowSolidPosition += aPosOfRooftopFloor;
            new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850ALV, "AdhesiveFT2850_A_PV",
                              aWhere, false, 0, OverlapCheck);

            nowSolidPosition =
                G4ThreeVector(0, dMFT2850BYDistanceFromCenter, detModuleFT2850Thickness / 2.);
            nowSolidPosition.rotateZ(120. * deg);
            nowSolidPosition += aPosOfRooftopFloor;
            nowAdhesiveFT2850RotMtx = new G4RotationMatrix();
            nowAdhesiveFT2850RotMtx->rotateZ(120. / 2. * deg);
            nowAdhesiveFT2850RotMtx->inverse();
            new G4PVPlacement(nowAdhesiveFT2850RotMtx, nowSolidPosition, adhesiveFT2850BLV,
                              "AdhesiveFT2850_B_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {0, copperFrameSizeY / 2. - detModuleFT2850CSizeY / 2.,
                                detModuleFT2850Thickness / 2.};
            nowSolidPosition += aPosOfRooftopFloor;
            new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850CLV, "AdhesiveFT2850_C_PV",
                              aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {0, copperFrameSizeY / 2. + detModuleTopPinConnectorSizeY / 2.,
                                detModulePinConnectorSizeZ / 2.};
            nowSolidPosition += aPosOfRooftopFloor;
            new G4PVPlacement(nullptr, nowSolidPosition, detModuleTopPinConnectorLV,
                              "PinConnectorTop_PV", aWhere, false, 0, OverlapCheck);
        } break;
        case 2: {
            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV", false);
            if (logicalVolumeChecker == nullptr) {
                photonDetectorPCBVertices.resize(7);
                photonDetectorPCBVertices[0] = {photDetPCBBaseBoxSizeX / 2.,
                                                -photDetPCBBaseBoxSizeY / 2.};
                photonDetectorPCBVertices[1] = {photDetPCBBaseBoxSizeX / 2.,
                                                photDetPCBBaseBoxSizeY / 2.};
                photonDetectorPCBVertices[5] = {-photDetPCBBaseBoxSizeX / 2.,
                                                photDetPCBBaseBoxSizeY / 2.};
                photonDetectorPCBVertices[6] = {-photDetPCBBaseBoxSizeX / 2.,
                                                -photDetPCBBaseBoxSizeY / 2.};

                nowVertexDisplacer = {0, photDetPCBTriangle1Length};
                nowVertexDisplacer.rotate(45. * deg);
                photonDetectorPCBVertices[2] = photonDetectorPCBVertices[1] + nowVertexDisplacer;

                nowVertexDisplacer = {0, photDetPCBTriangle2Length};
                nowVertexDisplacer.rotate(-45. * deg);
                photonDetectorPCBVertices[4] = photonDetectorPCBVertices[5] + nowVertexDisplacer;

                nowVertexDisplacer = {0, -photDetPCBTriangle0Length + photDetPCBTriangle2Length};
                nowVertexDisplacer.rotate(-45. * deg);
                photonDetectorPCBVertices[3] = photonDetectorPCBVertices[2] + nowVertexDisplacer;

                photonDetectorPCBSolid =
                    new G4ExtrudedSolid("PCBForPhotonDetectorSolid", photonDetectorPCBVertices,
                                        photDetPCBThick / 2., 0, 1., 0, 1.);
                photonDetectorPCBEpoxyBox =
                    new G4Box("PhotonDetectorPCBEpoxy_Box", photDetPCBEpoxyBoxSizeX / 2.,
                              photDetPCBEpoxyBoxSizeY / 2., photDetPCBEpoxyBoxSizeZ / 2.);
                photonDetectorPCBSolder1Box =
                    new G4Box("PhotonDetectorPCBSolder1_Box", photDetPCBSolderBox1SizeX / 2.,
                              photDetPCBSolderBox1SizeY / 2., photDetPCBThick / 2);
                photonDetectorPCBSolder2Box =
                    new G4Box("PhotonDetectorPCBSolder2_Box", photDetPCBSolderBox2SizeX / 2.,
                              photDetPCBSolderBox2SizeY / 2., photDetPCBThick / 2);

                photonDetectorPCBLV =
                    new G4LogicalVolume(photonDetectorPCBSolid, _copper, "PCBForPhotonDetector_LV");
                photonDetectorPCBEpoxyLV   = new G4LogicalVolume(photonDetectorPCBEpoxyBox, _solder,
                                                               "PhotonDetectorPCBEpoxy_LV");
                photonDetectorPCBSolder1LV = new G4LogicalVolume(
                    photonDetectorPCBSolder1Box, _solder, "PhotonDetectorPCBSolder1_LV");
                photonDetectorPCBSolder2LV = new G4LogicalVolume(
                    photonDetectorPCBSolder2Box, _solder, "PhotonDetectorPCBSolder2_LV");

                photonDetectorPCBLV->SetVisAttributes(fI_PCB_VisAttr);
                photonDetectorPCBEpoxyLV->SetVisAttributes(fI_LeadDefault_VisAttr);
                photonDetectorPCBSolder1LV->SetVisAttributes(fI_LeadDefault_VisAttr);
                photonDetectorPCBSolder2LV->SetVisAttributes(fI_LeadDefault_VisAttr);
            } else {
                photonDetectorPCBLV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV");
                photonDetectorPCBSolid =
                    static_cast<G4ExtrudedSolid *>(photonDetectorPCBLV->GetSolid());
                photonDetectorPCBVertices = photonDetectorPCBSolid->GetPolygon();
                photonDetectorPCBEpoxyLV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PhotonDetectorPCBEpoxy_LV");
                photonDetectorPCBSolder1LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PhotonDetectorPCBSolder1_LV");
                photonDetectorPCBSolder2LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PhotonDetectorPCBSolder2_LV");
            }

            photonDetectorEpoxyRotMtx = new G4RotationMatrix();
            photonDetectorEpoxyRotMtx->rotateZ(45 * deg);

            nowPCBPosition = {photDetPCBOriginDistX, photDetPCBOriginDistY,
                              photDetPCBEpoxyBoxSizeZ + photDetPCBThick / 2.};
            nowPCBPosition += aPosOfRooftopFloor;
            new G4PVPlacement(nullptr, nowPCBPosition, photonDetectorPCBLV,
                              "PCBForPhotonDetector_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {photonDetectorPCBVertices[2].x(), photDetPCBEpoxyBox1Dist,
                                -photDetPCBThick / 2. - photDetPCBEpoxyBoxSizeZ / 2.};
            new G4PVPlacement(photonDetectorEpoxyRotMtx, nowPCBPosition + nowSolidPosition,
                              photonDetectorPCBEpoxyLV, "PhotonDetectorPCBEpoxy_PV", aWhere, false,
                              0, OverlapCheck);

            nowSolidPosition = {photonDetectorPCBVertices[4].x(), photDetPCBEpoxyBox2Dist,
                                -photDetPCBThick / 2. - photDetPCBEpoxyBoxSizeZ / 2.};
            new G4PVPlacement(photonDetectorEpoxyRotMtx, nowPCBPosition + nowSolidPosition,
                              photonDetectorPCBEpoxyLV, "PhotonDetectorPCBEpoxy_PV", aWhere, false,
                              1, OverlapCheck);

            nowSolidPosition = {photonDetectorPCBVertices[2].x() + photDetPCBSolderBox1SizeX / 2.,
                                photDetPCBSolderBox1SizeY / 2., photDetPCBThick};
            new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition,
                              photonDetectorPCBSolder1LV, "PhotonDetectorPCBSolder1_PV", aWhere,
                              false, 0, OverlapCheck);

            nowSolidPosition = {
                (photonDetectorPCBVertices[3].x() + photonDetectorPCBVertices[4].x()) / 2., 0.,
                photDetPCBThick};
            new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition,
                              photonDetectorPCBSolder2LV, "PhotonDetectorPCBSolder2_PV", aWhere,
                              false, 0, OverlapCheck);
        } break;
        case 3: 
				{
					// Crystal PCB Bolts head -----------------------
					G4Tubs *crystalPCBBoltsHeadTubs = new G4Tubs("PCBBoltsHeadForCrystalTubs", 0,
							photDetPCBBoltsHoleRadius * 1.45, photDetPCBBoltsHeadHeight / 2., 0, 360. * deg);
					G4LogicalVolume *crystalPCBBoltsHeadLV = new G4LogicalVolume(
							crystalPCBBoltsHeadTubs, _brass, "PCBBoltsHeadForTop_LV");
					crystalPCBBoltsHeadLV->SetVisAttributes(fI_BrassBolt_VisAttr);
					nowSolidPosition = {type3_photDetPCBBaseBoxSizeX / 2. - photDetPCBBoltsHoleDistX,
															type3_photDetPCBBaseBoxSizeY / 2. - photDetPCBBoltsHoleDistY,
															+ photDetPCBBoltsHeadHeight /2.};
					nowSolidPosition += aPosOfRooftopFloor;

					new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
							"PCBBoltsForPhoton_PV", aWhere, false, 0, OverlapCheck);
					nowSolidPosition[1] += - type3_photDetPCBBaseBoxSizeY + photDetPCBBoltsHoleDistY*2;
					new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
							"PCBBoltsForPhoton_PV", aWhere, false, 1, OverlapCheck);
					nowSolidPosition[0] += -type3_photDetPCBBaseBoxSizeX + photDetPCBBoltsHoleDistX*2;
					new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
							"PCBBoltsForPhoton_PV", aWhere, false, 2, OverlapCheck);

					// photon-phonon bolts positioning
					G4Tubs *photonPhononBoltsHead2 = new G4Tubs("photonPhononBoltsHead2Tubs",0,
							type3_photonPhononBoltsHeadRadius, 
							(type3_photonPhononBoltsHeadThick - type3_photonFrameTopSpace) / 2., 0, 360. * deg);
					G4LogicalVolume *photonPhononBoltsHeadLV = new G4LogicalVolume(photonPhononBoltsHead2, 
							_brass, "photonPhononBoltsHeadForTop_LV");
					photonPhononBoltsHeadLV->SetVisAttributes(fI_BrassBolt_VisAttr);
					
					nowSolidPosition = {-type3_photonFrameSizeX / 2. + type3_photonFrameBaseWidthX / 2.,
						type3_photonFrameSizeY / 2. - type3_photonFrameBaseWidthY / 2.,
						(type3_photonPhononBoltsHeadThick - type3_photonFrameTopSpace) / 2. };
					nowSolidPosition += aPosOfRooftopFloor;

					new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsHeadLV,
							"photonPhononBolts_PV", aWhere, false, 0, OverlapCheck);
					nowSolidPosition[0] += type3_photonFrameSizeX - type3_photonFrameBaseWidthX;
					new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsHeadLV,
							"photonPhononBolts_PV", aWhere, false, 1, OverlapCheck);
					nowSolidPosition[1] += -type3_photonFrameSizeY + type3_photonFrameBaseWidthY;
					new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsHeadLV,
							"photonPhononBolts_PV", aWhere, false, 2, OverlapCheck);
					nowSolidPosition[0] -= type3_photonFrameSizeX - type3_photonFrameBaseWidthX;
					new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsHeadLV,
							"photonPhononBolts_PV", aWhere, false, 3, OverlapCheck);
				} break;
    }
}

// Modified by Jeewon at 2021. Mar
// case1 : Bolts
// case2 : PCBs
// case3 : Real geometry phonon PCBs
// case4 ; Real geometry bolts for phonon and photon detector
void AmoreDetectorConstruction::Place_I_GeometriesInDetModuleOf(G4LogicalVolume *aWhere,
                                                                G4int aType,
                                                                const CrystalModuleInfo &aCrystal,
                                                                G4bool aPlacePartsForJustBelow) {
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4LogicalVolume *logicalVolumeChecker;
    G4ThreeVector nowSolidPosition;

    G4ExtrudedSolid *crystalPCBSolid;
    G4LogicalVolume *crystalPCBLV;
    G4ThreeVector nowPCBPosition;
    G4TwoVector nowVertexDisplacer;
    std::vector<G4TwoVector> crystalPCBVertices;

		// For case 1
    G4Tubs *boltM4;
    G4Tubs *boltM3;
    G4Tubs *boltM3F;
    G4Tubs *boltAB;

    G4Box *detModuleTopPinConnectorBox;
    G4Box *detModuleB1PinConnectorBox;
    G4Box *detModuleB2PinConnectorBox;

    G4Box *adhesiveFT2850ABox;
    G4Box *adhesiveFT2850BBox;
    G4Box *adhesiveFT2850CBox;
    G4Box *adhesiveFT2850AForCrystalBox;
    G4Box *adhesiveFT2850BForCrystalBox;
    G4Box *adhesiveFT2850CForCrystalBox;

    G4LogicalVolume *boltM4LV;
    G4LogicalVolume *boltM3LV;
    G4LogicalVolume *boltM3FLV;
    G4LogicalVolume *boltABLV;

    G4LogicalVolume *detModuleTopPinConnectorLV;
    G4LogicalVolume *detModuleB1PinConnectorLV;
    G4LogicalVolume *detModuleB2PinConnectorLV;

    G4LogicalVolume *adhesiveFT2850ALV;
    G4LogicalVolume *adhesiveFT2850BLV;
    G4LogicalVolume *adhesiveFT2850CLV;
    G4LogicalVolume *adhesiveFT2850AForCrystalLV;
    G4LogicalVolume *adhesiveFT2850BForCrystalLV;
    G4LogicalVolume *adhesiveFT2850CForCrystalLV;

    G4double boltM3FYStartCoord = copperFrameSizeY / 2. + smallBlock1DistY - smallBlock1SizeY / 2.;
    G4double boltM3FYSpacing =
        (fabs(boltM3FYStartCoord) * 2. - 3. * 2. * detModuleBoltM3FRadius) / 4.;

    G4RotationMatrix *nowAdhesiveFT2850RotMtx;

		// For case 4
		G4LogicalVolume *crystalPCBBoltsBodyLV;
		G4LogicalVolume *crystalPCBBoltsHeadLV;
		G4LogicalVolume *crystalClampLV;
		G4LogicalVolume *photonPhononBoltsLV;
		G4LogicalVolume *clampBoltsLV;
			
		// For case 2
    G4Box *photonDetectorPCBEpoxyBox;
    G4Box *photonDetectorPCBSolder1Box;
    G4Box *photonDetectorPCBSolder2Box;
		G4ExtrudedSolid *photonDetectorPCBSolid;

    G4Box *crystalHeaterBox;
    G4Box *crystalHeaterElectrodeBox;
    G4Tubs *crystalPEEKTub;

    G4Box *crystalPCBSolder1Box;
    G4Box *crystalPCBSolder2Box;
    G4Box *crystalPCBEpoxy1Box;
    G4Box *crystalPCBEpoxy2Box;

    G4LogicalVolume *photonDetectorPCBLV;
    G4LogicalVolume *photonDetectorPCBEpoxyLV;
    G4LogicalVolume *photonDetectorPCBSolder1LV;
    G4LogicalVolume *photonDetectorPCBSolder2LV;

    G4LogicalVolume *crystalHeaterLV;
    G4LogicalVolume *crystalHeaterElectrodeLV;
    G4LogicalVolume *crystalPEEKLV;

    G4LogicalVolume *crystalPCBSolder1LV;
    G4LogicalVolume *crystalPCBSolder2LV;
    G4LogicalVolume *crystalPCBEpoxy1LV;
    G4LogicalVolume *crystalPCBEpoxy2LV;

    std::vector<G4TwoVector> photonDetectorPCBVertices;
    G4RotationMatrix *photonDetectorEpoxyRotMtx;

		// for type3
		G4double type3_moduleSizeZ = aCrystal.fCrystalHeight + type3_smallBlock1SizeZ + type3_smallBlock3SizeZ ;
		G4double clampSizeX = type3_copperFrameSizeX / 2. - type3_smallBlockSizeX - type3_photonFrameSizeX / 2.;
		G4ThreeVector moduleTopPosition = {0, 0, type3_moduleSizeZ/2.};

    switch (aType) {
        case 1: {
            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("BoltM3_LV", false);
            if (logicalVolumeChecker == nullptr) {
                boltM4  = new G4Tubs("BoltM4_Tub", 0, detModuleBoltM4Radius,
                                    detModuleBoltM4Length / 2., 0, 360 * deg);
                boltM3  = new G4Tubs("BoltM3_Tub", 0, detModuleBoltM3Radius,
                                    detModuleBoltM3Length / 2., 0, 360 * deg);
                boltM3F = new G4Tubs("BoltM3F_Tub", 0, detModuleBoltM3FRadius,
                                     detModuleBoltM3FLength / 2., 0, 360 * deg);
                boltAB  = new G4Tubs("BoltAB_Tub", 0, detModuleBoltABRadius,
                                    detModuleBoltABLength / 2., 0, 360 * deg);

                adhesiveFT2850ABox =
                    new G4Box("AdhesiveFT2850_A_Box", detModuleFT2850ASizeX / 2.,
                              detModuleFT2850ASizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850BBox =
                    new G4Box("AdhesiveFT2850_B_Box", detModuleFT2850BSizeX / 2.,
                              detModuleFT2850BSizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850CBox =
                    new G4Box("AdhesiveFT2850_C_Box", detModuleFT2850CSizeX / 2.,
                              detModuleFT2850CSizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850AForCrystalBox = new G4Box(
                    "AdhesiveFT2850_A_ForCrystal_Box", detModuleFT2850AForCrystalSizeX / 2.,
                    detModuleFT2850AForCrystalSizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850BForCrystalBox = new G4Box(
                    "AdhesiveFT2850_B_ForCrystal_Box", detModuleFT2850BForCrystalSizeX / 2.,
                    detModuleFT2850BForCrystalSizeY / 2., detModuleFT2850Thickness / 2.);
                adhesiveFT2850CForCrystalBox = new G4Box(
                    "AdhesiveFT2850_C_ForCrystal_Box", detModuleFT2850CForCrystalSizeX / 2.,
                    detModuleFT2850CForCrystalSizeY / 2., detModuleFT2850Thickness / 2.);

                boltM4LV  = new G4LogicalVolume(boltM4, _stainless, "BoltM4_LV");
                boltM3LV  = new G4LogicalVolume(boltM3, _stainless, "BoltM3_LV");
                boltM3FLV = new G4LogicalVolume(boltM3F, _stainless, "BoltM3F_LV");
                boltABLV  = new G4LogicalVolume(boltAB, _stainless, "BoltAB_LV");

                adhesiveFT2850ALV =
                    new G4LogicalVolume(adhesiveFT2850ABox, _stycast, "AdhesiveFT2850_A_LV");
                adhesiveFT2850BLV =
                    new G4LogicalVolume(adhesiveFT2850BBox, _stycast, "AdhesiveFT2850_B_LV");
                adhesiveFT2850CLV =
                    new G4LogicalVolume(adhesiveFT2850CBox, _stycast, "AdhesiveFT2850_C_LV");
                adhesiveFT2850AForCrystalLV = new G4LogicalVolume(
                    adhesiveFT2850AForCrystalBox, _stycast, "AdhesiveFT2850_A_ForCrystal_LV");
                adhesiveFT2850BForCrystalLV = new G4LogicalVolume(
                    adhesiveFT2850BForCrystalBox, _stycast, "AdhesiveFT2850_B_ForCrystal_LV");
                adhesiveFT2850CForCrystalLV = new G4LogicalVolume(
                    adhesiveFT2850CForCrystalBox, _stycast, "AdhesiveFT2850_C_ForCrystal_LV");

                boltM4LV->SetVisAttributes(fI_Bolt_VisAttr);
                boltM3LV->SetVisAttributes(fI_Bolt_VisAttr);
                boltM3FLV->SetVisAttributes(fI_Bolt_VisAttr);
                boltABLV->SetVisAttributes(fI_Bolt_VisAttr);

                adhesiveFT2850ALV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850BLV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850CLV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850AForCrystalLV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850BForCrystalLV->SetVisAttributes(fI_FT2850_VisAttr);
                adhesiveFT2850CForCrystalLV->SetVisAttributes(fI_FT2850_VisAttr);

                detModuleTopPinConnectorBox =
                    new G4Box("PinConnectorTop_Box", detModuleTopPinConnectorSizeX / 2.,
                              detModuleTopPinConnectorSizeY / 2., detModulePinConnectorSizeZ / 2.);
                detModuleB1PinConnectorBox =
                    new G4Box("PinConnectorB1_Box", detModuleB1PinConnectorSizeX / 2.,
                              detModuleB1PinConnectorSizeY / 2, detModulePinConnectorSizeZ / 2.);
                detModuleB2PinConnectorBox =
                    new G4Box("PinConnectorB2_Box", detModuleB2PinConnectorSizeX / 2.,
                              detModuleB2PinConnectorSizeY / 2., detModulePinConnectorSizeZ / 2.);

                detModuleTopPinConnectorLV =
                    new G4LogicalVolume(detModuleTopPinConnectorBox, _kevlar, "PinConnectorTop_LV");
                detModuleB1PinConnectorLV =
                    new G4LogicalVolume(detModuleB1PinConnectorBox, _kevlar, "PinConnectorB1_LV");
                detModuleB2PinConnectorLV =
                    new G4LogicalVolume(detModuleB2PinConnectorBox, _kevlar, "PinConnectorB2_LV");

            } else {
                G4LogicalVolumeStore *nowInstance = G4LogicalVolumeStore::GetInstance();

                boltM4LV  = nowInstance->GetVolume("BoltM4_LV", false);
                boltM3LV  = nowInstance->GetVolume("BoltM3_LV", false);
                boltM3FLV = nowInstance->GetVolume("BoltM3F_LV", false);
                boltABLV  = nowInstance->GetVolume("BoltAB_LV", false);

                adhesiveFT2850ALV = nowInstance->GetVolume("AdhesiveFT2850_A_LV", false);
                adhesiveFT2850BLV = nowInstance->GetVolume("AdhesiveFT2850_B_LV", false);
                adhesiveFT2850CLV = nowInstance->GetVolume("AdhesiveFT2850_C_LV", false);

                adhesiveFT2850AForCrystalLV =
                    nowInstance->GetVolume("AdhesiveFT2850_A_ForCrystal_LV", false);
                adhesiveFT2850BForCrystalLV =
                    nowInstance->GetVolume("AdhesiveFT2850_B_ForCrystal_LV", false);
                adhesiveFT2850CForCrystalLV =
                    nowInstance->GetVolume("AdhesiveFT2850_C_ForCrystal_LV", false);

                detModuleTopPinConnectorLV = nowInstance->GetVolume("PinConnectorTop_LV");
                detModuleB1PinConnectorLV  = nowInstance->GetVolume("PinConnectorB1_LV");
                detModuleB2PinConnectorLV  = nowInstance->GetVolume("PinConnectorB2_LV");
            }

            nowSolidPosition = G4ThreeVector(dMFT2850AForCrystalXDistanceFromCenter, 0,
                                             -aCrystal.fCrystalHeight / 2. - bottomSupportASizeZ -
                                                 detModuleFT2850Thickness / 2.);

            new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850AForCrystalLV,
                              "AdhesiveFT2850_A_ForCrystal_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = G4ThreeVector(dMFT2850BForCrystalXDistanceFromCenter,
                                             dMFT2850BForCrystalYDistanceFromCenter,
                                             -aCrystal.fCrystalHeight / 2. - bottomSupportASizeZ -
                                                 detModuleFT2850Thickness / 2.);

            new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850BForCrystalLV,
                              "AdhesiveFT2850_B_ForCrystal_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = G4ThreeVector(dMFT2850CForCrystalXDistanceFromCenter,
                                             dMFT2850CForCrystalYDistanceFromCenter,
                                             -aCrystal.fCrystalHeight / 2. - bottomSupportASizeZ -
                                                 detModuleFT2850Thickness / 2.);

            new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850CForCrystalLV,
                              "AdhesiveFT2850_C_ForCrystal_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {dMFT2850BForCrystalXDistanceFromCenter,
                                copperFrameSizeY / 2. + detModuleB1PinConnectorSizeY / 2.,
                                -aCrystal.fCrystalHeight / 2. - bottomSupportASizeZ -
                                    detModulePinConnectorSizeZ / 2.};
            new G4PVPlacement(nullptr, nowSolidPosition, detModuleB1PinConnectorLV,
                              "PinConnectorB1_LV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {dMFT2850CForCrystalXDistanceFromCenter,
                                -copperFrameSizeY / 2. - detModuleB2PinConnectorSizeY / 2.,
                                -aCrystal.fCrystalHeight / 2. - bottomSupportASizeZ -
                                    detModulePinConnectorSizeZ / 2.};

            new G4PVPlacement(nullptr, nowSolidPosition, detModuleB2PinConnectorLV,
                              "PinConnectorB2_LV", aWhere, false, 0, OverlapCheck);

            // SSCB BOLTS: M4-12
            nowSolidPosition =
                G4ThreeVector(-copperFrameSizeX / 2. + copperFrameHoleDistance,
                              -copperFrameSizeY / 2. + copperFrameHoleDistance,
                              aCrystal.fCrystalHeight / 2. + detModuleBoltM4Length / 2.);

            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                              0, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                              1, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            nowSolidPosition[1] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                              2, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                              3, OverlapCheck);
            nowSolidPosition *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                              0, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                              1, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            nowSolidPosition[1] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                              2, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                              3, OverlapCheck);

            // SSCB BOLTS on crystal support A in Bottom CopperFrame: M3-4
            nowSolidPosition = G4ThreeVector(
                bottomSupportASizeX + bottomSupportADistX - detModuleBoltM3FRadius -
                    detMBM3FAtSupportADistance,
                -copperFrameSizeY / 2. + copperFrameHoleDistance,
                -aCrystal.fCrystalHeight / 2. - detModuleBoltM3FLength / 2. - bottomSupportASizeZ);
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FBottom_AtSupportA_PV",
                              aWhere, false, 0, OverlapCheck);
            nowSolidPosition[1] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FBottom_AtSupportA_PV",
                              aWhere, false, 1, OverlapCheck);

            // SSCB BOLTS on crystal support B in Bottom CopperFrame: M3-4
            nowSolidPosition = G4ThreeVector(
                -bottomSupportBSizeX - bottomSupportBDistX + detModuleBoltM3FRadius +
                    detMBM3FAtSupportBDistance,
                -copperFrameSizeY / 2. + copperFrameHoleDistance,
                -aCrystal.fCrystalHeight / 2. - detModuleBoltM3FLength / 2. - bottomSupportBSizeZ);
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FBottom_AtSupportB_PV",
                              aWhere, false, 0, OverlapCheck);
            nowSolidPosition[1] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FBottom_AtSupportB_PV",
                              aWhere, false, 1, OverlapCheck);

            // SSCB BOLTS on Upper CopperFrame (along Y axis): M3-4
            nowSolidPosition =
                G4ThreeVector(-copperFrameSizeX / 2. - smallBlock1DistX,
                              boltM3FYStartCoord - boltM3FYSpacing - detModuleBoltM3FRadius,
                              aCrystal.fCrystalHeight / 2. + detModuleBoltM3FLength / 2.);

            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FUpper_AlongYAxis_PV",
                              aWhere, false, 0, OverlapCheck);
            nowSolidPosition[1] -= boltM3FYSpacing + 2. * detModuleBoltM3FRadius;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FUpper_AlongYAxis_PV",
                              aWhere, false, 1, OverlapCheck);
            nowSolidPosition[1] -= boltM3FYSpacing + 2. * detModuleBoltM3FRadius;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FUpper_AlongYAxis_PV",
                              aWhere, false, 2, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            nowSolidPosition[1] = boltM3FYStartCoord - boltM3FYSpacing - detModuleBoltM3FRadius;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FUpper_AlongYAxis_PV",
                              aWhere, false, 3, OverlapCheck);
            nowSolidPosition[1] -= boltM3FYSpacing + 2. * detModuleBoltM3FRadius;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FUpper_AlongYAxis_PV",
                              aWhere, false, 4, OverlapCheck);
            nowSolidPosition[1] -= boltM3FYSpacing + 2. * detModuleBoltM3FRadius;
            new G4PVPlacement(nullptr, nowSolidPosition, boltM3FLV, "BoltM3FUpper_AlongYAxis_PV",
                              aWhere, false, 5, OverlapCheck);

            // additional Bolts : AB
            // beside of pilar
            nowSolidPosition = G4ThreeVector(
                -copperFrameSizeX / 2. + copperFrameHoleDistance,
                copperFrameSizeY / 2. - dMBABDistanceXFromXEnd - detModuleBoltABRadius,
                -aCrystal.fCrystalHeight / 2. + copperFrameSizeZ + detModuleBoltABLength / 2.);

            new G4PVPlacement(nullptr, nowSolidPosition, boltABLV, "BoltAB_PV", aWhere, false, 0, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltABLV, "BoltAB_PV", aWhere, false, 1, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            nowSolidPosition[1] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltABLV, "BoltAB_PV", aWhere, false, 2, OverlapCheck);
            nowSolidPosition[0] *= -1.;
            new G4PVPlacement(nullptr, nowSolidPosition, boltABLV, "BoltAB_PV", aWhere, false, 3, OverlapCheck);

            // SSCB BOLTS on Upper CopperFrame to fix photon detector of just below module: M3-6
            if (aPlacePartsForJustBelow) {
                nowSolidPosition = G4ThreeVector(
                    -copperFrameSizeX / 2. + dMBM3DistanceXFromXEnd + detModuleBoltM3Radius,
                    copperFrameSizeY / 2. - copperFrameHoleDistance,
                    -aCrystal.fCrystalHeight / 2. - smallBlock1SizeZ + detModuleBoltM3Length / 2.);
                new G4PVPlacement(nullptr, nowSolidPosition, boltM3LV, "BoltM3_PV", aWhere, false,
                                  0, OverlapCheck);
                nowSolidPosition[0] *= -1.;
                new G4PVPlacement(nullptr, nowSolidPosition, boltM3LV, "BoltM3_PV", aWhere, false,
                                  1, OverlapCheck);
                nowSolidPosition[0] *= -1.;
                nowSolidPosition[1] *= -1.;
                new G4PVPlacement(nullptr, nowSolidPosition, boltM3LV, "BoltM3_PV", aWhere, false,
                                  2, OverlapCheck);
                nowSolidPosition[0] *= -1.;
                new G4PVPlacement(nullptr, nowSolidPosition, boltM3LV, "BoltM3_PV", aWhere, false,
                                  3, OverlapCheck);

                // Adhesive parts for optical detector cabling
                nowSolidPosition = G4ThreeVector(0, dMFT2850AYDistanceFromCenter,
                                                 -aCrystal.fCrystalHeight / 2. - smallBlock1SizeZ +
                                                     detModuleFT2850Thickness / 2.);

                new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850ALV,
                                  "AdhesiveFT2850_A_PV", aWhere, false, 0, OverlapCheck);

                nowSolidPosition = G4ThreeVector(0, dMFT2850BYDistanceFromCenter,
                                                 -aCrystal.fCrystalHeight / 2. - smallBlock1SizeZ +
                                                     detModuleFT2850Thickness / 2.);
                nowSolidPosition.rotateZ(120. * deg);
                nowAdhesiveFT2850RotMtx = new G4RotationMatrix();
                nowAdhesiveFT2850RotMtx->rotateZ(120. / 2. * deg);
                nowAdhesiveFT2850RotMtx->inverse();

                new G4PVPlacement(nowAdhesiveFT2850RotMtx, nowSolidPosition, adhesiveFT2850BLV,
                                  "AdhesiveFT2850_B_PV", aWhere, false, 0, OverlapCheck);

                nowSolidPosition =
                    G4ThreeVector(0, copperFrameSizeY / 2. - detModuleFT2850CSizeY / 2.,
                                  -aCrystal.fCrystalHeight / 2. - smallBlock1SizeZ +
                                      detModuleFT2850Thickness / 2.);

                new G4PVPlacement(nullptr, nowSolidPosition, adhesiveFT2850CLV,
                                  "AdhesiveFT2850_C_PV", aWhere, false, 0, OverlapCheck);

                nowSolidPosition = {0, copperFrameSizeY / 2. + detModuleTopPinConnectorSizeY / 2.,
                                    -aCrystal.fCrystalHeight / 2. - smallBlock1SizeZ +
                                        detModulePinConnectorSizeZ / 2.};
                new G4PVPlacement(nullptr, nowSolidPosition, detModuleTopPinConnectorLV,
                                  "PinConnectorTop_PV", aWhere, false, 0, OverlapCheck);
            }
            break;
        }
        case 2: {
            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV", false);
            if (logicalVolumeChecker == nullptr) {
                photonDetectorPCBVertices.resize(7);
                photonDetectorPCBVertices[0] = {photDetPCBBaseBoxSizeX / 2.,
                                                -photDetPCBBaseBoxSizeY / 2.};
                photonDetectorPCBVertices[1] = {photDetPCBBaseBoxSizeX / 2.,
                                                photDetPCBBaseBoxSizeY / 2.};
                photonDetectorPCBVertices[5] = {-photDetPCBBaseBoxSizeX / 2.,
                                                photDetPCBBaseBoxSizeY / 2.};
                photonDetectorPCBVertices[6] = {-photDetPCBBaseBoxSizeX / 2.,
                                                -photDetPCBBaseBoxSizeY / 2.};

                nowVertexDisplacer = {0, photDetPCBTriangle1Length};
                nowVertexDisplacer.rotate(45. * deg);
                photonDetectorPCBVertices[2] = photonDetectorPCBVertices[1] + nowVertexDisplacer;

                nowVertexDisplacer = {0, photDetPCBTriangle2Length};
                nowVertexDisplacer.rotate(-45. * deg);
                photonDetectorPCBVertices[4] = photonDetectorPCBVertices[5] + nowVertexDisplacer;

                nowVertexDisplacer = {0, -photDetPCBTriangle0Length + photDetPCBTriangle2Length};
                nowVertexDisplacer.rotate(-45. * deg);
                photonDetectorPCBVertices[3] = photonDetectorPCBVertices[2] + nowVertexDisplacer;

                photonDetectorPCBSolid =
                    new G4ExtrudedSolid("PCBForPhotonDetectorSolid", photonDetectorPCBVertices,
                                        photDetPCBThick / 2., 0, 1., 0, 1.);
                photonDetectorPCBEpoxyBox =
                    new G4Box("PhotonDetectorPCBEpoxy_Box", photDetPCBEpoxyBoxSizeX / 2.,
                              photDetPCBEpoxyBoxSizeY / 2., photDetPCBEpoxyBoxSizeZ / 2.);
                photonDetectorPCBSolder1Box =
                    new G4Box("PhotonDetectorPCBSolder1_Box", photDetPCBSolderBox1SizeX / 2.,
                              photDetPCBSolderBox1SizeY / 2., photDetPCBThick / 2);
                photonDetectorPCBSolder2Box =
                    new G4Box("PhotonDetectorPCBSolder2_Box", photDetPCBSolderBox2SizeX / 2.,
                              photDetPCBSolderBox2SizeY / 2., photDetPCBThick / 2);

                photonDetectorPCBLV =
                    new G4LogicalVolume(photonDetectorPCBSolid, _copper, "PCBForPhotonDetector_LV");
                photonDetectorPCBEpoxyLV   = new G4LogicalVolume(photonDetectorPCBEpoxyBox, _solder,
                                                               "PhotonDetectorPCBEpoxy_LV");
                photonDetectorPCBSolder1LV = new G4LogicalVolume(
                    photonDetectorPCBSolder1Box, _solder, "PhotonDetectorPCBSolder1_LV");
                photonDetectorPCBSolder2LV = new G4LogicalVolume(
                    photonDetectorPCBSolder2Box, _solder, "PhotonDetectorPCBSolder2_LV");

                photonDetectorPCBLV->SetVisAttributes(fI_PCB_VisAttr);
                photonDetectorPCBEpoxyLV->SetVisAttributes(fI_LeadDefault_VisAttr);
                photonDetectorPCBSolder1LV->SetVisAttributes(fI_LeadDefault_VisAttr);
                photonDetectorPCBSolder2LV->SetVisAttributes(fI_LeadDefault_VisAttr);

                crystalHeaterBox =
                    new G4Box("HeaterForCrystalBox", crystalHeaterBoxSizeX / 2.,
                              crystalHeaterBoxSizeY / 2., crystalHeaterBoxSizeZ / 2.);
                crystalHeaterElectrodeBox =
                    new G4Box("HeaterForCrystalElectrodeBox", crystalHeaterElectrodeSizeX / 2.,
                              crystalHeaterElectrodeSizeY / 2., crystalHeaterElectrodeSizeZ / 2.);

                crystalHeaterLV =
                    new G4LogicalVolume(crystalHeaterBox, _SiWafer, "HeaterForCrystal_LV");
                crystalHeaterElectrodeLV = new G4LogicalVolume(crystalHeaterElectrodeBox, _AuPd,
                                                               "HeaterForCrystal_Electrode_LV");

                nowSolidPosition = {crystalHeaterElectrodeStartX, 0,
                                    -crystalHeaterBoxSizeZ / 2. + crystalHeaterElectrodeSizeZ / 2.};
                for (G4int i = 0; i < 12; i++) {
                    new G4PVPlacement(nullptr, nowSolidPosition, crystalHeaterElectrodeLV,
                                      "HeaterForCrystal_Electrode_PV", crystalHeaterLV, false, i, OverlapCheck);
                    nowSolidPosition[0] += crystalHeaterElectrodeSpacing;
                }

                crystalHeaterLV->SetVisAttributes(fI_GeWafer_VisAttr);
                crystalHeaterElectrodeLV->SetVisAttributes(fI_GoldDefault_VisAttr);

                crystalPEEKTub = new G4Tubs("CrystalPEEK_Tub", 0, crystalPeekRadius,
                                            crystalPeekThick / 2., 0, 360. * deg);
                crystalPEEKLV  = new G4LogicalVolume(crystalPEEKTub, _peek1, "CrystalPEEK_LV");
                crystalPEEKLV->SetVisAttributes(fI_PEEK_VisAttr);

                crystalPCBVertices.resize(10);
                crystalPCBVertices[0] = {crystalPCBWidth / 2,
                                         -crystalPCBMiddleBoxLength / 2. + crystalPCBDownBoxLength};
                crystalPCBVertices[1] = {crystalPCBWidth / 2, crystalPCBMiddleBoxLength / 2.};
                crystalPCBVertices[6] = {-crystalPCBWidth / 2, crystalPCBMiddleBoxLength / 2.};
                crystalPCBVertices[7] = {-crystalPCBWidth / 2, -crystalPCBMiddleBoxLength / 2.};
                crystalPCBVertices[8] = {-crystalPCBWidth / 2 + crystalPCBDownBoxWidth,
                                         -crystalPCBMiddleBoxLength / 2.};
                crystalPCBVertices[9] = {-crystalPCBWidth / 2 + crystalPCBDownBoxWidth,
                                         -crystalPCBMiddleBoxLength / 2. + crystalPCBDownBoxLength};

                crystalPCBVertices[2] = {crystalPCBWidth / 2., -crystalPCBUpBoxLength / 2.};
                crystalPCBVertices[3] = {crystalPCBWidth / 2., crystalPCBUpBoxLength / 2.};
                crystalPCBVertices[4] = {-crystalPCBWidth / 2., crystalPCBUpBoxLength / 2.};
                crystalPCBVertices[5] = {-crystalPCBWidth / 2., -crystalPCBUpBoxLength / 2.};

                nowVertexDisplacer = {crystalPCBUpBoxXDistance, crystalPCBUpBoxYDistance};

                crystalPCBVertices[2] += nowVertexDisplacer;
                crystalPCBVertices[3] += nowVertexDisplacer;
                crystalPCBVertices[4] += nowVertexDisplacer;
                crystalPCBVertices[5] += nowVertexDisplacer;

                crystalPCBSolid = new G4ExtrudedSolid("PCBForCrystal_Solid", crystalPCBVertices,
                                                      crystalPCBThick / 2., 0, 1., 0, 1.);

                crystalPCBLV = new G4LogicalVolume(crystalPCBSolid, _copper, "PCBForCrystal_LV");
                crystalPCBLV->SetVisAttributes(fI_PCB_VisAttr);

                crystalPCBEpoxy1Box =
                    new G4Box("CrystalPCBEpoxy1_Box", crystalPCBEpoxyBox1SizeXY / 2.,
                              crystalPCBEpoxyBox1SizeXY / 2., crystalPCBEpoxyBox12SizeZ / 2.);
                crystalPCBEpoxy2Box =
                    new G4Box("CrystalPCBEpoxy2_Box", crystalPCBEpoxyBox2SizeX / 2.,
                              crystalPCBEpoxyBox2SizeY / 2., crystalPCBEpoxyBox12SizeZ / 2.);

                crystalPCBSolder1Box =
                    new G4Box("CrystalPCBSolder1_Box", crystalPCBSolderBox1SizeX / 2.,
                              crystalPCBSolderBox1SizeY / 2., crystalPCBThick / 2.);
                crystalPCBSolder2Box =
                    new G4Box("CrystalPCBSolder2_Box", crystalPCBSolderBox2SizeX / 2.,
                              crystalPCBSolderBox2SizeY / 2., crystalPCBThick / 2.);

                crystalPCBEpoxy1LV =
                    new G4LogicalVolume(crystalPCBEpoxy1Box, _solder, "CrystalPCBEpoxy1_LV");
                crystalPCBEpoxy2LV =
                    new G4LogicalVolume(crystalPCBEpoxy2Box, _solder, "CrystalPCBEpoxy2_LV");
                crystalPCBSolder1LV =
                    new G4LogicalVolume(crystalPCBSolder1Box, _solder, "CrystalPCBSolder1_LV");
                crystalPCBSolder2LV =
                    new G4LogicalVolume(crystalPCBSolder2Box, _solder, "CrystalPCBSolder2_LV");

                crystalPCBSolder1LV->SetVisAttributes(fI_LeadDefault_VisAttr);
                crystalPCBSolder2LV->SetVisAttributes(fI_LeadDefault_VisAttr);
                crystalPCBEpoxy1LV->SetVisAttributes(fI_Epoxy_VisAttr);
                crystalPCBEpoxy2LV->SetVisAttributes(fI_Epoxy_VisAttr);

            } else {
                photonDetectorPCBLV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForPhotonDetector_LV");
                photonDetectorPCBSolid =
                    static_cast<G4ExtrudedSolid *>(photonDetectorPCBLV->GetSolid());
                photonDetectorPCBVertices = photonDetectorPCBSolid->GetPolygon();
                photonDetectorPCBEpoxyLV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PhotonDetectorPCBEpoxy_LV");
                photonDetectorPCBSolder1LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PhotonDetectorPCBSolder1_LV");
                photonDetectorPCBSolder2LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("PhotonDetectorPCBSolder2_LV");
                crystalHeaterLV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("HeaterForCrystal_LV");

                crystalPEEKLV = G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalPEEK_LV");
                crystalPCBLV  = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForCrystal_LV");
                crystalPCBSolid    = static_cast<G4ExtrudedSolid *>(crystalPCBLV->GetSolid());
                crystalPCBVertices = crystalPCBSolid->GetPolygon();

                crystalPCBEpoxy1LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalPCBEpoxy1_LV");
                crystalPCBEpoxy2LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalPCBEpoxy2_LV");
                crystalPCBSolder1LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalPCBSolder1_LV");
                crystalPCBSolder2LV =
                    G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalPCBSolder2_LV");
            }

            new G4PVPlacement(nullptr,
                              {crystalHeaterXDist, crystalHeaterYDist,
                               -aCrystal.fCrystalHeight / 2. - crystalHeaterBoxSizeZ / 2.},
                              crystalHeaterLV, "HeaterForCrystal_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {crystalPeekDist2, -crystalPeekDist1,
                                -aCrystal.fCrystalHeight / 2. - crystalPeekThick / 2.};
            new G4PVPlacement(nullptr, nowSolidPosition, crystalPEEKLV, "CrystalPEEK_PV", aWhere,
                              false, 0, OverlapCheck);
            nowSolidPosition = {-crystalPeekDist1, 0,
                                -aCrystal.fCrystalHeight / 2. - crystalPeekThick / 2.};
            new G4PVPlacement(nullptr, nowSolidPosition, crystalPEEKLV, "CrystalPEEK_PV", aWhere,
                              false, 1, OverlapCheck);
            nowSolidPosition = {crystalPeekDist1, crystalPeekDist2,
                                -aCrystal.fCrystalHeight / 2. - crystalPeekThick / 2.};
            new G4PVPlacement(nullptr, nowSolidPosition, crystalPEEKLV, "CrystalPEEK_PV", aWhere,
                              false, 2, OverlapCheck);

            nowPCBPosition = {crystalPCBOriginDistX, crystalPCBOriginDistY,
                              -aCrystal.fCrystalHeight / 2. - bottomSupportASizeZ -
                                  crystalPCBEpoxyBox12SizeZ - crystalPCBThick / 2.};

            new G4PVPlacement(nullptr, nowPCBPosition, crystalPCBLV, "PCBForCrystal_PV", aWhere,
                              false, 0, OverlapCheck);

            nowSolidPosition = {crystalPCBEpoxyBox1XDist, crystalPCBEpoxyBox1YDist,
                                crystalPCBThick / 2. + crystalPCBEpoxyBox12SizeZ / 2.};
            new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition, crystalPCBEpoxy1LV,
                              "CrystalPCBEpoxy1_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {crystalPCBEpoxyBox2XDist, crystalPCBEpoxyBox2YDist,
                                crystalPCBThick / 2. + crystalPCBEpoxyBox12SizeZ / 2.};
            new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition, crystalPCBEpoxy2LV,
                              "CrystalPCBEpoxy2_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {crystalPCBSolderBox1XDist, crystalPCBSolderBox1YDist,
                                -crystalPCBThick};
            new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition, crystalPCBSolder1LV,
                              "CrystalPCBSolder1_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {crystalPCBSolderBox2XDist, crystalPCBSolderBox2YDist,
                                -crystalPCBThick};
            new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition, crystalPCBSolder2LV,
                              "CrystalPCBSolder2_PV", aWhere, false, 0, OverlapCheck);

            if (aPlacePartsForJustBelow) {
                photonDetectorEpoxyRotMtx = new G4RotationMatrix();
                photonDetectorEpoxyRotMtx->rotateZ(45 * deg);

                nowPCBPosition = {photDetPCBOriginDistX, photDetPCBOriginDistY,
                                  -aCrystal.fCrystalHeight / 2. - smallBlock1SizeZ +
                                      photDetPCBEpoxyBoxSizeZ + photDetPCBThick / 2.};
                new G4PVPlacement(nullptr, nowPCBPosition, photonDetectorPCBLV,
                                  "PCBForPhotonDetector_PV", aWhere, false, 0, OverlapCheck);

                nowSolidPosition = {
                    photonDetectorPCBVertices[2].x(),
                    (photonDetectorPCBVertices[1].y() + photonDetectorPCBVertices[2].y()) / 2.,
                    -photDetPCBThick / 2. - photDetPCBEpoxyBoxSizeZ / 2.};
                new G4PVPlacement(photonDetectorEpoxyRotMtx, nowPCBPosition + nowSolidPosition,
                                  photonDetectorPCBEpoxyLV, "PhotonDetectorPCBEpoxy_PV", aWhere,
                                  false, 0, OverlapCheck);

                nowSolidPosition = {
                    photonDetectorPCBVertices[4].x(),
                    (photonDetectorPCBVertices[5].y() + photonDetectorPCBVertices[4].y()) / 2.,
                    -photDetPCBThick / 2. - photDetPCBEpoxyBoxSizeZ / 2.};
                new G4PVPlacement(photonDetectorEpoxyRotMtx, nowPCBPosition + nowSolidPosition,
                                  photonDetectorPCBEpoxyLV, "PhotonDetectorPCBEpoxy_PV", aWhere,
                                  false, 1, OverlapCheck);

                nowSolidPosition = {photonDetectorPCBVertices[2].x() +
                                        photDetPCBSolderBox1SizeX / 2.,
                                    photDetPCBSolderBox1SizeY / 2., photDetPCBThick};
                new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition,
                                  photonDetectorPCBSolder1LV, "PhotonDetectorPCBSolder1_PV", aWhere,
                                  false, 0, OverlapCheck);

                nowSolidPosition = {
                    (photonDetectorPCBVertices[3].x() + photonDetectorPCBVertices[4].x()) / 2., 0.,
                    photDetPCBThick};
                new G4PVPlacement(nullptr, nowPCBPosition + nowSolidPosition,
                                  photonDetectorPCBSolder2LV, "PhotonDetectorPCBSolder2_PV", aWhere,
                                  false, 0, OverlapCheck);
            }break;
        }
				case 3: {// real configuration : phonon PCB
            logicalVolumeChecker = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForCrystal_LV", false);
            if (logicalVolumeChecker == nullptr) 
						{
							// PCB For crystal ---------------------
              crystalPCBVertices.resize(11);
              crystalPCBVertices[0] = {type3_crystalPCBWidth / 2., type3_crystalPCBLength / 2.};
							crystalPCBVertices[1] = {-type3_crystalPCBWidth / 2., type3_crystalPCBLength / 2.};
							crystalPCBVertices[4] = {-type3_crystalPCBWidth / 2. + type3_crystalPCBGrooveX,
									                    -type3_crystalPCBLength / 2.};

              nowVertexDisplacer = {0,-type3_crystalPCBShortY};
							crystalPCBVertices[10] = crystalPCBVertices[0] + nowVertexDisplacer;
							crystalPCBVertices[2]  = crystalPCBVertices[1] + nowVertexDisplacer; 

              nowVertexDisplacer = {type3_crystalPCBGrooveX,type3_crystalPCBGrooveY};
							crystalPCBVertices[9]  = crystalPCBVertices[10] - nowVertexDisplacer; 

							crystalPCBVertices[3] = {-type3_crystalPCBWidth / 2. + type3_crystalPCBGrooveX,
									                    -type3_crystalPCBLength / 2.+type3_crystalPCBBigY1};
							crystalPCBVertices[5] = {-type3_crystalPCBWidth / 2.+type3_crystalPCBGrooveX + type3_crystalPCBDownBoxX2,
									                    -type3_crystalPCBLength / 2.};
							crystalPCBVertices[6] = {-type3_crystalPCBWidth / 2.+type3_crystalPCBGrooveX + type3_crystalPCBDownBoxX2,
									                    -type3_crystalPCBLength / 2. + type3_crystalPCBDownBoxY};
							crystalPCBVertices[7] = {-type3_crystalPCBWidth / 2.+type3_crystalPCBGrooveX 
									                      +type3_crystalPCBDownBoxX1 + type3_crystalPCBDownBoxX2,
									                    -type3_crystalPCBLength / 2. + type3_crystalPCBDownBoxY};
							crystalPCBVertices[8] = {-type3_crystalPCBWidth / 2.+type3_crystalPCBGrooveX 
									                       +type3_crystalPCBDownBoxX1 + type3_crystalPCBDownBoxX2,
									                      -type3_crystalPCBLength / 2. + type3_crystalPCBDownBoxY + type3_crystalPCBBigY2};
              crystalPCBSolid = new G4ExtrudedSolid("PCBForCrystal_Solid", crystalPCBVertices,
                                                      crystalPCBThick / 2., 0, 1., 0, 1.);
              crystalPCBLV = new G4LogicalVolume(crystalPCBSolid, _solder, "PCBForCrystal_LV");
              crystalPCBLV->SetVisAttributes(fI_PCB_VisAttr);

							// Heater -----------------------------
              crystalHeaterBox = new G4Box("HeaterForCrystalBox", 
									crystalHeaterBoxSizeX / 2., crystalHeaterBoxSizeY / 2., crystalHeaterBoxSizeZ / 2.);
              crystalHeaterElectrodeBox = new G4Box("HeaterForCrystalElectrodeBox", 
									crystalHeaterElectrodeSizeX / 2., crystalHeaterElectrodeSizeY / 2., crystalHeaterElectrodeSizeZ / 2.);
              crystalHeaterLV = new G4LogicalVolume(crystalHeaterBox, _SiWafer, "HeaterForCrystal_LV");
              crystalHeaterElectrodeLV = new G4LogicalVolume( 
									crystalHeaterElectrodeBox, _AuPd, "HeaterForCrystal_Electrode_LV");
              crystalHeaterLV->SetVisAttributes(fI_GeWafer_VisAttr);
              crystalHeaterElectrodeLV->SetVisAttributes(fI_GoldDefault_VisAttr);

              nowSolidPosition = {crystalHeaterElectrodeStartX, 0,
								-crystalHeaterBoxSizeZ / 2. + crystalHeaterElectrodeSizeZ / 2.};
              for (G4int i = 0; i < 12; i++) 
							{
								new G4PVPlacement(nullptr, nowSolidPosition, crystalHeaterElectrodeLV,
										"HeaterForCrystal_Electrode_PV", crystalHeaterLV, false, i, OverlapCheck);
                nowSolidPosition[0] += crystalHeaterElectrodeSpacing;
              }

							// PEEK spacer ------------------------
							crystalPEEKTub = new G4Tubs("CrystalPEEK_Tub", 0, 
									crystalPeekRadius, crystalPeekThick / 2., 0, 360. * deg);
              crystalPEEKLV  = new G4LogicalVolume(crystalPEEKTub, _peek1, "CrystalPEEK_LV");
              crystalPEEKLV->SetVisAttributes(fI_PEEK_VisAttr);

            } else {
              crystalPCBLV  = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBForCrystal_LV");
              crystalHeaterLV = G4LogicalVolumeStore::GetInstance()->GetVolume("HeaterForCrystal_LV");
              crystalPEEKLV = G4LogicalVolumeStore::GetInstance()->GetVolume("CrystalPEEK_LV");
            }

						// Physical volume placement ///////
						// PCB positioning
						nowPCBPosition = {-type3_copperFrameSizeX / 2. + type3_smallBlockSizeX + type3_bottomSupportBDistX,
							type3_copperFrameSizeY / 2. - type3_crystalPCBLength / 2., 
							type3_smallBlock3SizeZ - type3_bottomSupportBSizeZ - crystalPCBThick / 2.};
						nowPCBPosition -= moduleTopPosition;

            new G4PVPlacement(nullptr, nowPCBPosition, crystalPCBLV, "PCBForCrystal_PV", aWhere,
                              false, 0, OverlapCheck);

						// Heater positioning 
						nowSolidPosition = { crystalHeaterXDist, crystalHeaterYDist,
							type3_smallBlock3SizeZ - crystalHeaterBoxSizeZ / 2.};
						nowSolidPosition -= moduleTopPosition;
            new G4PVPlacement(nullptr, nowSolidPosition,
								//{crystalHeaterXDist, crystalHeaterYDist, -aCrystal.fCrystalHeight / 2. - crystalHeaterBoxSizeZ / 2.},
								crystalHeaterLV, "HeaterForCrystal_PV", aWhere, false, 0, OverlapCheck);

						// PEEK positioning
            nowSolidPosition = {crystalPeekDist2, -crystalPeekDist1,
							type3_smallBlock3SizeZ - crystalPeekThick / 2.};
						nowSolidPosition -= moduleTopPosition;
            new G4PVPlacement(nullptr, nowSolidPosition, 
								crystalPEEKLV, "CrystalPEEK_PV", aWhere, false, 0, OverlapCheck);

            nowSolidPosition = {-crystalPeekDist1, 0,
							type3_smallBlock3SizeZ - crystalPeekThick / 2.};
						nowSolidPosition -= moduleTopPosition;
                                //-aCrystal.fCrystalHeight / 2. - crystalPeekThick / 2.};
            new G4PVPlacement(nullptr, nowSolidPosition, 
								crystalPEEKLV, "CrystalPEEK_PV", aWhere, false, 1, OverlapCheck);

            nowSolidPosition = {crystalPeekDist1, crystalPeekDist2,
							type3_smallBlock3SizeZ - crystalPeekThick / 2.};
						nowSolidPosition -= moduleTopPosition;
                                //-aCrystal.fCrystalHeight / 2. - crystalPeekThick / 2.};
            new G4PVPlacement(nullptr, nowSolidPosition, 
								crystalPEEKLV, "CrystalPEEK_PV", aWhere, false, 2, OverlapCheck);

        break;
			}
			case 4:{
        logicalVolumeChecker = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBBoltsForCrystal_LV", false);
        if (logicalVolumeChecker == nullptr) 
				{
					// Crystal PCB Bolts (3ea) ----------------------
					G4Tubs *crystalPCBBoltsBodyTubs = new G4Tubs("PCBBoltsBodyForCrystalTubs", 0,
							photDetPCBBoltsHoleRadius - solidBooleanTol, photDetPCBBoltsHoleHeight / 2, 0, 360. * deg);
					crystalPCBBoltsBodyLV = new G4LogicalVolume(
							crystalPCBBoltsBodyTubs, _brass, "PCBBoltsForCrystal_LV");
					crystalPCBBoltsBodyLV->SetVisAttributes(fI_BrassBolt_VisAttr);
					
					// Crystal PCB Bolts Head ----------------------------
					G4Tubs *crystalPCBBoltsHeadTubs = new G4Tubs("PCBBoltsHeadForCrystalTubs", 0,
							photDetPCBBoltsHoleRadius * 1.45, photDetPCBBoltsHeadHeight / 2., 0, 360. * deg);
					crystalPCBBoltsHeadLV = new G4LogicalVolume(
							crystalPCBBoltsHeadTubs, _brass, "PCBBoltsHeadForCrystal_LV");
					crystalPCBBoltsHeadLV->SetVisAttributes(fI_BrassBolt_VisAttr);

					// Clamp (4ea) ----------------------------------
					G4Tubs *crystalClampTubs = new G4Tubs("ClampTubs",0,
							type3_photonFrameClampRadius - solidBooleanTol/2., (type3_photonFrameClampThick-solidBooleanTol) /2., 90 * deg, 180. *deg);
					G4Box *crystalClampBox = new G4Box("ClampBox",
							clampSizeX-solidBooleanTol/2, type3_photonFrameClampRadius-solidBooleanTol/2, 
							(type3_photonFrameClampSizeZ-solidBooleanTol) / 2.);
					G4UnionSolid *crystalClampSolid = new G4UnionSolid("ClampSolid", 
							crystalClampTubs, crystalClampBox, nullptr, {clampSizeX/2., 0, 
							type3_photonFrameClampThick / 2. - type3_photonFrameClampSizeZ / 2.});

					crystalClampLV = new G4LogicalVolume( crystalClampSolid, _teflon, "Clamp_LV");
					crystalClampLV->SetVisAttributes(fI_Clamp_VisAttr);

					// Clamp bolts (4ea) -----------------------------
					G4Tubs *clampBoltsHeadTubs = new G4Tubs("clampBoltsHeadTubs", 0,
							type3_photonPhononBoltsHeadRadius - solidBooleanTol, 
							(type3_photonPhononBoltsHeadThick -solidBooleanTol)/ 2., 0, 360. * deg);

					clampBoltsLV = new G4LogicalVolume(clampBoltsHeadTubs, _brass, "clampBolts_LV");
					clampBoltsLV->SetVisAttributes(fI_BrassBolt_VisAttr);

					// photon-phonon connection Bolts (4ea) -----------
					G4Tubs *photonPhononBoltsBody = new G4Tubs("photonPhononBoltsTubs", 0,
							type3_photonPhononBoltsRadius, 
							(type3_photonFrameSizeZ)/ 2., 0, 360. * deg);
					G4Tubs *photonPhononBoltsHead1 = new G4Tubs("photonPhononBoltsHead1Tubs",0,
							type3_photonPhononBoltsHeadRadius, 
							type3_photonFrameTopSpace / 2., 0, 360. * deg);
					G4UnionSolid *photonPhononBoltsSolid = new G4UnionSolid("photonPhononBoltsSolid",
							photonPhononBoltsHead1, photonPhononBoltsBody, nullptr, {0, 0,
							-type3_photonFrameTopSpace / 2 - type3_photonFrameSizeZ / 2.});

					photonPhononBoltsLV = new G4LogicalVolume(photonPhononBoltsSolid, _brass, "photonPhononBolts_LV");
					photonPhononBoltsLV->SetVisAttributes(fI_BrassBolt_VisAttr);

				}else {
					crystalPCBBoltsBodyLV = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBBoltsForCrystal_LV");
					crystalPCBBoltsHeadLV = G4LogicalVolumeStore::GetInstance()->GetVolume("PCBBoltsHeadForCrystal_LV");
					crystalClampLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Clamp_LV");
					photonPhononBoltsLV = G4LogicalVolumeStore::GetInstance()->GetVolume("photonPhononBolts_LV");
					clampBoltsLV = G4LogicalVolumeStore::GetInstance()->GetVolume("clampBolts_LV");
				}

				// Physical Volume positioning ///////
				// Crystal PCB Bolts positioning
				nowSolidPosition = {type3_photDetPCBBaseBoxSizeX / 2. - photDetPCBBoltsHoleDistX,
					                 type3_photDetPCBBaseBoxSizeY / 2. - photDetPCBBoltsHoleDistY,
													 - photDetPCBBoltsHoleHeight /2.};
				nowSolidPosition += moduleTopPosition;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsBodyLV,
						"PCBBoltsForPhoton_PV", aWhere, false, 0, OverlapCheck);
				nowSolidPosition[1] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsBodyLV,
						"PCBBoltsForPhoton_PV", aWhere, false, 1, OverlapCheck);
				nowSolidPosition[0] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsBodyLV,
						"PCBBoltsForPhoton_PV", aWhere, false, 2, OverlapCheck);

				nowSolidPosition = {-type3_copperFrameSizeX / 2. + type3_smallBlockSizeX*3/2 + photDetPCBBoltsHoleRadius,
					                  type3_copperFrameSizeY / 2. - type3_crystalPCBShortY - type3_crystalPCBGrooveY,
					                  type3_smallBlock3SizeZ - type3_bottomSupportBSizeZ - photDetPCBBoltsHeadHeight / 2.};
				nowSolidPosition -= moduleTopPosition;

				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
						"PCBBoltsForCrystal_PV", aWhere, false, 0, OverlapCheck);
				nowSolidPosition[1] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
						"PCBBoltsForCrystal_PV", aWhere, false, 1, OverlapCheck);
				nowSolidPosition[0] = - type3_smallBlockSizeX/3;// + photDetPCBBoltsHoleRadius*1.45; 
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
						"PCBBoltsForCrystal_PV", aWhere, false, 2, OverlapCheck);
				nowSolidPosition[1] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
						"PCBBoltsForCrystal_PV", aWhere, false, 3, OverlapCheck);

				nowSolidPosition[0] = type3_copperFrameSizeX / 2. - type3_smallBlockSizeX - type3_bottomSupportADistX;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
						"PCBBoltsForCrystal_PV", aWhere, false, 4, OverlapCheck);
				nowSolidPosition[1] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
						"PCBBoltsForCrystal_PV", aWhere, false, 5, OverlapCheck);

				// Clamp positioning
				G4double clampPosX = type3_copperFrameSizeX / 2. - type3_smallBlockSizeX - clampSizeX * 2;
				G4double clampPosY = type3_copperFrameSizeY / 2. - type3_smallBlockSizeY - clampSizeX * 2;
				nowSolidPosition = G4ThreeVector(clampPosX, 0, 
					-type3_photonFrameTopSpace - type3_photonFrameHoleThick - type3_photonFramePeekThick 
						-type3_germaniumWaferThick - type3_photonFrameClampThick / 2.);
				nowSolidPosition += moduleTopPosition;

				G4RotationMatrix *nowClampRotMtx = new G4RotationMatrix();

				new G4PVPlacement(nullptr, nowSolidPosition, crystalClampLV,
						"ClampForCrystal_PV", aWhere, false, 0, OverlapCheck);

				nowSolidPosition[0] -= clampPosX;
				nowSolidPosition[1] += clampPosY;
				nowClampRotMtx->rotateZ(90.*deg);
				new G4PVPlacement(G4Transform3D(*nowClampRotMtx, nowSolidPosition), crystalClampLV,
						"ClampForCrystal_PV", aWhere, false, 1, OverlapCheck);

				nowSolidPosition[0] -= clampPosX;
				nowSolidPosition[1] -= clampPosY;
				nowClampRotMtx->rotateZ(90.*deg);
				new G4PVPlacement(G4Transform3D(*nowClampRotMtx, nowSolidPosition), crystalClampLV,
						"ClampForCrystal_PV", aWhere, false, 2, OverlapCheck);
				
				nowSolidPosition[0] += clampPosX;
				nowSolidPosition[1] -= clampPosY;
				nowClampRotMtx->rotateZ(90.*deg);
				new G4PVPlacement(G4Transform3D(*nowClampRotMtx, nowSolidPosition), crystalClampLV,
						"ClampForCrystal_PV", aWhere, false, 3, OverlapCheck);

				// clamp bolts positioning
				nowSolidPosition = {type3_copperFrameSizeX / 2 - type3_copperFrameXWidth / 2, 0, 
					-type3_photonFrameTopSpace - type3_photonFrameSizeZ + type3_photonPhononBoltsHeadThick / 2};
				nowSolidPosition += moduleTopPosition;

				new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsLV,
						"ClampBoltsForCrystal_PV", aWhere, false, 0, OverlapCheck);

				nowSolidPosition[0] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsLV,
						"ClampBoltsForCrystal_PV", aWhere, false, 1, OverlapCheck);

				nowSolidPosition[0] += (type3_copperFrameSizeX - type3_copperFrameXWidth)/ 2.;
				nowSolidPosition[1] += (type3_copperFrameSizeY - type3_copperFrameYWidth)/ 2.;
				new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsLV,
						"ClampBoltsForCrystal_PV", aWhere, false, 2, OverlapCheck);

				nowSolidPosition[1] *= -1;
				new G4PVPlacement(nullptr, nowSolidPosition, clampBoltsLV,
						"ClampBoltsForCrystal_PV", aWhere, false, 3, OverlapCheck);

				// photon-phonon bolts positioning
				nowSolidPosition = {-type3_photonFrameSizeX / 2. + type3_photonFrameBaseWidthX / 2.,
					type3_photonFrameSizeY / 2. - type3_photonFrameBaseWidthY / 2.,
					-type3_photonFrameTopSpace / 2. };
				nowSolidPosition += moduleTopPosition;

				for(int i = 0; i < 4; i++)
				{
					if(i==1 || i==3) nowSolidPosition[0] *= -1.;
					else if (i==2) nowSolidPosition[1] *= -1.;
					new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsLV,
						"photonPhononBolts_PV", aWhere, false, i, OverlapCheck);
				}

				// Mass
				//G4double mass_pp = photonPhononBoltsLV->GetMass()/g;
				//G4double mass_PCB = crystalPCBBoltsBodyLV->GetMass()/g;
				//G4double mass_clamp = crystalClampLV->GetMass()/g;
				//G4double mass_clampB = clampBoltsLV->GetMass()/g;

        if (aPlacePartsForJustBelow) {
					// Crystal PCB Bolts head -----------------------
					nowSolidPosition = {type3_photDetPCBBaseBoxSizeX / 2. - photDetPCBBoltsHoleDistX,
															type3_photDetPCBBaseBoxSizeY / 2. - photDetPCBBoltsHoleDistY,
															+ photDetPCBBoltsHeadHeight /2.};
					nowSolidPosition -= moduleTopPosition;

					new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
							"PCBBoltsForPhoton_PV", aWhere, false, 0, OverlapCheck);
					nowSolidPosition[1] *= -1;
					new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
							"PCBBoltsForPhoton_PV", aWhere, false, 1, OverlapCheck);
					nowSolidPosition[0] *= -1;
					new G4PVPlacement(nullptr, nowSolidPosition, crystalPCBBoltsHeadLV,
							"PCBBoltsForPhoton_PV", aWhere, false, 2, OverlapCheck);

					// photon-phonon bolts positioning
					G4Tubs *photonPhononBoltsHead2 = new G4Tubs("photonPhononBoltsHead2Tubs",0,
							type3_photonPhononBoltsHeadRadius, 
							(type3_photonPhononBoltsHeadThick - type3_photonFrameTopSpace) / 2., 0, 360. * deg);
					G4LogicalVolume *photonPhononBoltsHeadLV = new G4LogicalVolume(photonPhononBoltsHead2, 
							_brass, "photonPhononBoltsHead_LV");
					photonPhononBoltsHeadLV->SetVisAttributes(fI_BrassBolt_VisAttr);
					
					nowSolidPosition = {-type3_photonFrameSizeX / 2. + type3_photonFrameBaseWidthX / 2.,
						type3_photonFrameSizeY / 2. - type3_photonFrameBaseWidthY / 2.,
						(type3_photonPhononBoltsHeadThick - type3_photonFrameTopSpace) / 2. };
					nowSolidPosition -= moduleTopPosition;

					for(int i = 0; i < 4; i++)
					{
						if(i==1 || i==3) nowSolidPosition[0] *= -1.;
						else if (i==2) nowSolidPosition[1] *= -1.;
						new G4PVPlacement(nullptr, nowSolidPosition, photonPhononBoltsHeadLV,
							"photonPhononBolts_PV", aWhere, false, i, OverlapCheck);
					}

					// Mass print out
					//mass_pp += photonPhononBoltsHeadLV->GetMass() / g;
					//mass_PCB += crystalPCBBoltsHeadLV->GetMass() / g;
				}
				//G4cout << " ======== Brass Bolts Mass =========" << G4endl;
				//G4cout << " photon-phonon bolts : " << mass_pp << " g" << G4endl;
				//G4cout << " PCB botls           : " << mass_PCB << " g" << G4endl;
				//G4cout << " clamp               : " << mass_clamp << " g" << G4endl;
				//G4cout << " clamp bolts         : " << mass_clampB << " g" << G4endl;
				//G4cout << " ===================================" << G4endl;
				break;
			}
      case 5: {
				logicalVolumeChecker = G4LogicalVolumeStore::GetInstance()->GetVolume("BoltM4_LV", false);
        if (logicalVolumeChecker == nullptr) 
				{
          boltM4  = new G4Tubs("BoltM4_Tub", 0, detModuleBoltM4Radius,
                                    detModuleBoltM4Length / 2., 0, 360 * deg);
					boltM4LV  = new G4LogicalVolume(boltM4, _brass, "BoltM4_LV");
          boltM4LV->SetVisAttributes(fI_BrassBolt_VisAttr);
				}else 
				{
					G4LogicalVolumeStore *nowInstance = G4LogicalVolumeStore::GetInstance();
          boltM4LV  = nowInstance->GetVolume("BoltM4_LV", false);
				}
        // SSCB BOLTS: M4-12
        nowSolidPosition = G4ThreeVector(-type3_copperFrameSizeX / 2. + type3_copperFrameHoleDistance,
                            -type3_copperFrameSizeY / 2. + type3_copperFrameHoleDistance,
														-type3_smallBlock1SizeZ + detModuleBoltM4Length / 2.);
                            //aCrystal.fCrystalHeight / 2. + detModuleBoltM4Length / 2.);
				nowSolidPosition += moduleTopPosition;

        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                              0, OverlapCheck);
        nowSolidPosition[0] *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                            1, OverlapCheck);
        nowSolidPosition[0] *= -1.;
        nowSolidPosition[1] *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                            2, OverlapCheck);
        nowSolidPosition[0] *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Top_PV", aWhere, false,
                            3, OverlapCheck);

				nowSolidPosition[2] = type3_smallBlock3SizeZ - detModuleBoltM4Length / 2.;
				nowSolidPosition -= moduleTopPosition;
        //nowSolidPosition *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                            0, OverlapCheck);
        nowSolidPosition[0] *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                            1, OverlapCheck);
        nowSolidPosition[0] *= -1.;
        nowSolidPosition[1] *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                          2, OverlapCheck);
        nowSolidPosition[0] *= -1.;
        new G4PVPlacement(nullptr, nowSolidPosition, boltM4LV, "BoltM4Bottom_PV", aWhere, false,
                              3, OverlapCheck);

        break;
        }
    }
}

G4LogicalVolume *AmoreDetectorConstruction::Build_I_MuonVetoPMT(G4Material *aShieldMaterial,
                                                                G4Material *aGreaseMaterial) {
    G4LogicalVolume *resultLV;
    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    constexpr const char *pmtLVName = "PMTEnvelope_LV";

    G4double muonVetoPMTWindowSizeR     = AMoRE_unit_inch * 3. / 2.;
    G4double muonVetoPMTTopSizeRmax     = std::max(muonVetoPMTWindowSizeR, muonVetoPMTSizeR);
    G4double muonVetoPMTVacuSizeR       = muonVetoPMTSizeR - muonVetoPMTTubeThickness;
    G4double muonVetoPMTVacuSizeZ       = muonVetoPMTSizeZ - muonVetoPMTBacksideThickness;
    G4double muonVetoPMTVacuTopSizeRmax = muonVetoPMTWindowSizeR - muonVetoPMTTubeThickness;
    G4double muonVetoPMTGreaseSizeR     = muonVetoPMTWindowSizeR - muonVetoPMTGreaseShieldThick;
    G4double muonVetoPMTEnvelopeSizeZ =
        muonVetoPMTSizeZ + muonVetoPMTWindowThickness + muonVetoPMTGreaseThickness;
    G4double tempZCoord;

    G4VSolid *muonVetoPMTTube;
    G4VSolid *muonVetoPMTTop;
    G4VSolid *muonVetoPMTWindow;
    G4VSolid *muonVetoPMTNowSolid;
    G4VSolid *muonVetoPMTVacuTube;
    G4VSolid *muonVetoPMTVacuTop;
    G4VSolid *solidPMTVacu;
    G4VSolid *solidPMTGrease;
    G4VSolid *solidPMTGreaseShield;
    G4VSolid *solidPMTEnvelope;

    G4LogicalVolume *logicPMTGrease;
    G4LogicalVolume *logicPMTGS;

    G4VisAttributes *muonVetoPMTVisAtt;
    G4VisAttributes *muonVetoPMTVacuVisAtt;
    G4VisAttributes *muonVetoPMTEnvelopeVisAtt;
    G4VisAttributes *muonVetoPMTGreaseVisAtt;

    G4VPhysicalVolume *physiPMTVacu;
    G4VPhysicalVolume *physiPMT;

    resultLV = G4LogicalVolumeStore::GetInstance()->GetVolume(pmtLVName, false);

    if (resultLV == nullptr) {
        muonVetoPMTTube =
            new G4Tubs("PMTBody_Tub", 0, muonVetoPMTSizeR, muonVetoPMTSizeZ / 2., 0, 360. * deg);
        muonVetoPMTTop    = new G4Sphere("PMTTop_HalfSphere", muonVetoPMTTopSizeRmin,
                                      muonVetoPMTTopSizeRmax, 0, 360 * deg, 90 * deg, 180 * deg);
        muonVetoPMTWindow = new G4Tubs("PMTWindow_Tub", 0, muonVetoPMTTopSizeRmax,
                                       muonVetoPMTWindowThickness / 2., 0, 360. * deg);
        muonVetoPMTNowSolid =
            new G4UnionSolid("PMTSolidAssemStage1_Solid", muonVetoPMTTube, muonVetoPMTTop, 0,
                             G4ThreeVector(0, 0, muonVetoPMTSizeZ / 2));
        muonVetoPMTNowSolid = new G4UnionSolid(
            "PMTSolid", muonVetoPMTNowSolid, muonVetoPMTWindow, 0,
            G4ThreeVector(0, 0, muonVetoPMTSizeZ / 2 + muonVetoPMTWindowThickness / 2));

        muonVetoPMTVacuTube = new G4Tubs("PMTVacubody0_Tub", 0, muonVetoPMTVacuSizeR,
                                         muonVetoPMTVacuSizeZ / 2, 0, 360. * deg);
        muonVetoPMTVacuTop =
            new G4Sphere("PMTVacuTop", muonVetoPMTVacuTopSizeRmin, muonVetoPMTVacuTopSizeRmax, 0,
                         360 * deg, 90 * deg, 180 * deg);
        solidPMTVacu = new G4UnionSolid("PMTVacu_Solid", muonVetoPMTVacuTube, muonVetoPMTVacuTop, 0,
                                        G4ThreeVector(0, 0, muonVetoPMTVacuSizeZ / 2));
        solidPMTGrease       = new G4Tubs("PMTGrease_Tub", 0, muonVetoPMTGreaseSizeR,
                                    muonVetoPMTGreaseThickness / 2, 0, 360. * deg);
        solidPMTGreaseShield = new G4Tubs("PMTGreaseShield_Tub", 0, muonVetoPMTWindowSizeR,
                                          muonVetoPMTGreaseThickness / 2., 0, 360. * deg);
        logicPMTGrease       = new G4LogicalVolume(solidPMTGrease, aGreaseMaterial, "PMTGrease_LV");
        logicPMTGS =
            new G4LogicalVolume(solidPMTGreaseShield, aShieldMaterial, "PMTGreaseShield_LV");

        fI_logicPMT     = new G4LogicalVolume(muonVetoPMTNowSolid, _glass, "PMT_LV");
        fI_logicPMTVacu = new G4LogicalVolume(solidPMTVacu, _vacuum, "PMTVacu_LV");

        solidPMTEnvelope = new G4Tubs("PMTEnvelope", 0, muonVetoPMTTopSizeRmax,
                                      muonVetoPMTEnvelopeSizeZ / 2., 0, 360. * deg);
        resultLV         = new G4LogicalVolume(solidPMTEnvelope, _N2_Gas, pmtLVName);

        tempZCoord = -muonVetoPMTEnvelopeSizeZ / 2. + muonVetoPMTSizeZ / 2.;

        physiPMT = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, tempZCoord), fI_logicPMT,
                                     "PMT_PV", resultLV, false, 0, OverlapCheck);

        tempZCoord = -muonVetoPMTBacksideThickness / 2;

        physiPMTVacu = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, tempZCoord), fI_logicPMTVacu,
                                         "phys_pmtvacu", fI_logicPMT, false, 0, OverlapCheck);

        tempZCoord = muonVetoPMTEnvelopeSizeZ / 2 - muonVetoPMTGreaseThickness / 2;

        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, tempZCoord), logicPMTGS,
                          "PMTGreaseShield_PV", resultLV, false, 0, OverlapCheck);
        new G4PVPlacement(nullptr, {}, logicPMTGrease, "PMTGrease_PV", logicPMTGS, false, 0, OverlapCheck);

        // for visual
        muonVetoPMTVisAtt         = new G4VisAttributes(G4Colour::Magenta());
        muonVetoPMTVacuVisAtt     = new G4VisAttributes(G4Colour::Brown());
        muonVetoPMTEnvelopeVisAtt = new G4VisAttributes(G4Colour::Gray());
        muonVetoPMTGreaseVisAtt   = new G4VisAttributes(G4Colour::Gray());

        muonVetoPMTVisAtt->SetForceWireframe(true);
        muonVetoPMTVacuVisAtt->SetForceWireframe(true);
        muonVetoPMTGreaseVisAtt->SetForceWireframe(true);

        fI_logicPMT->SetVisAttributes(muonVetoPMTVisAtt);
        fI_logicPMTVacu->SetVisAttributes(muonVetoPMTVacuVisAtt);
        resultLV->SetVisAttributes(muonVetoPMTEnvelopeVisAtt);
        logicPMTGrease->SetVisAttributes(muonVetoPMTGreaseVisAtt);

        new G4LogicalBorderSurface("MLCS_photocathode_logsurf",
                                   physiPMT, // exiting glass into vac.
                                   physiPMTVacu, Photocathode_opsurf);

        G4Region *PmtRegion = new G4Region("MLCS");
        PmtRegion->AddRootLogicalVolume(fI_logicPMT);
        new CupPMTOpticalModel("MLCS_optical_model", physiPMT);
    }

    return resultLV;
}

void AmoreDetectorConstruction::Place_I_RealNeutronShieldingAt(G4LogicalVolume *aWhere,
                                                               G4ThreeVector aTlate, G4int aType) {

    using namespace AmoreDetectorStaticInfo;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;

    G4LogicalVolume *logicalVolumeChecker;
    // Polyethylene geometry
    G4Box *polyEthylene_Top_Box;
    G4Box *polyEthylene_Bottom_Box;
    G4Box *polyEthylene_SideFB_Box;
    G4Box *polyEthylene_SideLR1_Box;
    G4Box *polyEthylene_SideLR2_Box;
    G4Box *polyEthylene_MufflerLR_Box;
    G4Box *polyEthylene_MufflerFB_Box;
    G4Box *bottomScintDummy_Box;

    G4LogicalVolume *polyEthylene_Top_BoxLogical;
    G4LogicalVolume *polyEthylene_Bottom_BoxLogical;
    G4LogicalVolume *polyEthylene_SideFB_BoxLogical;
    G4LogicalVolume *polyEthylene_SideLR1_BoxLogical;
    G4LogicalVolume *polyEthylene_SideLR2_BoxLogical;
    G4LogicalVolume *polyEthylene_MufflerLR_BoxLogical;
    G4LogicalVolume *polyEthylene_MufflerFB_BoxLogical;
    G4LogicalVolume *bottomScintDummy_BoxLogical;

    // Borated PE geometry
    G4Box *boratedPE_Top_Box;
    G4Box *boratedPE_Bottom_Box;
    G4Box *boratedPE_SideFB_Box;
    G4Box *boratedPE_SideLR_Box;
    G4Box *boratedPE_SideLR1_Box;
    G4Box *boratedPE_SideLR2_Box;
    G4Box *boratedPE_MufflerLR_Box;
    G4Box *boratedPE_MufflerFB_Box;

    G4LogicalVolume *boratedPE_Top_BoxLogical;
    G4LogicalVolume *boratedPE_Bottom_BoxLogical;
    G4LogicalVolume *boratedPE_SideFB_BoxLogical;
    G4LogicalVolume *boratedPE_SideLR_BoxLogical;
    G4LogicalVolume *boratedPE_SideLR1_BoxLogical;
    G4LogicalVolume *boratedPE_SideLR2_BoxLogical;
    G4LogicalVolume *boratedPE_MufflerLR_BoxLogical;
    G4LogicalVolume *boratedPE_MufflerFB_BoxLogical;

		// Boric acid rubber
		G4VSolid *BoricAcid_Rubber;
		G4Box *BoricAcid_Bottom;

		G4LogicalVolume *BoricAcid_RubberLogical;
		G4LogicalVolume *BoricAcid_BottomLogical;

    // Boric acid panel geometry
    G4Box *sideLRBoricAcid_HousingBox; // X-Axis: FB | Y-Axis: LR
    G4Box *sideLRBoricAcid_Block;
    G4Box *sideFBBoricAcid_HousingBox;
    G4Box *sideFBBoricAcid_Block;
    G4Box *bottomBoricAcid_HousingBox;
    G4Box *bottomBoricAcid_Block;
    G4Box *topPbAboveRubberB4C_Box;
    G4Box *detectorRubberB4C_Box;

    G4LogicalVolume *sideLRBoricAcid_MotherLogical;
    G4LogicalVolume *sideLRBoricAcid_BlockLogical;
    G4LogicalVolume *sideFBBoricAcid_MotherLogical;
    G4LogicalVolume *sideFBBoricAcid_BlockLogical;
    G4LogicalVolume *bottomBoricAcid_MotherLogical;
    G4LogicalVolume *bottomBoricAcid_BlockLogical;
    G4LogicalVolume *topPbAboveRubberB4C_BoxLogical;
    G4LogicalVolume *detectorRubberB4C_BoxLogical;

    G4ThreeVector sideLRBoricAcid_Pos[2];
    G4ThreeVector sideFBBoricAcid_Pos[2];
    G4ThreeVector bottomBoricAcid_Pos;
    G4ThreeVector nowSolidDisplacer;
    G4ThreeVector FB_Trans;
    G4ThreeVector LR_Trans;

    G4double sideLRScint_YPos;
    G4double leadBoxsize;
    G4double leadBoxzsize;
    G4double length_Bottom;
    G4double width_Bottom;
    G4double height_LR;
    G4double height_FB;
    G4double width_LR;
    G4double width_FB;
    G4double totalThickness_BoricAcidPanel;
    G4double topScint_boxThick_half;

    switch (aType) {
        case 1: { // For AMoRE-I (From J.S. Lee)
									// Side BPE, PE and bottom BPE,PE,PS

            leadBoxsize  = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
            leadBoxzsize = (leadBoxZinnerlength + leadBoxThickness) / 2.0;

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("PolyEthylene_FB_LV", false);
            if (logicalVolumeChecker == nullptr) {
                polyEthylene_SideFB_Box =
                    new G4Box("PolyEthylene_FB_Box", polyEthylene_commonThickness / 2.,
                              2. * sideFBScint_boxWidth_half, sideFBScint_boxHeight_half);
                polyEthylene_SideFB_BoxLogical = new G4LogicalVolume(
                    polyEthylene_SideFB_Box, _polyethylene, "PolyEthylene_FB_LV");

                polyEthylene_SideLR1_Box = new G4Box(
                    "PolyEthylene_LR1_Box", 2. * sideLRScint_boxWidth_half + sideLRScint_trapZ_half,
                    polyEthylene_commonThickness / 2.,
                    sideLRScint_boxHeight_half - polyEthyleneLR_cuttingLength / 2.);
                polyEthylene_SideLR1_BoxLogical = new G4LogicalVolume(
                    polyEthylene_SideLR1_Box, _polyethylene, "PolyEthylene_LR1_LV");

                polyEthylene_Bottom_Box = new G4Box(
                    "PolyEthylene_Bottom_Box", (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    polyEthylene_commonThickness / 2.);
                polyEthylene_Bottom_BoxLogical = new G4LogicalVolume(
                    polyEthylene_Bottom_Box, _polyethylene, "PolyEthylene_Bottom_LV");

                bottomScintDummy_Box = new G4Box(
                    "BottomScintillator_Dummy_Box",
                    (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2., scint_bottom_Thick_half);
                bottomScintDummy_BoxLogical = new G4LogicalVolume(bottomScintDummy_Box, _vinylt,
                                                                  "BottomScintillator_Dummy_LV");

                polyEthylene_SideFB_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                polyEthylene_SideLR1_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                polyEthylene_Bottom_BoxLogical->SetVisAttributes(fI_PE_VisAttr);

                boratedPE_SideFB_Box =
                    new G4Box("BoratedPE_SideFB_Box", boratedPE_thickness / 2.,
                              2. * sideFBScint_boxWidth_half, sideFBScint_boxHeight_half);
                boratedPE_SideFB_BoxLogical = new G4LogicalVolume(
                    boratedPE_SideFB_Box, _BoratedPE_5perCent, "BoratedPE_SideFB_LV");

                boratedPE_SideLR_Box = new G4Box(
                    "BoratedPE_SideLR_Box", 2. * sideLRScint_boxWidth_half + sideLRScint_trapZ_half,
                    boratedPE_thickness / 2., sideLRScint_boxHeight_half);
                boratedPE_SideLR_BoxLogical = new G4LogicalVolume(
                    boratedPE_SideLR_Box, _BoratedPE_5perCent, "BoratedPE_SideLR_LV");

                boratedPE_Bottom_Box = new G4Box(
                    "BoratedPE_Bottom_Box", (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2., boratedPE_thickness / 2.);
                boratedPE_Bottom_BoxLogical = new G4LogicalVolume(
                    boratedPE_Bottom_Box, _BoratedPE_5perCent, "BoratedPE_Bottom_LV");

                boratedPE_SideFB_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
                boratedPE_SideLR_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
                boratedPE_Bottom_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
            } else {
                G4LogicalVolumeStore *instance  = G4LogicalVolumeStore::GetInstance();
                boratedPE_SideFB_BoxLogical     = instance->GetVolume("BoratedPE_SideFB_LV");
                boratedPE_SideLR_BoxLogical     = instance->GetVolume("BoratedPE_SideLR_LV");
                boratedPE_Bottom_BoxLogical     = instance->GetVolume("BoratedPE_Bottom_LV");
                polyEthylene_SideFB_BoxLogical  = instance->GetVolume("PolyEthylene_FB_LV");
                polyEthylene_SideLR1_BoxLogical = instance->GetVolume("PolyEthylene_LR1_LV");
                polyEthylene_Bottom_BoxLogical  = instance->GetVolume("PolyEthylene_Bottom_LV");
                bottomScintDummy_BoxLogical = instance->GetVolume("BottomScintillator_Dummy_LV");
            }

            FB_Trans = {0, 0, -leadBoxzsize + (sideFBScint_boxHeight_half + scint_reflector_thick)};
            LR_Trans = {0, 0, -leadBoxzsize + (sideLRScint_boxHeight_half + scint_reflector_thick)};
            sideLRScint_YPos =
                std::fmax(2 * (sideFBScint_boxWidth_half + scint_reflector_thick), leadBoxsize);

            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans +
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      boratedPE_thickness / 2.,
                                  0, 0),
                boratedPE_SideFB_BoxLogical, "BoratedPE_SideFB_PV", aWhere, false, 0, OverlapCheck);
            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans -
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      boratedPE_thickness / 2.,
                                  0, 0),
                boratedPE_SideFB_BoxLogical, "BoratedPE_SideFB_PV", aWhere, false, 1, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans +
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      boratedPE_thickness + polyEthylene_commonThickness / 2.,
                                  0, 0),
                polyEthylene_SideFB_BoxLogical, "PolyEthylene_SideFB_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans -
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      boratedPE_thickness + polyEthylene_commonThickness / 2.,
                                  0, 0),
                polyEthylene_SideFB_BoxLogical, "PolyEthylene_SideFB_PV", aWhere, false, 1, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans +
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      boratedPE_thickness / 2.,
                                  0),
                boratedPE_SideLR_BoxLogical, "BoratedPE_SideLR_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans -
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      boratedPE_thickness / 2.,
                                  0),
                boratedPE_SideLR_BoxLogical, "BoratedPE_SideLR_PV", aWhere, false, 1, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans -
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      boratedPE_thickness - polyEthylene_commonThickness / 2.,
                                  0) +
                    G4ThreeVector(0, 0, -polyEthyleneLR_cuttingLength / 2.),
                polyEthylene_SideLR1_BoxLogical, "PolyEthylene_SideLR_PV", aWhere, false, 1, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans +
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      boratedPE_thickness - polyEthylene_commonThickness / 2.,
                                  0) +
                    G4ThreeVector(0, 0, -polyEthyleneLR_cuttingLength / 2.),
                polyEthylene_SideLR1_BoxLogical, "PolyEthylene_SideLR_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate - G4ThreeVector(0, 0,
                                       (leadBoxZinnerlength + leadBoxThickness) / 2. +
                                           scint_bottom_Thick_half + bottomNeutronShieldGap),
                bottomScintDummy_BoxLogical, "BottomScintillator_Dummy_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate - G4ThreeVector(0, 0,
                                       (leadBoxZinnerlength + leadBoxThickness) / 2. +
                                           scint_bottom_Thick_half * 2. + bottomNeutronShieldGap +
                                           polyEthylene_commonThickness / 2.),
                polyEthylene_Bottom_BoxLogical, "PolyEthylene_Bottom_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate - G4ThreeVector(0, 0,
                                       (leadBoxZinnerlength + leadBoxThickness) / 2. +
                                           scint_bottom_Thick_half * 2. + boratedPE_thickness / 2. +
                                           bottomNeutronShieldGap + polyEthylene_commonThickness),
                boratedPE_Bottom_BoxLogical, "BoratedPE_Bottom_PV", aWhere, false, 0, OverlapCheck);
        } break;
        case 2: { // Top BPE, PE
            topScint_boxThick_half = scint_Thick_half;

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("PolyEthylene_Top_LV", false);
            if (logicalVolumeChecker == nullptr) {
                polyEthylene_Top_Box = new G4Box(
                    "PolyEthylene_Top_Box", (topScint_boxWidth_half + scint_reflector_thick) * 2.,
                    topScint_boxHeight_half + scint_reflector_thick,
                    polyEthylene_topThickness / 2.);
                boratedPE_Top_Box = new G4Box(
                    "BoratedPE_Top_Box", (topScint_boxWidth_half + scint_reflector_thick) * 2.,
                    topScint_boxHeight_half + scint_reflector_thick, boratedPE_thickness / 2.);

                polyEthylene_Top_BoxLogical =
                    new G4LogicalVolume(polyEthylene_Top_Box, _polyethylene, "PolyEthylene_Top_LV");
                polyEthylene_Top_BoxLogical->SetVisAttributes(fI_PE_VisAttr);

                boratedPE_Top_BoxLogical =
                    new G4LogicalVolume(boratedPE_Top_Box, _BoratedPE_5perCent, "BoratedPE_Top_LV");
                boratedPE_Top_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);

            } else {
                G4LogicalVolumeStore *instance = G4LogicalVolumeStore::GetInstance();
                polyEthylene_Top_BoxLogical    = instance->GetVolume("PolyEthylene_Top_LV");
                boratedPE_Top_BoxLogical       = instance->GetVolume("BoratedPE_Top_LV");
            }

            new G4PVPlacement(
                nullptr,
                aTlate + G4ThreeVector(0, 0,
                                       pbtopbox_zsize +
                                           (topScint_boxThick_half + scint_reflector_thick) * 2. +
                                           boratedPE_thickness / 2.),
                boratedPE_Top_BoxLogical, "BoratedPE_Top_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + G4ThreeVector(0, 0,
                                       pbtopbox_zsize +
                                           (topScint_boxThick_half + scint_reflector_thick) * 2. +
                                           boratedPE_thickness + polyEthylene_topThickness / 2.),
                polyEthylene_Top_BoxLogical, "PolyEthylene_Top_PV", aWhere, false, 0, OverlapCheck);
        } break;
        case 3: { // BoricAcid LR, FB, Bottom
            totalThickness_BoricAcidPanel = thickness_BoricAcid + thickness_BoricAcidHousing * 2.;
            leadBoxsize                   = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
            leadBoxzsize                  = (leadBoxZinnerlength + leadBoxThickness) / 2.0;
            height_LR                     = leadBoxzsize - leadBoxThickness / 2.;
            width_LR                      = leadBoxsize - leadBoxThickness;
            height_FB                     = leadBoxzsize - leadBoxThickness / 2.;
            width_FB      = leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel;
            length_Bottom = width_LR - totalThickness_BoricAcidPanel;
            width_Bottom  = width_FB - totalThickness_BoricAcidPanel;

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("BoricAcid_LR_Mother_LV", false);
            if (logicalVolumeChecker == nullptr) {
                sideLRBoricAcid_HousingBox =
                    new G4Box("BoricAcid_LR_Mother_Box", width_LR,
                              totalThickness_BoricAcidPanel / 2., height_LR);
                sideFBBoricAcid_HousingBox =
                    new G4Box("BoricAcid_FB_Mother_Box", totalThickness_BoricAcidPanel / 2.,
                              width_FB, height_FB);
                bottomBoricAcid_HousingBox =
                    new G4Box("BoricAcid_Bottom_Mother_Box", length_Bottom, width_Bottom,
                              totalThickness_BoricAcidPanel / 2.);

                sideLRBoricAcid_Block =
                    new G4Box("BoricAcid_LR_Block", width_LR - thickness_BoricAcidSpacing,
                              thickness_BoricAcid / 2., height_LR - thickness_BoricAcidSpacing);
                sideFBBoricAcid_Block = new G4Box("BoricAcid_FBBlock_Box", thickness_BoricAcid / 2.,
                                                  width_FB - thickness_BoricAcidSpacing,
                                                  height_FB - thickness_BoricAcidSpacing);
                bottomBoricAcid_Block = new G4Box(
                    "BoricAcid_BottomBlock_Box", length_Bottom - thickness_BoricAcidSpacing,
                    width_Bottom - thickness_BoricAcidSpacing, thickness_BoricAcid / 2.);

                sideLRBoricAcid_MotherLogical = new G4LogicalVolume(
                    sideLRBoricAcid_HousingBox, _polyethylene, "BoricAcid_LR_Mother_LV");
                sideFBBoricAcid_MotherLogical = new G4LogicalVolume(
                    sideFBBoricAcid_HousingBox, _polyethylene, "BoricAcid_FB_Mother_LV");
                bottomBoricAcid_MotherLogical = new G4LogicalVolume(
                    bottomBoricAcid_HousingBox, _polyethylene, "BoricAcid_Bottom_Mother_LV");

                sideLRBoricAcid_BlockLogical = new G4LogicalVolume(
                    sideLRBoricAcid_Block, _BoricAcidPowder, "BoricAcid_LR_Block_LV");
                sideFBBoricAcid_BlockLogical = new G4LogicalVolume(
                    sideFBBoricAcid_Block, _BoricAcidPowder, "BoricAcid_FB_Block_LV");
                bottomBoricAcid_BlockLogical = new G4LogicalVolume(
                    bottomBoricAcid_Block, _BoricAcidPowder, "BoricAcid_Bottom_Block_LV");

                new G4PVPlacement(nullptr, {}, sideLRBoricAcid_BlockLogical,
                                  "BoricAcid_LR_Block_PV", sideLRBoricAcid_MotherLogical, false, 0, OverlapCheck);
                new G4PVPlacement(nullptr, {}, sideFBBoricAcid_BlockLogical,
                                  "BoricAcid_FB_Block_PV", sideFBBoricAcid_MotherLogical, false, 0, OverlapCheck);
                new G4PVPlacement(nullptr, {}, bottomBoricAcid_BlockLogical,
                                  "BoricAcid_Bottom_Block_PV", bottomBoricAcid_MotherLogical, false,
                                  0, OverlapCheck);

                bottomBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
                bottomBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);

                sideLRBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
                sideLRBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);

                sideFBBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
                sideFBBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);
            } else {
                G4LogicalVolumeStore *instance = G4LogicalVolumeStore::GetInstance();
                sideLRBoricAcid_MotherLogical  = instance->GetVolume("BoricAcid_LR_Mother_LV");
                sideFBBoricAcid_MotherLogical  = instance->GetVolume("BoricAcid_FB_Mother_LV");
                bottomBoricAcid_MotherLogical  = instance->GetVolume("BoricAcid_Bottom_Mother_LV");
            }

            sideLRBoricAcid_Pos[0] = {
                0, leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel / 2.,
                -leadBoxzsize + leadBoxThickness + height_LR};
            sideLRBoricAcid_Pos[1] = {
                0, -leadBoxsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2.,
                -leadBoxzsize + leadBoxThickness + height_LR};
            sideFBBoricAcid_Pos[0] = {leadBoxsize - leadBoxThickness -
                                          totalThickness_BoricAcidPanel / 2.,
                                      0, -leadBoxzsize + leadBoxThickness + height_FB};
            sideFBBoricAcid_Pos[1] = {-leadBoxsize + leadBoxThickness +
                                          totalThickness_BoricAcidPanel / 2.,
                                      0, -leadBoxzsize + leadBoxThickness + height_FB};

            bottomBoricAcid_Pos = {
                0, 0, -leadBoxzsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2.};

            new G4PVPlacement(nullptr, aTlate + sideLRBoricAcid_Pos[0],
                              sideLRBoricAcid_MotherLogical, "BoricAcid_LR_Mother_PV", aWhere,
                              false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, aTlate + sideLRBoricAcid_Pos[1],
                              sideLRBoricAcid_MotherLogical, "BoricAcid_LR_Mother_PV", aWhere,
                              false, 1, OverlapCheck);

            new G4PVPlacement(nullptr, aTlate + sideFBBoricAcid_Pos[0],
                              sideFBBoricAcid_MotherLogical, "BoricAcid_FB_Mother_PV", aWhere,
                              false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, aTlate + sideFBBoricAcid_Pos[1],
                              sideFBBoricAcid_MotherLogical, "BoricAcid_FB_Mother_PV", aWhere,
                              false, 1, OverlapCheck);
            new G4PVPlacement(nullptr, aTlate + bottomBoricAcid_Pos, bottomBoricAcid_MotherLogical,
                              "BoricAcid_Bottom_Mother_PV", aWhere, false, 0, OverlapCheck);
        } break;
        case 4: { // Muffler PE, BPE
            G4LogicalVolumeStore *instance = G4LogicalVolumeStore::GetInstance();

            logicalVolumeChecker = instance->GetVolume("PolyEthylene_MufflerLR_LV", false);
            if (logicalVolumeChecker == nullptr) {
                polyEthylene_MufflerLR_Box =
                    new G4Box("PolyEthylene_MufflerLR_Box", nShield_mufflerLR_SizeX / 2.,
                              polyEthylene_commonThickness / 2., nShield_mufflerLR_SizeY / 2.);
                polyEthylene_MufflerLR_BoxLogical = new G4LogicalVolume(
                    polyEthylene_MufflerLR_Box, _polyethylene, "PolyEthylene_MufflerLR_LV");

                polyEthylene_MufflerFB_Box =
                    new G4Box("PolyEthylene_MufflerFB_Box", polyEthylene_commonThickness / 2.,
                              nShield_mufflerFB_SizeX / 2., nShield_mufflerFB_SizeY / 2.);
                polyEthylene_MufflerFB_BoxLogical = new G4LogicalVolume(
                    polyEthylene_MufflerFB_Box, _polyethylene, "PolyEthylene_MufflerFB_LV");

                boratedPE_MufflerLR_Box =
                    new G4Box("BoratedPE_MufflerLR_Box", nShield_mufflerLR_SizeX / 2.,
                              boratedPE_thickness / 2., nShield_mufflerLR_SizeY / 2.);
                boratedPE_MufflerLR_BoxLogical = new G4LogicalVolume(
                    boratedPE_MufflerLR_Box, _BoratedPE_5perCent, "BoratedPE_MufflerLR_LV");

                boratedPE_MufflerFB_Box =
                    new G4Box("BoratedPE_MufflerFB_Box", boratedPE_thickness / 2.,
                              nShield_mufflerFB_SizeX / 2., nShield_mufflerFB_SizeY / 2.);
                boratedPE_MufflerFB_BoxLogical = new G4LogicalVolume(
                    boratedPE_MufflerFB_Box, _BoratedPE_5perCent, "BoratedPE_MufflerFB_LV");

                polyEthylene_MufflerLR_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                polyEthylene_MufflerFB_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                boratedPE_MufflerLR_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
                boratedPE_MufflerFB_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
            } else {
                polyEthylene_MufflerLR_BoxLogical =
                    instance->GetVolume("PolyEthylene_MufflerLR_LV");
                polyEthylene_MufflerFB_BoxLogical =
                    instance->GetVolume("PolyEthylene_MufflerFB_LV");
                boratedPE_MufflerLR_BoxLogical = instance->GetVolume("BoratedPE_MufflerLR_LV");
                boratedPE_MufflerFB_BoxLogical = instance->GetVolume("BoratedPE_MufflerFB_LV");
            }

            leadBoxsize  = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
            leadBoxzsize = (leadBoxZinnerlength + leadBoxThickness) / 2.0;
            FB_Trans = {0, 0, -leadBoxzsize + (sideFBScint_boxHeight_half + scint_reflector_thick)};
            LR_Trans = {0, 0, -leadBoxzsize + (sideLRScint_boxHeight_half + scint_reflector_thick)};
            sideLRScint_YPos =
                std::fmax(2 * (sideFBScint_boxWidth_half + scint_reflector_thick), leadBoxsize);
            nowSolidDisplacer = {
                0,
                sideLRScint_YPos + (scint_Thick_half + scint_reflector_thick) * 2. +
                    boratedPE_thickness * 2. + polyEthylene_commonThickness * 3 / 2. +
                    nShield_mufflerLR_DistY,
                nShield_mufflerLR_DistX};

            new G4PVPlacement(nullptr, aTlate + LR_Trans + nowSolidDisplacer,
                              polyEthylene_MufflerLR_BoxLogical, "PolyEthylene_MufflerLR_PV",
                              aWhere, false, 0, OverlapCheck);
            nowSolidDisplacer[1] *= -1.;
            new G4PVPlacement(nullptr, aTlate + LR_Trans + nowSolidDisplacer,
                              polyEthylene_MufflerLR_BoxLogical, "PolyEthylene_MufflerLR_PV",
                              aWhere, false, 1, OverlapCheck);

            nowSolidDisplacer = {0,
                                 sideLRScint_YPos +
                                     (scint_Thick_half + scint_reflector_thick) * 2. +
                                     boratedPE_thickness * 3 / 2. + polyEthylene_commonThickness +
                                     nShield_mufflerLR_DistY,
                                 nShield_mufflerLR_DistX};
            new G4PVPlacement(nullptr, aTlate + LR_Trans + nowSolidDisplacer,
                              boratedPE_MufflerLR_BoxLogical, "BoratedPE_MufflerLR_PV", aWhere,
                              false, 0, OverlapCheck);
            nowSolidDisplacer[1] *= -1.;
            new G4PVPlacement(nullptr, aTlate + LR_Trans + nowSolidDisplacer,
                              boratedPE_MufflerLR_BoxLogical, "BoratedPE_MufflerLR_PV", aWhere,
                              false, 1, OverlapCheck);

            nowSolidDisplacer = {leadBoxsize + scint_muffler_Thick_half * 2. +
                                     scint_reflector_thick * 2. + boratedPE_thickness +
                                     polyEthylene_commonThickness / 2. +
                                     mufflerFBScint_YDistanceFromLead,
                                 0,
                                 sideFBScint_boxHeight_half + mufflerFBScint_boxHeight_half * 2. -
                                     nShield_mufflerFB_SizeY / 2.};
            new G4PVPlacement(nullptr, aTlate + FB_Trans + nowSolidDisplacer,
                              polyEthylene_MufflerFB_BoxLogical, "PolyEthylene_MufflerFB_PV",
                              aWhere, false, 0, OverlapCheck);
            nowSolidDisplacer[0] *= -1;
            new G4PVPlacement(nullptr, aTlate + FB_Trans + nowSolidDisplacer,
                              polyEthylene_MufflerFB_BoxLogical, "PolyEthylene_MufflerFB_PV",
                              aWhere, false, 1, OverlapCheck);

            nowSolidDisplacer = {leadBoxsize + scint_muffler_Thick_half * 2. +
                                     scint_reflector_thick * 2. + boratedPE_thickness / 2. +
                                     mufflerFBScint_YDistanceFromLead,
                                 0,
                                 sideFBScint_boxHeight_half + mufflerFBScint_boxHeight_half * 2. -
                                     nShield_mufflerFB_SizeY / 2.};
            new G4PVPlacement(nullptr, aTlate + FB_Trans + nowSolidDisplacer,
                              boratedPE_MufflerFB_BoxLogical, "BoratedPE_MufflerFB_PV", aWhere,
                              false, 0, OverlapCheck);
            nowSolidDisplacer[0] *= -1;
            new G4PVPlacement(nullptr, aTlate + FB_Trans + nowSolidDisplacer,
                              boratedPE_MufflerFB_BoxLogical, "BoratedPE_MufflerFB_PV", aWhere,
                              false, 1, OverlapCheck);

        } break;
				case 5: {
						// BoricAcidRubber place at OVC surface with tube shape  instead of case 3(J.S.Lee's version) or case 13(Pilot Run7 version)

            totalThickness_BoricAcidPanel = thickness_BoricAcid;
            leadBoxsize                   = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
            leadBoxzsize                  = (leadBoxZinnerlength + leadBoxThickness) / 2.0;
            height_LR                     = leadBoxzsize - leadBoxThickness / 2.;
            width_LR                      = leadBoxsize - leadBoxThickness;
            height_FB                     = leadBoxzsize - leadBoxThickness / 2.;
            width_FB      = leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel;
            length_Bottom = width_LR - totalThickness_BoricAcidPanel;
            width_Bottom  = width_FB - totalThickness_BoricAcidPanel;
						bottomBoricAcid_Pos = {0,0, - ss_height_half - (leadBoxZinnerlength-ss_height_half*2 - ssb_zsize*2)+ totalThickness_BoricAcidPanel/2};
            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("BoricAcid_Rubber_LV", false);
            if (logicalVolumeChecker == nullptr) {
							BoricAcid_Rubber = new G4Tubs("BoricAcid_Rubber",ss_radius + ss_thick + solidBooleanTol, ss_radius + ss_thick + totalThickness_BoricAcidPanel, ss_height_half+ssb_zsize,0,360.*deg);
							BoricAcid_Bottom = new G4Box("BoricAcid_Bottom", length_Bottom, width_Bottom, totalThickness_BoricAcidPanel/2.);

							BoricAcid_RubberLogical = new G4LogicalVolume(BoricAcid_Rubber, _BoricAcidRubber, "BoricAcid_Rubber_LV");
							BoricAcid_BottomLogical = new G4LogicalVolume(BoricAcid_Bottom, _BoricAcidRubber, "BoricAcid_Bottom_LV");
							BoricAcid_RubberLogical -> SetVisAttributes(fI_BoricAcid_VisAttr);
							BoricAcid_BottomLogical -> SetVisAttributes(fI_BoricAcid_VisAttr);
						} else{
                G4LogicalVolumeStore *instance = G4LogicalVolumeStore::GetInstance();
                BoricAcid_RubberLogical  = instance->GetVolume("BoricAcid_Rubber_LV");
								BoricAcid_BottomLogical  = instance->GetVolume("BoricAcid_Bottom_LV");
						}

						new G4PVPlacement(nullptr, aTlate, 
								BoricAcid_RubberLogical, "BoricAcid_Rubber_PV", aWhere, false, 0, OverlapCheck);
						new G4PVPlacement(nullptr, aTlate + bottomBoricAcid_Pos,
								BoricAcid_BottomLogical, "BoricAcid_Bottom_PV", aWhere, false, 0, OverlapCheck);
				} break;
        case 11: { // For mimicking configuration for Pilot (Run 7 configuration)
            namespace PilotValues = AmoreDetectorStaticInfo::AMoRE_Pilot_Run7Save;

            leadBoxsize  = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
            leadBoxzsize = (leadBoxZinnerlength + leadBoxThickness) / 2.0;

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("PolyEthylene_FB_LV", false);
            if (logicalVolumeChecker == nullptr) {
                polyEthylene_SideFB_Box =
                    new G4Box("PolyEthylene_FB_Box", PilotValues::polyEthylene_thickness1 / 2.,
                              2. * sideFBScint_boxWidth_half, sideFBScint_boxHeight_half);
                polyEthylene_SideFB_BoxLogical = new G4LogicalVolume(
                    polyEthylene_SideFB_Box, _polyethylene, "PolyEthylene_FB_LV");

                polyEthylene_SideLR1_Box = new G4Box(
                    "PolyEthylene_LR1_Box", 2. * sideLRScint_boxWidth_half,
                    PilotValues::polyEthylene_thickness1 / 2.,
                    sideLRScint_boxHeight_half - PilotValues::polyEthylene_LRCuttingLength1 / 2.);
                polyEthylene_SideLR1_BoxLogical = new G4LogicalVolume(
                    polyEthylene_SideLR1_Box, _polyethylene, "PolyEthylene_LR1_LV");

                polyEthylene_SideLR2_Box = new G4Box(
                    "PolyEthylene_LR2_Box", 2. * sideLRScint_boxWidth_half,
                    PilotValues::polyEthylene_thickness2 / 2.,
                    sideLRScint_boxHeight_half - PilotValues::polyEthylene_LRCuttingLength2 / 2.);
                polyEthylene_SideLR2_BoxLogical = new G4LogicalVolume(
                    polyEthylene_SideLR2_Box, _polyethylene, "PolyEthylene_LR2_LV");

                polyEthylene_Bottom_Box = new G4Box(
                    "PolyEthylene_Bottom_Box", (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    PilotValues::polyEthylene_thickness1 / 2.);
                polyEthylene_Bottom_BoxLogical = new G4LogicalVolume(
                    polyEthylene_Bottom_Box, _polyethylene, "PolyEthylene_Bottom_LV");

                polyEthylene_SideFB_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                polyEthylene_SideLR1_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                polyEthylene_SideLR2_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
                polyEthylene_Bottom_BoxLogical->SetVisAttributes(fI_PE_VisAttr);

                boratedPE_SideFB_Box =
                    new G4Box("BoratedPE_SideFB_Box", PilotValues::boratedPE_thickness / 2.,
                              2. * sideFBScint_boxWidth_half, sideFBScint_boxHeight_half);
                boratedPE_SideFB_BoxLogical = new G4LogicalVolume(
                    boratedPE_SideFB_Box, _BoratedPE_5perCent, "BoratedPE_SideFB_LV");

                boratedPE_SideLR1_Box =
                    new G4Box("BoratedPE_SideLR1_Box", 2. * sideLRScint_boxWidth_half,
                              PilotValues::boratedPE_thickness / 2., sideLRScint_boxHeight_half);
                boratedPE_SideLR1_BoxLogical = new G4LogicalVolume(
                    boratedPE_SideLR1_Box, _BoratedPE_5perCent, "BoratedPE_SideLR1_LV");

                boratedPE_SideLR2_Box =
                    new G4Box("BoratedPE_SideLR2_Box", 2. * sideLRScint_boxWidth_half,
                              PilotValues::boratedPE_thickness / 2., 
															sideLRScint_boxHeight_half- PilotValues::polyEthylene_LRCuttingLength1 / 2.
															+ PilotValues::polyEthylene_LRCuttingLength2 / 2.);
                boratedPE_SideLR2_BoxLogical = new G4LogicalVolume(
                    boratedPE_SideLR2_Box, _BoratedPE_5perCent, "BoratedPE_SideLR2_LV");

                boratedPE_Bottom_Box = new G4Box(
                    "BoratedPE_Bottom_Box", (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    (leadBoxThickness * 2. + leadBoxXYinnerlength) / 2.,
                    PilotValues::boratedPE_thickness / 2.);
                boratedPE_Bottom_BoxLogical = new G4LogicalVolume(
                    boratedPE_Bottom_Box, _BoratedPE_5perCent, "BoratedPE_Bottom_LV");

                boratedPE_SideFB_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
                boratedPE_SideLR1_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
                boratedPE_SideLR2_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);
                boratedPE_Bottom_BoxLogical->SetVisAttributes(fI_BoratedPE_VisAttr);

                detectorRubberB4C_Box = new G4Box(
                    "DetectorRubberB4C_Box", PilotValues::boronCarbideRubber_UnderPb_XYsize / 2.,
                    PilotValues::boronCarbideRubber_UnderPb_XYsize / 2.,
                    PilotValues::boronCarbideRubber_thickness / 2.);
                detectorRubberB4C_BoxLogical = new G4LogicalVolume(
                    detectorRubberB4C_Box, _B4CRubber24perCent, "DetectorRubberB4C_LV");
                detectorRubberB4C_BoxLogical->SetVisAttributes(fI_B4CRubber_VisAttr);
            } else {
                G4LogicalVolumeStore *instance  = G4LogicalVolumeStore::GetInstance();
                boratedPE_SideFB_BoxLogical     = instance->GetVolume("BoratedPE_SideFB_LV");
                boratedPE_SideLR1_BoxLogical     = instance->GetVolume("BoratedPE_SideLR1_LV");
                boratedPE_SideLR2_BoxLogical     = instance->GetVolume("BoratedPE_SideLR2_LV");
                boratedPE_Bottom_BoxLogical     = instance->GetVolume("BoratedPE_Bottom_LV");
                polyEthylene_SideFB_BoxLogical  = instance->GetVolume("PolyEthylene_FB_LV");
                polyEthylene_SideLR1_BoxLogical = instance->GetVolume("PolyEthylene_LR1_LV");
                polyEthylene_SideLR2_BoxLogical = instance->GetVolume("PolyEthylene_LR2_LV");
                polyEthylene_Bottom_BoxLogical  = instance->GetVolume("PolyEthylene_Bottom_LV");
                detectorRubberB4C_BoxLogical    = instance->GetVolume("DetectorRubberB4C_LV");
            }

            FB_Trans = {0, 0, -leadBoxzsize + (sideFBScint_boxHeight_half + scint_reflector_thick)};
            LR_Trans = {0, 0, -leadBoxzsize + (sideLRScint_boxHeight_half + scint_reflector_thick)};
            sideLRScint_YPos =
                std::fmax(2 * (sideFBScint_boxWidth_half + scint_reflector_thick), leadBoxsize);

            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans +
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      PilotValues::boratedPE_thickness / 2.,
                                  0, 0),
                boratedPE_SideFB_BoxLogical, "BoratedPE_SideFB_PV", aWhere, false, 0, OverlapCheck);
            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans -
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      PilotValues::boratedPE_thickness / 2.,
                                  0, 0),
                boratedPE_SideFB_BoxLogical, "BoratedPE_SideFB_PV", aWhere, false, 1, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans +
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      PilotValues::boratedPE_thickness +
                                      PilotValues::polyEthylene_thickness1 / 2.,
                                  0, 0),
                polyEthylene_SideFB_BoxLogical, "PolyEthylene_SideFB_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + FB_Trans -
                    G4ThreeVector(leadBoxsize + (scint_Thick_half + scint_reflector_thick) * 2. +
                                      PilotValues::boratedPE_thickness +
                                      PilotValues::polyEthylene_thickness1 / 2.,
                                  0, 0),
                polyEthylene_SideFB_BoxLogical, "PolyEthylene_SideFB_PV", aWhere, false, 1, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans +
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      PilotValues::boratedPE_thickness / 2.,
                                  0) +
									  G4ThreeVector(0, 0, -PilotValues::polyEthylene_LRCuttingLength1 / 2.
											                  +PilotValues::polyEthylene_LRCuttingLength2 / 2.)	,
                boratedPE_SideLR2_BoxLogical, "BoratedPE_SideLR2_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans -
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      PilotValues::boratedPE_thickness / 2.,
                                  0),
                boratedPE_SideLR1_BoxLogical, "BoratedPE_SideLR1_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans -
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      PilotValues::boratedPE_thickness -
                                      PilotValues::polyEthylene_thickness2 / 2.,
                                  0) +
                    G4ThreeVector(0, 0, -PilotValues::polyEthylene_LRCuttingLength2 / 2.),
                polyEthylene_SideLR2_BoxLogical, "PolyEthylene_SideLR2_PV", aWhere, false, 0, OverlapCheck);

            new G4PVPlacement(
                nullptr,
                aTlate + LR_Trans +
                    G4ThreeVector(0.,
                                  -sideLRScint_YPos -
                                      (scint_Thick_half + scint_reflector_thick) * 2. -
                                      PilotValues::boratedPE_thickness -
                                      PilotValues::polyEthylene_thickness1 / 2.,
                                  0) +
                    G4ThreeVector(0, 0, -PilotValues::polyEthylene_LRCuttingLength1 / 2.),
                polyEthylene_SideLR1_BoxLogical, "PolyEthylene_SideLR1_PV", aWhere, false, 0, OverlapCheck);

						/*
            new G4PVPlacement(nullptr,
                              aTlate - G4ThreeVector(0, 0,
                                                     (leadBoxZinnerlength + leadBoxThickness) / 2. +
                                                         PilotValues::boratedPE_thickness / 2. +
                                                         PilotValues::boratedPE_BottomGap),
                              boratedPE_Bottom_BoxLogical, "BoratedPE_Bottom_PV", aWhere, false, 0, OverlapCheck);
*/
            new G4PVPlacement(nullptr,
                              aTlate - G4ThreeVector(0, 0,
                                                     (leadBoxZinnerlength + leadBoxThickness) / 2. +
																												 (scint_muffler_Thick_half+scint_reflector_thick)*2 +
                                                         //PilotValues::boratedPE_thickness +
                                                         PilotValues::boratedPE_BottomGap +
                                                         PilotValues::polyEthylene_thickness1 / 2.),
                              polyEthylene_Bottom_BoxLogical, "PolyEthylene_Bottom_PV", aWhere,
                              false, 0, OverlapCheck);

            new G4PVPlacement(nullptr,
                              aTlate - G4ThreeVector(0, 0,
                                                     (leadBoxZinnerlength + leadBoxThickness) / 2. +
																												 (scint_muffler_Thick_half+scint_reflector_thick)*2 +
                                                         PilotValues::boratedPE_thickness /2. +
                                                         //PilotValues::boratedPE_thickness * 1.5 +
                                                         PilotValues::boratedPE_BottomGap +
                                                         PilotValues::polyEthylene_thickness1),
                              boratedPE_Bottom_BoxLogical, "BoratedPE_Bottom_PV", aWhere, false, 0, OverlapCheck);

						/*
            new G4PVPlacement(
                nullptr,
                aTlate + G4ThreeVector(0, 0,
                                       leadBoxzsize + sst_zsize_half * 2. +
                                           PilotValues::boronCarbideRubber_GapFromSSOVC),
                detectorRubberB4C_BoxLogical, "B4CRubber_Detector_PV", aWhere, false, 0, OverlapCheck);
								*/
        } break;
        case 12: {
            namespace PilotValues  = AmoreDetectorStaticInfo::AMoRE_Pilot_Run7Save;
            topScint_boxThick_half = scint_Thick_half;

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("TopPbAboveRubberB4C_LV", false);
            if (logicalVolumeChecker == nullptr) {
                topPbAboveRubberB4C_Box =
                    new G4Box("TopPbAboveRubberB4C_Box", pbtopbox_size, pbtopbox_size,
                              PilotValues::boronCarbideRubber_thickness / 2.);
                polyEthylene_Top_Box =
                    new G4Box("PolyEthylene_Top_Box", pbtopbox_size, pbtopbox_size,
                              PilotValues::polyEthylene_thickness1 / 2.);

                topPbAboveRubberB4C_BoxLogical = new G4LogicalVolume(
                    topPbAboveRubberB4C_Box, _B4CRubber24perCent, "TopPbAboveRubberB4C_LV");
                topPbAboveRubberB4C_BoxLogical->SetVisAttributes(fI_B4CRubber_VisAttr);

                polyEthylene_Top_BoxLogical =
                    new G4LogicalVolume(polyEthylene_Top_Box, _polyethylene, "PolyEthylene_Top_LV");
                polyEthylene_Top_BoxLogical->SetVisAttributes(fI_PE_VisAttr);
            } else {
                G4LogicalVolumeStore *instance = G4LogicalVolumeStore::GetInstance();
                topPbAboveRubberB4C_BoxLogical = instance->GetVolume("TopPbAboveRubberB4C_LV");
                polyEthylene_Top_BoxLogical    = instance->GetVolume("PolyEthylene_Top_LV");
            }

						/*
            new G4PVPlacement(
                nullptr,
                aTlate + G4ThreeVector(0, 0,
                                       pbtopbox_zsize +
                                           (topScint_boxThick_half + scint_reflector_thick) * 2. +
                                           PilotValues::boronCarbideRubber_thickness / 2.),
                topPbAboveRubberB4C_BoxLogical, "B4CRubber_TopPbAbove_PV", aWhere, false, 0, OverlapCheck);
								*/
            new G4PVPlacement(
                nullptr,
                aTlate + G4ThreeVector(0, 0,
                                       pbtopbox_zsize +
                                           (topScint_boxThick_half + scint_reflector_thick) * 2. +
                                           //PilotValues::boronCarbideRubber_thickness +
                                           PilotValues::polyEthylene_thickness1 / 2.),
                polyEthylene_Top_BoxLogical, "PolyEthylene_Top_PV", aWhere, false, 0, OverlapCheck);
        } break;
        case 13: {
            namespace PilotValues = AmoreDetectorStaticInfo::AMoRE_Pilot_Run7Save;

            totalThickness_BoricAcidPanel =
                PilotValues::thickness_BoricAcid + PilotValues::thickness_BoricAcidHousing * 2.;
            leadBoxsize   = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
            leadBoxzsize  = (leadBoxZinnerlength + leadBoxThickness) / 2.0;
            height_LR     = leadBoxzsize - leadBoxThickness / 2.;
            width_LR      = leadBoxsize - leadBoxThickness;
            height_FB     = leadBoxzsize - leadBoxThickness / 2.;
            width_FB      = leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel;
            length_Bottom = width_LR - totalThickness_BoricAcidPanel;
            width_Bottom  = width_FB - totalThickness_BoricAcidPanel;

            logicalVolumeChecker =
                G4LogicalVolumeStore::GetInstance()->GetVolume("BoricAcid_LR_Mother_LV", false);
            if (logicalVolumeChecker == nullptr) {
                sideLRBoricAcid_HousingBox =
                    new G4Box("BoricAcid_LR_Mother_Box", width_LR,
                              totalThickness_BoricAcidPanel / 2., height_LR);
                sideFBBoricAcid_HousingBox =
                    new G4Box("BoricAcid_FB_Mother_Box", totalThickness_BoricAcidPanel / 2.,
                              width_FB, height_FB);
                bottomBoricAcid_HousingBox =
                    new G4Box("BoricAcid_Bottom_Mother_Box", length_Bottom, width_Bottom,
                              totalThickness_BoricAcidPanel / 2.);

                sideLRBoricAcid_Block = new G4Box(
                    "BoricAcid_LR_Block", width_LR - PilotValues::thickness_BoricAcidSpacing,
                    PilotValues::thickness_BoricAcid / 2.,
                    height_LR - PilotValues::thickness_BoricAcidSpacing);
                sideFBBoricAcid_Block =
                    new G4Box("BoricAcid_FBBlock_Box", PilotValues::thickness_BoricAcid / 2.,
                              width_FB - PilotValues::thickness_BoricAcidSpacing,
                              height_FB - PilotValues::thickness_BoricAcidSpacing);
                bottomBoricAcid_Block =
                    new G4Box("BoricAcid_BottomBlock_Box",
                              length_Bottom - PilotValues::thickness_BoricAcidSpacing,
                              width_Bottom - PilotValues::thickness_BoricAcidSpacing,
                              PilotValues::thickness_BoricAcid / 2.);

                sideLRBoricAcid_MotherLogical = new G4LogicalVolume(
                    sideLRBoricAcid_HousingBox, _polyethylene, "BoricAcid_LR_Mother_LV");
                sideFBBoricAcid_MotherLogical = new G4LogicalVolume(
                    sideFBBoricAcid_HousingBox, _polyethylene, "BoricAcid_FB_Mother_LV");
                bottomBoricAcid_MotherLogical = new G4LogicalVolume(
                    bottomBoricAcid_HousingBox, _polyethylene, "BoricAcid_Bottom_Mother_LV");

                sideLRBoricAcid_BlockLogical = new G4LogicalVolume(
                    sideLRBoricAcid_Block, _BoricAcidPowder, "BoricAcid_LR_Block_LV");
                sideFBBoricAcid_BlockLogical = new G4LogicalVolume(
                    sideFBBoricAcid_Block, _BoricAcidPowder, "BoricAcid_FB_Block_LV");
                bottomBoricAcid_BlockLogical = new G4LogicalVolume(
                    bottomBoricAcid_Block, _BoricAcidPowder, "BoricAcid_Bottom_Block_LV");

                new G4PVPlacement(nullptr, {}, sideLRBoricAcid_BlockLogical,
                                  "BoricAcid_LR_Block_PV", sideLRBoricAcid_MotherLogical, false, 0, OverlapCheck);
                new G4PVPlacement(nullptr, {}, sideFBBoricAcid_BlockLogical,
                                  "BoricAcid_FB_Block_PV", sideFBBoricAcid_MotherLogical, false, 0, OverlapCheck);
                new G4PVPlacement(nullptr, {}, bottomBoricAcid_BlockLogical,
                                  "BoricAcid_Bottom_Block_PV", bottomBoricAcid_MotherLogical, false,
                                  0, OverlapCheck);

                bottomBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
                bottomBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);

                sideLRBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
                sideLRBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);

                sideFBBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
                sideFBBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);
            } else {
                G4LogicalVolumeStore *instance = G4LogicalVolumeStore::GetInstance();
                sideLRBoricAcid_MotherLogical  = instance->GetVolume("BoricAcid_LR_Mother_LV");
                sideFBBoricAcid_MotherLogical  = instance->GetVolume("BoricAcid_FB_Mother_LV");
                bottomBoricAcid_MotherLogical  = instance->GetVolume("BoricAcid_Bottom_Mother_LV");
            }

            sideLRBoricAcid_Pos[0] = {
                0, leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel / 2.,
                -leadBoxzsize + leadBoxThickness + height_LR};
            sideLRBoricAcid_Pos[1] = {
                0, -leadBoxsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2.,
                -leadBoxzsize + leadBoxThickness + height_LR};
            sideFBBoricAcid_Pos[0] = {leadBoxsize - leadBoxThickness -
                                          totalThickness_BoricAcidPanel / 2.,
                                      0, -leadBoxzsize + leadBoxThickness + height_FB};
            sideFBBoricAcid_Pos[1] = {-leadBoxsize + leadBoxThickness +
                                          totalThickness_BoricAcidPanel / 2.,
                                      0, -leadBoxzsize + leadBoxThickness + height_FB};

            bottomBoricAcid_Pos = {
                0, 0, -leadBoxzsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2.};

            new G4PVPlacement(nullptr, aTlate + sideLRBoricAcid_Pos[0],
                              sideLRBoricAcid_MotherLogical, "BoricAcid_LR_Mother_PV", aWhere,
                              false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, aTlate + sideLRBoricAcid_Pos[1],
                              sideLRBoricAcid_MotherLogical, "BoricAcid_LR_Mother_PV", aWhere,
                              false, 1, OverlapCheck);

            new G4PVPlacement(nullptr, aTlate + sideFBBoricAcid_Pos[0],
                              sideFBBoricAcid_MotherLogical, "BoricAcid_FB_Mother_PV", aWhere,
                              false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, aTlate + sideFBBoricAcid_Pos[1],
                              sideFBBoricAcid_MotherLogical, "BoricAcid_FB_Mother_PV", aWhere,
                              false, 1, OverlapCheck);
            new G4PVPlacement(nullptr, aTlate + bottomBoricAcid_Pos, bottomBoricAcid_MotherLogical,
                              "BoricAcid_Bottom_Mother_PV", aWhere, false, 0, OverlapCheck);
        } break;
        default:
            G4Exception(__PRETTY_FUNCTION__, "AMORE_I_TYPEERROR",
                        G4ExceptionSeverity::FatalErrorInArgument,
                        "Wrong type for neutron shielding configuration");
            break;
    }
}

void AmoreDetectorConstruction::ConstructAMoRE_I() {
    G4int flagInvisible       = 0;
    G4bool neutronShieldBuilt = false;
    using namespace AmoreDetectorStaticInfo::AMoRE_I;
    using namespace AmoreDetectorStaticInfo::ColorTable;
    namespace PilotValues = AmoreDetectorStaticInfo::AMoRE_Pilot_Run7Save;

    ///////////////////////////////////////////////////////
    // -- Make visualization attributes

    fI_Invisible_VisAttr        = new G4VisAttributes(false);
    fI_CopperDefault_VisAttr    = new G4VisAttributes(brown);
    fI_LeadDefault_VisAttr      = new G4VisAttributes(grey);
    fI_GoldDefault_VisAttr      = new G4VisAttributes(G4Colour(1, 0.87, 0));
    fI_CopperFrame_VisAttr      = new G4VisAttributes(brownl);
    fI_OpticalFrame_VisAttr     = new G4VisAttributes(greenl);
    fI_GeWafer_VisAttr          = new G4VisAttributes(bluel);
    fI_LeadShield_VisAttr       = new G4VisAttributes(greyl);
    fI_Reflector_VisAttr        = new G4VisAttributes(G4Colour(0.75, 0.0, 0.75, 0.3));
    fI_Crystal_VisAttr          = new G4VisAttributes(yellow);
    fI_Bolt_VisAttr             = new G4VisAttributes(grey);
    fI_FT2850_VisAttr           = new G4VisAttributes(red);
    fI_Vikuiti_VisAttr          = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75, 0.3));
    fI_PEEK_VisAttr             = new G4VisAttributes(yellow);
    fI_Clamp_VisAttr            = new G4VisAttributes(white);
    fI_ClampBolt_VisAttr        = new G4VisAttributes(brown);
    fI_PCB_VisAttr              = new G4VisAttributes(white);
    fI_Epoxy_VisAttr            = new G4VisAttributes(yellow);
    fI_BoratedPE_VisAttr        = new G4VisAttributes(redl);
    fI_PE_VisAttr               = new G4VisAttributes(lgreenl);
    fI_BoricAcid_VisAttr        = new G4VisAttributes(whitel);
    fI_BoricAcidHousing_VisAttr = new G4VisAttributes(G4Colour(1, 1, 0.9, 0.4));
    fI_B4CRubber_VisAttr        = new G4VisAttributes(greyl);
		fI_BrassBolt_VisAttr        = new G4VisAttributes(greyl);

    fI_CopperDefault_VisAttr->SetForceSolid(true);
    fI_LeadDefault_VisAttr->SetForceSolid(true);
    fI_GoldDefault_VisAttr->SetForceSolid(true);
    fI_CopperFrame_VisAttr->SetForceSolid(true);
    fI_OpticalFrame_VisAttr->SetForceSolid(false);
    fI_GeWafer_VisAttr->SetForceSolid(true);
    fI_LeadShield_VisAttr->SetForceSolid(true);
    fI_Reflector_VisAttr->SetForceSolid(true);
    fI_Crystal_VisAttr->SetForceSolid(true);
    fI_Bolt_VisAttr->SetForceSolid(true);
    fI_FT2850_VisAttr->SetForceSolid(true);
    fI_Vikuiti_VisAttr->SetForceSolid(true);
    fI_PEEK_VisAttr->SetForceSolid(true);
    fI_Clamp_VisAttr->SetForceSolid(true);
    fI_ClampBolt_VisAttr->SetForceSolid(true);
    fI_PCB_VisAttr->SetForceSolid(true);
    fI_Epoxy_VisAttr->SetForceSolid(true);
    fI_BoratedPE_VisAttr->SetForceSolid(true);
    fI_PE_VisAttr->SetForceSolid(true);
    fI_BoricAcid_VisAttr->SetForceSolid(true);
    fI_BoricAcidHousing_VisAttr->SetForceSolid(true);
    fI_B4CRubber_VisAttr->SetForceSolid(true);
		fI_BrassBolt_VisAttr->SetForceSolid(true);

    ///////////////////////////////////////////////
    // Primitive values
    ///////////////////////////////////////////////

    // H Beam points
    std::vector<G4TwoVector> large_points; // Large beam points.
    large_points.push_back(G4TwoVector(100.0, 100.0));
    large_points.push_back(G4TwoVector(100.0, -100.0));
    large_points.push_back(G4TwoVector(88.0, -100.0));
    large_points.push_back(G4TwoVector(88.0, -4.0));
    large_points.push_back(G4TwoVector(-88.0, -4.0));
    large_points.push_back(G4TwoVector(-88.0, -100.0));
    large_points.push_back(G4TwoVector(-100.0, -100.0));
    large_points.push_back(G4TwoVector(-100.0, 100.0));
    large_points.push_back(G4TwoVector(-88.0, 100.0));
    large_points.push_back(G4TwoVector(-88.0, 4.0));
    large_points.push_back(G4TwoVector(88.0, 4.0));
    large_points.push_back(G4TwoVector(88.0, 100.0));

    std::vector<G4TwoVector> small_points; // Small beam points.
    small_points.push_back(G4TwoVector(75.0, 75.0));
    small_points.push_back(G4TwoVector(75.0, -75.0));
    small_points.push_back(G4TwoVector(65.0, -75.0));
    small_points.push_back(G4TwoVector(65.0, -3.5));
    small_points.push_back(G4TwoVector(-65.0, -3.5));
    small_points.push_back(G4TwoVector(-65.0, -75.0));
    small_points.push_back(G4TwoVector(-75.0, -75.0));
    small_points.push_back(G4TwoVector(-75.0, 75.0));
    small_points.push_back(G4TwoVector(-65.0, 75.0));
    small_points.push_back(G4TwoVector(-65.0, 3.5));
    small_points.push_back(G4TwoVector(65.0, 3.5));
    small_points.push_back(G4TwoVector(65.0, 75.0));

    ///////////////////////////////////////////////
    // Composite values
    ///////////////////////////////////////////////

    //  Mounting structure of H-beam
    G4double posZ_Bplate1    = Hbeam1_length + Bplate_thick_half;
    G4double posZ_Bplate2    = Hbeam2_length + Bplate_thick_half;
    G4double posMountOrigin1 = -workboxH + Hbeam1_length + 2. * Bplate_thick_half;
    G4double posMountOrigin2 = -workboxH + Hbeam2_length + 2. * Bplate_thick_half - 1 * mm;
    G4ThreeVector posBplate1 = G4ThreeVector(0., 0., posZ_Bplate1);
    G4ThreeVector posBplate2 = G4ThreeVector(0., 0., posZ_Bplate2 - 1 * mm);

    /// Vertical H beam Rotatrion Matrix
    G4RotationMatrix *VeBeamRotation = new G4RotationMatrix();
    VeBeamRotation->rotateX(0. * deg);
    VeBeamRotation->rotateY(0. * deg);

    /// Horizontal H beam Rotation Matrix X
    G4RotationMatrix *HoBXRotation = new G4RotationMatrix();
    HoBXRotation->rotateX(90. * deg);
    HoBXRotation->rotateZ(90. * deg);
    HoBXRotation->rotateX(90. * deg);

    /// Horizontal H beam Rotation Matrix Y
    G4RotationMatrix *HoBYRotation = new G4RotationMatrix();
    HoBYRotation->rotateY(90. * deg);
    HoBYRotation->rotateX(90. * deg);

    // PMT
    G4double muonVetoPMTEnvelopeSizeZ =
        muonVetoPMTSizeZ + muonVetoPMTWindowThickness + muonVetoPMTGreaseThickness;

    // Scintillator

    G4RotationMatrix *scintMother_alignMtx_LRSide1 = new G4RotationMatrix();
    G4RotationMatrix *scintMother_alignMtx_LRSide2 = new G4RotationMatrix();

    scintMother_alignMtx_LRSide1->rotateX(-90 * deg);

    scintMother_alignMtx_LRSide2->rotateX(-90 * deg);
    scintMother_alignMtx_LRSide2->rotateY(-180 * deg);

    //  Pb Top Box
    G4double CenterPbTZ = Hbeam1_length + posMountOrigin1 + pbtopbox_zsize + 1 + lifting_PbTopBox;

    //  Pb  Box
    G4double leadBoxsize      = (leadBoxXYinnerlength + 2 * leadBoxThickness) / 2.0;
    G4double leadBoxzsize     = (leadBoxZinnerlength + leadBoxThickness) / 2.0;
    G4double verticalGapToTop = -pbtopbox_zsize + verticalGapToMuVeto;
    G4double PbboxTopZ        = CenterPbTZ;
    G4double CenterPbZ        = PbboxTopZ - leadBoxzsize - verticalGapToTop;
    G4ThreeVector PosPbBox    = G4ThreeVector(0., 0., CenterPbZ);
    G4ThreeVector PosInnerDet = PosPbBox;

    // SSOVC
    G4ThreeVector SSTopShift = G4ThreeVector(0., 0., ss_height_half + ssb_zsize + sst_zsize_half);

    // Layer4_ 50K-SHIELD Copper
    G4double cu4_total_height_half = cu4_inner_height + cu4b_zsize + cu4t_zsize;
    G4ThreeVector PosCu4Outer =
        G4ThreeVector(0., 0., ss_height_half - cu4_total_height_half - cu4_gap_fromTop);
    G4ThreeVector PosCu4Inner = G4ThreeVector(0., 0., cu4b_zsize - cu4t_zsize);

    // Layer3_Cu IVC
    G4double cu3_total_height_half = cu3_inner_height + cu3b_zsize + cu3t_zsize;
    G4ThreeVector PosCu3Outer(0., 0., cu4_inner_height - cu3_total_height_half - cu3_gap_fromTop);
    G4ThreeVector PosCu3Inner(0, 0, cu3b_zsize - cu3t_zsize);

    // Layer2_Cu SHIELD-STILL
    G4double cu2_total_height_half = cu2_inner_height + cu2b_zsize + cu2t_zsize;
    G4ThreeVector PosCu2Outer =
        G4ThreeVector(0., 0., cu3_inner_height - cu2_total_height_half - cu2_gap_fromTop);
    G4ThreeVector PosCu2Inner = G4ThreeVector(0., 0., cu2b_zsize - cu2t_zsize);

    // Layer1_Cu 50mK-SHIELD
    G4double cu1_total_height_half = cu1_inner_height + cu1b_zsize + cu1t_zsize;
    G4ThreeVector PosCu1Outer(0, 0, cu2_inner_height - cu1_total_height_half - cu1_gap_fromTop);
    G4ThreeVector PosCu1Inner(0, 0, cu1b_zsize - cu1t_zsize);

    // Cu Plate: Mixing Chamber SHIELD
    G4double CenterCuMCPZ  = cu1_inner_height - copperMixingChamberPlatThickness - plateGap;
    G4ThreeVector PosCuMCP = G4ThreeVector(0., 0., CenterCuMCPZ);

    // Plate1_Cu
    G4double CenterCuP1Z =
        CenterCuMCPZ - copperMixingChamberPlatThickness - copperPlate1ZSize_half - plateGap2;
    G4ThreeVector PosCuP1 = G4ThreeVector(0., 0., CenterCuP1Z);

    // Plate2_Pb Shield
    G4double CenterPbP2Z  = CenterCuP1Z - copperPlate1ZSize_half - leadPlate2ZSize_half;
    G4ThreeVector PosPbP2 = G4ThreeVector(0., 0., CenterPbP2Z);

    // Plate3_Cu
    G4double CenterCuP3Z  = CenterPbP2Z - leadPlate2ZSize_half - cup3_zsize;
    G4ThreeVector PosCuP3 = G4ThreeVector(0., 0., CenterCuP3Z);

    // Plate4_Cu
    G4double CenterCuP4Z  = CenterCuP3Z - cup3_zsize - copperPlate1ZSize_half - plateGap2;
    G4ThreeVector PosCuP4 = G4ThreeVector(0., 0., CenterCuP4Z);

    // Plate5_Pb Shield
    G4double CenterPbP5Z  = CenterCuP4Z - copperPlate1ZSize_half - leadPlate2ZSize_half;
    G4ThreeVector PosPbP5 = G4ThreeVector(0., 0., CenterPbP5Z);

    // Plate6_Cu
    G4double CenterCuP6Z  = CenterPbP5Z - leadPlate2ZSize_half - cup3_zsize;
    G4ThreeVector PosCuP6 = G4ThreeVector(0., 0., CenterCuP6Z);

    // Boric acid panel
    G4double totalThickness_BoricAcidPanel = thickness_BoricAcid + thickness_BoricAcidHousing * 2.;

    G4double height_LR = leadBoxzsize - leadBoxThickness / 2.;
    G4double width_LR  = leadBoxsize - leadBoxThickness;

    G4double height_FB = leadBoxzsize - leadBoxThickness / 2.;
    G4double width_FB  = leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel;

    G4double length_Top = width_LR;
    G4double width_Top  = width_FB;

    G4ThreeVector sideLRBoricAcid_Pos[2];
    sideLRBoricAcid_Pos[0] =
        G4ThreeVector(0, leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel / 2.,
                      -leadBoxzsize + leadBoxThickness + height_LR);
    sideLRBoricAcid_Pos[1] =
        G4ThreeVector(0, -leadBoxsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2.,
                      -leadBoxzsize + leadBoxThickness + height_LR);

    G4ThreeVector sideFBBoricAcid_Pos[2];
    sideFBBoricAcid_Pos[0] =
        G4ThreeVector(leadBoxsize - leadBoxThickness - totalThickness_BoricAcidPanel / 2., 0,
                      -leadBoxzsize + leadBoxThickness + height_FB);
    sideFBBoricAcid_Pos[1] =
        G4ThreeVector(-leadBoxsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2., 0,
                      -leadBoxzsize + leadBoxThickness + height_FB);

    G4ThreeVector bottomBoricAcid_Pos =
        G4ThreeVector(0, 0, -leadBoxzsize + leadBoxThickness + totalThickness_BoricAcidPanel / 2.);
    G4ThreeVector topBoricAcid_Pos = G4ThreeVector(
        0, 0, +leadBoxzsize + sst_zsize_half * 2. + totalThickness_BoricAcidPanel / 2.);

    // G10 material holders between SS and 50K-SHIELD Copper
    G4double Cuholder_r = Cuholder_R - 1.9038 * mm;
    G4double CuholderH  = plateGap; // height

    G4double xG10[6][6];
    G4double yG10[6][6];
    G4double CenterPlateGap4VolZ[5] = {
        cu1_inner_height - CuholderH / 2, cu2_inner_height - CuholderH / 2,
        cu3_inner_height - CuholderH / 2, cu4_inner_height - CuholderH / 2,
        ss_height_half - CuholderH / 2};

    G4double PlateGapR[5] = {cu1_radius - cu1_thick, cu2_radius - cu2_thick, cu3_radius - cu3_thick,
                             cu4_radius - cu4_thick, ss_radius - ss_thick};
    ///////////////////////////////////////////////
    // Variables
    ///////////////////////////////////////////////

    // * Define a world physical volume.
    // * Set the class variable "world_phys" to point to the world phys volume.
    G4Box *boxHall;
    G4LogicalVolume *logiHall;
    G4VisAttributes *logiHallVis;
    G4VPhysicalVolume *physHall;

    /////////////////////////////////////////////////////////////////
    // Cavity : Air Room 
    /////////////////////////////////////////////////////////////////
    G4VSolid *RockCavitySolid;
    G4LogicalVolume *logiRockCavity;
    G4VisAttributes *logiRockCavityVis;

		/////////////////////////////////////////
		// ROCK
		/////////////////////////////////////////
		G4Sphere *RockSolid;
		G4LogicalVolume *logiRock;
		G4VisAttributes *logiRockVis;

    /////////////////////////////////////////
    // Rock Disk solid
    /////////////////////////////////////////
    G4Tubs *RockDiskSolid;
    G4LogicalVolume *logiRockDisk;
    G4VisAttributes *logiRockDiskVis;

    /////////////////////////////////////////
    // Work Area  
    /////////////////////////////////////////
    G4Box *WorkArea;
    G4LogicalVolume *logiWorkArea;
    G4VisAttributes *logiWorkAreaVis;

    ////////////////////////////////////////////////////////////
    //  Mounting structure of H-beam properties & variables
    ////////////////////////////////////////////////////////////
    G4Box *Bplate;
    G4UnionSolid *UnionHbeam1;
    G4UnionSolid *UnionHbeam2;
    G4ExtrudedSolid *Hbeam1;
    G4ExtrudedSolid *Hbeam2;
    G4ExtrudedSolid *Ho_beamX;
    G4ExtrudedSolid *Ho_beamXLar;
    G4ExtrudedSolid *Ho_beamY;
    G4ExtrudedSolid *Ho_beamYLar;

    G4LogicalVolume *logicHbeam1;
    G4LogicalVolume *logicHbeam2;
    G4LogicalVolume *logicHo_beamX;
    G4LogicalVolume *logicHo_beamXLar;
    G4LogicalVolume *logicHo_beamY;
    G4LogicalVolume *logicHo_beamYLar;

    G4VisAttributes *Hb1VisAtt;
    G4VisAttributes *Hb2VisAtt;
    G4VisAttributes *HoBXVisAtt;
    G4VisAttributes *HoBXLVisAtt;
    G4VisAttributes *HoBYLVisAtt;
    G4VisAttributes *HoBYVisAtt;

    /////////////////////////////////////////
		// CHECK!!!!!!!!!
    // G4double PbTBoxVSize = 175.45 * cm; //  114.05 cm (center to top of Stainless Steel cover)
    //   + 46.4 cm (gap)
    //   + 15.0 cm (Lead top thickmess)

    //  Pb Top Box
    G4Box *TopPbBox;
    G4Box *TopPbBlock;
    G4LogicalVolume *logiTopPbBox;
    G4VPhysicalVolume *physTopPbBox1 = nullptr;
    G4LogicalVolume *logiTopPbBlock;

    //  Pb  Box
    G4Box *PbBoxOut;
    G4Box *PbBoxIn;
    G4Box *PbBlockOut;
    G4Box *PbBlockIn;
    G4SubtractionSolid *PbBox;
    G4SubtractionSolid *PbBlock;

    G4LogicalVolume *logiPbBox;
    G4LogicalVolume *logiPbBlock;
    G4VisAttributes *logiPbBoxVis;

    ////////////////
    // Inner detector
    ////////////////

    G4Box *InnerDetectorBox;
    G4Box *SSTopBoxForInnerDet;
    G4VSolid *InnerDetector;
    G4LogicalVolume *logiInnerDetector;

    G4ThreeVector PosSS;

    G4Tubs *SSCylinder;
    G4Box *SSTop;
    G4VSolid *SSOVC0;
    G4Tubs *SSOVCInnerSolid;
    G4LogicalVolume *logiSSOVCOuter;
    G4LogicalVolume *logiSSOVCInner;
    G4VisAttributes *logiSSOVCVis;

    /////////////////////////////////////////////////////////////////
    // Layer4_ 50K-SHIELD Copper

    G4Tubs *Cu4OuterCylinder;
    G4Tubs *Cu4InnerCylinder;
    G4LogicalVolume *logiCu4Outer;
    G4LogicalVolume *logiCu4Inner;
    G4VisAttributes *logiCu4Vis;

    ///////////////////////////////////////////////////////
    // Layer3_Cu IVC

    G4Tubs *Cu3OuterCylinder;
    G4Tubs *Cu3InnerCylinder;
    G4LogicalVolume *logiCu3Outer;
    G4LogicalVolume *logiCu3Inner;
    G4VisAttributes *logiCu3Vis;

    ///////////////////////////////////////////////////////
    // Layer2_Cu SHIELD-STILL

    G4Tubs *Cu2OuterCylinder;
    G4Tubs *Cu2InnerCylinder;
    G4LogicalVolume *logiCu2Outer;
    G4LogicalVolume *logiCu2Inner;
    G4VisAttributes *logiCu2Vis;

    ///////////////////////////////////////////////////////
    // Layer1_Cu 50mK-SHIELD

    G4Tubs *Cu1OuterCylinder;
    G4Tubs *Cu1InnerCylinder;
    G4LogicalVolume *logiCu1Outer;
    G4LogicalVolume *logiCu1Inner;
    G4VisAttributes *logiCu1Vis;

    ///////////////////////////////////////////////////////
    // Cu Plate: Mixing Chamber SHIELD

    G4Tubs *CuMCPlate;
    G4LogicalVolume *logiCuMCP;
    G4VisAttributes *logiCuMCPVis;

    ///////////////////////////////////////////////////////
    // Plate1_Cu

    G4Tubs *CuPlate1;
    G4LogicalVolume *logiCuP1;
    G4VisAttributes *logiCuP1Vis;

    ///////////////////////////////////////////////////////
    // Plate2_Pb Shield

    G4Tubs *PbPlate2;
    G4LogicalVolume *logiPbP2;
    G4VisAttributes *logiPbP2Vis;

    ///////////////////////////////////////////////////////
    // Plate3_Cu

    G4Tubs *CuPlate3;
    G4LogicalVolume *logiCuP3;
    G4VisAttributes *logiCuP3Vis;

    ///////////////////////////////////////////////////////
    // G10
    ///////////////////////////////////////////////////////
    G4ThreeVector PosPlateGap4Vol[5];
    G4LogicalVolume *PlateGap4VolLogi[5];

    G4ThreeVector PosG10holder;

    G4Box *toyRubberB4C_Box;
    G4LogicalVolume *toyRubberB4C_BoxLogical;

    ///////////////////////////////////////////////////////
    // Polyethylene geometry
    ///////////////////////////////////////////////////////
    G4Box *toyPE_Box;
    G4LogicalVolume *toyPE_BoxLogical;

    G4double toyPE_thickness = 30 * cm;

    ///////////////////////////////////////////////////////
    // Boric acid panel geometry
    ///////////////////////////////////////////////////////
    G4Box *topBoricAcid_HousingBox;
    G4Box *topBoricAcid_Block;

    G4LogicalVolume *topBoricAcid_MotherLogical;
    G4LogicalVolume *topBoricAcid_BlockLogical;

    // Superconducting magnetic shielding
    G4LogicalVolume *superMagnetShield_LV;
    G4ThreeVector superMagnetTlate;

    // Detector array
    G4LogicalVolume *detectorArray_LV;
    G4ThreeVector detectorArrayTlate;

    ////////////////////////////////////////////////////////////////
    // Physics Hall  30m x 30m x 30m
    ////////////////////////////////////////////////////////////////

    boxHall  = new G4Box("hallbox", bounding_size, bounding_size, bounding_size);
    logiHall = new G4LogicalVolume(boxHall, _air, "logiHall");
    physHall = new G4PVPlacement(nullptr, {}, // translation
                                 logiHall,    // associated logical vol
                                 "physHall",  // name
                                 NULL,        // parent
                                 false,       // no "Many"
                                 0);          // copy number

    // the class variable "world_phys" must be set to the world phys volume.
    world_phys = physHall;

    /////////////////////////////////////////////////////////////////
    // Cavity and rock
		// For neutron mode : Spher shape cavity , r= 3.5 m + 5 mm
		// For normal mode  : Half Tube shape cavity and Tube shape bottom rock, r = 5 m (cavity), 8 m (rock) 
    ////////////////////////////////////////////////////////////////
    G4RotationMatrix *tube_alignMtx = new G4RotationMatrix();
    tube_alignMtx->rotateX(-90 * deg);
    if (fNeutronMode) 
		{
        //RockSolid->SetDeltaThetaAngle(180. * deg);
        RockCavitySolid =
            new G4Sphere("Inner Cavity", 0, NEUT_CavityRadius, 0, 360 * deg, 0, 180 * deg);
    } else 
		{
        RockCavitySolid =
            new G4Tubs("Inner tunnel", 0, RockCavityRadius, RockCavityRadius, 0, 180 * deg);
    }
    logiRockCavity    = new G4LogicalVolume(RockCavitySolid, _air, "Rock_Cavity_LV");
    logiRockCavityVis = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));

		RockSolid = new G4Sphere("Rock_Solid", 0, RockRadius, 0, 360 * deg, 0, 90 * deg);
		logiRock  = new G4LogicalVolume(RockSolid, _rock, "Rock_LV");

    RockDiskSolid = new G4Tubs("Rock Disk", 0, RockRadius, RockDiskThick_half, 0 * deg, 360 * deg);
    logiRockDisk  = new G4LogicalVolume(RockDiskSolid, _rock, "Rock_Disk_LV");

    ///////////////////////////////////////////////////////
    // Work Area   3.7 m x 3.6 m x 4.5 m
    ///////////////////////////////////////////////////////
    WorkArea = new G4Box("WorkArea", workboxX, workboxY, workboxH);

    G4RotationMatrix *WA_alignMtx = new G4RotationMatrix();
    WA_alignMtx->rotateX(90 * deg);
    logiWorkArea = new G4LogicalVolume(WorkArea, _air, "logiWorkArea");

    ///////////////////////////////////////////////////////
    // Layer6_Pb Top Box 1500x1500x150 mm 
    TopPbBox   = new G4Box("TopPbBox", pbtopbox_size, pbtopbox_size / 2, pbtopbox_zsize);
    TopPbBlock = new G4Box("TopPbBlock", pbtopbox_size - pbtopbox_housing_thickness,
                           (pbtopbox_size / 2.) - pbtopbox_housing_thickness,
                           pbtopbox_zsize - pbtopbox_housing_thickness);

    logiTopPbBox = new G4LogicalVolume(TopPbBox, _stainless, "logiTopPbBox");
    logiTopPbBlock = new G4LogicalVolume(TopPbBlock, _lead, "logiTopPbBlock");

    new G4PVPlacement(nullptr, {}, logiTopPbBlock, "physTopPbBlock", logiTopPbBox, false, 0, OverlapCheck);

    ///////////////////////////////////////////////////////
    // Layer6_Pb body of Lead Box ( Thickness of lead : 150-3 mm )
		// Outer volume : 1100 x 1100 x 1771 mm 
		// Inner volume : 800 x 800 x 1621 mm
    PbBoxOut  = new G4Box("PbBoxOut", leadBoxsize / 2., leadBoxsize, leadBoxzsize);
    PbBoxIn   = new G4Box("PbBoxIn", (leadBoxsize - leadBoxThickness) / 2.,
                        leadBoxsize - leadBoxThickness, leadBoxzsize - leadBoxThickness / 2.);
    PbBox     = new G4SubtractionSolid("PbBox", PbBoxOut, PbBoxIn, 0,
                                   G4ThreeVector(leadBoxThickness / 2., 0, leadBoxThickness / 2.));
    logiPbBox = new G4LogicalVolume(PbBox, _stainless, "logiPbBox");

    PbBlockOut =
        new G4Box("PbBlockOut", leadBoxsize / 2. - leadBoxHousingThickness,
                  leadBoxsize - leadBoxHousingThickness, leadBoxzsize - leadBoxHousingThickness);

    PbBlockIn =
        new G4Box("PbBlockIn", (leadBoxsize - leadBoxThickness + leadBoxHousingThickness) / 2.,
                  leadBoxsize - leadBoxThickness + leadBoxHousingThickness,
                  leadBoxzsize - leadBoxThickness / 2. + leadBoxHousingThickness / 2.);

    PbBlock =
        new G4SubtractionSolid("PbBlock", PbBlockOut, PbBlockIn, nullptr,
                               G4ThreeVector((leadBoxThickness - leadBoxHousingThickness) / 2., 0,
                                             (leadBoxThickness - leadBoxHousingThickness) / 2.));

    logiPbBlock = new G4LogicalVolume(PbBlock, _lead, "logiPbBlock");
    new G4PVPlacement(nullptr, {}, logiPbBlock, "physPbBlock", logiPbBox, false, 0, OverlapCheck);
		//G4cout << " top Pb shield mass : " << logiTopPbBlock->GetMass() / kg << G4endl;

    ///////////////////////////////////////////////////////
    // InnerDetector Space
    ///////////////////////////////////////////////////////
    InnerDetectorBox = new G4Box("InnerDetectorBox", leadBoxsize, leadBoxsize, leadBoxzsize);
    SSTopBoxForInnerDet =
        new G4Box("SSTopBoxForInnerDet", sstopbox_xsize, sstopbox_ysize, sst_zsize_half);
    InnerDetector = new G4UnionSolid("InnerDetectorSolid", InnerDetectorBox, SSTopBoxForInnerDet,
                                     nullptr, G4ThreeVector(0, 0, leadBoxzsize + sst_zsize_half));

    logiInnerDetector = new G4LogicalVolume(InnerDetector, _air, "logiInnerDetector");
    logiInnerDetector->SetVisAttributes(fI_Invisible_VisAttr);

    ///////////////////////////////////////////////////////
    // Layer5_OVC Stainless Steel (SS)
    // top : 748 x 748 x 30 mm
		// cylinder : r = 320 mm, height = 1538 + 20  mm  = 1558 mm
    SSCylinder = new G4Tubs("SSCylinder", 0., ss_radius, ss_height_half + ssb_zsize, 0, 360. * deg);
    SSTop      = new G4Box("SSTop", sstopbox_xsize, sstopbox_ysize, sst_zsize_half);
    SSOVC0     = new G4UnionSolid("SSOVC0", SSCylinder, SSTop, 0, SSTopShift);
    logiSSOVCOuter = new G4LogicalVolume(SSOVC0, _stainless, "logiSSOVCOuterLV");

    SSOVCInnerSolid =
        new G4Tubs("SSOVCInnerCylinder", 0, ss_radius - ss_thick, ss_height_half, 0, 360 * deg);
    logiSSOVCInner = new G4LogicalVolume(SSOVCInnerSolid, _vacuum, "SSOVCInnerLV");
    logiSSOVCInner->SetVisAttributes(fI_Invisible_VisAttr);

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, ssb_zsize), logiSSOVCInner, "physSSOVCInner",
                      logiSSOVCOuter, false, 0, OverlapCheck);

    //////////////////////////////////////////////////////////
    // OVC welding parts

    G4double OVCwe_thickness = 3 * mm / 2.;
    G4ThreeVector OVCwePos   = G4ThreeVector(0, 0, leadBoxzsize - ss_height_half * 2.);
    G4Tubs *OVCwRing =
        new G4Tubs("SSRing", ss_radius, ss_radius + ss_thick, OVCwe_thickness, 0, 360 * deg);

    G4Box *OVCwBox = new G4Box("SSweldingbox", OVCwe_thickness, ss_thick / 2., ss_height_half);
    G4VSolid *OVCwel =
        new G4UnionSolid("SS_Welding", OVCwRing, OVCwBox, 0,
                         G4ThreeVector(0, ss_radius + ss_thick / 2., ss_height_half));

    G4LogicalVolume *lv_OVCwe = new G4LogicalVolume(OVCwel, _aluminum, "OVCwelding");

    new G4PVPlacement(nullptr, OVCwePos, lv_OVCwe, "physOVCwelding", logiInnerDetector, false, 0, OverlapCheck);

    ///////////////////////////////////////////////////////
    // Layer4_ 50K-SHIELD Copper

    Cu4OuterCylinder = new G4Tubs("Cu4OuterCylinder", 0., cu4_radius,
                                  cu4_inner_height + cu4b_zsize + cu4t_zsize, 0, 360. * deg);
    Cu4InnerCylinder =
        new G4Tubs("Cu4InnerCylinder", 0, cu4_radius - cu4_thick, cu4_inner_height, 0., 360. * deg);
    logiCu4Outer = new G4LogicalVolume(Cu4OuterCylinder, _copper, "logiCu4");
    logiCu4Inner = new G4LogicalVolume(Cu4InnerCylinder, _vacuum, "logiCu4Inner");
    logiCu4Inner->SetVisAttributes(fI_Invisible_VisAttr);
    logiCu4Inner->SetVisAttributes(fI_Invisible_VisAttr);

    new G4PVPlacement(nullptr, PosCu4Inner, logiCu4Inner, "physCu4Inner", logiCu4Outer, false, 0,
                      OverlapCheck);

    ///////////////////////////////////////////////////////
    // Layer3_Cu IVC

    Cu3OuterCylinder = new G4Tubs("Cu3OuterCylinder", 0, cu3_radius,
                                  cu3t_zsize + cu3_inner_height + cu3b_zsize, 0, 360 * deg);
    Cu3InnerCylinder =
        new G4Tubs("Cu3InnerCylinder", 0, cu3_radius - cu3_thick, cu3_inner_height, 0, 360 * deg);
    logiCu3Outer = new G4LogicalVolume(Cu3OuterCylinder, _stainless, "logiCu3");
    logiCu3Inner = new G4LogicalVolume(Cu3InnerCylinder, _vacuum, "logiCu3Inner");
    logiCu3Inner->SetVisAttributes(fI_Invisible_VisAttr);

    new G4PVPlacement(nullptr, PosCu3Outer, logiCu3Outer, "physIVCOuter", logiCu4Inner, false, 0,OverlapCheck);
    new G4PVPlacement(nullptr, PosCu3Inner, logiCu3Inner, "physIVCInner", logiCu3Outer, false, 0,OverlapCheck);

    ////////////////////////////////////////////////////////
    // IVC welding part

    G4double IVCwe_thickness = 3 * mm / 2.;
    G4ThreeVector IVCwePos   = G4ThreeVector(
        0, 0, cu4_inner_height - cu3_gap_fromTop - (cu3_inner_height + cu3t_zsize) * 2.);
    G4Tubs *IVCwRing =
        new G4Tubs("Cu3Ring", cu3_radius, cu3_radius + cu3_thick, IVCwe_thickness, 0, 360 * deg);

    G4Box *IVCwBox = new G4Box("Cu3weldingbox", IVCwe_thickness, cu3_thick / 2., cu3_inner_height);
    G4VSolid *IVCwel =
        new G4UnionSolid("IVC_Welding", IVCwRing, IVCwBox, 0,
                         G4ThreeVector(0, cu3_radius + cu3_thick / 2., cu3_inner_height));

    G4LogicalVolume *lv_IVCwe = new G4LogicalVolume(IVCwel, _copper, "IVCwelding");

    new G4PVPlacement(nullptr, IVCwePos, lv_IVCwe, "physIVCwelding", logiCu4Inner, false, 0, OverlapCheck);

    ///////////////////////////////////////////////////////
    // Layer2_Cu3 SHIELD-STILL

    Cu2OuterCylinder = new G4Tubs("Cu2OuterCylinder", 0, cu2_radius,
                                  cu2_inner_height + cu2b_zsize + cu2t_zsize, 0, 360. * deg);
    Cu2InnerCylinder =
        new G4Tubs("Cu2InnerCylinder", 0, cu2_radius - cu2_thick, cu2_inner_height, 0, 360. * deg);
    logiCu2Outer = new G4LogicalVolume(Cu2OuterCylinder, _copper, "logiCu2Outer");
    logiCu2Inner = new G4LogicalVolume(Cu2InnerCylinder, _vacuum, "logiCu2Inner");
    logiCu2Inner->SetVisAttributes(fI_Invisible_VisAttr);

    new G4PVPlacement(nullptr, PosCu2Outer, logiCu2Outer, "physCu2Outer", logiCu3Inner, false, 0,OverlapCheck);
    new G4PVPlacement(nullptr, PosCu2Inner, logiCu2Inner, "physCu2Inner", logiCu2Outer, false, 0,OverlapCheck);

    ///////////////////////////////////////////////////////
    // Layer1_Cu1 50mk-SHIELD

    Cu1OuterCylinder = new G4Tubs("Cu1OuterCylinder", 0., cu1_radius,
                                  cu1_inner_height + cu1b_zsize + cu1t_zsize, 0., 360. * deg);
    Cu1InnerCylinder = new G4Tubs("Cu1InnerCylinder", 0., cu1_radius - cu1_thick, cu1_inner_height,
                                  0., 360. * deg);

    logiCu1Outer = new G4LogicalVolume(Cu1OuterCylinder, _copper, "logiCu1Outer");
    logiCu1Inner = new G4LogicalVolume(Cu1InnerCylinder, _vacuum, "logiCu1Inner");
    logiCu1Inner->SetVisAttributes(fI_Invisible_VisAttr);

    new G4PVPlacement(nullptr, PosCu1Inner, logiCu1Inner, "physCu1Inner", logiCu1Outer, false, 0, OverlapCheck);
    new G4PVPlacement(nullptr, PosCu1Outer, logiCu1Outer, "physCu1Outer", logiCu2Inner, false, 0, OverlapCheck);

    ///////////////////////////////////////////////////////
    // Mixing Chamber Cu Plate

    CuMCPlate = new G4Tubs("CuMCPlate", 0, copperMixingChamberPlateRadius,
                           copperMixingChamberPlatThickness, 0, 360. * deg);
    logiCuMCP = new G4LogicalVolume(CuMCPlate, _copper, "logiCuMCP");

    ///////////////////////////////////////////////////////
    // Plate1_Cu
    CuPlate1 = new G4Tubs("CuPlate1", 0, copperMixingChamberPlateRadius, copperPlate1ZSize_half, 0,
                          360. * deg);
    logiCuP1 = new G4LogicalVolume(CuPlate1, _copper, "logiCuP1");

    ///////////////////////////////////////////////////////
    // Plate2_Pb Shield
    PbPlate2 = new G4Tubs("PbPlate2", 0, copperMixingChamberPlateRadius, leadPlate2ZSize_half, 0,
                          360. * deg);
    logiPbP2 = new G4LogicalVolume(PbPlate2, _lead, "logiPbP2");

    ///////////////////////////////////////////////////////
    CuPlate3 = new G4Tubs("CuPlate3", 0, copperMixingChamberPlateRadius, cup3_zsize, 0, 360. * deg);
    logiCuP3 = new G4LogicalVolume(CuPlate3, _copper, "logiCuP3");

    /////////////////////////////////////////////////////////
    // G10 material holders between SS and 50K-SHIELD Copper

    for (int g10Lyr = 0; g10Lyr < 5; g10Lyr++) {
        G4Tubs *PlateGap4Vol = new G4Tubs("PlateGap4Vol", 0, PlateGapR[g10Lyr] - 3. * mm,
                                          CuholderH / 2., 0. * deg, 360. * deg);
        PlateGap4VolLogi[g10Lyr] =
            new G4LogicalVolume(PlateGap4Vol, _vacuum, "PlateGap4VolLogical");
        PlateGap4VolLogi[g10Lyr]->SetVisAttributes(fI_Invisible_VisAttr);

        G4Tubs *G10holderSolid =
            new G4Tubs("G10holderSolid", // name
                       Cuholder_r, Cuholder_R, CuholderH / 2, 0. * deg, 360. * deg);
        G4LogicalVolume *G10holderLogi =
            new G4LogicalVolume(G10holderSolid, g10material, "G10holderLogical");

        G4int copyNo = 0;
        if (g10Lyr < 4) {
            for (int NinLyr = 0; NinLyr < 3; NinLyr++) {
                copyNo               = NinLyr;
                xG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * cos((NinLyr * 120.) * deg));
                yG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * sin((NinLyr * 120.) * deg));
                PosG10holder         = G4ThreeVector(xG10[NinLyr][g10Lyr], yG10[NinLyr][g10Lyr], 0);
                new G4PVPlacement(nullptr, PosG10holder, G10holderLogi, "G10holderPhys",
                                  PlateGap4VolLogi[g10Lyr], false, copyNo, OverlapCheck);
            }
        } else {
            for (int NinLyr = 0; NinLyr < 6; NinLyr++) {
                copyNo               = NinLyr;
                xG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * cos((NinLyr * 60.) * deg));
                yG10[NinLyr][g10Lyr] = (G10holderXY[g10Lyr] * sin((NinLyr * 60.) * deg));
                PosG10holder         = G4ThreeVector(xG10[NinLyr][g10Lyr], yG10[NinLyr][g10Lyr], 0);
                new G4PVPlacement(nullptr, PosG10holder, G10holderLogi, "G10holderPhys",
                                  PlateGap4VolLogi[g10Lyr], false, copyNo, OverlapCheck);
            }
        }
        G4VisAttributes *G10holderVisAtt = new G4VisAttributes(lgreen);
        G10holderVisAtt->SetVisibility(true);
        G10holderVisAtt->SetForceSolid(true);
        G10holderLogi->SetVisAttributes(G10holderVisAtt);
				/*
        G4cout << "###   " << G10holderLogi->GetName() << "  " << G10holderLogi->GetMass() / kg
               << G4endl;
        G4cout << "###   " << PlateGap4VolLogi[g10Lyr]->GetName() << "  "
               << PlateGap4VolLogi[g10Lyr]->GetMass() / kg << G4endl;
							 */
    }

    // -----------------------------------
    // H beam construction
    // -----------------------------------
    //
    /// Vertical H beam
    //

    Bplate = new G4Box("Bplate", Bplate_size, Bplate_size, Bplate_thick_half);
    Hbeam1 = new G4ExtrudedSolid("Hbeam1", large_points, Hbeam1_length, 0.0, 1.0, 0.0, 1.0);
    Hbeam2 = new G4ExtrudedSolid("Hbeam2", large_points, Hbeam2_length, 1.0, 1.0, 1.0, 1.0);

    UnionHbeam1 = new G4UnionSolid("UnionHbeam1", Hbeam1, Bplate, 0, -posBplate1);
    UnionHbeam2 = new G4UnionSolid("UnionHbeam2", Hbeam2, Bplate, 0, -posBplate2);

    logicHbeam1 = new G4LogicalVolume(UnionHbeam1, _steel, "logicHbeam1");
    logicHbeam2 = new G4LogicalVolume(UnionHbeam2, _steel, "logicHbeam2");

    Hb1VisAtt = new G4VisAttributes(yellow);
    Hb2VisAtt = new G4VisAttributes(yellow);

    logicHbeam1->SetVisAttributes(Hb1VisAtt);
    logicHbeam2->SetVisAttributes(Hb2VisAtt);

    /// Horizontal H -beam X
    Ho_beamX =
        new G4ExtrudedSolid("Ho_beamX", small_points, Horizontal_beamX_length, 1.0, 1.0, 1.0, 1.0);

    logicHo_beamX = new G4LogicalVolume(Ho_beamX, _steel, "logicHo_beamX");

    HoBXVisAtt = new G4VisAttributes(green);
    logicHo_beamX->SetVisAttributes(HoBXVisAtt);

    Ho_beamXLar = new G4ExtrudedSolid("Ho_beamXLar", large_points, Horizontal_beamXLarge_length,
                                      1.0, 1.0, 1.0, 1.0);
    logicHo_beamXLar = new G4LogicalVolume(Ho_beamXLar, _steel, "logicHo_beamXLar");

    HoBXLVisAtt = new G4VisAttributes(yellow);
    logicHo_beamXLar->SetVisAttributes(HoBXLVisAtt);

    G4ThreeVector Gantry_Top_Trans(0, 0, -Hbeam1_length - 20. + HoLayer[2] - 75. + posMountOrigin1);

    /// Horizontal H -beam Y

    Ho_beamY =
        new G4ExtrudedSolid("Ho_beamY", small_points, Horizontal_beamY_length, 1.0, 1.0, 1.0, 1.0);

    logicHo_beamY = new G4LogicalVolume(Ho_beamY, _steel, "logicHo_beamY");

    HoBYVisAtt = new G4VisAttributes(green);
    logicHo_beamY->SetVisAttributes(HoBYVisAtt);

    Ho_beamYLar = new G4ExtrudedSolid("Ho_beamYLar", large_points, Horizontal_beamYLarge_length,
                                      1.0, 1.0, 1.0, 1.0);

    logicHo_beamYLar = new G4LogicalVolume(Ho_beamYLar, _steel, "logicHo_beamYLar");

    HoBYLVisAtt = new G4VisAttributes(yellow);
    logicHo_beamYLar->SetVisAttributes(HoBYLVisAtt);

		///////////////////////////////////
		// Neutron Mode Geometry Test Volumes-------------------------------------
    // B4C Rubber (24%)
    toyRubberB4C_Box =
        new G4Box("ToyRubberB4C_Box", workboxX + boronCarbideRubber_thickness,
                  workboxY + boronCarbideRubber_thickness, workboxH + boronCarbideRubber_thickness);
    toyRubberB4C_BoxLogical =
        new G4LogicalVolume(toyRubberB4C_Box, _B4CRubber24perCent, "ToyRubberB4C_LV");

    // Polyethelene
    toyPE_Box =
        new G4Box("ToyPolyEthylene_Box", workboxX + boronCarbideRubber_thickness + toyPE_thickness,
                  workboxY + boronCarbideRubber_thickness + toyPE_thickness,
                  workboxH + boronCarbideRubber_thickness + toyPE_thickness);
    toyPE_BoxLogical = new G4LogicalVolume(toyPE_Box, _polyethylene, "PolyEthylene_LV");

    toyPE_BoxLogical->SetVisAttributes(fI_PE_VisAttr);

    // Boric acid panel
    topBoricAcid_HousingBox = new G4Box("BoricAcid_Top_Mother_Box", length_Top, width_Top,
                                        totalThickness_BoricAcidPanel / 2.);
    topBoricAcid_Block =
        new G4Box("BoricAcid_Top_Block", length_Top - thickness_BoricAcidSpacing,
                  width_Top - thickness_BoricAcidSpacing, thickness_BoricAcid / 2.);

    topBoricAcid_MotherLogical =
        new G4LogicalVolume(topBoricAcid_HousingBox, _polyethylene, "BoricAcid_Top_Mother_LV");

    topBoricAcid_BlockLogical =
        new G4LogicalVolume(topBoricAcid_Block, _BoricAcidPowder, "BoricAcid_Top_Block_LV");

    new G4PVPlacement(nullptr, {}, topBoricAcid_BlockLogical, "BoricAcid_Top_Block_PV",
                      topBoricAcid_MotherLogical, false, 0, OverlapCheck);

    topBoricAcid_MotherLogical->SetVisAttributes(fI_BoricAcidHousing_VisAttr);
    topBoricAcid_BlockLogical->SetVisAttributes(fI_BoricAcid_VisAttr);
		// --------------------------------------------------------------------------
		///////////////////////////////////

    // Calculate position of superconducting magnet shield and detector array
    superMagnetTlate = G4ThreeVector(
        0, 0, CenterCuP6Z - cup3_zsize - scShieldSizeZ / 2. - scsCuFrameTopDiskThickness);

    detectorArrayTlate = G4ThreeVector(0, 0,
                                       scShieldSizeZ / 2. - scShieldLeadThickness / 2. -
                                           arrayHeight / 2. - scsCuFrameTopDiskThickness / 2.);
    //////////////////////////////
    // Place all the geometries at Work Area
    //////////////////////////////

    // Calculation of position parameters
    G4ThreeVector FB_Trans =
        G4ThreeVector(0, 0, -leadBoxzsize + (sideFBScint_boxHeight_half + scint_reflector_thick));
    G4ThreeVector LR_Trans =
        G4ThreeVector(0, 0, -leadBoxzsize + (sideLRScint_boxHeight_half + scint_reflector_thick));
    G4double sideLRScint_YPos =
        std::fmax(2 * (sideFBScint_boxWidth_half + scint_reflector_thick), leadBoxsize);

    if (fEnable_OriginalGeom) {
        G4double topScint_boxThick_half    = scint_Thick_half;
        G4double sideLRScint_boxThick_half = scint_Thick_half;
        G4double sideFBScint_boxThick_half = scint_Thick_half;
        physTopPbBox1 =
            new G4PVPlacement(nullptr, G4ThreeVector(0., -pbtopbox_size / 2., CenterPbTZ),
                              logiTopPbBox, "physTopPbBox1", logiWorkArea, false, 0,OverlapCheck);
        new G4PVPlacement(nullptr, G4ThreeVector(0., pbtopbox_size / 2., CenterPbTZ), logiTopPbBox,
                          "physTopPbBox2", logiWorkArea, false, 0,OverlapCheck);
        new G4PVPlacement(nullptr, PosInnerDet, logiInnerDetector, "physInnerDetector",
                          logiWorkArea, false, 0,OverlapCheck);

        if (fEnable_InnerDetector) {
            G4ThreeVector PbBoxXHalfTrans1(-leadBoxsize / 2., 0, 0);
            G4ThreeVector PbBoxXHalfTrans2(leadBoxsize / 2., 0, 0);
            G4ThreeVector ShiftPbBoxInInnerDet(0, 0, 0.);
            G4RotationMatrix *PbBox2_alignMtx = new G4RotationMatrix();
            PbBox2_alignMtx->rotateZ(180 * deg);
            new G4PVPlacement(nullptr, PbBoxXHalfTrans1 + ShiftPbBoxInInnerDet, logiPbBox,
                              "physPbBox1", logiInnerDetector, false, 0, OverlapCheck);
            new G4PVPlacement(PbBox2_alignMtx, PbBoxXHalfTrans2 + ShiftPbBoxInInnerDet, logiPbBox,
                              "physPbBox2", logiInnerDetector, false, 0, OverlapCheck);
            PosSS = G4ThreeVector(0, 0, leadBoxzsize - ss_height_half - ssb_zsize);
            new G4PVPlacement(nullptr, PosSS, logiSSOVCOuter, "physSSOVC", logiInnerDetector, false,
                              0, OverlapCheck);
            new G4PVPlacement(nullptr, PosCu4Outer, logiCu4Outer, "physCu4Outer", logiSSOVCInner,
                              false, 0, OverlapCheck);
        }

        if (fEnable_Innermost) {
            new G4PVPlacement(nullptr, PosCuMCP, logiCuMCP, "physCuMCP", logiCu1Inner, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, PosCuP1, logiCuP1, "physCuP1", logiCu1Inner, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, PosPbP2, logiPbP2, "physPbP2", logiCu1Inner, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, PosCuP3, logiCuP3, "physCuP3", logiCu1Inner, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, PosCuP4, logiCuP1, "physCuP4", logiCu1Inner, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, PosPbP5, logiPbP2, "physPbP5", logiCu1Inner, false, 0, OverlapCheck);
            new G4PVPlacement(nullptr, PosCuP6, logiCuP3, "physCuP6", logiCu1Inner, false, 0, OverlapCheck);

            G4LogicalVolume *G10MotherLV[5] = {logiCu1Inner, logiCu2Inner, logiCu3Inner,
                                               logiCu4Inner, logiSSOVCInner};

            const char G10PVName[5][30] = {"PlateGap4VolPhys_Cu1", "PlateGap4VolPhys_Cu2",
                                           "PlateGap4VolPhys_Cu3", "PlateGap4VolPhys_Cu4",
                                           "PlateGap4VolPhys_SSOVC"};
            for (G4int i = 0; i < 5; i++) {
                new G4PVPlacement(nullptr, G4ThreeVector(0, 0, CenterPlateGap4VolZ[i]),
                                  PlateGap4VolLogi[i], G10PVName[i], G10MotherLV[i], false, 0, OverlapCheck);
            }

            if (fI_Enable_SuperConductingShield)
                superMagnetShield_LV = Build_I_SCMagnetShieldAt(logiCu1Inner, superMagnetTlate, 1,
                                                                _copper, _vacuum, _lead);
            else
                superMagnetShield_LV = Build_I_SCMagnetShieldAt(logiCu1Inner, superMagnetTlate, 2,
                                                                _copper, _vacuum, _lead);
						
            if (fI_Enable_CrystalArray) {
                fI_DetectorModuleRegion = new G4Region("DetectorModuleRegion");
								fI_crystalsRegion = new G4Region("crystals");

                detectorArray_LV = Build_I_DetectorArray(3, _copper, _vacuum);
                detectorArray_LV->SetVisAttributes(fI_Invisible_VisAttr);

                new G4PVPlacement(nullptr, detectorArrayTlate, detectorArray_LV, "DetectorArray_PV",
                                  superMagnetShield_LV, false, 0, OverlapCheck);
            }

						if(fDbgMsgOn){
							std::cout << "Cu1 inner center position     : " << PosSS+G4ThreeVector(0,0,ssb_zsize)+PosCu4Outer+PosCu4Inner+PosCu3Outer+PosCu3Inner+PosCu2Outer+PosCu2Inner+PosCu1Outer+PosCu1Inner << std::endl;
							std::cout << "SC magnet center position     : " << superMagnetTlate << std::endl;
							std::cout << "Work area height              : " << workboxH << std::endl;
							std::cout << "Inner detector position       : " << PosInnerDet << std::endl;
							std::cout << "Detector array center position: " << detectorArrayTlate << std::endl;
						}
        }

        if (fEnable_Gantry) {
            static constexpr char beamX_PVname[][50] = {"HoriXBeamLyr1PV", "HoriXBeamLyr2PV",
                                                        "HoriXBeamLyr3PV"};
            static constexpr char beamY_PVname[][50] = {"HoriYBeamLyr1PV", "HoriYBeamLyr2PV",
                                                        "HoriYBeamLyr3PV"};
            for (int i = 0; i < 4; i++) {
                new G4PVPlacement(
                    VeBeamRotation,
                    G4ThreeVector(x_beam[i] * 750., y_beam[i] * 1600., posMountOrigin1)
                        .rotateZ(90. * deg),
                    logicHbeam1, "VerticalBeamPV", logiWorkArea, false, i, OverlapCheck);
            }

            for (int j = 0; j < 2; j++) {
                new G4PVPlacement(
                    nullptr,
                    G4ThreeVector(HoLineY[j] * 750., 0, posMountOrigin2).rotateZ(90. * deg),
                    logicHbeam2, "VerticalBottomBeamPV", logiWorkArea, false, j, OverlapCheck);
                for (int i = 0; i < 3; i++) {
                    new G4PVPlacement(HoBXRotation,
                                      G4ThreeVector(HoLine[j] * 750., 0.,
                                                    -Hbeam1_length - 20. + HoLayer[i] -
                                                        ((i == 0) ? 80. : 0.) + posMountOrigin1)
                                          .rotateZ(90. * deg),
                                      (i == 0) ? logicHo_beamXLar : logicHo_beamX, beamX_PVname[i],
                                      logiWorkArea, false, j, OverlapCheck);
                }
            }
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 3; i++) {
                    new G4PVPlacement(HoBYRotation,
                                      G4ThreeVector(0., HoLineY[j] * 1600.,
                                                    -Hbeam1_length - 20. + HoLayerY[i] -
                                                        ((i == 0) ? 80. : 0.) + posMountOrigin1)
                                          .rotateZ(90. * deg),
                                      (i == 0) ? logicHo_beamYLar : logicHo_beamY, beamY_PVname[i],
                                      logiWorkArea, false, j, OverlapCheck);
                }
            }
        }
        if (fEnable_Scintillator) {
            G4LogicalVolume *logicPMTEnvelope = Build_I_MuonVetoPMT(_aluminum, _grease);
            G4ThreeVector Top_Trans =
                G4ThreeVector(0, 0,
                              physTopPbBox1->GetTranslation().z() + TopPbBox->GetZHalfLength() +
                                  (topScint_boxThick_half + scint_reflector_thick));
            G4ThreeVector Top_PMTTrans = G4ThreeVector(
                0,
                topScint_boxHeight_half + topScint_trapZ_half * 2. + muonVetoPMTEnvelopeSizeZ / 2.,
                0);

            G4RotationMatrix *PMT_alignMtx_Top1 = new G4RotationMatrix();
            PMT_alignMtx_Top1->rotateX(-90 * deg);
            G4RotationMatrix *PMT_alignMtx_Top2 = new G4RotationMatrix(*PMT_alignMtx_Top1);
            PMT_alignMtx_Top2->rotateX(180 * deg);
            G4LogicalVolume *topScintLV = Build_I_TopScintillator(_vinylt, _tyvek, Tyvek_opsurf);

            G4VPhysicalVolume *T1_Phys = new G4PVPlacement(
                nullptr,
                G4ThreeVector(-(topScint_boxWidth_half + scint_reflector_thick), 0, 0) + Top_Trans,
                topScintLV, "Envelope_MuonTopScintillator_PV", logiWorkArea, false, 0, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T1_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 0, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T1_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 1, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T1_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 2, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T1_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 3, OverlapCheck);

            G4VPhysicalVolume *T2_Phys = new G4PVPlacement(
                0,
                G4ThreeVector((topScint_boxWidth_half + scint_reflector_thick), 0, 0) + Top_Trans,
                topScintLV, "Envelope_MuonTopScintillator_PV", logiWorkArea, false, 1, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T2_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 4, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top2,
                              T2_Phys->GetTranslation() - Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 5, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T2_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(-topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 6, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_Top1,
                              T2_Phys->GetTranslation() + Top_PMTTrans +
                                  G4ThreeVector(+topScint_trapX2_half, 0, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 7, OverlapCheck);

            G4ThreeVector LR_PMTTrans =
                G4ThreeVector(sideLRScint_trapZ_half * 2 + sideLRScint_boxWidth_half +
                                  muonVetoPMTEnvelopeSizeZ / 2.,
                              0, 0);
            G4ThreeVector LR_Mod_PMTTrans =
                G4ThreeVector(0, 0, sideLRScint_boxHeight_half - sideLRScint_trapX2_half * 3);
            G4RotationMatrix *PMT_alignMtx_LRSide1 = new G4RotationMatrix();
            PMT_alignMtx_LRSide1->rotateY(90 * deg);
            G4RotationMatrix *PMT_alignMtx_LRSide2 = new G4RotationMatrix();
            PMT_alignMtx_LRSide2->rotateY(-90 * deg);

            G4LogicalVolume *sideLRScintLV =
                Build_I_Side_LRScintillator(_vinylt, _tyvek, Tyvek_opsurf);
            G4VPhysicalVolume *L1_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide2,
                PosPbBox +
                    G4ThreeVector((sideLRScint_boxWidth_half + scint_reflector_thick),
                                  -sideLRScint_YPos -
                                      (sideLRScint_boxThick_half + scint_reflector_thick),
                                  0) +
                    LR_Trans,
                sideLRScintLV, "Envelope_MuonSideLRScintillator_PV", logiWorkArea, false, 2, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              L1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 8, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              L1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 9, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              L1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 10, OverlapCheck);

            G4VPhysicalVolume *L2_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide1,
                PosPbBox +
                    G4ThreeVector(-(sideLRScint_boxWidth_half + scint_reflector_thick),
                                  -sideLRScint_YPos -
                                      (sideLRScint_boxThick_half + scint_reflector_thick),
                                  0) +
                    LR_Trans,
                sideLRScintLV, "Envelope_MuonSideLRScintillator_PV", logiWorkArea, false, 3, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              L2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 11, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              L2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 12, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              L2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 13, OverlapCheck);

            G4VPhysicalVolume *R1_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide2,
                PosPbBox +
                    G4ThreeVector(
                        +(sideLRScint_boxWidth_half + scint_reflector_thick),
                        sideLRScint_YPos + (sideLRScint_boxThick_half + scint_reflector_thick), 0) +
                    LR_Trans,
                sideLRScintLV, "Envelope_MuonSideLRScintillator_PV", logiWorkArea, false, 4, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              R1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 14, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              R1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 15, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide1,
                              R1_Phys->GetTranslation() + LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 16, OverlapCheck);

            G4VPhysicalVolume *R2_Phys = new G4PVPlacement(
                scintMother_alignMtx_LRSide1,
                PosPbBox +
                    G4ThreeVector(
                        -(sideLRScint_boxWidth_half + scint_reflector_thick),
                        sideLRScint_YPos + (sideLRScint_boxThick_half + scint_reflector_thick), 0) +
                    LR_Trans,
                sideLRScintLV, "Envelope_MuonSideLRScintillator_PV", logiWorkArea, false, 5, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              R2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans +
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 17, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              R2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 18, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_LRSide2,
                              R2_Phys->GetTranslation() - LR_PMTTrans - LR_Mod_PMTTrans -
                                  G4ThreeVector(0, 0, sideLRScint_trapX2_half * 2.),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 19, OverlapCheck);

            G4RotationMatrix *scintMother_alignMtx_FBSide = new G4RotationMatrix;
            scintMother_alignMtx_FBSide->rotateY(-90 * deg);
            scintMother_alignMtx_FBSide->rotateZ(-90 * deg);
            G4ThreeVector FB_PMTTrans =
                G4ThreeVector(0, 0,
                              sideFBScint_boxHeight_half + sideFBScint_trapZ_half * 2. +
                                  muonVetoPMTEnvelopeSizeZ / 2.);
            G4RotationMatrix *PMT_alignMtx_FBSide = new G4RotationMatrix();
            PMT_alignMtx_FBSide->rotateX(180 * deg);
            G4LogicalVolume *sideFBScintLV =
                Build_I_Side_FBScintillator(_vinylt, _tyvek, Tyvek_opsurf);

            G4VPhysicalVolume *B1_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(-leadBoxsize -
                                      (sideFBScint_boxThick_half + scint_reflector_thick),
                                  -(sideFBScint_boxWidth_half + scint_reflector_thick), 0) +
                    FB_Trans,
                sideFBScintLV, "Envelope_MuonSideFBScintillator_PV", logiWorkArea, false, 6, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 20, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 21, OverlapCheck);

            G4VPhysicalVolume *B2_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(-leadBoxsize -
                                      (sideFBScint_boxThick_half + scint_reflector_thick),
                                  (sideFBScint_boxWidth_half + scint_reflector_thick), 0) +
                    FB_Trans,
                sideFBScintLV, "Envelope_MuonSideFBScintillator_PV", logiWorkArea, false, 7, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 22, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              B2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 23, OverlapCheck);

            G4VPhysicalVolume *F1_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(leadBoxsize + (sideFBScint_boxThick_half + scint_reflector_thick),
                                  -(sideFBScint_boxWidth_half + scint_reflector_thick), 0) +
                    FB_Trans,
                sideFBScintLV, "Envelope_MuonSideFBScintillator_PV", logiWorkArea, false, 8, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 24, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F1_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 25, OverlapCheck);

            G4VPhysicalVolume *F2_Phys = new G4PVPlacement(
                scintMother_alignMtx_FBSide,
                PosPbBox +
                    G4ThreeVector(leadBoxsize + (sideFBScint_boxThick_half + scint_reflector_thick),
                                  (sideFBScint_boxWidth_half + scint_reflector_thick), 0) +
                    FB_Trans,
                sideFBScintLV, "Envelope_MuonSideFBScintillator_PV", logiWorkArea, false, 9, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, -sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 26, OverlapCheck);
            new G4PVPlacement(PMT_alignMtx_FBSide,
                              F2_Phys->GetTranslation() + FB_PMTTrans +
                                  G4ThreeVector(0, +sideFBScint_trapX2_half, 0),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 27, OverlapCheck);

						// Additional Plastic Scintillator for Muon Veto (FB)
            G4RotationMatrix *scintMother_alignMtx_mufflerFBSide = new G4RotationMatrix();
            scintMother_alignMtx_mufflerFBSide->rotateZ(90. * deg);
            scintMother_alignMtx_mufflerFBSide->rotateX(90. * deg);

            G4ThreeVector mufflerFBScintDisplacer = {
                leadBoxsize + scint_muffler_Thick_half + scint_reflector_thick +
                    mufflerFBScint_YDistanceFromLead,
                (mufflerFBScint_boxWidth_half + scint_reflector_thick),
                sideFBScint_boxHeight_half + mufflerFBScint_boxHeight_half};
            G4ThreeVector mufflerFBPMTDisplacer = {0, 0,
                                                   -mufflerFBScint_boxHeight_half -
                                                       mufflerFBScint_trapZ_half * 2. -
                                                       muonVetoPMTEnvelopeSizeZ / 2.};

            G4LogicalVolume *mufflerFBScintLV = Build_I_Muffler_FBScintillator(
                _vinylt, _tyvek, _vm2000, Tyvek_opsurf, Vikuiti_opsurf);

            G4VPhysicalVolume *mufflerFB1ScintPV = new G4PVPlacement(
                scintMother_alignMtx_mufflerFBSide, PosPbBox + mufflerFBScintDisplacer + FB_Trans,
                mufflerFBScintLV, "Envelope_MuonMufflerFBScintillator_PV", logiWorkArea, false, 10, OverlapCheck);
            new G4PVPlacement(nullptr, mufflerFB1ScintPV->GetTranslation() + mufflerFBPMTDisplacer,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 28, OverlapCheck);
            mufflerFBScintDisplacer[1] *= -1.;
            G4VPhysicalVolume *mufflerFB2ScintPV = new G4PVPlacement(
                scintMother_alignMtx_mufflerFBSide, PosPbBox + mufflerFBScintDisplacer + FB_Trans,
                mufflerFBScintLV, "Envelope_MuonMufflerFBScintillator_PV", logiWorkArea, false, 11, OverlapCheck);
            new G4PVPlacement(nullptr, mufflerFB2ScintPV->GetTranslation() + mufflerFBPMTDisplacer,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 29, OverlapCheck);
            mufflerFBScintDisplacer[0] *= -1.;
            mufflerFBScintDisplacer[1] *= -1.;
            G4VPhysicalVolume *mufflerFB3ScintPV = new G4PVPlacement(
                scintMother_alignMtx_mufflerFBSide, PosPbBox + mufflerFBScintDisplacer + FB_Trans,
                mufflerFBScintLV, "Envelope_MuonMufflerFBScintillator_PV", logiWorkArea, false, 12, OverlapCheck);
            new G4PVPlacement(nullptr, mufflerFB3ScintPV->GetTranslation() + mufflerFBPMTDisplacer,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 30, OverlapCheck);
            mufflerFBScintDisplacer[1] *= -1.;
            G4VPhysicalVolume *mufflerFB4ScintPV = new G4PVPlacement(
                scintMother_alignMtx_mufflerFBSide, PosPbBox + mufflerFBScintDisplacer + FB_Trans,
                mufflerFBScintLV, "Envelope_MuonMufflerFBScintillator_PV", logiWorkArea, false, 13, OverlapCheck);
            new G4PVPlacement(nullptr, mufflerFB4ScintPV->GetTranslation() + mufflerFBPMTDisplacer,
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 31, OverlapCheck);

						// Additional Plastic Scintillator for Muon Veto (LR)
            G4RotationMatrix *scintMother_alignMtx_mufflerLRSide = new G4RotationMatrix();
            scintMother_alignMtx_mufflerLRSide->rotateZ(90. * deg);
            scintMother_alignMtx_mufflerLRSide->rotateX(90. * deg);

            G4RotationMatrix *PMTMother_alignMtx_mufflerLRSide = new G4RotationMatrix();
            PMTMother_alignMtx_mufflerLRSide->rotateY(90. * deg);

						G4ThreeVector mufflerLRScintDisplacer = {0,
							                                    -leadBoxsize - scint_muffler_Thick_half 
																									- scint_reflector_thick 
																									- mufflerLRScint_YDistanceFromLead*1.6,
							                                    leadBoxzsize + mufflerLRScint_boxWidth_half*2 };

						G4ThreeVector mufflerLRPMTDisplacer = { -mufflerLRScint_boxHeight_half 
							                                    - mufflerLRScint_trapZ_half * 2. 
																									- muonVetoPMTEnvelopeSizeZ / 2.,
																									0, 0};

            G4LogicalVolume *mufflerLRScintLV = Build_I_Muffler_LRScintillator(
                _vinylt, _tyvek, _vm2000, Tyvek_opsurf, Vikuiti_opsurf);

						// ___ 1st additional LR placement
            G4VPhysicalVolume *mufflerLR1ScintPV = new G4PVPlacement(
                G4Transform3D(*scintMother_alignMtx_mufflerLRSide, 
									PosPbBox + mufflerLRScintDisplacer + LR_Trans),
                mufflerLRScintLV, "Envelope_MuonMufflerLRScintillator_PV", 
								logiWorkArea, false, 14, OverlapCheck);
            
						new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR1ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 32, OverlapCheck);
						
						mufflerLRPMTDisplacer[0] *= -1.;
						PMTMother_alignMtx_mufflerLRSide->rotateY(180. * deg);
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR1ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 33, OverlapCheck);

						// ___ 2nd additional LR placement
						mufflerLRScintDisplacer[2] -= (mufflerLRScint_boxWidth_half + scint_reflector_thick)*2;
            G4VPhysicalVolume *mufflerLR2ScintPV = new G4PVPlacement(
                G4Transform3D(*scintMother_alignMtx_mufflerLRSide, 
								 PosPbBox + mufflerLRScintDisplacer + LR_Trans),
                mufflerLRScintLV, "Envelope_MuonMufflerLRScintillator_PV", 
								logiWorkArea, false, 15, OverlapCheck);

            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR2ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
                  logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 34, OverlapCheck);

						mufflerLRPMTDisplacer[0] *= -1.;
						PMTMother_alignMtx_mufflerLRSide->rotateY(180. * deg);
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR2ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 35, OverlapCheck);

						// ___ 3rd additional LR placement
						mufflerLRScintDisplacer[1] += mufflerLRScint_YDistanceFromLead*0.6;
            mufflerLRScintDisplacer[1] *= -1.;
						mufflerLRScintDisplacer[2] += (mufflerLRScint_boxWidth_half + scint_reflector_thick)*2;
            scintMother_alignMtx_mufflerLRSide->rotateX(30. * deg);

            G4VPhysicalVolume *mufflerLR3ScintPV = new G4PVPlacement(
                G4Transform3D(*scintMother_alignMtx_mufflerLRSide, 
									PosPbBox + mufflerLRScintDisplacer + LR_Trans),
                mufflerLRScintLV, "Envelope_MuonMufflerLRScintillator_PV", 
								logiWorkArea, false, 16, OverlapCheck);
							 
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR3ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 36, OverlapCheck);

						mufflerLRPMTDisplacer[0] *= -1.;
						PMTMother_alignMtx_mufflerLRSide->rotateY(180. * deg);
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR3ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 37, OverlapCheck);

						// ___ 4th additional LR placement
						mufflerLRScintDisplacer[1] += mufflerLRScint_boxWidth_half + scint_reflector_thick;
						mufflerLRScintDisplacer[2] -= sqrt(3)*(mufflerLRScint_boxWidth_half + scint_reflector_thick);
            G4VPhysicalVolume *mufflerLR4ScintPV = new G4PVPlacement(
                G4Transform3D(*scintMother_alignMtx_mufflerLRSide, 
									PosPbBox + mufflerLRScintDisplacer + LR_Trans),
                mufflerLRScintLV, "Envelope_MuonMufflerLRScintillator_PV", 
								logiWorkArea, false, 17, OverlapCheck);

            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR4ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 38, OverlapCheck);

						mufflerLRPMTDisplacer[0] *= -1.;
						PMTMother_alignMtx_mufflerLRSide->rotateY(180. * deg);
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_mufflerLRSide, 
									mufflerLR4ScintPV->GetTranslation() + mufflerLRPMTDisplacer),
                              logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 39, OverlapCheck);

						// Ground Plastic Scintillator for Muon Veto
            G4RotationMatrix *scintMother_alignMtx_groundSide = new G4RotationMatrix();
            scintMother_alignMtx_groundSide->rotateZ(90. * deg);

            G4RotationMatrix *PMTMother_alignMtx_groundSide = new G4RotationMatrix();
            PMTMother_alignMtx_groundSide->rotateY(90. * deg);

						G4ThreeVector groundScintDisplacer = {0, 
							                                    mufflerLRScint_boxWidth_half + scint_reflector_thick,
							                                    (leadBoxZinnerlength + leadBoxThickness) / 2. +
																										scint_muffler_Thick_half + scint_reflector_thick 
																										+ PilotValues::boratedPE_BottomGap};
																										//scint_bottom_Thick_half + bottomNeutronShieldGap};

						G4ThreeVector groundPMTDisplacer = { -mufflerLRScint_boxHeight_half 
							                                    - mufflerLRScint_trapZ_half * 2. 
																									- muonVetoPMTEnvelopeSizeZ / 2.,
																									0, 0};

            G4LogicalVolume *groundScintLV = Build_I_Muffler_LRScintillator(
                _vinylt, _tyvek, _vm2000, Tyvek_opsurf, Vikuiti_opsurf);

						// ___ 1st ground PS placement
            G4VPhysicalVolume *ground1ScintPV = new G4PVPlacement(
                G4Transform3D(*scintMother_alignMtx_groundSide, 
									PosPbBox - groundScintDisplacer),
                groundScintLV, "Envelope_MuonGroundScintillator_PV", 
								logiWorkArea, false, 18, OverlapCheck);
            
						new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_groundSide, 
									ground1ScintPV->GetTranslation() + groundPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 40, OverlapCheck);
						
						groundPMTDisplacer[0] *= -1.;
						PMTMother_alignMtx_groundSide->rotateY(180. * deg);
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_groundSide, 
									ground1ScintPV->GetTranslation() + groundPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 41, OverlapCheck);

						// ___ 2nd ground PS placement
						groundScintDisplacer[1] -= (mufflerLRScint_boxWidth_half + scint_reflector_thick) * 2; 
            G4VPhysicalVolume *ground2ScintPV = new G4PVPlacement(
                G4Transform3D(*scintMother_alignMtx_groundSide, 
									PosPbBox - groundScintDisplacer),
                groundScintLV, "Envelope_MuonGroundScintillator_PV", 
								logiWorkArea, false, 19, OverlapCheck);
            
						new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_groundSide, 
									ground2ScintPV->GetTranslation() + groundPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 42, OverlapCheck);
						
						groundPMTDisplacer[0] *= -1.;
						PMTMother_alignMtx_groundSide->rotateY(180. * deg);
            new G4PVPlacement(G4Transform3D(*PMTMother_alignMtx_groundSide, 
									ground2ScintPV->GetTranslation() + groundPMTDisplacer),
									logicPMTEnvelope, "MuVetoPMT_PV", logiWorkArea, false, 43, OverlapCheck);

        }
        if (fEnable_NeutronShield) {
            Place_I_RealNeutronShieldingAt(logiWorkArea, PosPbBox, 11);
            Place_I_RealNeutronShieldingAt(logiWorkArea,
                                           {0, 0, physTopPbBox1->GetTranslation().z()}, 12);
						Place_I_RealNeutronShieldingAt(logiInnerDetector, PosSS, 5);
						
            neutronShieldBuilt = true;
        }
    } else { // if EnableOrigGeom is false;
			 // Muon plastic scintillator test
        G4RotationMatrix *rotMtx = new G4RotationMatrix;
				rotMtx->rotateZ(90 * deg);
        rotMtx->rotateY(90 * deg);
        G4LogicalVolume *tempMuffLV =
            Build_I_Muffler_LRScintillator(_vinylt, _tyvek, _vm2000, Tyvek_opsurf, Vikuiti_opsurf);
        new G4PVPlacement(rotMtx, {}, tempMuffLV, "physTest", logiWorkArea, false, 0, OverlapCheck);
    }
    if (fNeutronMode) {
        // Neutron simulation section
        logiRockCavityVis->SetVisibility(true);

        switch (whichNShieldingConf) {
            default:
                G4cerr << "Wrong neutron shielding configuration! I will fuck up this simulation!"
                       << G4endl;
                exit(-1);
                return;
                break;
            case kNS_RealConf:
                G4cout << "Real configuration has been selected." << G4endl;
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea", logiRockCavity, false,
                                  0, OverlapCheck);
                if (!neutronShieldBuilt) {
                    Place_I_RealNeutronShieldingAt(logiWorkArea, PosPbBox, 1);
                    Place_I_RealNeutronShieldingAt(logiWorkArea,
                                                   {0, 0, physTopPbBox1->GetTranslation().z()}, 2);
                    Place_I_RealNeutronShieldingAt(logiWorkArea, PosPbBox, 4);
                    Place_I_RealNeutronShieldingAt(logiInnerDetector, {}, 3);
                    neutronShieldBuilt = true;
                }
                break;
            case kNS_Naked:
                G4cout << "Naked shielding configuration has been selected." << G4endl;
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea", logiRockCavity, false,
                                  0, OverlapCheck);
                break;
            case kNS_B4C:
                G4cout << "B4C rubber shielding configuration has been selected." << G4endl;
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea",
                                  toyRubberB4C_BoxLogical, false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {}, toyRubberB4C_BoxLogical, "physB4CRubber",
                                  logiRockCavity, false, 0, OverlapCheck);
                break;
            case kNS_B4CnPE:
                G4cout << "B4C rubber + PE 20cm shielding configuration has been selected."
                       << G4endl;
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea",
                                  toyRubberB4C_BoxLogical, false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {}, toyRubberB4C_BoxLogical, "physB4CRubber",
                                  toyPE_BoxLogical, false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {}, toyPE_BoxLogical, "physPolyEthylene", logiRockCavity,
                                  false, 0, OverlapCheck);
                break;
            case kNS_PE20:
                G4cout << "PE 20cm shielding configuration has been selected." << G4endl;
                toyPE_thickness = 20 * cm;
                toyPE_Box->SetXHalfLength(workboxX + toyPE_thickness);
                toyPE_Box->SetYHalfLength(workboxY + toyPE_thickness);
                toyPE_Box->SetZHalfLength(workboxH + toyPE_thickness);
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea", toyPE_BoxLogical,
                                  false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {}, toyPE_BoxLogical, "physPolyEthylene", logiRockCavity,
                                  false, 0, OverlapCheck);
                break;
            case kNS_PE10:
                G4cout << "PE 10cm shielding configuration has been selected." << G4endl;
                toyPE_thickness = 10 * cm;
                toyPE_Box->SetXHalfLength(workboxX + toyPE_thickness);
                toyPE_Box->SetYHalfLength(workboxY + toyPE_thickness);
                toyPE_Box->SetZHalfLength(workboxH + toyPE_thickness);
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea", toyPE_BoxLogical,
                                  false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {}, toyPE_BoxLogical, "physPolyEthylene", logiRockCavity,
                                  false, 0, OverlapCheck);
                break;
            case kNS_B4CnBAcid:
                G4cout << "B4C rubber + BoricAcid shielding configuration has been selected."
                       << G4endl;
                new G4PVPlacement(nullptr, {}, logiWorkArea, "physWorkArea",
                                  toyRubberB4C_BoxLogical, false, 0, OverlapCheck);

                if (!neutronShieldBuilt) {
                    Place_I_RealNeutronShieldingAt(logiInnerDetector, {}, 3);
                    neutronShieldBuilt = true;
                }

                new G4PVPlacement(nullptr, PosInnerDet + topBoricAcid_Pos,
                                  topBoricAcid_MotherLogical, "BoricAcid_Top_MotherPV",
                                  logiWorkArea, false, 0, OverlapCheck);

                new G4PVPlacement(nullptr, {}, toyRubberB4C_BoxLogical, "physB4CRubber",
                                  logiRockCavity, false, 0, OverlapCheck);
                break;
        }

        fCavernPhysical =
            new G4PVPlacement(nullptr, {}, logiRockCavity, "physRockCavity", logiHall, false, 0, OverlapCheck);
    } else { // if NeutronMod is false;
        fCavernPhysical = new G4PVPlacement(tube_alignMtx,    // rotation
                                            {},               // translation
                                            logiRockCavity,   // associated logical vol
                                            "physRockCavity", // name
                                            //logiHall,         // parent
                                            logiRock,         // parent
                                            false,            // no "Many"
                                            0,
																						OverlapCheck);               // copy number
				new G4PVPlacement(nullptr, {0,0,0}, logiRock, "physRock",
						              logiHall, false, 0, OverlapCheck);
        new G4PVPlacement(nullptr, {0, 0, -RockDiskThick_half}, logiRockDisk, "physRockDisk",
                          logiHall, false, 0, OverlapCheck);
        new G4PVPlacement(WA_alignMtx, {0, workboxH, 0}, logiWorkArea, "physWorkArea",
                          logiRockCavity, false, 0, OverlapCheck);
    }

    //////////////////////////////
    // Set Attributes
    //////////////////////////////
    G4cout << "Set Geometry Attributes...\n";
    logiHallVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.1));
    logiHall->SetVisAttributes(logiHallVis);
    //logiHall->SetVisAttributes(G4VisAttributes::Invisible);

    //	  IVC_Welding//////////

    G4VisAttributes *va_IVCwe = new G4VisAttributes(red);
    va_IVCwe->SetForceSolid(true);
    lv_IVCwe->SetVisAttributes(va_IVCwe);

    G4VisAttributes *va_OVCwe = new G4VisAttributes(red);
    va_OVCwe->SetForceSolid(true);
    lv_OVCwe->SetVisAttributes(va_OVCwe);

    logiRockVis = new G4VisAttributes(G4Colour(1, 0, 0, 0.5));
    logiRockVis->SetForceSolid(true);
    logiRockDiskVis = new G4VisAttributes(red);
		logiRockDiskVis->SetForceSolid(true);
    logiWorkAreaVis = new G4VisAttributes(white);
    //logiWorkAreaVis->SetVisibility(false);
    logiPbBoxVis = new G4VisAttributes(grey);
    logiPbBoxVis = new G4VisAttributes(greyl);
    logiSSOVCVis = new G4VisAttributes(cyanl);
    logiCu4Vis   = new G4VisAttributes(bluel);
    logiCu3Vis   = new G4VisAttributes(green);
    logiCu2Vis   = new G4VisAttributes(cyan);
    logiCu1Vis   = new G4VisAttributes(brownl);
    logiCuMCPVis = new G4VisAttributes(lblue);
    logiCuMCPVis->SetForceSolid(true);
    logiCuP1Vis = new G4VisAttributes(lgreen);
    logiPbP2Vis = new G4VisAttributes(grey);
    logiCuP3Vis = new G4VisAttributes(red);

    if (flagInvisible) {
        logiRock->SetVisAttributes(G4VisAttributes::Invisible);
        logiWorkArea->SetVisAttributes(G4VisAttributes::Invisible);
        logiTopPbBox->SetVisAttributes(G4VisAttributes::Invisible);
        logiPbBox->SetVisAttributes(G4VisAttributes::Invisible);
        logiSSOVCOuter->SetVisAttributes(G4VisAttributes::Invisible);
        logiCu4Outer->SetVisAttributes(G4VisAttributes::Invisible);
        logiCu3Outer->SetVisAttributes(G4VisAttributes::Invisible);
        logiCu2Outer->SetVisAttributes(G4VisAttributes::Invisible);
        logiCu1Outer->SetVisAttributes(G4VisAttributes::Invisible);

        if (flagInvisible == 999) {
            logiCuMCP->SetVisAttributes(G4VisAttributes::Invisible);
            logiCuP1->SetVisAttributes(G4VisAttributes::Invisible);
            logiPbP2->SetVisAttributes(G4VisAttributes::Invisible);
            logiCuP3->SetVisAttributes(G4VisAttributes::Invisible);
        }
    } else {
        logiCuP1Vis->SetForceSolid(true);
        logiPbP2Vis->SetForceSolid(true);
        logiCuP3Vis->SetForceSolid(true);

        logiRock->SetVisAttributes(logiRockVis);
        logiRockCavity->SetVisAttributes(logiRockCavityVis);
        logiRockDisk->SetVisAttributes(logiRockDiskVis);
        logiWorkArea->SetVisAttributes(logiWorkAreaVis);

        logiTopPbBox->SetVisAttributes(logiPbBoxVis);
        logiPbBox->SetVisAttributes(logiPbBoxVis);
        logiSSOVCVis->SetForceWireframe(true);
        logiSSOVCOuter->SetVisAttributes(logiSSOVCVis);

        logiCu4Vis->SetForceWireframe(true);
        logiCu4Outer->SetVisAttributes(logiCu4Vis);
        logiCu3Vis->SetForceWireframe(true);
        logiCu4Outer->SetVisAttributes(logiCu4Vis);
        logiCu3Vis->SetForceWireframe(true);
        logiCu3Outer->SetVisAttributes(logiCu3Vis);
        logiCu2Outer->SetVisAttributes(logiCu2Vis);
        logiCu1Outer->SetVisAttributes(logiCu1Vis);
        logiCuMCP->SetVisAttributes(logiCuMCPVis);
        logiCuP1->SetVisAttributes(logiCuP1Vis);
        logiPbP2->SetVisAttributes(logiPbP2Vis);
        logiCuP3->SetVisAttributes(logiCuP3Vis);
    }

#if G4VERSION_NUMBER < 1000
    Construct_I_SDandField();
#endif

}

void AmoreDetectorConstruction::Construct_I_SDandField() {
    //////////////////////////////
    // --- TG sensitive detector
    //////////////////////////////

    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    AmoreModuleSD *moduleSD;
    if (fEnable_Innermost && fModuleSDInfos.size() > 0) {
        moduleSD = new AmoreModuleSD("/CupDet/MDSD", fModuleSDInfos);
        SDman->AddNewDetector(moduleSD);
        for (auto nowLV : fI_CrystalLVs)
            nowLV->SetSensitiveDetector(moduleSD);
        for (auto nowLV : fI_GeWaferLVs)
            nowLV->SetSensitiveDetector(moduleSD);
        for (auto nowLV : fI_CrystalGoldFilmLVs)
            nowLV->SetSensitiveDetector(moduleSD);
        for (auto nowLV : fI_GeWaferGoldFilmLVs)
            nowLV->SetSensitiveDetector(moduleSD);
    }

    CupVetoSD *MuonScintSD = nullptr;
    auto checkAndCreateSD  = [&]() -> CupVetoSD * {
        if (MuonScintSD == nullptr) {
            MuonScintSD = new CupVetoSD("/CupDet/MuonVetoSD",
                                        2 + 4 + 4 + 4 + 4 + 2); // Top 2, SideX 4, SideY 4, MufflerSideX 4, MufflerLR 4, ground 2
            SDman->AddNewDetector(MuonScintSD);
        }
        return MuonScintSD;
    };
    if (fI_TopScint_BoxLogical != nullptr) {
        checkAndCreateSD();
        fI_TopScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
        fI_TopScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
        fI_TopScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);
    }

    if (fI_SideFBScint_BoxLogical != nullptr) {
        checkAndCreateSD();
        fI_SideFBScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
        fI_SideFBScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
        fI_SideFBScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);
    }

    if (fI_SideLRScint_BoxLogical != nullptr) {
        checkAndCreateSD();
        fI_SideLRScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
        fI_SideLRScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
        fI_SideLRScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);
    }

    if (fI_MufflerFBScint_BoxLogical != nullptr) {
        checkAndCreateSD();
        fI_MufflerFBScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
        fI_MufflerFBScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
        if (fI_MufflerFBScint_PMTTrapLogical != nullptr)
				{fI_MufflerFBScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);}
    }

    if (fI_MufflerLRScint_BoxLogical != nullptr) {
        checkAndCreateSD();
        fI_MufflerLRScint_BoxLogical->SetSensitiveDetector(MuonScintSD);
        fI_MufflerLRScint_FlatTrapLogical->SetSensitiveDetector(MuonScintSD);
        if (fI_MufflerLRScint_PMTTrapLogical != nullptr)
				{fI_MufflerLRScint_PMTTrapLogical->SetSensitiveDetector(MuonScintSD);}
    }
    if (fEnable_Scintillator && MuonScintSD == nullptr) {
        G4Exception(__PRETTY_FUNCTION__, "AmoreSD_ERROR", JustWarning,
                    "Scintillator flag was set to be enabled but one of LV for scintillators is "
                    "has not been set. SD was not configured.");
    }

    CupPMTSD *pmtSD;
    pmtSD = new CupPMTSD("/cupdet/pmt/MLCS", 44); 
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD);
    if (fI_logicPMT != nullptr) {
        fI_logicPMT->SetSensitiveDetector(pmtSD);
        fI_logicPMTVacu->SetSensitiveDetector(pmtSD);
    }
}
