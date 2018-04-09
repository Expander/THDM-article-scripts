Get["models/HSSUSY/HSSUSY_librarylink.m"];
Get["addons/SplitMSSMTower/SplitMSSMTower_librarylink.m"];
Get["models/SplitMSSM/SplitMSSM_librarylink.m"];
Get["models/THDMIIMSSMBCFull/THDMIIMSSMBCFull_librarylink.m"];
Get["models/HGTHDMIIMSSMBCFull/HGTHDMIIMSSMBCFull_librarylink.m"];
Get["addons/SplitTHDMTHDMTower/SplitTHDMTHDMTower_librarylink.m"];
Get["addons/SplitTHDMSplitTower/SplitTHDMSplitTower_librarylink.m"];
Get["models/MSSMEFTHiggs/MSSMEFTHiggs_librarylink.m"];

invalid;
Mtpole = 173.21;
MZpole = 91.1876;
GFInput = 1.16637*^-5;
mbAtmb = 4.18;
alphaSAtMZ = 0.1181;
alphaEmAtMZ = 1/127.950;

sigmaAlphaS = 0.0006;
sigmaMt = 0.98;

SMParameters = {
    alphaEmMZ -> alphaEmAtMZ, (* SMINPUTS[1] *)
    GF -> GFInput,            (* SMINPUTS[2] *)
    alphaSMZ -> alphaSAtMZ,   (* SMINPUTS[3] *)
    MZ -> MZpole,             (* SMINPUTS[4] *)
    mbmb -> mbAtmb,           (* SMINPUTS[5] *)
    Mt -> Mtpole,             (* SMINPUTS[6] *)
    Mtau -> 1.777,            (* SMINPUTS[7] *)
    Mv3 -> 0,                 (* SMINPUTS[8] *)
    MW -> 80.384,             (* SMINPUTS[9] *)
    Me -> 0.000510998902,     (* SMINPUTS[11] *)
    Mv1 -> 0,                 (* SMINPUTS[12] *)
    Mm -> 0.1056583715,       (* SMINPUTS[13] *)
    Mv2 -> 0,                 (* SMINPUTS[14] *)
    md2GeV -> 0.00475,        (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,         (* SMINPUTS[22] *)
    ms2GeV -> 0.104,          (* SMINPUTS[23] *)
    mcmc -> 1.27,             (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

LinearRange[start_, stop_, steps_] :=
    Table[start + i/steps (stop - start), {i, 0, steps}];

LogRange[start_, stop_, steps_] :=
    Module[{i, result = {}},
           For[i = 0, i <= steps, i++,
               result = AppendTo[result, Exp[Log[start] + (Log[stop] - Log[start]) i / steps]];
              ];
           result
          ];

RunHSSUSY[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ,
          ytLoops_:2, Qpole_:0, QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}},
           handle = FSHSSUSYOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   MSUSY   -> MS,
                   M1Input -> MS,
                   M2Input -> MS,
                   M3Input -> MS,
                   MuInput -> MS,
                   mAInput -> MS,
                   MEWSB   -> Mtpole,
                   AtInput -> MS/TB + Xt MS,
                   TanBeta -> TB,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   msq2 -> MS^2 IdentityMatrix[3],
                   msu2 -> {{ MS^2, 0   , 0    },
                            { 0   , MS^2, 0    },
                            { 0   , 0   , MS^2 }},
                   msd2 -> MS^2 IdentityMatrix[3],
                   msl2 -> MS^2 IdentityMatrix[3],
                   mse2 -> MS^2 IdentityMatrix[3]
               }
           ];
           Print["using handle = ", handle];
           spectrum = FSHSSUSYCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSHSSUSYCloseHandle[handle];
              $Failed
              ,
              spectrum = HSSUSY /. spectrum;
              If[calcUncerts,
                 uncerts = FSHSSUSYCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSHSSUSYCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = HSSUSY /. uncerts;
                ];
              FSHSSUSYCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunHSSUSYDeg[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ,
             MQ_, MQ3_, MU_, MU3_, MD_, ML_, ME_,
             Mi_, M3_, Mu_, mA_] :=
    Module[{handle, spectrum},
           handle = FSHSSUSYOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> 0,                (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> 0          (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   MSUSY   -> MS,
                   M1Input -> Mi,
                   M2Input -> Mi,
                   M3Input -> M3,
                   MuInput -> Mu,
                   mAInput -> mA,
                   MEWSB   -> Mtpole,
                   AtInput -> Mu/TB + Xt Sqrt[MQ3 MU3],
                   TanBeta -> TB,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> 0,
                   DeltaYt -> 0,
                   msq2 -> {{ MQ^2, 0   , 0    },
                            { 0   , MQ^2, 0    },
                            { 0   , 0   , MQ3^2 }},
                   msu2 -> {{ MU^2, 0   , 0    },
                            { 0   , MU^2, 0    },
                            { 0   , 0   , MU3^2 }},
                   msd2 -> MD^2 IdentityMatrix[3],
                   msl2 -> ML^2 IdentityMatrix[3],
                   mse2 -> ME^2 IdentityMatrix[3]
               }
           ];
           spectrum = FSHSSUSYCalculateSpectrum[handle];
           FSHSSUSYCloseHandle[handle];
           If[spectrum === $Failed,
              $Failed,
              HSSUSY /. spectrum
             ]
          ];

RunSplitMSSMTower[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ, M3_?NumericQ,
                  ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0,
                  QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}, Qi = Qin, Qm = Qmat},
           If[Qi === 0, Qi = Mi];
           If[Qm === 0, Qm = Mi];
           handle = FSSplitMSSMTowerOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   MSUSY   -> MS,
                   M1Input -> Mi,
                   M2Input -> Mi,
                   M3Input -> M3,
                   MuInput -> Mi,
                   mAInput -> MS,
                   MEWSB   -> Mtpole,
                   AtInput -> MS/TB + Xt MS,
                   Qinput -> Qi,
                   Qmatch -> Qm,
                   TanBeta -> TB,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   msq2 -> MS^2 IdentityMatrix[3],
                   msu2 -> MS^2 IdentityMatrix[3],
                   msd2 -> MS^2 IdentityMatrix[3],
                   msl2 -> MS^2 IdentityMatrix[3],
                   mse2 -> MS^2 IdentityMatrix[3]
               }
           ];
           spectrum = FSSplitMSSMTowerCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSSplitMSSMTowerCloseHandle[handle];
              $Failed
              ,
              spectrum = SplitMSSM /. spectrum;
              If[calcUncerts,
                 uncerts = FSSplitMSSMTowerCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSSplitMSSMTowerCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = SplitMSSM /. uncerts;
                ];
              FSSplitMSSMTowerCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunSplitMSSM[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ, M3_?NumericQ,
                  ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0,
                  QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}, Qi = Qin, Qm = Qmat},
           If[Qi === 0, Qi = Mi];
           If[Qm === 0, Qm = Mi];
           handle = FSSplitMSSMOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   MSUSY   -> MS,
                   M1Input -> Mi,
                   M2Input -> Mi,
                   M3Input -> M3,
                   MuInput -> Mi,
                   mAInput -> MS,
                   MEWSB   -> Mtpole,
                   AtInput -> MS/TB + Xt MS,
                   Qinput -> Qi,
                   Qmatch -> Qm,
                   TanBeta -> TB,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   msq2 -> MS^2 IdentityMatrix[3],
                   msu2 -> MS^2 IdentityMatrix[3],
                   msd2 -> MS^2 IdentityMatrix[3],
                   msl2 -> MS^2 IdentityMatrix[3],
                   mse2 -> MS^2 IdentityMatrix[3]
               }
           ];
           spectrum = FSSplitMSSMCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSSplitMSSMCloseHandle[handle];
              $Failed
              ,
              spectrum = SplitMSSM /. spectrum;
              If[calcUncerts,
                 uncerts = FSSplitMSSMCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSSplitMSSMCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = SplitMSSM /. uncerts;
                ];
              FSSplitMSSMCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunSplitMSSMDeg[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ,
                MQ_, MQ3_, MU_, MU3_, MD_, ML_, ME_, mA_] :=
    Module[{handle, spectrum, Mu = Mi, M3 = Mi},
           handle = FSSplitMSSMOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> 0,                (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> 0          (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   MSUSY   -> MS,
                   M1Input -> Mi,
                   M2Input -> Mi,
                   M3Input -> M3,
                   MuInput -> Mu,
                   mAInput -> mA,
                   MEWSB   -> Mtpole,
                   AtInput -> Mu/TB + Xt Sqrt[MQ3 MU3],
                   Qinput -> 0,
                   Qmatch -> 0,
                   TanBeta -> TB,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> 0,
                   DeltaYt -> 0,
                   msq2 -> {{ MQ^2, 0   , 0    },
                            { 0   , MQ^2, 0    },
                            { 0   , 0   , MQ3^2 }},
                   msu2 -> {{ MU^2, 0   , 0    },
                            { 0   , MU^2, 0    },
                            { 0   , 0   , MU3^2 }},
                   msd2 -> MD^2 IdentityMatrix[3],
                   msl2 -> ML^2 IdentityMatrix[3],
                   mse2 -> ME^2 IdentityMatrix[3]
               }
           ];
           spectrum = FSSplitMSSMCalculateSpectrum[handle];
           FSSplitMSSMCloseHandle[handle];
           If[spectrum === $Failed,
              $Failed,
              SplitMSSM /. spectrum
             ]
          ];

RunTHDMIIMSSMBCFull[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                    ytLoops_:2, Qpole_:0, QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}},
           handle = FSTHDMIIMSSMBCFullOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   TanBeta -> TB,
                   MSUSY -> MS,
                   MEWSB -> Mtpole,
                   MuInput -> MS,
                   M1Input -> MS,
                   M2Input -> MS,
                   MAInput -> MA,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   AeInput -> 0 MS TB IdentityMatrix[3],
                   AdInput -> 0 MS TB IdentityMatrix[3],
                   AuInput -> {
                       {0, 0, 0            },
                       {0, 0, 0            },
                       {0, 0, MS/TB + Xt MS}
                   },
                   msqInput -> MS {1,1,1},
                   msuInput -> MS {1,1,1},
                   msdInput -> MS {1,1,1},
                   mslInput -> MS {1,1,1},
                   mseInput -> MS {1,1,1}
               }
           ];
           spectrum = FSTHDMIIMSSMBCFullCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSTHDMIIMSSMBCFullCloseHandle[handle];
              $Failed
              ,
              spectrum = THDMIIMSSMBCFull /. spectrum;
              If[calcUncerts,
                 uncerts = FSTHDMIIMSSMBCFullCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSTHDMIIMSSMBCFullCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = THDMIIMSSMBCFull /. uncerts;
                ];
              FSTHDMIIMSSMBCFullCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunHGTHDMIIMSSMBCFull[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                      Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                      ytLoops_:2, Qpole_:0, QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}},
           handle = FSHGTHDMIIMSSMBCFullOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   TanBeta -> TB,
                   MSUSY -> MS,
                   MEWSB -> Mtpole,
                   MuInput -> Mu,
                   M1Input -> M12,
                   M2Input -> M12,
                   M3Input -> M3,
                   MAInput -> MA,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   AeInput -> 0 MS TB IdentityMatrix[3],
                   AdInput -> 0 MS TB IdentityMatrix[3],
                   AuInput -> {
                       {0, 0, 0            },
                       {0, 0, 0            },
                       {0, 0, Mu/TB + Xt MS}
                   },
                   msqInput -> MS {1,1,1},
                   msuInput -> MS {1,1,1},
                   msdInput -> MS {1,1,1},
                   mslInput -> MS {1,1,1},
                   mseInput -> MS {1,1,1}
               }
           ];
           spectrum = FSHGTHDMIIMSSMBCFullCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSHGTHDMIIMSSMBCFullCloseHandle[handle];
              $Failed
              ,
              spectrum = HGTHDMIIMSSMBCFull /. spectrum;
              If[calcUncerts,
                 uncerts = FSHGTHDMIIMSSMBCFullCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSHGTHDMIIMSSMBCFullCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = HGTHDMIIMSSMBCFull /. uncerts;
                ];
              FSHGTHDMIIMSSMBCFullCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunSplitTHDMTHDMTower[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                      Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                      ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}, Qi = Qin, Qm = Qmat},
           If[Qi === 0, Qi = M12];
           If[Qm === 0, Qm = M12];
           handle = FSSplitTHDMTHDMTowerOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   TanBeta -> TB,
                   MSUSY -> MS,
                   MEWSB -> Mtpole,
                   MuInput -> Mu,
                   M1Input -> M12,
                   M2Input -> M12,
                   M3Input -> M3,
                   MAInput -> MA,
                   Qinput -> Qi,
                   Qmatch -> Qm,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   AeInput -> 0 MS TB IdentityMatrix[3],
                   AdInput -> 0 MS TB IdentityMatrix[3],
                   AuInput -> {
                       {0, 0, 0            },
                       {0, 0, 0            },
                       {0, 0, Mu/TB + Xt MS}
                   },
                   msqInput -> MS {1,1,1},
                   msuInput -> MS {1,1,1},
                   msdInput -> MS {1,1,1},
                   mslInput -> MS {1,1,1},
                   mseInput -> MS {1,1,1}
               }
           ];
           spectrum = FSSplitTHDMTHDMTowerCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSSplitTHDMTHDMTowerCloseHandle[handle];
              $Failed
              ,
              spectrum = HGTHDMIIMSSMBCFull /. spectrum;
              If[calcUncerts,
                 uncerts = FSSplitTHDMTHDMTowerCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSSplitTHDMTHDMTowerCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = HGTHDMIIMSSMBCFull /. uncerts;
                ];
              FSSplitTHDMTHDMTowerCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunSplitTHDMSplitTower[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                       Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                       ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}, Qi = Qin, Qm = Qmat},
           If[Qi === 0, Qi = MA];
           If[Qm === 0, Qm = MA];
           handle = FSSplitTHDMSplitTowerOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   TanBeta -> TB,
                   MSUSY -> MS,
                   MEWSB -> Mtpole,
                   MuInput -> Mu,
                   M1Input -> M12,
                   M2Input -> M12,
                   M3Input -> M3,
                   MAInput -> MA,
                   Qinput -> Qi,
                   Qmatch -> Qm,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   AeInput -> 0 MS TB IdentityMatrix[3],
                   AdInput -> 0 MS TB IdentityMatrix[3],
                   AuInput -> {
                       {0, 0, 0            },
                       {0, 0, 0            },
                       {0, 0, Mu/TB + Xt MS}
                   },
                   msqInput -> MS {1,1,1},
                   msuInput -> MS {1,1,1},
                   msdInput -> MS {1,1,1},
                   mslInput -> MS {1,1,1},
                   mseInput -> MS {1,1,1}
               }
           ];
           spectrum = FSSplitTHDMSplitTowerCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSSplitTHDMSplitTowerCloseHandle[handle];
              $Failed
              ,
              spectrum = HGTHDMIIMSSMBCFull /. spectrum;
              If[calcUncerts,
                 uncerts = FSSplitTHDMSplitTowerCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSSplitTHDMSplitTowerCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = HGTHDMIIMSSMBCFull /. uncerts;
                ];
              FSSplitTHDMSplitTowerCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

RunMSSMEFTHiggs[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ,
             ytLoops_:2, Qpole_:0, QDR_:0, calcUncerts_:False, yt_:2] :=
    Module[{handle, spectrum, uncerts = {}},
           handle = FSMSSMEFTHiggsOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> 0,                (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> Qpole,         (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> yt,      (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   MSUSY -> MS,
                   M1Input -> MS,
                   M2Input -> MS,
                   M3Input -> MS,
                   MuInput -> MS,
                   mAInput -> MS,
                   TanBeta -> TB,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   AeInput -> 0 MS TB IdentityMatrix[3],
                   AdInput -> 0 MS TB IdentityMatrix[3],
                   AuInput -> {
                       {0, 0, 0            },
                       {0, 0, 0            },
                       {0, 0, MS/TB + Xt MS}
                   },
                   mq2Input -> MS^2 IdentityMatrix[3],
                   mu2Input -> MS^2 IdentityMatrix[3],
                   md2Input -> MS^2 IdentityMatrix[3],
                   ml2Input -> MS^2 IdentityMatrix[3],
                   me2Input -> MS^2 IdentityMatrix[3]
               }
           ];
           spectrum = FSMSSMEFTHiggsCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSMSSMEFTHiggsCloseHandle[handle];
              $Failed
              ,
              spectrum = MSSMEFTHiggs /. spectrum;
              If[calcUncerts,
                 uncerts = FSMSSMEFTHiggsCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSMSSMEFTHiggsCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = MSSMEFTHiggs /. uncerts;
                ];
              FSMSSMEFTHiggsCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

GetPar[spec_, par__] :=
    GetPar[spec, #]& /@ {par};

GetPar[spec_, par_] :=
    If[spec =!= $Failed, (par /. spec), invalid];

GetPar[spec_, par_[n__?IntegerQ]] :=
    If[spec =!= $Failed, (par /. spec)[[n]], invalid];

RunHSSUSYMh[args__] :=
    Module[{spec = RunHSSUSY[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]]]
             ]
          ];

RunHSSUSYDegMh[args__] :=
    Module[{spec = RunHSSUSYDeg[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]]]
             ]
          ];

RunSplitMSSMTowerMh[args__] :=
    Module[{spec = RunSplitMSSMTower[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]]]
             ]
          ];

RunSplitMSSMMh[args__] :=
    Module[{spec = RunSplitMSSM[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]]]
             ]
          ];

RunSplitMSSMDegMh[args__] :=
    Module[{spec = RunSplitMSSMDeg[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]]]
             ]
          ];

RunTHDMIIMSSMBCFullMh[args__] :=
    Module[{spec = RunTHDMIIMSSMBCFull[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunHGTHDMIIMSSMBCFullMh[args__] :=
    Module[{spec = RunHGTHDMIIMSSMBCFull[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunSplitTHDMTHDMTowerMh[args__] :=
    Module[{spec = RunSplitTHDMTHDMTower[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunSplitTHDMSplitTowerMh[args__] :=
    Module[{spec = RunSplitTHDMSplitTower[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunMSSMEFTHiggsMh[args__] :=
    Module[{spec = RunMSSMEFTHiggs[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunMSSMEFTHiggsDegMh[args__] :=
    Module[{spec = RunMSSMEFTHiggsDeg[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunHSSUSYUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ,
                       ytLoops_:2, Qpole_:0, QDR_:0] :=
    Module[{uncerts, MhEFT},
           uncerts = RunHSSUSY[MS, TB, Xt, ytLoops, Qpole, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunHSSUSY[MS, TB, Xt, ytLoops, Qpole, QDR, False, 1, 0], Pole[M[hh]]];
           MhYt  = GetPar[RunHSSUSY[MS, TB, Xt, ytLoops, Qpole, QDR, False, 0, 1], Pole[M[hh]]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)),
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunSplitMSSMTowerUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ, M3_?NumericQ,
                               ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0] :=
    Module[{uncerts, MhEFT},
           uncerts = RunSplitMSSMTower[MS, TB, Xt, Mi, M3, ytLoops, Qpole, Qin, Qmat, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunSplitMSSMTower[MS, TB, Xt, Mi, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 1, 0], Pole[M[hh]]];
           MhYt  = GetPar[RunSplitMSSMTower[MS, TB, Xt, Mi, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 0, 1], Pole[M[hh]]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)),
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunSplitMSSMUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ, M3_?NumericQ,
                               ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0] :=
    Module[{uncerts, MhEFT},
           uncerts = RunSplitMSSM[MS, TB, Xt, Mi, M3, ytLoops, Qpole, Qin, Qmat, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunSplitMSSM[MS, TB, Xt, Mi, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 1, 0], Pole[M[hh]]];
           MhYt  = GetPar[RunSplitMSSM[MS, TB, Xt, Mi, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 0, 1], Pole[M[hh]]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)),
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunTHDMIIMSSMBCFullUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                                 ytLoops_:2, Qpole_:0, QDR_:0] :=
    Module[{uncerts, MhEFT, MhYt},
           uncerts = RunTHDMIIMSSMBCFull[MS, TB, Xt, MA, ytLoops, Qpole, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunTHDMIIMSSMBCFull[MS, TB, Xt, MA, ytLoops, Qpole, QDR, False, 1, 0], Pole[M[hh]][1]];
           MhYt  = GetPar[RunTHDMIIMSSMBCFull[MS, TB, Xt, MA, ytLoops, Qpole, QDR, False, 0, 1], Pole[M[hh]][1]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  (Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)))[[1]],
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunHGTHDMIIMSSMBCFullUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                                   Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                                   ytLoops_:2, Qpole_:0, QDR_:0] :=
    Module[{uncerts, MhEFT, MhYt},
           uncerts = RunHGTHDMIIMSSMBCFull[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunHGTHDMIIMSSMBCFull[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, QDR, False, 1, 0], Pole[M[hh]][1]];
           MhYt  = GetPar[RunHGTHDMIIMSSMBCFull[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, QDR, False, 0, 1], Pole[M[hh]][1]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  (Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)))[[1]],
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunSplitTHDMTHDMTowerUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                                   Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                                   ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0] :=
    Module[{uncerts, MhEFT, MhYt},
           uncerts = RunSplitTHDMTHDMTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunSplitTHDMTHDMTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 1, 0], Pole[M[hh]][1]];
           MhYt  = GetPar[RunSplitTHDMTHDMTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 0, 1], Pole[M[hh]][1]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  (Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)))[[1]],
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunSplitTHDMSplitTowerUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                                   Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                                   ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0] :=
    Module[{uncerts, MhEFT, MhYt},
           uncerts = RunSplitTHDMSplitTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunSplitTHDMSplitTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 1, 0], Pole[M[hh]][1]];
           MhYt  = GetPar[RunSplitTHDMSplitTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 0, 1], Pole[M[hh]][1]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  (Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)))[[1]],
                  MhEFT,
                  MhYt
              }
             ]
          ];

RunMSSMEFTHiggsUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ,
                          ytLoops_:2, Qpole_:0, QDR_:0] :=
    Module[{uncerts, MhEFT},
           uncerts = RunMSSMEFTHiggs[MS, TB, Xt, ytLoops, Qpole, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhYt  = GetPar[RunMSSMEFTHiggs[MS, TB, Xt, ytLoops, Qpole, QDR, False, 0], Pole[M[hh]][1]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  (Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)))[[1]],
                  0,
                  MhYt
              }
             ]
          ];

IsValid[ex_]          := FreeQ[ex, invalid];
RemoveInvalid[l_List] := Select[l, IsValid];
MaxDiff[l_List]       := Max[Abs[RemoveInvalid[l]]] - Min[Abs[RemoveInvalid[l]]];
MaxDiff[{}]           := 0;
MinMax[l_List]        := { Min[Abs[RemoveInvalid[l]]], Max[Abs[RemoveInvalid[l]]] };
MinMax[{}]            := { invalid, invalid };

RunEFT[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunHSSUSYMh[MS, TB, Xt, 1, 0, MS];
           Mhyt2L = RunHSSUSYMh[MS, TB, Xt, 2, 0, MS];
           Mhyt3L = RunHSSUSYMh[MS, TB, Xt, 3, 0, MS];
           DMh = RunHSSUSYUncertainties[MS, TB, Xt, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunSplitTower[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ, M3_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunSplitMSSMTowerMh[MS, TB, Xt, Mi, M3, 1];
           Mhyt2L = RunSplitMSSMTowerMh[MS, TB, Xt, Mi, M3, 2];
           Mhyt3L = RunSplitMSSMTowerMh[MS, TB, Xt, Mi, M3, 3];
           DMh = RunSplitMSSMTowerUncertainties[MS, TB, Xt, Mi, M3, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunSplit[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, Mi_?NumericQ, M3_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunSplitMSSMMh[MS, TB, Xt, Mi, M3, 1];
           Mhyt2L = RunSplitMSSMMh[MS, TB, Xt, Mi, M3, 2];
           Mhyt3L = RunSplitMSSMMh[MS, TB, Xt, Mi, M3, 3];
           DMh = RunSplitMSSMUncertainties[MS, TB, Xt, Mi, M3, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunTHDM[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunTHDMIIMSSMBCFullMh[MS, TB, Xt, MA, 1, 0, MS];
           Mhyt2L = RunTHDMIIMSSMBCFullMh[MS, TB, Xt, MA, 2, 0, MS];
           Mhyt3L = RunTHDMIIMSSMBCFullMh[MS, TB, Xt, MA, 3, 0, MS];
           DMh = RunTHDMIIMSSMBCFullUncertainties[MS, TB, Xt, MA, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunHGTHDM[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
          Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunHGTHDMIIMSSMBCFullMh[MS, TB, Xt, MA, Mu, M12, M3, 1, 0, MS];
           Mhyt2L = RunHGTHDMIIMSSMBCFullMh[MS, TB, Xt, MA, Mu, M12, M3, 2, 0, MS];
           Mhyt3L = RunHGTHDMIIMSSMBCFullMh[MS, TB, Xt, MA, Mu, M12, M3, 3, 0, MS];
           DMh = RunHGTHDMIIMSSMBCFullUncertainties[MS, TB, Xt, MA, Mu, M12, M3, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunTHDMSplitTHDM[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                 Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunSplitTHDMTHDMTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 1];
           Mhyt2L = RunSplitTHDMTHDMTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 2];
           Mhyt3L = RunSplitTHDMTHDMTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 3];
           DMh = RunSplitTHDMTHDMTowerUncertainties[MS, TB, Xt, MA, Mu, M12, M3, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunTHDMSplitMSSM[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                 Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunSplitTHDMSplitTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 1];
           Mhyt2L = RunSplitTHDMSplitTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 2];
           Mhyt3L = RunSplitTHDMSplitTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 3];
           DMh = RunSplitTHDMSplitTowerUncertainties[MS, TB, Xt, MA, Mu, M12, M3, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

RunFSEFTHiggs[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunMSSMEFTHiggsMh[MS, TB, Xt, 1, 0, MS];
           Mhyt2L = RunMSSMEFTHiggsMh[MS, TB, Xt, 2, 0, MS];
           Mhyt3L = RunMSSMEFTHiggsMh[MS, TB, Xt, 3, 0, MS];
           DMh = RunMSSMEFTHiggsUncertainties[MS, TB, Xt, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

steps = 60;

MSstart = 500;
MSstop  = 1.0 10^16;
Xtstart = -3.5;
Xtstop  = 3.5;

(********** THDM uncertainty plot ***********)

CombineTHDMUncertainties[{_, Mh2L_?NumericQ, Mh3L_?NumericQ,
                          QpoleMin_?NumericQ, QpoleMax_?NumericQ,
                          ASMin_?NumericQ, ASMax_?NumericQ,
                          MTMin_?NumericQ, MTMax_?NumericQ,
                          MhEFT_?NumericQ, MhYt_?NumericQ}] :=
    (
        + Abs[Mh2L - Mh3L]
        + Max[Abs[Mh2L - QpoleMin], Abs[Mh2L - QpoleMax]]
        + Abs[ASMin - ASMax]
        + Abs[MTMin - MTMax]
        + Abs[Mh2L - MhYt]
        + Abs[Mh2L - MhEFT]
    );

RunTHDMIIMSSMBCFullMhDMh[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ] :=
    Module[{spec, mh, dmh},
           spec = RunTHDM[MS, TB, Xt, MA];
           If[spec === invalid, Return[invalid]];
           mh = spec[[2]];
           dmh = CombineTHDMUncertainties[spec];
           { mh, dmh }
          ];

FindTHDMPointForMh[Mh_, MS_, TBfix_, Xtfix_, MA_] :=
    Module[{point, xt = Xtfix, tb = TBfix, digits = 4, prec = 0.1, mh, dmh, F, DF, excl = 0},
           F[xMS_?NumericQ, xTB_?NumericQ, xXt_?NumericQ, xMA_?NumericQ] :=
              Module[{xres},
                     xres = RunTHDMIIMSSMBCFullMhDMh[xMS, xTB, xXt, xMA];
                     If[xres === invalid, 1000, First[xres]]
                    ];
           DF[xMS_?NumericQ, xTB_?NumericQ, xXt_?NumericQ, xMA_?NumericQ] :=
              Module[{xres},
                     xres = RunTHDMIIMSSMBCFullMhDMh[xMS, xTB, xXt, xMA];
                     If[xres === invalid, 0, Last[xres]]
                    ];
           If[!FSUncertaintyFirstRunDone,
              point = FindRoot[F[MS, TBfix, x, MA] == Mh, {x, 1}, AccuracyGoal -> digits, PrecisionGoal -> digits];
              xt = x /. point;
              If[point === {} || Abs[F[MS, TBfix, xt, MA] - Mh] > prec,
                 mh  =  F[MS, TBfix, xt, MA];
                 dmh = DF[MS, TBfix, xt, MA];
                 If[mh < Mh && mh + dmh > Mh,
                    excl = 1,
                    FSUncertaintyFirstRunDone = True
                   ];
                 (* Print["Problem: Mh = ", F[MS, TBfix, xt, MA], " for xt = ", xt, " TB = ", TBfix]; *)
                 Return[{ TB -> tb, Xt -> xt, excluded -> excl }];
                ];
             ];
           If[FSUncertaintyFirstRunDone,
              point = FindRoot[F[MS, t, Xtfix, MA] == Mh, {t, 2}, AccuracyGoal -> digits, PrecisionGoal -> digits];
              tb = t /. point;
              If[point === {} || Abs[F[MS, tb, Xtfix, MA] - Mh] > prec,
                 (* Print["   Problem: Mh = ", F[MS, tb, Xtfix, MA], " for xt = ", Xtfix, " TB = ", tb]; *)
                 tb = invalid
                ]
             ];
           { TB -> tb, Xt -> xt, excluded -> excl }
          ];

RunTHDMOverContourMh[Mh_, MS_, TBfix_, Xtfix_, MA_] :=
    Module[{point},
           point = FindTHDMPointForMh[Mh, MS, TBfix, Xtfix, MA];
           Print[{N[MS], point}];
           { MS, TB /. point, Xt /. point, MA, Sequence @@ RunTHDM[MS, TB /. point, Xt /. point, MA], excluded /. point }
          ];

ScanTHDMUncertainty[TBfix_, Xtfix_, MA_, start_:200, stop_:1.0 10^16, steps_:500] :=
    Module[{res},
           Off[FindRoot::lstol];
           Off[FSTHDMIIMSSMBCFullCalculateSpectrum::error];
           FSUncertaintyFirstRunDone = False;
           res =  RunTHDMOverContourMh[125, N[#], TBfix, Xtfix, MA]& /@ LogRange[start, stop, steps];
           Export["THDMIIMSSMBCFull_uncertainty_MS_MA-" <> ToString[MA] <> ".dat", res, "Table"];
           res
          ];

(* ScanTHDMUncertainty[20, 0, 800] *)

(********** HSSUSY scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

ScanHSSUSYMS[TB_, Xt_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = {#, TB, Xt, Sequence @@ RunEFT[#, TB, Xt]}& /@ LogRange[start, stop, steps];
           Export["HSSUSY_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <> ".dat", res, "Table"];
           res
          ];

(* ScanHSSUSYMS[#, N@Sqrt[6]]& /@ {2, 10, 20, 50} *)
(* ScanHSSUSYMS[#, N@Sqrt[6], 500, 1.0 10^4]& /@ {2, 10, 20, 50} *)

(********** HSSUSY scenario 2: TB = 2, MS = 2 TeV **********)

ScanHSSUSYXt[TB_, MS_, start_:-4, stop_:4, steps_:60] :=
    Module[{res},
           res = {MS, TB, N[#], Sequence @@ RunEFT[MS, TB, #]}& /@ LinearRange[start, stop, steps];
           Export["HSSUSY_Xt_TB-" <> ToString[TB] <> "_MS-" <> ToString[MS] <> ".dat", res, "Table"];
           res
          ];

(* ScanHSSUSYXt[20, 2000]; *)

(********** FlexibleEFTHiggs scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

ScanFlexibleEFTHiggsMS[TB_, Xt_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = {#, TB, Xt, Sequence @@ RunFSEFTHiggs[#, TB, Xt]}& /@ LogRange[start, stop, steps];
           Export["FlexibleEFTHiggs_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <> ".dat", res, "Table"];
           res
          ];

(* ScanFlexibleEFTHiggsMS[#, N@Sqrt[6]]& /@ {2, 10, 20, 50} *)

(********** FlexibleEFTHiggs scenario 2: TB = 2, MS = 2 TeV **********)

ScanFlexibleEFTHiggsXt[TB_, MS_, start_:-4, stop_:4, steps_:60] :=
    Module[{res},
           res = {MS, TB, N[#], Sequence @@ RunFSEFTHiggs[MS, TB, #]}& /@ LinearRange[start, stop, steps];
           Export["FlexibleEFTHiggs_Xt_TB-" <> ToString[TB] <> "_MS-" <> ToString[MS] <> ".dat", res, "Table"];
           res
          ];

(* ScanFlexibleEFTHiggsXt[20, 2000]; *)

(********** SplitMSSMTower scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

ScanSplitMSSMTowerMS[TB_, Xt_, Mi_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = {#, Mi, Mi, Mi, TB, Xt, Sequence @@ RunSplitTower[#, TB, Xt, Mi, Mi]}& /@ LogRange[start, stop, steps];
           Export["SplitMSSMTower_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <>
                  "_Mi-" <> ToString[Mi] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitMSSMTowerMS[#, N@Sqrt[6], 2000]& /@ {2, 10, 20, 50} *)
(* ScanSplitMSSMTowerMS[#, N@Sqrt[6], 2000, 500, 1.0 10^4]& /@ {2, 10, 20, 50} *)

(********** SplitMSSMTower scenario 2: TB = 10, MS = 5 TeV **********)

ScanSplitMSSMTowerXt[TB_, MS_, Mi_, start_:-4, stop_:4, steps_:60] :=
    Module[{res},
           res = {MS, Mi, Mi, Mi, TB, N[#], Sequence @@ RunSplitTower[MS, TB, #, Mi, Mi]}& /@ LinearRange[start, stop, steps];
           Export["SplitMSSMTower_Xt_TB-" <> ToString[TB] <> "_MS-" <> ToString[MS] <>
                  "_Mi-" <> ToString[Mi] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitMSSMTowerXt[10, 5000, 2000]; *)

(********** SplitMSSM scenario 1: TB = 2, 10, 20, 50, Xt = 0, MS = 10 TeV, 15 TeV **********)

ScanSplitMSSMMi[TB_, Xt_, MS_, M3fac_, start_:500, stop_:3000, steps_:60] :=
    Module[{res},
           res = {MS, N[#], N[#], N[#] M3fac, TB, Xt, Sequence @@ RunSplit[MS, TB, Xt, N[#], N[#] M3fac]}& /@ LogRange[start, stop, steps];
           Export["SplitMSSM_Mi_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <>
                  "_MS-" <> ToString[MS] <> "_M3-MiX" <> ToString[M3fac] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitMSSMMi[#, 0, 10000, 1]& /@ {2, 10, 20, 50} *)
(* ScanSplitMSSMMi[#, 0, 10000, 2]& /@ {2, 10, 20, 50} *)
(* ScanSplitMSSMMi[#, 0, 15000, 1]& /@ {2, 10, 20, 50} *)
(* ScanSplitMSSMMi[#, 0, 15000, 2]& /@ {2, 10, 20, 50} *)

(********** HSSUSY non-degenerate masses: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

HSSUSYDegVary[MS_, TB_, Xt_] :=
    Module[{data, tuples},
           tuples = Tuples[{MS/2, 2 MS}, 11];
           data = RunHSSUSYDegMh[MS, TB, Xt, Sequence @@ #]& /@ tuples;
           {MS, TB, Xt, Sequence @@ MinMax[data]}
          ];

ScanHSSUSYDeg[TB_, Xt_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = HSSUSYDegVary[#, TB, Xt]& /@ LogRange[start, stop, steps];
           Export["HSSUSY_degenerate_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <> ".dat", res, "Table"];
           res
          ];

(* ScanHSSUSYDeg[#, N@Sqrt[6]]& /@ {2, 10, 20, 50}; *)

(********** SplitMSSM non-degenerate masses: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

SplitMSSMDegVary[MS_, TB_, Xt_, Mlow_] :=
    Module[{data, tuples},
           tuples = Tuples[{
               {Mlow/2, 2 Mlow}, (* Mi = M3 = Mu *)
               {MS/2, 2 MS}, (* MQ *)
               {MS/2, 2 MS}, (* MQ3 *)
               {MS/2, 2 MS}, (* MU *)
               {MS/2, 2 MS}, (* MU3 *)
               {MS/2, 2 MS}, (* MD *)
               {MS/2, 2 MS}, (* ML *)
               {MS/2, 2 MS}, (* ME *)
               {MS/2, 2 MS}  (* mA *)
           }];
           data = RunSplitMSSMDegMh[MS, TB, Xt, Sequence @@ #]& /@ tuples;
           {MS, TB, Xt, Sequence @@ MinMax[data]}
          ];

ScanSplitMSSMDeg[TB_, Xt_, Mlow_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = SplitMSSMDegVary[#, TB, Xt, Mlow]& /@ LogRange[start, stop, steps];
           Export["SplitMSSM_degenerate_MS_TB-" <> ToString[TB] <> "_Xt-" <>
                  ToString[Xt] <> "_Mlow-" <> ToString[Mlow] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitMSSMDeg[#, N@Sqrt[6], 1500]& /@ {2, 10, 20, 50}; *)

(********** THDM degenerate masses: TB = [2,50], MS = [1000, 10^16], Xt = ? **********)

ScanTHDMIIMSSMBCFullMSTB[Xt_, MA_, MSstart_:1000, MSstop_:1.0 10^16, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LogRange[MSstart, MSstop, steps], LinearRange[TBstart, TBstop, steps]}];
           res = {Sequence @@ N[#], MA, Xt, Sequence @@ RunTHDM[Sequence @@ #, Xt, MA]}& /@ tuples;
           Export["THDMIIMSSMBCFull_TB_MS_Xt-" <> ToString[Xt] <> "_MA-" <> ToString[MA] <> ".dat", res, "Table"];
           res
          ];

(* ScanTHDMIIMSSMBCFullMSTB[0, 400]; *)
(* ScanTHDMIIMSSMBCFullMSTB[0, 800]; *)
(* ScanTHDMIIMSSMBCFullMSTB[N@Sqrt[6], 400]; *)
(* ScanTHDMIIMSSMBCFullMSTB[N@Sqrt[6], 800]; *)

(********** THDM degenerate masses: TB = [2,50], MA = [100, 500], Xt = ? **********)

ScanTHDMIIMSSMBCFullTBMA[Xt_, MS_, MAstart_:100, MAstop_:500, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[TBstart, TBstop, steps], LinearRange[MAstart, MAstop, steps]}];
           res = {MS, Sequence @@ N[#], Xt, Sequence @@ RunTHDM[MS, #[[1]], Xt, #[[2]]]}& /@ tuples;
           Export["THDMIIMSSMBCFull_TB_MA_Xt-" <> ToString[Xt] <> "_MS-" <> ToString[MS] <> ".dat", res, "Table"];
           res
          ];

(* ScanTHDMIIMSSMBCFullTBMA[0, 5000]; *)
(* ScanTHDMIIMSSMBCFullTBMA[0, 10^4]; *)
(* ScanTHDMIIMSSMBCFullTBMA[0, 5 10^4]; *)
(* ScanTHDMIIMSSMBCFullTBMA[N@Sqrt[6], 5000]; *)
(* ScanTHDMIIMSSMBCFullTBMA[N@Sqrt[6], 10^4]; *)
(* ScanTHDMIIMSSMBCFullTBMA[N@Sqrt[6], 5 10^4]; *)

(********** THDM degenerate masses: MA = [100, 500], MS = [1000, 10^16], Xt = ? **********)

ScanTHDMIIMSSMBCFullMSMA[Xt_, TB_, MSstart_:1000, MSstop_:1.0 10^16, MAstart_:100, MAstop_:500, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], TB, N[#[[1]]], Xt, Sequence @@ RunTHDM[#[[2]], TB, Xt, #[[1]]]}& /@ tuples;
           Export["THDMIIMSSMBCFull_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <> ".dat", res, "Table"];
           res
          ];

(* ScanTHDMIIMSSMBCFullMSMA[0, #]& /@ {2, 10, 20, 50}; *)
(* ScanTHDMIIMSSMBCFullMSMA[N@Sqrt[6], #]& /@ {2, 10, 20, 50}; *)

(********** HGTHDM degenerate masses: TB = [2,50], MS = [1000, 10^16], Xt = ? **********)

ScanHGTHDMIIMSSMBCFullMSTB[Xt_, MA_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LogRange[MSstart, MSstop, steps], LinearRange[TBstart, TBstop, steps]}];
           res = {N[#[[1]]], Mui, Mui, Mui, N[#[[2]]], MA, Xt,
                  Sequence @@ RunHGTHDM[Sequence @@ #, Xt, MA, Mui, Mui, Mui]}& /@ tuples;
           Export["HGTHDMIIMSSMBCFull_TB_MS_Xt-" <> ToString[Xt] <> "_MA-" <> ToString[MA] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanHGTHDMIIMSSMBCFullMSTB[0, 400, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullMSTB[0, 800, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullMSTB[N@Sqrt[6], 400, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullMSTB[N@Sqrt[6], 800, 2000]; *)

(********** HGTHDM degenerate masses: TB = [2,50], MA = [100, 500], Xt = ? **********)

ScanHGTHDMIIMSSMBCFullTBMA[Xt_, MS_, Mui_, MAstart_:100, MAstop_:500, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[TBstart, TBstop, steps], LinearRange[MAstart, MAstop, steps]}];
           res = {MS, Mui, Mui, Mui, Sequence @@ N[#], Xt,
                  Sequence @@ RunHGTHDM[MS, #[[1]], Xt, #[[2]], Mui, Mui, Mui]}& /@ tuples;
           Export["HGTHDMIIMSSMBCFull_TB_MA_Xt-" <> ToString[Xt] <> "_MS-" <> ToString[MS] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanHGTHDMIIMSSMBCFullTBMA[0, 5000, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullTBMA[0, 10^4, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullTBMA[0, 5 10^4, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullTBMA[N@Sqrt[6], 5000, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullTBMA[N@Sqrt[6], 10^4, 2000]; *)
(* ScanHGTHDMIIMSSMBCFullTBMA[N@Sqrt[6], 5 10^4, 2000]; *)

(********** HGTHDM degenerate masses: MA = [100, 500], MS = [1000, 10^16], Xt = ? **********)

ScanHGTHDMIIMSSMBCFullMSMA[Xt_, TB_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, MAstart_:100, MAstop_:500, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], Mui, Mui, Mui, TB, N[#[[1]]], Xt,
                  Sequence @@ RunHGTHDM[#[[2]], TB, Xt, #[[1]], Mui, Mui, Mui]}& /@ tuples;
           Export["HGTHDMIIMSSMBCFull_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanHGTHDMIIMSSMBCFullMSMA[0, #, 2000]& /@ {2, 10, 20, 50}; *)
(* ScanHGTHDMIIMSSMBCFullMSMA[N@Sqrt[6], #, 2000]& /@ {2, 10, 20, 50}; *)

(********** THDM+split -> THDM tower degenerate masses: TB = [2,50], MS = [1000, 10^16], Xt = ? **********)

ScanSplitTHDMTHDMTowerMSTB[Xt_, MA_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LogRange[MSstart, MSstop, steps], LinearRange[TBstart, TBstop, steps]}];
           res = {N[#[[1]]], Mui, Mui, Mui, N[#[[2]]], MA, Xt,
                  Sequence @@ RunTHDMSplitTHDM[Sequence @@ #, Xt, MA, Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMTHDMTower_TB_MS_Xt-" <> ToString[Xt] <> "_MA-" <> ToString[MA] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitTHDMTHDMTowerMSTB[0, 400, 2000]; *)
(* ScanSplitTHDMTHDMTowerMSTB[0, 800, 2000]; *)
(* ScanSplitTHDMTHDMTowerMSTB[N@Sqrt[6], 400, 2000]; *)
(* ScanSplitTHDMTHDMTowerMSTB[N@Sqrt[6], 800, 2000]; *)

(********** THDM+split -> THDM degenerate masses: TB = [2,50], MA = [100, 500], Xt = ? **********)

ScanSplitTHDMTHDMTowerTBMA[Xt_, MS_, Mui_, MAstart_:100, MAstop_:500, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[TBstart, TBstop, steps], LinearRange[MAstart, MAstop, steps]}];
           res = {MS, Mui, Mui, Mui, Sequence @@ N[#], Xt,
                  Sequence @@ RunTHDMSplitTHDM[MS, #[[1]], Xt, #[[2]], Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMTHDMTower_TB_MA_Xt-" <> ToString[Xt] <> "_MS-" <> ToString[MS] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitTHDMTHDMTowerTBMA[0, 5000, 2000]; *)
(* ScanSplitTHDMTHDMTowerTBMA[0, 10^4, 2000]; *)
(* ScanSplitTHDMTHDMTowerTBMA[0, 5 10^4, 2000]; *)
(* ScanSplitTHDMTHDMTowerTBMA[N@Sqrt[6], 5000, 2000]; *)
(* ScanSplitTHDMTHDMTowerTBMA[N@Sqrt[6], 10^4, 2000]; *)
(* ScanSplitTHDMTHDMTowerTBMA[N@Sqrt[6], 5 10^4, 2000]; *)

(********** THDM+split -> THDM tower degenerate masses: MA = [100, 500], MS = [1000, 10^16], Xt = ? **********)

ScanSplitTHDMTHDMTowerMSMA[Xt_, TB_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, MAstart_:100, MAstop_:500, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], Mui, Mui, Mui, TB, N[#[[1]]], Xt,
                  Sequence @@ RunTHDMSplitTHDM[#[[2]], TB, Xt, #[[1]], Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMTHDMTower_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitTHDMTHDMTowerMSMA[0, #, 2000]& /@ {2, 10, 20, 50}; *)
(* ScanSplitTHDMTHDMTowerMSMA[N@Sqrt[6], #, 2000]& /@ {2, 10, 20, 50}; *)

(********** THDM+split -> SM+split tower degenerate masses: TB = [2,50], MS = [1000, 10^16], Xt = ? **********)

ScanSplitTHDMSplitTowerMSTB[Xt_, MA_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LogRange[MSstart, MSstop, steps], LinearRange[TBstart, TBstop, steps]}];
           res = {N[#[[1]]], Mui, Mui, Mui, N[#[[2]]], MA, Xt,
                  Sequence @@ RunTHDMSplitMSSM[Sequence @@ #, Xt, MA, Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMSplitTower_TB_MS_Xt-" <> ToString[Xt] <> "_MA-" <> ToString[MA] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitTHDMSplitTowerMSTB[0, 2000, 1000]; *)
(* ScanSplitTHDMSplitTowerMSTB[0, 5000, 1000]; *)
(* ScanSplitTHDMSplitTowerMSTB[N@Sqrt[6], 2000, 1000]; *)
(* ScanSplitTHDMSplitTowerMSTB[N@Sqrt[6], 5000, 1000]; *)

(********** THDM+split -> SM+split tower degenerate masses: TB = [2,50], MA = [2000, 10000], Xt = ? **********)

ScanSplitTHDMSplitTowerTBMA[Xt_, MS_, Mui_, MAstart_:2000, MAstop_:10000, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[TBstart, TBstop, steps], LinearRange[MAstart, MAstop, steps]}];
           res = {MS, Mui, Mui, Mui, Sequence @@ N[#], Xt,
                  Sequence @@ RunTHDMSplitMSSM[MS, #[[1]], Xt, #[[2]], Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMSplitTower_TB_MA_Xt-" <> ToString[Xt] <> "_MS-" <> ToString[MS] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitTHDMSplitTowerTBMA[0, 5000, 1000]; *)
(* ScanSplitTHDMSplitTowerTBMA[0, 10^4, 1000]; *)
(* ScanSplitTHDMSplitTowerTBMA[0, 5 10^4, 1000]; *)
(* ScanSplitTHDMSplitTowerTBMA[N@Sqrt[6], 5000, 1000]; *)
(* ScanSplitTHDMSplitTowerTBMA[N@Sqrt[6], 10^4, 1000]; *)
(* ScanSplitTHDMSplitTowerTBMA[N@Sqrt[6], 5 10^4, 1000]; *)

(********** THDM+split -> SM+split tower degenerate masses: MA = [2000, 10000], MS = [1000, 10^16], Xt = ? **********)

ScanSplitTHDMSplitTowerMSMA[Xt_, TB_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, MAstart_:2000, MAstop_:10000, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], Mui, Mui, Mui, TB, N[#[[1]]], Xt,
                  Sequence @@ RunTHDMSplitMSSM[#[[2]], TB, Xt, #[[1]], Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMSplitTower_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitTHDMSplitTowerMSMA[0, #, 1000]& /@ {2, 10, 20, 50}; *)
(* ScanSplitTHDMSplitTowerMSMA[N@Sqrt[6], #, 1000]& /@ {2, 10, 20, 50}; *)
