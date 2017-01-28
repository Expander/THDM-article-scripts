Get["models/HSSUSY/HSSUSY_librarylink.m"];

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
                   precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
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
           spectrum = FSHSSUSYCalculateSpectrum[handle];
           If[calcUncerts,
              uncerts = FSHSSUSYCalculateUncertainties[handle];
             ];
           FSHSSUSYCloseHandle[handle];
           If[calcUncerts, uncerts, spectrum]
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

RunHSSUSYQDR[args__] :=
    Module[{spec},
           spec = RunHSSUSY[args];
           If[spec === $Failed,
              { invalid, invalid, invalid, invalid }
              ,
              {
                  GetPar[spec, Pole[M[hh]]],
                 -GetPar[spec, Yu[3,3]],   (* yt^SM(QDR)     *)
                  GetPar[spec, \[Lambda]], (* lambda^SM(QDR) *)
                  GetPar[spec, g3]         (* g3^SM(QDR)     *)
              }
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
                  Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)),
                  Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)),
                  Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)),
                  MhEFT,
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
           MhEFT  = RunHSSUSYMh[MS, TB, Xt, 2, 0, MS];
           DMh = RunHSSUSYUncertainties[MS, TB, Xt, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

steps = 60;

MSstart = 500;
MSstop  = 1.0 10^18;
Xtstart = -3.5;
Xtstop  = 3.5;

(********** HSSUSY scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

ScanHSSUSYMS[TB_, Xt_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = {#, TB, Xt, Sequence @@ RunEFT[#, TB, Xt]}& /@ LogRange[start, stop, steps];
           Export["HSSUSY_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <> ".dat", res, "Table"];
           res
          ];

ScanHSSUSYMS[#, N@Sqrt[6]]& /@ {2, 10, 20, 50}
