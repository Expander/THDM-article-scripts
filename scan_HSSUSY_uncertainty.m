Get["models/HSSUSY/HSSUSY_librarylink.m"];
Get["model_files/HSSUSY/HSSUSY_uncertainty_estimate.m"];

(* To run this script, please first apply the following patch:

diff --git a/model_files/HSSUSY/HSSUSY_uncertainty_estimate.m b/model_files/HSSUSY/HSSUSY_uncertainty_estimate.m
index bca32d3d..fd127971 100644
--- a/model_files/HSSUSY/HSSUSY_uncertainty_estimate.m
+++ b/model_files/HSSUSY/HSSUSY_uncertainty_estimate.m
@@ -174,5 +174,5 @@ CalcHSSUSYDMh[args__] :=
            DMhSUSY = Max[Abs[Max[varyQmatch] - Mh],
                          Abs[Min[varyQmatch] - Mh]] +
                      Abs[Mh - MhYtMSSM];
-           { Mh0, DMhSM + DMhEFT + DMhSUSY }
+           { Mh0, DMhSM, DMhEFT, DMhSUSY }
           ];

*)

Mtpole = 173.34;
AS = 0.1184;

CalcDMh[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MG_?NumericQ, mQ3_?NumericQ, mU3_?NumericQ] :=
    CalcHSSUSYDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            poleMassLoopOrder -> 2,
            ewsbLoopOrder -> 2,
            betaFunctionLoopOrder -> 3,
            thresholdCorrectionsLoopOrder -> 2,
            thresholdCorrections -> 122111121
        },
        fsModelParameters -> {
            TanBeta -> TB,
            MEWSB -> Mtpole,
            MSUSY -> MS,
            M1Input -> MS,
            M2Input -> MS,
            M3Input -> MG,
            MuInput -> MS,
            mAInput -> MS,
            AtInput -> (Xt + 1/TB) * MS,
            msq2 -> DiagonalMatrix[{MS^2, MS^2, mQ3^2}],
            msu2 -> DiagonalMatrix[{MS^2, MS^2, mU3^2}],
            msd2 -> MS^2 IdentityMatrix[3],
            msl2 -> MS^2 IdentityMatrix[3],
            mse2 -> MS^2 IdentityMatrix[3],
            LambdaLoopOrder -> 2,
            TwoLoopAtAs -> 1,
            TwoLoopAbAs -> 1,
            TwoLoopAtAb -> 1,
            TwoLoopAtauAtau -> 1,
            TwoLoopAtAt -> 1
        },
        fsSMParameters -> {
            alphaSMZ -> AS,
            Mt -> Mtpole
        }
   ];

CalcMh[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MG_?NumericQ, mQ3_?NumericQ, mU3_?NumericQ, Mtp_:Mtpole, as_:AS] :=
    Module[{handle, spec},
           handle = FSHSSUSYOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,
                   calculateStandardModelMasses -> 1,
                   poleMassLoopOrder -> 2,
                   ewsbLoopOrder -> 2,
                   betaFunctionLoopOrder -> 3,
                   thresholdCorrectionsLoopOrder -> 2,
                   thresholdCorrections -> 122111121
               },
               fsModelParameters -> {
                   TanBeta -> TB,
                   MEWSB -> Mtpole,
                   MSUSY -> MS,
                   M1Input -> MS,
                   M2Input -> MS,
                   M3Input -> MG,
                   MuInput -> MS,
                   mAInput -> MS,
                   AtInput -> (Xt + 1/TB) * MS,
                   msq2 -> DiagonalMatrix[{MS^2, MS^2, mQ3^2}],
                   msu2 -> DiagonalMatrix[{MS^2, MS^2, mU3^2}],
                   msd2 -> MS^2 IdentityMatrix[3],
                   msl2 -> MS^2 IdentityMatrix[3],
                   mse2 -> MS^2 IdentityMatrix[3],
                   LambdaLoopOrder -> 2,
                   TwoLoopAtAs -> 1,
                   TwoLoopAbAs -> 1,
                   TwoLoopAtAb -> 1,
                   TwoLoopAtauAtau -> 1,
                   TwoLoopAtAt -> 1
               },
               fsSMParameters -> {
                   alphaSMZ -> as,
                   Mt -> Mtp
               }
           ];
           spec = FSHSSUSYCalculateSpectrum[handle];
           FSHSSUSYCloseHandle[handle];
           If[spec === $Failed, $Failed,
              Pole[M[hh]] /. (HSSUSY /. spec)]
    ];

CalcDMhParam[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MG_?NumericQ, mQ3_?NumericQ, mU3_?NumericQ] :=
    Module[{Mh, MhMtp, MhMtm, MhASp, MhASm},
           Mh    = CalcMh[MS, TB, Xt, MG, mQ3, mU3];
           MhMtp = CalcMh[MS, TB, Xt, MG, mQ3, mU3, Mtpole + 0.98];
           MhMtm = CalcMh[MS, TB, Xt, MG, mQ3, mU3, Mtpole - 0.98];
           MhASp = CalcMh[MS, TB, Xt, MG, mQ3, mU3, Mtpole, AS + 6 10^-4];
           MhASm = CalcMh[MS, TB, Xt, MG, mQ3, mU3, Mtpole, AS - 6 10^-4];
           (
               Abs[Max[{Mh, MhMtp, MhMtm}] - Min[{Mh, MhMtp, MhMtm}]] +
               Abs[Max[{Mh, MhASp, MhASm}] - Min[{Mh, MhASp, MhASm}]]
           ) / 2
          ];

CalcDMhFixedTB[MS_, TB_, Xt_, MG_, mQ3_, mU3_, Mh_] :=
    Module[{xt, point},
           point = FindRoot[CalcMh[MS,TB,xt,MG,mQ3,mU3] == Mh, {xt, 1, 0, Xt},
                            AccuracyGoal -> 4, PrecisionGoal -> 4];
           xt = If[point =!= {}, xt /. point, Xt];
           xt = If[NumericQ[xt], xt, Xt];
           { MS, TB, xt, Sequence @@ CalcDMh[MS, TB, xt, MG, mQ3, mU3], CalcDMhParam[MS, TB, xt, MG, mQ3, mU3] }
          ];

CalcDMhFixedXt[MS_, TB_, Xt_, MG_, mQ3_, mU3_, Mh_] :=
    Module[{tb, point},
           point = FindRoot[CalcMh[MS,tb,Xt,MG,mQ3,mU3] == Mh, {tb, 2, 1, TB},
                            AccuracyGoal -> 4, PrecisionGoal -> 4];
           tb = If[point =!= {}, tb /. point, TB];
           tb = If[NumericQ[tb], tb, TB];
           { MS, tb, Xt, Sequence @@ CalcDMh[MS, tb, Xt, MG, mQ3, mU3], CalcDMhParam[MS, tb, Xt, MG, mQ3, mU3]}
          ];

LaunchKernels[];
DistributeDefinitions[CalcDMhFixedTB, CalcDMhFixedXt];

(* mU3 = mQ3 = MG = MS *)
dataXt = ParallelMap[CalcDMhFixedTB[N[#], 20, N[Sqrt[6]], N[#], N[#], N[#], 125]&, LogRange[1000  , 2 10^4, 60]];
dataTB = ParallelMap[CalcDMhFixedXt[N[#], 20,          0, N[#], N[#], N[#], 125]&, LogRange[1 10^4, 10^10 , 60]];

Export["HSSUSY_uncertainty_plot_vary_Xt.dat", dataXt];
Export["HSSUSY_uncertainty_plot_vary_TB.dat", dataTB];
Export["HSSUSY_uncertainty_plot.dat", Join[dataXt, dataTB]];

(* mU3 = mQ3 = MG/2 = MS *)
dataXt = ParallelMap[CalcDMhFixedTB[N[#], 20, N[Sqrt[6]], N[2 #], N[#], N[#], 125]&, LogRange[1000  , 2 10^4, 60]];
dataTB = ParallelMap[CalcDMhFixedXt[N[#], 20,          0, N[2 #], N[#], N[#], 125]&, LogRange[1 10^4, 10^10 , 60]];

Export["HSSUSY_uncertainty_plot_vary_Xt_MG-2MS.dat", dataXt];
Export["HSSUSY_uncertainty_plot_vary_TB_MG-2MS.dat", dataTB];
Export["HSSUSY_uncertainty_plot_MG-2MS.dat", Join[dataXt, dataTB]];

(* mU3 = mQ3/1.5 = MG = MS *)
dataXt = ParallelMap[CalcDMhFixedTB[N[#], 20, N[Sqrt[6]], N[#], N[1.5 #], N[#], 125]&, LogRange[1000  , 2 10^4, 60]];
dataTB = ParallelMap[CalcDMhFixedXt[N[#], 20,          0, N[#], N[1.5 #], N[#], 125]&, LogRange[1 10^4, 10^10 , 60]];

Export["HSSUSY_uncertainty_plot_vary_Xt_mQ3-1.5MS.dat", dataXt];
Export["HSSUSY_uncertainty_plot_vary_TB_mQ3-1.5MS.dat", dataTB];
Export["HSSUSY_uncertainty_plot_mQ3-1.5MS.dat", Join[dataXt, dataTB]];

(* mU3/1.5 = mQ3 = MG = MS *)
dataXt = ParallelMap[CalcDMhFixedTB[N[#], 20, N[Sqrt[6]], N[#], N[#], N[1.5 #], 125]&, LogRange[1000  , 2 10^4, 60]];
dataTB = ParallelMap[CalcDMhFixedXt[N[#], 20,          0, N[#], N[#], N[1.5 #], 125]&, LogRange[1 10^4, 10^10 , 60]];

Export["HSSUSY_uncertainty_plot_vary_Xt_mU3-1.5MS.dat", dataXt];
Export["HSSUSY_uncertainty_plot_vary_TB_mU3-1.5MS.dat", dataTB];
Export["HSSUSY_uncertainty_plot_mU3-1.5MS.dat", Join[dataXt, dataTB]];

(* mU3/0.8 = mQ3/1.2 = MG = MS *)
dataXt = ParallelMap[CalcDMhFixedTB[N[#], 20, N[Sqrt[6]], N[#], N[1.2 #], N[0.8 #], 125]&, LogRange[1000  , 2 10^4, 60]];
dataTB = ParallelMap[CalcDMhFixedXt[N[#], 20,          0, N[#], N[1.2 #], N[0.8 #], 125]&, LogRange[1 10^4, 10^10 , 60]];

Export["HSSUSY_uncertainty_plot_vary_Xt_mQ3-1.2MS_mU3-0.8MS.dat", dataXt];
Export["HSSUSY_uncertainty_plot_vary_TB_mQ3-1.2MS_mU3-0.8MS.dat", dataTB];
Export["HSSUSY_uncertainty_plot_mQ3-1.2MS_mU3-0.8MS.dat", Join[dataXt, dataTB]];

(* mU3/0.8 = mQ3/1.2 = MG/1.8 = MS *)
dataXt = ParallelMap[CalcDMhFixedTB[N[#], 20, N[Sqrt[6]], N[1.8 #], N[1.2 #], N[0.8 #], 125]&, LogRange[1000  , 2 10^4, 60]];
dataTB = ParallelMap[CalcDMhFixedXt[N[#], 20,          0, N[1.8 #], N[1.2 #], N[0.8 #], 125]&, LogRange[1 10^4, 10^10 , 60]];

Export["HSSUSY_uncertainty_plot_vary_Xt_MG-1.8MS_mQ3-1.2MS_mU3-0.8MS.dat", dataXt];
Export["HSSUSY_uncertainty_plot_vary_TB_MG-1.8MS_mQ3-1.2MS_mU3-0.8MS.dat", dataTB];
Export["HSSUSY_uncertainty_plot_MG-1.8MS_mQ3-1.2MS_mU3-0.8MS.dat", Join[dataXt, dataTB]];
