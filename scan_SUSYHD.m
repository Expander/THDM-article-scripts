Needs["SUSYHD`"];

Mtpole = 173.21;
alphaSAtMZ = 0.1181;

LinearRange[start_, stop_, steps_] :=
    Table[start + i/steps (stop - start), {i, 0, steps}];

LogRange[start_, stop_, steps_] :=
    Module[{i, result = {}},
           For[i = 0, i <= steps, i++,
               result = AppendTo[result, Exp[Log[start] + (Log[stop] - Log[start]) i / steps]];
              ];
           result
          ];

RunSUSYHD[MS_, TB_, Xt_] :=
    Module[{gauginoM1 = MS, gauginoM2 = MS, gauginoM3 = MS, higgsMu = MS, higgsAt, mq3=MS,
            mu3=MS, md3=MS, mq2=MS, mu2=MS, md2=MS, mq1=MS, mu1=MS, md1=MS,
            ml3=MS, me3=MS, ml2=MS, me2=MS, ml1=MS, me1=MS, mA=MS,
            mass, dmass },

           SetSMparameters[Mtpole, alphaSAtMZ];
           higgsAt=(higgsMu/TB + Xt MS);

           mass = MHiggs[{TB, gauginoM1, gauginoM2, gauginoM3,
                          higgsMu, higgsAt, mq3, mu3, md3, mq2, mu2,
                          md2, mq1, mu1, md1, ml3, me3, ml2, me2, ml1,
                          me1, mA}, Rscale->MS, scheme->"DRbar",
                          hiOrd->{1,1,1,0}, numerical->True,
                          split->False];

           dmass= \[CapitalDelta]MHiggs[{TB, gauginoM1, gauginoM2,
                                         gauginoM3, higgsMu, higgsAt,
                                         mq3, mu3, md3, mq2, mu2, md2,
                                         mq1, mu1, md1, ml3, me3, ml2,
                                         me2, ml1, me1, mA},
                                         Rscale->MS, scheme->"DRbar",
                                         sources->{1,1,1},
                                         numerical->False];

           {mass, dmass}
    ];

RunSUSYHD[MS_, TB_, Xt_, Mi_, M3_] :=
    Module[{gauginoM1 = Mi, gauginoM2 = Mi, gauginoM3 = M3, higgsMu = Mi, higgsAt, mq3=MS,
            mu3=MS, md3=MS, mq2=MS, mu2=MS, md2=MS, mq1=MS, mu1=MS, md1=MS,
            ml3=MS, me3=MS, ml2=MS, me2=MS, ml1=MS, me1=MS, mA=MS,
            mass, dmass },

           SetSMparameters[Mtpole, alphaSAtMZ];
           higgsAt=(higgsMu/TB + Xt MS);

           mass = MHiggs[{TB, gauginoM1, gauginoM2, gauginoM3,
                          higgsMu, higgsAt, mq3, mu3, md3, mq2, mu2,
                          md2, mq1, mu1, md1, ml3, me3, ml2, me2, ml1,
                          me1, mA}, Rscale->MS, scheme->"DRbar",
                          hiOrd->{1,1,1,0}, numerical->True,
                          split->True];

           dmass= \[CapitalDelta]MHiggs[{TB, gauginoM1, gauginoM2,
                                         gauginoM3, higgsMu, higgsAt,
                                         mq3, mu3, md3, mq2, mu2, md2,
                                         mq1, mu1, md1, ml3, me3, ml2,
                                         me2, ml1, me1, mA},
                                         Rscale->MS, scheme->"DRbar",
                                         sources->{1,1,1},
                                         numerical->False];

           {mass, dmass}
    ];

(********** HSSUSY scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

ScanMS[TB_, Xt_] :=
    Module[{data},
           data = {N[#], N[TB], N[Xt], Sequence @@ RunSUSYHD[N[#],TB,Xt]}& /@ LogRange[500, 10^16, 60];
           Export["SUSYHD_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[N[Xt]] <> ".dat", data, "Table"];
          ];

(* ScanMS[#, Sqrt[6]]& /@ {2, 10, 20, 50}; *)

(********** HSSUSY scenario 2: TB = 2, MS = 2 TeV **********)

ScanXt[TB_, MS_, start_:-4, stop_:4, steps_:60] :=
    Module[{res},
           res = {MS, TB, N[#], Sequence @@ RunSUSYHD[MS, TB, #]}& /@ LinearRange[start, stop, steps];
           Export["SUSYHD_Xt_TB-" <> ToString[TB] <> "_MS-" <> ToString[MS] <> ".dat", res, "Table"];
           res
          ];

(* ScanXt[20, 2000]; *)

(********** split-MSSM scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********)

ScanSplitMS[TB_, Xt_, Mi_, start_:500, stop_:1.0 10^16, steps_:60] :=
    Module[{res},
           res = {#, Mi, Mi, Mi, TB, Xt, Sequence @@ RunSUSYHD[#, TB, Xt, Mi, Mi]}& /@ LogRange[start, stop, steps];
           Export["SUSYHD_MS_TB-" <> ToString[TB] <> "_Xt-" <> ToString[Xt] <>
                  "_Mi-" <> ToString[Mi] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitMS[#, N@Sqrt[6], 2000]& /@ {2, 10, 20, 50}; *)

(********** split-MSSM scenario 2: TB = 10, MS = 5 TeV **********)

ScanSplitXt[TB_, MS_, Mi_, start_:-4, stop_:4, steps_:60] :=
    Module[{res},
           res = {MS, Mi, Mi, Mi, TB, N[#], Sequence @@ RunSUSYHD[MS, TB, #, Mi, Mi]}& /@ LinearRange[start, stop, steps];
           Export["SUSYHD_Xt_TB-" <> ToString[TB] <> "_MS-" <> ToString[MS] <>
                  "_Mi-" <> ToString[Mi] <> ".dat", res, "Table"];
           res
          ];

(* ScanSplitXt[10, 5000, 2000]; *)
