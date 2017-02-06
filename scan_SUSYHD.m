Needs["SUSYHD`"];

Mtpole = 173.21;
alphaSAtMZ = 0.1181;

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
                          hiOrd->{1,1,1,1}, numerical->True,
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

ScanMS[TB_, Xt_] :=
    Module[{data},
           data = {N[#], N[TB], N[Xt], Sequence @@ RunSUSYHD[N[#],TB,Xt]}& /@ LogRange[500, 10^16, 60];
           Export["SUSYHD_MS_TB-" <> ToString[N[TB]] <> "_Xt-" <> ToString[N[Xt]] <> ".dat", data, "Table"];
          ];

ScanMS[#, Sqrt[6]]& /@ {2, 10, 20, 50}
