Get["MhEFT.m"];

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

RunMhEFT[MS_, TB_, Xt_, mA_, Mu_, Mi_, M3_] :=
    Module[{gauginoM1 = Mi, gauginoM2 = Mi, gauginoM3 = M3,
            Ab = 0, Atau = 0, verbose = 0, mass2L, mass3L, dmass = 0 },

           SMobs = {Mt -> Mtpole, MW -> 80.384, MZ -> 91.1876, 
                    Mh -> 125.15, \[Alpha]3MZ -> alphaSAtMZ, 
                    V -> 246.21971, \[Delta]ytth -> 0, \[Delta]\[Lambda]th -> 0};

           mass2L = mh /. evalMh[{MS,mA,TB,Xt,Ab,Atau,Mu,gauginoM1,gauginoM2,2,verbose}][[2]];
           mass3L = mh /. evalMh[{MS,mA,TB,Xt,Ab,Atau,Mu,gauginoM1,gauginoM2,3,verbose}][[2]];
           dmass = Abs[mass2L - mass3L];

           {mass2L, dmass}
    ];

RunMhEFT[MS_, TB_, Xt_, MA_] := RunMhEFT[MS, TB, Xt, MA, MS, MS, MS];

(********** THDM degenerate masses: TB = [2,50], MS = [1000, 10^16], Xt = ? **********)

ScanMhEFTMSTB[Xt_, MA_, MSstart_:1000, MSstop_:1.0 10^16, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LogRange[MSstart, MSstop, steps], LinearRange[TBstart, TBstop, steps]}];
           res = {Sequence @@ N[#], MA, Xt, Sequence @@ RunMhEFT[Sequence @@ #, Xt, MA]}& /@ tuples;
           Export["MhEFT_TB_MS_Xt-" <> ToString[Xt] <> "_MA-" <> ToString[MA] <> ".dat", res, "Table"];
           res
          ];

ScanMhEFTMSTB[0, 400];
ScanMhEFTMSTB[0, 800];

(********** THDM degenerate masses: TB = [2,50], MA = [100, 500], Xt = ? **********)

ScanMhEFTTBMA[Xt_, MS_, MAstart_:100, MAstop_:500, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[TBstart, TBstop, steps], LinearRange[MAstart, MAstop, steps]}];
           res = {MS, Sequence @@ N[#], Xt, Sequence @@ RunMhEFT[MS, #[[1]], Xt, #[[2]]]}& /@ tuples;
           Export["MhEFT_TB_MA_Xt-" <> ToString[Xt] <> "_MS-" <> ToString[MS] <> ".dat", res, "Table"];
           res
          ];

ScanMhEFTTBMA[0, 5000];
ScanMhEFTTBMA[0, 10^4];
ScanMhEFTTBMA[0, 5 10^4];

(********** THDM degenerate masses: MA = [100, 500], MS = [1000, 10^16], Xt = ? **********)

ScanMhEFTMSMA[Xt_, TB_, MSstart_:1000, MSstop_:1.0 10^16, MAstart_:100, MAstop_:500, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], TB, N[#[[1]]], Xt, Sequence @@ RunMhEFT[#[[2]], TB, Xt, #[[1]]]}& /@ tuples;
           Export["MhEFT_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <> ".dat", res, "Table"];
           res
          ];

ScanMhEFTMSMA[0, #]& /@ {2, 10, 20, 50};

(********** HGTHDM degenerate masses: TB = [2,50], MS = [1000, 10^16], Xt = ? **********)

ScanHGMhEFTMSTB[Xt_, MA_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LogRange[MSstart, MSstop, steps], LinearRange[TBstart, TBstop, steps]}];
           res = {N[#[[1]]], Mui, Mui, Mui, N[#[[2]]], MA, Xt,
                  Sequence @@ RunMhEFT[Sequence @@ #, Xt, MA, Mui, Mui, Mui]}& /@ tuples;
           Export["MhEFT_TB_MS_Xt-" <> ToString[Xt] <> "_MA-" <> ToString[MA] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

ScanHGMhEFTMSTB[0, 400, 2000];
ScanHGMhEFTMSTB[0, 800, 2000];

(********** HGTHDM degenerate masses: TB = [2,50], MA = [100, 500], Xt = ? **********)

ScanHGMhEFTTBMA[Xt_, MS_, Mui_, MAstart_:100, MAstop_:500, TBstart_:2, TBstop_:50, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[TBstart, TBstop, steps], LinearRange[MAstart, MAstop, steps]}];
           res = {MS, Mui, Mui, Mui, Sequence @@ N[#], Xt,
                  Sequence @@ RunMhEFT[MS, #[[1]], Xt, #[[2]], Mui, Mui, Mui]}& /@ tuples;
           Export["MhEFT_TB_MA_Xt-" <> ToString[Xt] <> "_MS-" <> ToString[MS] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

ScanHGMhEFTTBMA[0, 5000, 2000];
ScanHGMhEFTTBMA[0, 10^4, 2000];
ScanHGMhEFTTBMA[0, 5 10^4, 2000];

(********** HGTHDM degenerate masses: MA = [100, 500], MS = [1000, 10^16], Xt = ? **********)

ScanHGMhEFTMSMA[Xt_, TB_, Mui_, MSstart_:1000, MSstop_:1.0 10^16, MAstart_:100, MAstop_:500, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], Mui, Mui, Mui, TB, N[#[[1]]], Xt,
                  Sequence @@ RunMhEFT[#[[2]], TB, Xt, #[[1]], Mui, Mui, Mui]}& /@ tuples;
           Export["MhEFT_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

ScanHGMhEFTMSMA[0, #, 2000]& /@ {2, 10, 20, 50};
