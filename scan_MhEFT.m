Get["MhEFT.m"];

Mtpole = 173.21;
alphaSAtMZ = 0.1181;

RunMhEFT[MS_, TB_, Xt_, Mi_, M3_] :=
    Module[{gauginoM1 = Mi, gauginoM2 = Mi, gauginoM3 = M3, higgsMu = Mi, mq3=MS,
            mu3=MS, md3=MS, mq2=MS, mu2=MS, md2=MS, mq1=MS, mu1=MS, md1=MS,
            ml3=MS, me3=MS, ml2=MS, me2=MS, ml1=MS, me1=MS, mA=MS,
            Ab = 0, Atau = 0, ytLo = 3, verbose = 0, mass, dmass = 0 },

           SMobs = {Mt -> Mtpole, MW -> 80.384, MZ -> 91.1876, 
                    Mh -> 125.15, \[Alpha]3MZ -> alphaSAtMZ, 
                    V -> 246.21971, \[Delta]ytth -> 0, \[Delta]\[Lambda]th -> 0};

           mass = mh /. evalMh[{MS,mA,TB,Xt,Ab,Atau,higgsMu,gauginoM1,gauginoM2,ytLo,verbose}][[2]];

           {mass, dmass}
    ];

RunMhEFT[MS_, TB_, Xt_] := RunMhEFT[MS, TB, Xt, MS, MS];

Print @ RunMhEFT[2000, 5, 0]
