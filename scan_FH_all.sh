#!/bin/sh

BASEDIR=$(dirname $0)

#********** HSSUSY scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********

for TB in 2 10 20 50; do
    Xt=2.44949
    ${BASEDIR}/scan_FH.sh --parameter=MS --TB=${TB} --Xt=${Xt} --start=500 --stop=5000 \
        | awk "{ print \$1 \"    \" ${TB} \"    \" ${Xt} \"    \" \$2 \"    \" \$3 }" \
        | tee FH-2.12.2_MS_TB-${TB}_Xt-${Xt}.dat
done

Xt=2.44949
TB=20
MS=2000

#********** HSSUSY scenario 2: TB = 2, MS = 2 TeV **********

${BASEDIR}/scan_FH.sh --parameter=Xt --MS=${MS} --TB=${TB} --start=-4 --stop=4 --step-size=linear \
    | awk "{ print ${MS} \"    \" ${TB} \"    \" \$1 \"    \" \$2 \"    \" \$3 }" \
    | tee FH-2.12.2_Xt_TB-${TB}_MS-${MS}.dat

#********** split-MSSM scenario 1: TB = 2, 10, 20, 50, Xt = Sqrt[6] **********

for TB in 2 10 20 50; do
    Xt=2.44949
    Mi=2000
    ${BASEDIR}/scan_FH.sh --parameter=MS --TB=${TB} --Xt=${Xt} --Mi=${Mi} --start=500 --stop=5000 \
        | awk "{ print \$1 \"    \" ${Mi} \"    \" ${Mi} \"    \" ${Mi} \"    \" ${TB} \"    \" ${Xt} \"    \" \$2 \"    \" \$3 }" \
        | tee FH-2.12.2_MS_TB-${TB}_Xt-${Xt}_Mi-${Mi}.dat
done

#********** split-MSSM scenario 2: TB = 10, MS = 5 TeV **********

Xt=2.44949
TB=10
MS=5000
Mi=2000

${BASEDIR}/scan_FH.sh --parameter=Xt --MS=${MS} --TB=${TB} --Mi=${Mi} --start=-4 --stop=4 --step-size=linear \
    | awk "{ print ${MS} \"    \" ${Mi} \"    \" ${Mi} \"    \" ${Mi} \"    \" ${TB} \"    \" \$1 \"    \" \$2 \"    \" \$3 }" \
    | tee FH-2.12.2_Xt_TB-${TB}_MS-${MS}_Mi-${Mi}.dat
