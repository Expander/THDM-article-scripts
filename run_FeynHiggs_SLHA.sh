#!/bin/sh

if [ $# -ne 6 ] ; then
    echo "Error: 6 Arguments required!"
    echo "  $0 <FH-dir> <MS> <tan(beta)> <Xt> <Mi> <looplevel>"
    exit 1
fi

fh_dir="$1"
shift
MS="$1"
shift
TB="$1"
shift
Xt="$1"
shift
Mi="$1"
shift
looplevel="$1"

fh="${fh_dir}/FeynHiggs"
fh_in=fh.in
fh_out="${fh_in}.fh-001"

# calculate At for Mu = MS
At=$(echo "scale=16; $MS * $Xt + $Mi / $TB" | bc -l)
# calculate Ab for Mu = MS and Xb = 0
Ab=$(echo "scale=16; $MS * $TB" | bc -l)

print_slha_block_awk='
BEGIN {
   is_block = 0
   if (block == "") {
      print "Error: block name not defined"
      print "   Please define the block name with -v block=<block-name>"
      exit 1
   }
}
{
   pattern     = "^block[[:blank:]]*" tolower(block) "([^[:graph:]].*)?$"
   not_pattern = "^block[[:blank:]]*.*$"

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   }

   if (is_block)
      print $0
}
'

print_slha_block_entry_awk='
BEGIN {
   if (entries == "") {
      print "Error: entries not defined"
      print "   Please define the block entries with -v entries=<entries> ."
      print "   Multiple entries can be given using : as separator."
      exit 1
   }

   len = split(entries,k,":");
}
{
  matches = 1;

  for (i in k) {
     if ($(i) != k[i])
        matches = 0
  }

  if (matches == 1)
     print $(len + 1)
}
'

cat <<EOF > "${fh_in}"
Block MINPAR                 # Input parameters
    4   1.000000000e+00      # sign(mu)
Block SMINPUTS               # Standard Model inputs
    1   1.279500000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.181000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.732100000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.384               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
BLOCK EXTPAR
         0     ${MS}   # Q
         1     ${Mi}   # M1
         2     ${Mi}   # M2
         3     ${Mi}   # M3
        11     ${At}   # At
        12     ${Ab}   # Ab
        13     ${Ab}   # Atau
        23     ${Mi}   # MUE
        26     ${MS}   # MA0
        25     ${TB}   # TB at Q
        31     ${MS}   # MSL(1)
        32     ${MS}   # MSL(2)
        33     ${MS}   # MSL(3)
        34     ${MS}   # MSE(1)
        35     ${MS}   # MSE(2)
        36     ${MS}   # MSE(3)
        41     ${MS}   # MSQ(1)
        42     ${MS}   # MSQ(2)
        43     ${MS}   # MSQ(3)
        44     ${MS}   # MSU(1)
        45     ${MS}   # MSU(2)
        46     ${MS}   # MSU(3)
        47     ${MS}   # MSD(1)
        48     ${MS}   # MSD(2)
        49     ${MS}   # MSD(3)
EOF

rm -f "${fh_out}"
Mh=
DMh=

${fh} "${fh_in}" 40020${looplevel}3110 >/dev/null 2>&1

if [ -e "${fh_out}" ] ; then
    Mh=$(awk -v block=MASS "${print_slha_block_awk}" "${fh_out}" | \
         awk -v entries=25 "${print_slha_block_entry_awk}")

    DMh=$(awk -v block=DMASS "${print_slha_block_awk}" "${fh_out}" | \
          awk -v entries=25 "${print_slha_block_entry_awk}")
fi

[ "x$Mh" = "x" ] && Mh=-
[ "x$DMh" = "x" ] && DMh=-

echo "$Mh   $DMh"

rm -f "${fh_in}" "${fh_out}"
