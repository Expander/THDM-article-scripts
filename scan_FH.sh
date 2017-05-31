#!/bin/sh

start=91.1876
stop=100000
steps=60
step_size=log

Xt=0
TB=5
MS=2000
UsedMi=MS

parameter=MS

BASEDIR=$(dirname $0)
FHdir="/undefined/FeynHiggs/directory"
FHlooplevel=2

run_fh() {
    local fh_dir="$1"
    ${BASEDIR}/run_FeynHiggs_SLHA.sh "$fh_dir" \
              "${MS}" "${TB}" "${Xt}" "${Mi}" "${FHlooplevel}"
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --FH-dir=*)              FHdir=$optarg ;;
            --parameter=*)           parameter=$optarg ;;
            --start=*)               start=$optarg ;;
            --stop=*)                stop=$optarg ;;
            --steps=*)               steps=$optarg ;;
            --step-size=*)           step_size=$optarg ;;
            --MS=*)                  MS=$optarg ;;
            --TB=*)                  TB=$optarg ;;
            --Xt=*)                  Xt=$optarg ;;
            --Mi=*)                  UsedMi=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

# printf "# MS = ${MS}, TanBeta = ${TB}, Xt = ${Xt}\n"
# printf "# %14s %16s %16s\n" "$parameter" "FeynHiggs" "DeltaFeyhnh"

for i in `seq 0 $steps`; do
    # calculate current value for the scanned variable
    case "$step_size" in
        linear)
            value=$(cat <<EOF | bc
scale=10
$start + ($stop - $start)*${i} / $steps
EOF
                 ) ;;
        log)
            value=$(cat <<EOF | bc -l
scale=10
e(l($start) + (l($stop) - l($start))*${i} / $steps)
EOF
                 ) ;;
        *) echo "Error: unknown step size: $step_size"
           exit 1 ;;
    esac

    eval "${parameter}=${value}"

    if [ "x$UsedMi" = xMS ] ; then
        Mi="$MS"
    else
        Mi="$UsedMi"
    fi

    FHlooplevel=2
    FHout=$(run_fh "$FHdir")
    MhFH=$(echo "$FHout" | awk '{ print $1 }')
    DeltaMhFH=$(echo "$FHout" | awk '{ print $2 }')

    FHlooplevel=0
    FHout=$(run_fh "$FHdir")
    MhFHEFT=$(echo "$FHout" | awk '{ print $1 }')
    DeltaMhFHEFT=$(echo "$FHout" | awk '{ print $2 }')

    printf "%16s %16s %16s %16s %16s\n" "$value" "$MhFH" "$DeltaMhFH" "$MhFHEFT" "$DeltaMhFHEFT"
done
