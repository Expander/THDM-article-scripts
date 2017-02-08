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

run_fh() {
    local fh_dir="$1"
    ${BASEDIR}/run_FeynHiggs_2.12.0_SLHA.sh "$fh_dir" "${MS}" "${TB}" "${Xt}" "${Mi}"
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
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

    FHout=$(run_fh "${HOME}/packages/FeynHiggs-2.12.2/build")
    MhFH=$(echo "$FHout" | awk '{ print $1 }')
    DeltaMhFH=$(echo "$FHout" | awk '{ print $2 }')

    printf "%16s %16s %16s\n" "$value" "$MhFH" "$DeltaMhFH"
done
