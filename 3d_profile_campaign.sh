#!/bin/bash
# profile_3d_campaign.sh

mkdir -p results_3d

OPTS="O1 O2 O3"
LIKWID_GROUPS="CACHE FLOPS_DP L2CACHE L3CACHE MEM_DP"
SIZES="30 50 80 100"
SWEEPS=500
CORE=2

for opt in $OPTS; do
    for N in $SIZES; do
        for grp in $LIKWID_GROUPS; do
            outfile="results_3d/${opt}_N${N}_${grp}.txt"
            echo "Running: $opt N=$N group=$grp"
            
            # warm-up run (discarded) to prime caches and resolve page faults
            likwid-perfctr -C $CORE -g $grp -m \
                ./poisson_benchmark_${opt} profile $N $SWEEPS \
                > /dev/null 2>&1
            
            # measurement run
            likwid-perfctr -C $CORE -g $grp -m \
                ./poisson_benchmark_${opt} profile $N $SWEEPS \
                > "$outfile" 2>&1
        done
    done
done

echo "Campaign complete. Results in results_3d/"
