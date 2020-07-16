#!/bin/bash
#First step for running simulation. 
#This is relativley quick, will run kmc and then make a histogram
#After this is finished, look at the histo file and determine L and U values 
#Then use the L and U values for run_sim_pt2.sh

gs=$1
cov=$2
het=$3
rs=10000
rand_het=1
folder=$4
rc="rc"
k=21
gc=$5
err=$6
rep=$7
#KMC k=24
# simulation.gs1000000.cov40.het1.5.rs10000.random_het1.gc50.rc.err1.repeats0.fasta
/seq/schatz/tbenavi/software/KMC/bin/kmc -k${k} -t6 -m32 -ci1 -cs100000 -fm ${folder}/simulation.chr20.gs$gs.cov$cov.het$het.rs$rs.random_het${rand_het}.gc${gc}.${rc}.err${err}.repeats${rep}.fasta \
        ${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep} tmp/

#kmc_tools transform histogram
kmc_tools transform ${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep} \
        histogram ${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_k${k}_err${err}_repeats${rep}.histo -cx100000

#After this, look at the histograms and determine L and U values to send to the next script, run_sim_pt2.sh