#!/bin/bash
#Part two of running smudge pairs 
#Need the L and U values 
#will reduce the kmc database
#Then will run smudge pairs 
#Then remove het
#then run kmc again and make a histogram 
gs=$1
cov=$2
het=$3
rs=10000
rand_het=1
ci=$4
cx=$5
folder=$6
rc="rc"
k=21
gc=$7
err=$8
rep=$9

echo $k

#KMC k=24
# simulation.gs1000000.cov40.het1.5.rs10000.random_het1.gc50.rc.err1.repeats0.fasta
#/seq/schatz/tbenavi/software/KMC/bin/kmc -k24 -t6 -m32 -ci1 -cs100000 -fm "coverage/simulation.gs$gs.cov$cov.het$het.rs$rs.random_het$rand_het.fasta" \
#       coverage/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het} tmp/

#kmc_tools transform histogram
#/seq/schatz/tbenavi/software/KMC/bin/kmc_tools transform coverage/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het} \
#       histogram coverage/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}.histo -cx100000


#kmc_tools transform reduce
/seq/schatz/kjenike/ctenophores/smooth/tests/chr20/KMC-middle_one_away/bin/kmc_tools transform \
	${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep} -ci$ci -cx$cx reduce \
	${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}_reduce

#smudge_pair
/seq/schatz/kjenike/ctenophores/smooth/tests/chr20/KMC-middle_two_away/bin/smudge_pairs \
	${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}_reduce \
	${folder}/coverages_kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}.tsv \
	${folder}/pairs_kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}.tsv > ${folder}/familysizes_kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}.tsv

#remove_het_my_copy.py
time python3 /seq/schatz/kjenike/ctenophores/smooth/tests/chr20/smooth_het-two_away/remove_het.py \
      ${folder}/pairs_kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}.tsv \
      ${folder}/simulation.chr20.gs$gs.cov$cov.het$het.rs$rs.random_het${rand_het}.gc${gc}.${rc}.err${err}.repeats${rep}.fasta \
      ${folder}/simulation.het_removed.gs$gs.cov$cov.het$het.rs$rs.random_het${rand_het}.gc${gc}.${rc}.err${err}.repeats${rep}.fasta \
      ${folder}/num_replacements_simulation.gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_${rc}_k${k}_gc${gc}_err${err}_repeats${rep}.txt \
      ${folder}/loc_replacements_simulation.gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_${rc}_k${k}_gc${gc}_err${err}_repeats${rep}.txt > ${folder}/output.gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}.msg 2>&1 

#kmc -k24
/seq/schatz/kjenike/ctenophores/smooth/tests/chr20/KMC-middle_one_away/bin/kmc -k${k} -t6 -m32 -ci1 -cs100000 -fm \
	${folder}/simulation.het_removed.gs$gs.cov$cov.het$het.rs$rs.random_het${rand_het}.gc${gc}.${rc}.err${err}.repeats${rep}.fasta \
	${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}_het_rmvd tmp/

#kmc_tools transform histogram
/seq/schatz/kjenike/ctenophores/smooth/tests/chr20/KMC-middle_one_away/bin/kmc_tools transform \
	${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}_het_rmvd histogram \
	${folder}/kmer_counts_gs${gs}_cov${cov}_het${het}_rs${rs}_random_het${rand_het}_gc${gc}_${rc}_err${err}_repeats${rep}_k${k}_het_rmvd.histo