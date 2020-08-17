#!/usr/bin/env python
#python compare_pairs.py true_pairs.tsv smudge_pairs.tsv kmer2snp_pairs.txt ebwt2indel_pairs.fasta discosnp_pairs.fasta results.tsv
import sys

def reverse_complement(string):
  """
  Assumes characters of string are upper case.
  """
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  return "".join(complement[base] for base in reversed(string))

def standardize_kmer_pair(kmer1, kmer2):
    kmer1_rc = reverse_complement(kmer1)
    kmer2_rc = reverse_complement(kmer2)
    kmer_to_pair = dict()
    num_diff = sum(1 for a, b in zip(kmer1, kmer2) if a!=b)
    if num_diff <= 2:
        kmer_to_pair[kmer1] = kmer2
        kmer_to_pair[kmer2] = kmer1
        kmer_to_pair[kmer1_rc] = kmer2_rc
        kmer_to_pair[kmer2_rc] = kmer1_rc
    else:
        kmer_to_pair[kmer1] = kmer2_rc
        kmer_to_pair[kmer2] = kmer1_rc
        kmer_to_pair[kmer1_rc] = kmer2
        kmer_to_pair[kmer2_rc] = kmer1
    first = min(kmer1, kmer2, kmer1_rc, kmer2_rc)
    return (first, kmer_to_pair[first])
i=0
with open(sys.argv[1], "r") as in_true_pairs:
    true_pairs = set()
    for line in in_true_pairs:
        kmer1, kmer2 = line.strip().split("\t")
        kmer1, kmer2 = standardize_kmer_pair(kmer1, kmer2)
        #if (kmer1, kmer2) in true_pairs:
        #    i+=1
        #    print(i)
        true_pairs.add((kmer1, kmer2))

with open(sys.argv[2], "r") as in_smudge_pairs:
    smudge_pairs = set()
    for line in in_smudge_pairs:
        kmer1, kmer2 = line.strip().split("\t")
        kmer1, kmer2 = standardize_kmer_pair(kmer1, kmer2)
        smudge_pairs.add((kmer1, kmer2))

with open(sys.argv[3], "r") as in_kmer2snp_pairs:
    kmer2snp_pairs = set()
    for line in in_kmer2snp_pairs:
        kmer1, kmer2, _ = line.strip().split(" ")
        kmer1, kmer2 = standardize_kmer_pair(kmer1, kmer2)
        kmer2snp_pairs.add((kmer1, kmer2))

with open(sys.argv[4], "r") as in_ebwt2indel_pairs:
    ebwt2indel_pairs = set()
    for i, line in enumerate(in_ebwt2indel_pairs):
        if i%4==1:
            kmer1 = line.strip()
        if i%4==3:
            kmer2 = line.strip()
            kmer1, kmer2 = standardize_kmer_pair(kmer1, kmer2)
            #if (kmer1, kmer2) in ebwt2indel_pairs:
            #    print(i)
            #    print(kmer1)
            #    print(kmer2)
            ebwt2indel_pairs.add((kmer1, kmer2))

with open(sys.argv[5], "r") as in_discosnp_pairs:
    discosnp_pairs = set()
    for i, line in enumerate(in_discosnp_pairs):
        if i%4==0:
            snp_locs = []
            num_snps = int(line.split('|')[3].split('_')[2])
            for j in range(num_snps):
                snp_loc = int(line.split('|')[1].split(',')[j].split(':')[1].split('_')[0])
                snp_locs.append(snp_loc)
        if i%4==1:
            kmer1s = []
            for j in range(num_snps):
                snp_loc = snp_locs[j]
                kmer_start = snp_loc-10 #assume k=21
                kmer_end = snp_loc+10
                kmer1 = line[kmer_start:kmer_end+1]
                kmer1s.append(kmer1)
        if i%4==3:
            for j in range(num_snps):
                snp_loc = snp_locs[j] #0-indexed
                kmer_start = snp_loc-10 #assume k=21
                kmer_end = snp_loc+10
                kmer1 = kmer1s[j]
                kmer2 = line[kmer_start:kmer_end+1]
                kmer1, kmer2 = standardize_kmer_pair(kmer1, kmer2)
                #print(kmer1)
                #print(kmer2)
                discosnp_pairs.add((kmer1, kmer2))
#print('length of discosnp_pairs')
#print(len(discosnp_pairs))

found_pairs = smudge_pairs | kmer2snp_pairs | ebwt2indel_pairs | discosnp_pairs

with open(sys.argv[6], "w") as out_results:
    for true_pair in true_pairs:
        kmer1, kmer2 = true_pair
        smudge = str(int(true_pair in smudge_pairs))
        kmer2snp = str(int(true_pair in kmer2snp_pairs))
        ebwt2indel = str(int(true_pair in ebwt2indel_pairs))
        discosnp = str(int(true_pair in discosnp_pairs)) 
        out_results.write("\t".join([kmer1, kmer2, smudge, kmer2snp, ebwt2indel, discosnp])+'\n')
    #for found_pair in found_pairs:
    #    kmer1, kmer2 = found_pair
    #    smudge = str(int(found_pair in smudge_pairs))
    #    kmer2snp = str(int(found_pair in kmer2snp_pairs))
    #    ebwt2indel = str(int(found_pair in ebwt2indel_pairs))
    #    discosnp = str(int(found_pair in discosnp_pairs))
    #    out_results.write("\t".join([kmer1, kmer2, smudge, kmer2snp, ebwt2indel, discosnp])+'\n')

tp_smudge_pairs = true_pairs & smudge_pairs
fn_smudge_pairs = true_pairs - smudge_pairs
fp_smudge_pairs = smudge_pairs - true_pairs
tp_kmer2snp_pairs = true_pairs & kmer2snp_pairs
fn_kmer2snp_pairs = true_pairs - kmer2snp_pairs
fp_kmer2snp_pairs = kmer2snp_pairs - true_pairs
tp_ebwt2indel_pairs = true_pairs & ebwt2indel_pairs
fn_ebwt2indel_pairs = true_pairs - ebwt2indel_pairs
fp_ebwt2indel_pairs = ebwt2indel_pairs - true_pairs
tp_discosnp_pairs = true_pairs & discosnp_pairs
fn_discosnp_pairs = true_pairs - discosnp_pairs
fp_discosnp_pairs = discosnp_pairs - true_pairs

num_tp_smudge = len(tp_smudge_pairs)
num_fn_smudge = len(fn_smudge_pairs)
num_fp_smudge = len(fp_smudge_pairs)
num_tp_kmer2snp = len(tp_kmer2snp_pairs)
num_fn_kmer2snp = len(fn_kmer2snp_pairs)
num_fp_kmer2snp = len(fp_kmer2snp_pairs)
num_tp_ebwt2indel = len(tp_ebwt2indel_pairs)
num_fn_ebwt2indel = len(fn_ebwt2indel_pairs)
num_fp_ebwt2indel = len(fp_ebwt2indel_pairs)
num_tp_discosnp = len(tp_discosnp_pairs)
num_fn_discosnp = len(fn_discosnp_pairs)
num_fp_discosnp = len(fp_discosnp_pairs)

#print(num_tp_smudge)
#print(num_fn_smudge)
#print(num_fp_smudge)
#print(num_tp_kmer2snp)
#print(num_fn_kmer2snp)
#print(num_fp_kmer2snp)
#print(num_tp_ebwt2indel)
#print(num_fn_ebwt2indel)
#print(num_fp_ebwt2indel)
print("\t".join([str(x) for x in [num_tp_smudge, num_fn_smudge, num_fp_smudge, num_tp_kmer2snp, num_fn_kmer2snp, num_fp_kmer2snp, num_tp_ebwt2indel, num_fn_ebwt2indel, num_fp_ebwt2indel, num_tp_discosnp, num_fn_discosnp, num_fp_discosnp]]))
#print(num_tp_discosnp)
#print(num_fn_discosnp)
#print(num_fp_discosnp)
