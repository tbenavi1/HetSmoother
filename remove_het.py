#!/usr/bin/env python
#remove_het.py takes an unzipped fasta (single line per sequence) or fastq file and a list of kmer pairs (with up to two SNPs different) as input. It reads in the list of kmer pairs first and prioritizes how to edit the SNPs. Then, it reads the fasta/fastq file line by line and edits the SNPs accordingly. It also produces a text file which records how many changes were made for each read.
#python remove_het.py kmer_pairs.tsv input_reads.fastq edited_reads.fastq replacements.txt
import sys
from collections import defaultdict

def reverse_complement(kmer):
	"""
	Assumes kmers are upper case.
	"""
	complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	return "".join(complement[base] for base in reversed(kmer))

def get_snp_locations(kmer1, kmer2, k):
	"""
	Assumes kmers have up to two SNPS different. Also accounts for the
	possibility that the kmers may pair only if one them is noncanonical.
	"""
	locations = []
	already_found_first = False
	for i in range(k):
		if kmer1[i] != kmer2[i]:
			locations.append(i)
			if already_found_first:
				if kmer1[i+1:] == kmer2[i+1:]:
					return kmer2, locations
				else:
					return get_snp_locations(kmer1, reverse_complement(kmer2), k)
			already_found_first = True
	return kmer2, locations

def prioritize_snp(i, kmer1, kmer2):
	window = -1
	while True:
		window += 1
		flipped = False
		A = kmer1[i-window:i+window+1]
		B = kmer2[i-window:i+window+1]
		if len(A) != 2*window+1:
			return
		if A > B:
			flipped = True
			A, B = B, A
		A_rc = reverse_complement(A)
		B_rc = reverse_complement(B)
		if (A < B < B_rc < A_rc) or (A < B_rc < B < A_rc):
			if flipped:
				return (kmer2, kmer1[i])
			else:
				return (kmer1, kmer2[i])
		if (B_rc < A < A_rc < B) or (B_rc < A_rc < A < B):
			if flipped:
				return (kmer1, kmer2[i])
			else:
				return (kmer2, kmer1[i])
		if (A == B_rc < B == A_rc):
			continue
		if (B_rc < A < B < A_rc):
			if flipped:
				return (kmer2, kmer1[i])
			else:
				return (kmer1, kmer2[i])
		if flipped:
			return (kmer1, kmer2[i])
		else:
			return (kmer2, kmer1[i])

def prioritize_kmers(kmer1, kmer2, k):
	prioritize_kmers_results = []
	kmer2, locations = get_snp_locations(kmer1, kmer2, k)
	for i in locations:
		prioritize_snp_results = prioritize_snp(i, kmer1, kmer2)
		if prioritize_snp_results:
			kmer, snp = prioritize_snp_results
			prioritize_kmers_results.append((kmer, i, snp))
			prioritize_kmers_results.insert(0, (reverse_complement(kmer), k-1-i, reverse_complement(snp)))
	return prioritize_kmers_results

#k=24 for beroe data
k=21

#set number of lines per read for input read file
if sys.argv[2].endswith(('.fq','.fastq')):
	num_lines_per_read = 4
	print('Input read file is a fastq file.')
elif sys.argv[2].endswith(('.fa', '.fasta')):
	num_lines_per_read = 2
	print('Input read file is a fasta file.')
else:
	print('Input read file must be an unzipped fastq or fasta file. The file name must end with .fa, .fasta, .fq, or .fastq')
	sys.exit()

#how often to update user about progress (1000 for long reads, 10000 for short reads)
reads_per_set = 1000

#make dictionary of replacements
rep = defaultdict(list)

print("Loading kmer pairs")
with open(sys.argv[1], 'r') as file_kmer_pairs:
	for line in file_kmer_pairs:
		kmer1, kmer2 = line.strip().split('\t')
		results = prioritize_kmers(kmer1, kmer2, k)
		for kmer, location, snp in results:
			rep[kmer].append((location, snp))
print("Kmer pairs loaded")

print("replacing reads")
with open(sys.argv[2], 'r') as file_input_reads, open(sys.argv[3], 'w') as file_output_reads, open(sys.argv[4], 'w') as file_output_numrep:
	total_replacements = 0
	total_reads_with_replacements = 0
	set_replacements = 0
	set_reads_with_replacements = 0
	for i,line in enumerate(file_input_reads):
		#If we've completed reads_per_set reads (fasta files need to have each sequence on exactly one line)
		if i%(num_lines_per_read*reads_per_set)==0:
			print(str(int(i/num_lines_per_read))+ " reads have been analyzed")
			print("this set: " + str(set_replacements) + " replacements have been made across " + str(set_reads_with_replacements) + " reads")
			print("running total: "+ str(total_replacements) + " replacements have been made across " + str(total_reads_with_replacements) + " reads")
			set_replacements = 0
			set_reads_with_replacements = 0
		#If this is the line with the read sequence
		if i%num_lines_per_read==1:
			read_len = len(line.strip())
			#n is the number of replacements made for this read
			n = 0
			#j iterates over the read
			j = 0
			#m is the final position of the edited read so far
			m = -1
			edited_read = ""
			read_replaced = False
			while j < read_len - k + 1:
				kmer1 = line[j:j+k]
				#if need to make a replacement
				if kmer1 in rep:
					read_replaced = True
					for location, snp in rep[kmer1]:
						#If the location of the snp on the read is past where we have already edited
						if location + j > m:
							replacement = line[m+1:location+j]+snp
							edited_read += replacement
							n += 1
							set_replacements += 1
							total_replacements += 1
							m += len(replacement)
					j += 1
				else:
					if j > m:
						edited_read += line[j]
						m+=1
					j += 1
			if read_replaced:
				set_reads_with_replacements += 1
				total_reads_with_replacements += 1
			file_output_numrep.write(str(n)+'\n')
			line = edited_read + line[m+1:]
		file_output_reads.write(line)
