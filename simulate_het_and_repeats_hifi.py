#Makes simulated reads 
#Arg 1: genome size (ex. 2000000)
#Arg 2: Coverage (ex. 30 for 30x)
#Arg 3: Read size (ex. 10000)
#Arg 4: Percentage of heterozygosity, must be between 0-100 (ex. 1)
#Arg 5: Whether the het should be randomly distributed or not (1 == randomly distributed)
#Arg 6: Percent GC, must be between 0-100 (ex. 50) 
#Arg 7: Error rate, must be between 0-100 (ex. 1 for 1% error in the reads)
#Arg 8: Repeat rate, must be between 0-100 (ex. 2 for 2% repeats)
#Arg 9: Path to a template genome (must be single line fasta file with no "N" or "n")

import csv
import sys
import random
import statistics
import numpy as np
import scipy.stats as stats

genome_size = int(sys.argv[1]) #1000000
coverage = int(sys.argv[2]) #40
read_size = int(sys.argv[3]) #10000
het = float(sys.argv[4]) #0 #This needs to be a percentage (so 1 would mean 1%)
random_het = int(sys.argv[5])
rc = "rc" 
gc = float(sys.argv[6])
err_rate = float(sys.argv[7]) #This also needs to be a percentage (so "1" would mean 1%)
spread=1
repeats = float(sys.argv[8]) #This means there will be X % of the genome that is repetative 
copy_cat_genome = sys.argv[9]
out = open("simulation.chr20.gs" + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] + "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".fasta", "w")
out_mat = open("mat_genome_chr20.gs" + sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] +  "."  + "rc"  + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".txt", "w")
out_pat = open("pat_genome_chr20.gs"+ sys.argv[1] + ".cov" + sys.argv[2] + ".het" + sys.argv[4] + ".rs" + sys.argv[3] + ".random_het" + sys.argv[5] + ".gc" + sys.argv[6] +  "." + "rc" + ".err" + sys.argv[7] + ".repeats" + sys.argv[8] + ".txt", "w")

def find(ch1, ch2,string1):
    pos = []
    for i in range(len(string1)):
        if ch1 == string1[i] or ch2 == string1[i]:
            pos.append(i)
    return pos        

def copycat(copy_cat_genome, gs, gc):
	#Use a provided genome as a template 
	copy_cat = ""
	with open(copy_cat_genome, encoding="utf8", errors='ignore') as f: 
		print("Opening template genome")
		Lines = f.readlines()
		for line in Lines:
			line = line.strip()
			if (line.count(">") == 0):
				copy_cat = copy_cat + line
	#Take the first X base pairs from the template to use 	
	mat_genome = copy_cat[:gs]
	
	#Option to change the amount of GC 
	#The main issue with this is it could break up repeats within the template genome. 
	if int(gc) != -1:
		print("changing gc") 
		#What is the current gc content? 
		gc_content = (mat_genome.count("C")+mat_genome.count("G"))/(len(mat_genome))

		if (gc_content*100) > gc: #The mat genome has too many gcs 
			while ((gc_content*100) > gc):
				num_bp_to_chng = int(((((gc_content*100) - gc)/100)*len(mat_genome))+1)
				nums = random.sample(find("G", "C", mat_genome), num_bp_to_chng)
				counter = 0 #Just for our progress report. 
				for i in nums:
					counter = counter + 1
					if (i%2) == 0:
						mat_genome = mat_genome[:(i)] + "A" + mat_genome[(i+1):]
					else:
						mat_genome = mat_genome[:(i)] + "T" + mat_genome[(i+1):]
					if counter%100000 == 0:
						gc_content = (mat_genome.count("C") + mat_genome.count("G"))/(len(mat_genome))
						print(str((gc_content*100)))
				gc_content = (mat_genome.count("C")+mat_genome.count("G"))/(len(mat_genome))

		elif (gc_content*100) < gc: #Here, our mat genome doesn't have enough gc! 		
			while ((gc_content*100) < gc):
				num_bp_to_chng = int((((gc - (gc_content*100))/100)*len(mat_genome))+1)
				nums = random.sample(find("A", "T", mat_genome), num_bp_to_chng)
				counter = 0 #Just for our progress report. 
				for i in nums:
					counter = counter + 1
					if (i%2) == 0:
						mat_genome = mat_genome[:(i)] + "G" + mat_genome[(i+1):]
					else:
						mat_genome = mat_genome[:(i)] + "C" + mat_genome[(i+1):]
					if counter%100000 == 0:
						gc_content = (mat_genome.count("C") + mat_genome.count("G"))/(len(mat_genome))
						print(str((gc_content*100)))
				gc_content = (mat_genome.count("C")+mat_genome.count("G"))/(len(mat_genome))

	print("Done with GC addition")
	return mat_genome

def add_errors(string, err_rate): #Will want to add errors to the reads. This will be on a per-read basis
	read_w_errors=string
	nums = random.sample(range(0,len(string)), int((err_rate/100)*len(string)))
	for i in nums:
		read_w_errors = read_w_errors[:(i-1)] + random_string_exclude(read_w_errors[i-1]) + read_w_errors[i:]
	
	return read_w_errors

def find_het(het_is_here, read_length, starting_position, rc_or_c):
	#There totally must be an easier way to do this. But this is sufficient for now 
	het_loc_for_read = ""
	y = [x for x in het_is_here if starting_position <= x < (starting_position+read_length)]
	
	if rc_or_c == 1:
		z = [(read_length-1 - (x-starting_position)) for x in y]
	else:
		z = [x-starting_position for x in y]
	
	for x in z:
		het_loc_for_read = het_loc_for_read  + str(x) + ","
	
	return het_loc_for_read

def reverse_comp(string):
	reverse = string[::-1]
	rc = ""
	for i in reverse:
		if i == "A":
			rc = rc + "T"
		elif i == "T":
			rc = rc + "A"
		elif i == "C":
			rc = rc + "G"
		elif i == "G":
			rc = rc + "C"
		elif i == "N":
			rc = rc + "N"
	return rc

def random_string(len):
	#Self explainatory 
	string = ""
	for i in range (0, len):
		nuc = random.choice("ACTG")
		string = string + nuc
	return string

def weighted_string(len, gc):
	#Like random strin, but will randomly choose from a string that contains a specified amount
	#of Gs and Cs
	ref = "C"*int(gc*100) + "G"*int(gc*100) + "A"*(100-int(gc*100)) + "T"*(100-int(gc*100))
	string = ""
	for i in range (0, len):
		nuc = random.choice(ref)
		string = string + nuc
	return string

def random_string_exclude(char):
	#When we want to pick a random letter (ACTG) but it can't be the same character as provided
	new_char = char
	while new_char == char:
		new_char = random.choice("ACTG")
	return new_char

def reverse_complement(kmer):
  """
  Assumes kmers are upper case.
  """
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  return "".join(complement[base] for base in reversed(kmer))

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
				return 1
			else:
				return 0
		if (B_rc < A < A_rc < B) or (B_rc < A_rc < A < B):
			if flipped:
				return 0
			else:
				return 1
		if (A == B_rc < B == A_rc):
			continue
		if (B_rc < A < B < A_rc):
			if flipped:
				return 1
			else:
				return 0
		if flipped:
			return 0
		else:
			return 1

def add_het(first_genome, het_level): 
	#Adds the het randomly  
	nums = random.sample(range(0,len(first_genome)), int((het_level/100)*len(first_genome)) )
	nums.sort()
	mat_nums = []
	pat_nums = []
	new_genome = first_genome
#	x = 0
	for i in nums:
#		x = x + 1
#		if x%10000==0:
#			print(str(x/len(nums)*100) + "%")
		new_genome = new_genome[:i] + random_string_exclude(new_genome[i]) + new_genome[i+1:]
		kmer1 = first_genome[i-10:i+11]
		kmer2 = new_genome[i-10:i+11]
		priority = prioritize_snp(10, kmer1, kmer2)
		if priority == 0:
			mat_nums.append(i)
		elif priority == 1:
			pat_nums.append(i)
		else:
			print("het change at position " + str(i) + " is ambigious.")
	print("added het")
	return new_genome, mat_nums, pat_nums
	
def add_het_evenly(first_genome, het_level):
	#This time we aren't going to distribute the het randomly. This time it will be even. 
	new_genome = first_genome
	nums = list(range(0,(len(first_genome))))[(int(100 / het_level)-1)::int(100 / het_level)]
	for i in nums:
		new_genome = new_genome[:i] + random_string_exclude(new_genome[i]) + new_genome[i+1:]
	return new_genome, nums

def add_repeats(sequence, repeat_rate, spread):
	#There must be a better name than "spread", but basically this is how many groups 
	#of repeats to add 
	spread = int(len(sequence)/1000) 
	if spread == 0:
		spread = 1
	print(spread)
	repeat_seq = sequence
	#First, calculate how many bases will be effected.
	effected_bases = len(sequence) * (repeat_rate / 100) 
	#Next calculate the size of the repeat "chunks" 
	chunk_size = effected_bases / spread
	#Then determine where in the genome we want to place the repeated regions
	nums = random.sample(range(0,len(sequence)), spread) 
	#Now, let's go through and replace them
	for i in nums:
		
		mer_len=1 #Only adding homopolymers for now 
		# Make the mer
		mer = repeat_seq[i:(i+mer_len)]
		# Make the entire repeat seq
		add_this = ""	
		while len(add_this) <= chunk_size :
			add_this = add_this + mer
		add_this = add_this[:int(chunk_size)]
		# Add the repeats to the genome 
		repeat_seq = repeat_seq[:i] + add_this + repeat_seq[(i+len(add_this)):]
	print("Added repeats")
	return repeat_seq[:len(sequence)] 

#####################################
#Main chunk
if sys.argv[9] == "0": #We aren't given a fasta file to use as a template
	mat_genome = weighted_string(genome_size,gc)
else: #We are given a copycat genome! 
	print("Using this as a template: "+copy_cat_genome)
	mat_genome = copycat(copy_cat_genome, genome_size, gc)

#Add repeats
if repeats != 0.0:
	mat_genome = add_repeats(mat_genome, repeats, spread)
mat_genome = mat_genome[:genome_size]

#Add het
if (het != 0.0) and (random_het == 1): #We want to randomly dist the het
	pat_genome, mat_het_is_here, pat_het_is_here = add_het(mat_genome, het)
elif (het != 0.0) and (random_het == 0): #We want to evenly dist the het
	pat_genome, het_is_here = add_het_evenly(mat_genome, het)
else: 
	pat_genome = mat_genome
	het_is_here = []
	mat_het_is_here = []
	pat_het_is_here = []

#We want to have a copy of the mat and pat actual genome (as a reference)
out_mat.write(mat_genome)
out_pat.write(pat_genome)

#Calculate the number of reads that we'll need to make 
number_reads = (genome_size * coverage) / read_size

for i in range(0,int(number_reads)): #This is where we make the reads. 
	#Pick a random place to start 
	random_number = random.randint(0,genome_size)
	#read_length = 
	read_het_locations = "" #Keeping track of where we add the het. 
	if (i%2) == 0: #Use mat 
		this_read = mat_genome[random_number:(random_number + read_size)]
		which_genome="mat"
	else: # Use pat 
		this_read = pat_genome[random_number:(random_number + read_size)]
		which_genome = "pat"
	
	#Decide if we will use the reverse compliment or not 
	rc_or_c = random.randint(0,1)
	rc="compliment"
	if rc_or_c == 1:
		rc = "reverse_compliment"
		this_read = reverse_comp(this_read)
	
	#Find where the het was added. There is totally an esier way, but for now this works well enough.
	if which_genome == "mat":
		read_het_locations = find_het(mat_het_is_here, len(this_read), random_number, rc_or_c)
	else:
		read_het_locations = find_het(pat_het_is_here, len(this_read), random_number, rc_or_c)

	#Finally, add error to each read. 
	this_read = add_errors(this_read, err_rate) 

	#write this in a fasta file
	#The name for the read will include the read #, whether it is from the mat or pat,
	#Whether it is the complimemt or reverse, and the het locations
	out.write(">read_"+str(i) + "|" + which_genome + "|" + rc +  "|" + str(read_het_locations) + "|" + "\n")
	out.write(this_read + "\n")

out.close()
out_mat.close()
out_pat.close()





