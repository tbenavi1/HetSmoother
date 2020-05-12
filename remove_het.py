#remove_het.py takes an unzipped fastq file and a list of kmer pairs (k1 replaced as k2) as input. It reads in the list of kmers first. Then, it reads the fastq file line by line, if it finds k1 it replaces it with k2.
#remove_het.py kmer_pairs.tsv input_reads.fastq output_reads.fastq
import re
import sys

def reverse_complement(kmer):
  """
  Assumes kmers are upper case, no Ns.
  """
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
  return "".join(complement[base] for base in reversed(kmer))

def prioritize_kmers(kmer1, kmer2):
  """
  Assumes kmers are the same size and are 1 SNP different.
  """
  for i in range(len(kmer1)):
    if kmer1[i] != kmer2[i]:
      assert kmer1[:i]+kmer1[i+1:]==kmer2[:i]+kmer2[i+1:]
      break
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
    if (A < B < B_rc < A_rc) or (A < B_rc < B < A_rc) or (A < B_rc == B < A_rc):
      if flipped:
        return (kmer1, kmer2)
      else:
        return (kmer2, kmer1)
    if (B_rc < A < A_rc < B) or (B_rc < A_rc < A < B) or (B_rc < A == A_rc < B):
      if flipped:
        return (kmer2, kmer1)
      else:
        return (kmer1, kmer2)
    assert A == B_rc < B == A_rc 

#k=24 for beroe data
k=21

#make dictionary of replacements
rep = dict()

print("Loading kmer pairs")
with open(sys.argv[1], 'r') as file_kmer_pairs:
  for line in file_kmer_pairs:
    kmer1, kmer2 = line.strip().split('\t')
    results = prioritize_kmers(kmer1, kmer2)
    if results:
      kmer1, kmer2 = results
      rep[kmer1] = kmer2
print("Kmer pairs loaded")

#print(rep)

#print("making pattern")
##make a pattern that will be matched if any kmer1 is matched
#pattern = re.compile("|".join(rep.keys()))
#print("pattern made")

print("replacing reads")
with open(sys.argv[2], 'r') as file_input_reads, open(sys.argv[3], 'w') as file_output_reads, open(sys.argv[4], 'w') as file_output_numrep:
  total_replacements = 0
  total_reads_with_replacements = 0
  set_replacements = 0
  set_reads_with_replacements = 0
  for i,line in enumerate(file_input_reads):
    #If we've completed 1,000 reads (fasta file where each sequence is exactly on one line)
    if i%2000==0: #i%4000==0 for 1,000 reads of fastq file
      print(str(int(i/2))+ " reads have been analyzed") #str(int(i/4)) for fastq file
      print("this set: " + str(set_replacements) + " replacements have been made across " + str(set_reads_with_replacements) + " reads")
      print("running total: "+ str(total_replacements) + " replacements have been made across " + str(total_reads_with_replacements) + " reads")
      set_replacements = 0
      set_reads_with_replacements = 0
    #If this is the line with the read sequence
    if i%2==1: #i%4==1 for fastq file
      #line,n = pattern.subn(lambda m: rep[m.group(0)], line)
      read_len = len(line.strip())
      #n is the nubmber of replacements made for this read
      n = 0
      #j iterates over the read
      j = 0
      edited_read = ""
      read_replaced = False
      while j < read_len - k + 1:
        kmer1 = line[j:j+k]
        #print(kmer1)
        #if need to make a replacement
        if kmer1 in rep:
          read_replaced = True
          edited_read += rep[kmer1]
          n += 1
          set_replacements += 1
          total_replacements += 1
          j += k
        else:
          edited_read += line[j]
          j += 1
      if read_replaced:
        set_reads_with_replacements += 1
        total_reads_with_replacements += 1
      file_output_numrep.write(str(n)+'\n')
      #if n>0:
      #  print(str(n) + " replacements made")
      line = edited_read + line[j:]
    file_output_reads.write(line)
