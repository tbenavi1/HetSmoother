#remove_het.py takes an unzipped fasta (single line per sequence) or fastq file and a list of kmer pairs as input. It reads in the list of kmers first and prioritizes which kmer to change. Then, it reads the fasta/fastq file line by line, if it finds first kmer it replaces it with second kmer. It also produces a text file which records how many changes were made for each read.
#python remove_het.py kmer_pairs.tsv input_reads.fastq edited_reads.fastq replacements.txt
import re
import sys

def reverse_complement(kmer):
  """
  Assumes kmers are upper case.
  """
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
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
k=24

#change to 2 for single line multi fasta file, change to 4 for fastq file
num_lines_per_read = 4

#how often to update user about progress (1000 for long reads, 10000 for short reads)
reads_per_set = 10000

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
    #If we've completed reads_per_set reads (fasta files need to have each sequence on exactly one line)
    if i%(num_lines_per_read*reads_per_set)==0:
      print(str(int(i/num_lines_per_read))+ " reads have been analyzed")
      print("this set: " + str(set_replacements) + " replacements have been made across " + str(set_reads_with_replacements) + " reads")
      print("running total: "+ str(total_replacements) + " replacements have been made across " + str(total_reads_with_replacements) + " reads")
      set_replacements = 0
      set_reads_with_replacements = 0
    #If this is the line with the read sequence
    if i%num_lines_per_read==1:
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
        elif reverse_complement(kmer1) in rep:
          read_replaced = True
          edited_read += reverse_complement(rep[reverse_complement(kmer1)])
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
