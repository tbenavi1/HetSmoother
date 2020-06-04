#python squish_reads.py input_reads.fasta squished_reads.fasta

from itertools import groupby
import sys

#2 for single line multi fasta file, 4 for fastq file
num_lines_per_read = 2

def squish(read):
  return "".join(x[0] for x in groupby(read))

with open(sys.argv[1], 'r') as file_input_reads, open(sys.argv[2], 'w') as file_output_reads:
  for i, line in enumerate(file_input_reads):
    if i%num_lines_per_read==1:
      line = squish(line)
    file_output_reads.write(line)
