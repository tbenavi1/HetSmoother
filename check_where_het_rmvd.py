import sys
#python check_where_het_rmvd.py file.fasta loc_replacements.txt track2.mode.bed num_ue_mode.txt num_uu_mode.txt num_ce_mode.txt num_ie_mode.txt num_iu_mode.txt

#load the genome het locations that are expected to be editable by mode (i.e. not too dense)
with open(sys.argv[3], "r") as in_genome_ee_het_spots:
	genome_ee_het_spots = set()
	for line in in_genome_ee_het_spots:
		genome_ee_het_spot = int(line.strip().split("\t")[1])
		genome_ee_het_spots.add(genome_ee_het_spot)
print('Finished loading genome het locations that are expected to be editable.')

#open fasta input file, loc_replacements input file, num_ue output file, and num_uu output file
with open(sys.argv[1], "r") as in_fasta, open(sys.argv[2], "r") as in_loc_replacements, open(sys.argv[4], "w") as out_num_ue, open(sys.argv[5], "w") as out_num_uu, open(sys.argv[6], "w") as out_num_ce, open(sys.argv[7], "w") as out_num_ie, open(sys.argv[8], "w") as out_num_iu:
	for line in in_fasta:
		if (line.count(">") > 0):
			line1 = line.strip()
			line2 = in_loc_replacements.readline().strip()
			if line2:
				line2 += ","
			read_het_spots=line1.split("|")[4].split(",")[:-1]
			#read_err_spots=line1.split("|")[5].split(",")[:-1]
			edited_read_het_spots=line2.split(",")[:-1]
			read_start = int(line1.split("|")[1])
			rc = line1.split("|")[3]
			#translate read_het_spots and edited_read_het_spots to genomic coordinate system
			if (rc == "original") or (rc == "compliment"): #backwards compatibility
				read_het_spots = set(read_start + int(x) for x in read_het_spots)
				edited_read_het_spots = set(read_start + int(x) for x in edited_read_het_spots)
				#read_err_spots = set(read_start + int(x) for x in read_err_spots)
			else:
				read_size = len(in_fasta.readline().strip()) #added as temporary fix
				read_het_spots = set(read_start + read_size - 1 - int(x) for x in read_het_spots)
				edited_read_het_spots = set(read_start + read_size - 1 - int(x) for x in edited_read_het_spots)
				#read_err_spots = set(read_start + read_size - 1 - int(x) for x in read_err_spots)
			correctly_edited = edited_read_het_spots & read_het_spots
			incorrectly_edited = edited_read_het_spots - read_het_spots
			incorrectly_unedited = read_het_spots - edited_read_het_spots
			editable_read_het_spots = read_het_spots & genome_ee_het_spots
			read_ue_het_spots = edited_read_het_spots - editable_read_het_spots #- read_err_spots
			read_uu_het_spots = editable_read_het_spots - edited_read_het_spots
			out_num_ce.write(str(len(correctly_edited))+'\n')
			out_num_ie.write(str(len(incorrectly_edited))+'\n')
			out_num_iu.write(str(len(incorrectly_unedited))+'\n')
			out_num_ue.write(str(len(read_ue_het_spots))+'\n')
			out_num_uu.write(str(len(read_uu_het_spots))+'\n')
