import sys


#Ex.
# python ./check_where_het_rmvd.py het_rmvd_loc_test.txt ~/Desktop/smooth_het/simulation.chr20.gs28100000.cov30.het4.rs10000.random_het1.gc-1.rc.err1.repeats0.noErrRecorded.fasta ~/Downloads/loc_replacements_simulation.gs28100000_cov30_het4_rs10000_random_het1_rc_k21_gc-1_err1_repeats0.txt
# python ./check_where_het_rmvd.py het_rmvd_loc_test.txt ~/Desktop/smooth_het/p_arm_ch20/simulation.chr20.gs28100000.cov30.het1.rs10000.random_het1.gc-1.rc.err1.repeats0.fasta ~/Desktop/smooth_het/loc_replacements_simulation.gs28100000_cov30_het1_rs10000_random_het1_rc_k21_gc-1_err1_repeats0.txt
# python ./check_where_het_rmvd.py het_rmvd_loc_test.txt ~/Desktop/smooth_het/simulation.chr20.gs27981280.cov30.het0.rs10000.random_het1.gc-1.rc.err1.repeats0.fasta ~/Desktop/smooth_het/loc_files/loc_replacements_simulation.gs27981280_cov30_het0_rs10000_random_het1_rc_k21_gc-1_err1_repeats0.txt

out = open(sys.argv[1], "w")
out_all_het = open(sys.argv[2], "w") # bed file with the list of all het sites 
out_missed_het = open(sys.argv[3], "w") #bed file with the list of the missed 

with open(sys.argv[4]) as f: #Fasta file
	all_het_spots = []
	all_err_spots = []
	starting_spots = []
	Lines = f.readlines()
	for line in Lines:
		line = line.strip()
		if (line.count(">") > 0):
			het_spots=line.split("|")[4].split(",")
			all_het_spots.append(het_spots)
			err_spots=line.split("|")[5].split(",")
			starting_spots.append(int(line.split("|")[1]))
			for ele in het_spots:
				if ele != "":
					out_all_het.write("chr_20\t" + str(int(line.split("|")[1]) + int(ele)) + "\t" + str(int(line.split("|")[1]) + int(ele) + 1) + "\n")
			all_err_spots.append(err_spots)

with open(sys.argv[5]) as f:
	i = 0
	Lines = f.readlines()
	for line in Lines:
		line = line.strip()
		rmvd_het = line.split(",")
		if rmvd_het[0] == '':
			union_het = 0
		else:
			union_het = len(set(rmvd_het) & set(all_het_spots[i]))
			all_het_union = set(rmvd_het) & set(all_het_spots[i])
			missed_het = set(all_het_spots[i]).difference(rmvd_het)
			for each in all_het_union:
				out.write("chr_20\t" + str(starting_spots[i]+int(each)) + "\t" + str(starting_spots[i]+int(each) + 1) + "\n")
			for each in missed_het:
				if each != '':
					out_missed_het.write("chr_20\t" + str(starting_spots[i]+int(each)) + "\t" + str(starting_spots[i]+int(each) + 1) + "\n")
		union_err = len(set(rmvd_het) & set(all_err_spots[i]))
		missed_het = len(all_het_spots[i]) - union_het
		false_rmvd = len(rmvd_het) - union_het
		if all_het_spots[i] == '':
			length_het_spots = 0
		else:
			length_het_spots = len(all_het_spots[i])
		print("read_" + str(i) + "\t" + str(union_het/length_het_spots*100) + "\t" + str(union_err))
		i = i + 1

out_missed_het.close()
out_all_het.close()
out.close()




