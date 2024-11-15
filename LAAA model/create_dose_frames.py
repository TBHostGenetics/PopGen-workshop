#!/usr/bin/env python

import sys

allele_file_name=sys.argv[1]
vit_file_name=sys.argv[2]
snp_info_file_name=sys.argv[3]
sample_id_file_name=sys.argv[4]
begin_pos=sys.argv[5]
end_pos=sys.argv[6]
outfile_prefix=sys.argv[7]

#Find index of lines of interest
snp_info_file=open(snp_info_file_name)
i=0
for line in snp_info_file:
    el=line.strip().split()
    if el[0] == begin_pos:
        begin_pos_index=i
    if el[0] == end_pos:
        end_pos_index=i
        break
    i=i+1
snp_info_file.close()

#Open output files and write header line by processing sample ID file
allele_dose_file=open(outfile_prefix + "allele_dose.txt", "w")
afr_dose_file=open(outfile_prefix + "afr_dose.txt", "w")
allele_afr_dose_file=open(outfile_prefix + "allele_afr_dose.txt", "w")
allele_dose_file.write("position\tref\talt")
afr_dose_file.write("position\tref\talt")
allele_afr_dose_file.write("position\tref\talt")
sample_id_file=open(sample_id_file_name)
for line in sample_id_file:
    allele_dose_file.write("\t" + line.strip())
    afr_dose_file.write("\t" + line.strip())
    allele_afr_dose_file.write("\t" + line.strip())
sample_id_file.close()
allele_dose_file.write("\n")
afr_dose_file.write("\n")
allele_afr_dose_file.write("\n")

#Process the phased allele and local ancestry files
snp_info_file=open(snp_info_file_name)
allele_file=open(allele_file_name)
vit_file=open(vit_file_name)
i=0
while i <= end_pos_index:
    snp_line=snp_info_file.readline()
    allele_line=allele_file.readline()
    vit_line=vit_file.readline()
    if i >= begin_pos_index:
        #Write the snp information to all the output dose files
        snp_el=snp_line.strip().split()
        allele_dose_file.write(snp_el[0] + "\t" + snp_el[1] + "\t" + snp_el[2])
        afr_dose_file.write(snp_el[0] + "\t" + snp_el[1] + "\t" + snp_el[2])
        allele_afr_dose_file.write(snp_el[0] + "\t" + snp_el[1] + "\t" + snp_el[2])
        allele_str=allele_line.strip()
        vit_el=vit_line.strip().split() 
        if len(allele_str) != len(vit_el):
            print("ERROR!! Number of alleles and local ancestry does not match at line " + str(i))
            break
        for c in range(0, len(allele_str), 2):
            #Write allele dose - simply sum the alleles (or or 1) for the current sample
            dose=str(int(allele_str[c]) + int(allele_str[c+1]))
            allele_dose_file.write("\t" + dose)        
            #Write afr dose - subtract 1 from the local ancestry call to get 0 or 1 copies of african ancestry, 
            #and then sum the counts for the current sample
            dose=str(int(vit_el[c]) + int(vit_el[c+1]))
            afr_dose_file.write("\t" + dose)        
            #Write allele afr dose - both allele dose and afr dose on same haplotype = 1 (0 otherwise), 
            #sum the 2 results
	    dose=str((int(allele_str[c]) & (int(vit_el[c]))) + (int(allele_str[c+1]) & (int(vit_el[c+1]))))
            allele_afr_dose_file.write("\t" + dose)        
        #Write new line
        allele_dose_file.write("\n")
        afr_dose_file.write("\n")
        allele_afr_dose_file.write("\n")
    i=i+1
snp_info_file.close()
allele_file.close()
vit_file.close()
allele_dose_file.close()
afr_dose_file.close()
allele_afr_dose_file.close()

