# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:08:35 2018

@author: Viktoria
"""

from __future__ import print_function

import argparse

id2gt = {} # Dictionary/hash table for genotypes with variant id as key
samples2index = {} #Hash tables for sample indexes, indices as key
samples = [] #List of samples
four_samples = [] #List of samples in four file
id2info = {} # Hash table for info (AD)

#Argument parser for input
parser = argparse.ArgumentParser()
parser.add_argument('vcf', type=argparse.FileType('r'), help="VCF file with variants")
parser.add_argument('four', type=argparse.FileType('r'), help="Four file with imputed genotypes")
parser.add_argument('--min_gq', type=int, default=0, help="Minimum GQ value")
parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
parser.add_argument('outfile', type=argparse.FileType('w'), help="Output file")
#args = parser.parse_args("test.vcf test.four --min_gq 30 output.txt".split(" "))
args = parser.parse_args()

for line in args.vcf:
    #Skip all lines that start with ##
    if line.startswith("##"):
        continue
    #Extract the header
    elif line.startswith("#"):
        header = line.rstrip().split("\t")
        for s in range(9, len(header)):
            samples2index[header[s]] = s-9 #Index of each sample found and put into hash table
            samples.append(header[s]) #Sample names appended to the list samples
        continue
        
        
    columns = line.split("\t") #Seperate columns
    variant_id = columns[2] #Current variant id
    id2gt[variant_id] = [-1]*(len(columns) - 9) #Initialize the id2gt values
    gq_index = columns[8].split(":").index("GQ") #Find index of GQ value
    id2info[variant_id] = [-1]*(len(columns) - 9) #initialize id2info
    ad_index = columns[8].split(":").index("AD") #Find index of AD value
        
    for s in range(9, len(columns)):
        gt = columns[s].split(":")[0] #Current genotype
        gq = columns[s].split(":")[gq_index] #Current GQ value
        temp_ad = columns[s].split(":")[ad_index] #Current AD value
        ad = temp_ad.split(",") #Split REF AD and ALT AD into two columns
        id2info[variant_id][s-9] = ad #Update values of ad for each variant to current AD value
        #Only want to check samples that have a high enough GQ value
        if args.min_gq > int(gq):
            continue
        else:
            id2gt[variant_id][s-9] = gt.count("1") #Update values of gt for each varian id to 0, 1 or 2 (how many changes)

#Initilize arrays for counting of genotypes zeros, ones and twos 
zeros = [0, 0, 0]
ones = [0, 0, 0]
twos = [0, 0, 0]
#Initilize variables for counting of correct, incorrect and no data variants
correct_count = 0
incorrect_count = 0
nodata_count = 0
prev_variant_id = ""
print("Variant", "Correct", "Incorrect", "No Data", "Incorrect sample with highest GQ") #Print header for table
mismatch = {} #Hash table for mismatched variant
summary = {} #Hsh table for correct/incorrect/no data information for variant

args.outfile.write("GT") #Skrifa inn header fyrir outfile
args.outfile.write("\t")
args.outfile.write("REF_AD")
args.outfile.write("\t")
args.outfile.write("ALT_AD")
args.outfile.write("\t")
args.outfile.write("Variant")
args.outfile.write("\t")
args.outfile.write("Sample")
args.outfile.write("\n")
    
for line in args.four:
    four_columns = line.rstrip().split("\t") #Find columns in four file
    sample_name = four_columns[0] #Sample names from four file
    variant_id = four_columns[1] #New current variant id
    gt1 = four_columns[2] #Genotype 1
    gt2 = four_columns[3] #Genotype 2
    
    #Start counting for each new variant_id
    if variant_id != prev_variant_id:
        #Update prev-variant_id to the first variant_id
        if prev_variant_id == "":
            prev_variant_id = variant_id
            mismatch[variant_id] = [False] * len(samples2index) #Initialize all mismatches as false
        else:
            summary[prev_variant_id] = [correct_count, incorrect_count, nodata_count] #Add counts for each variant to hash table
            prev_variant_id = variant_id #Update to next variant id
            #Initialize counts and mismatch for next variant id
            correct_count = 0
            incorrect_count = 0
            nodata_count = 0
            mismatch[variant_id] = [False] * len(samples2index)
        
    #Skip sample if it's not in vcf file  
    if sample_name not in samples2index:
        continue
    
    #Skip variant if itðs not in vcf file    
    if variant_id not in id2gt:
        continue
        
    sample_index = samples2index[sample_name] #Extract relevant sample index from vcf file for current sample in four file
    variant_gt = id2gt[variant_id] #Extract relevant genotype from vcf file for the current variant in four file
    
    #Only want to compare and count valid genotypes in four file (0 or 1)    
    if variant_gt[sample_index] != -1 and ((gt1 == '0') or (gt1 == '1')) and ((gt2 == '0') or (gt2 == '1')):        
        #Change genotypes from string to int
        gt1 = int(gt1)
        gt2 = int(gt2)
        #Count zeros, ones and twos and add to count arrays
        if (gt1 + gt2) == 0:
            zeros[variant_gt[sample_index]] +=1
        elif (gt1 + gt2) == 1:
            ones[variant_gt[sample_index]] +=1
        else:
            twos[variant_gt[sample_index]] +=1
            
        gt = gt1 + gt2
        ad = id2info[variant_id][sample_index]
        args.outfile.write(str(gt)) #Skrifa inn genatýpu
        args.outfile.write("\t")
        args.outfile.write(ad[0]) #Skrifa inn AD REF
        args.outfile.write("\t")
        args.outfile.write(ad[1]) #Skrifa inn AD ALT
        args.outfile.write("\t")
        args.outfile.write(variant_id) #Skrifa inn variant
        args.outfile.write("\t")
        args.outfile.write(sample_name) #Skrifa inn sample
        args.outfile.write("\n")
        
        #Check if the genotype (0, 1 or 2) is the same or not  in vcf and four file and count instances    
        if variant_gt[sample_index] == (gt1 + gt2):
            correct_count +=1
            if args.verbose:
                print(sample_name, "is the same") #Optional print for match
        else:
            incorrect_count +=1
            #For mismatched samples: change mismatch value at that sample index to true
            mismatch[variant_id][sample_index] = True      
            if args.verbose:
                print("Mismatch on", variant_id, "in", sample_name, sample_index) #Optional print for mismatch
    
    #Count where data is missing or is not valid    
    else:
        nodata_count +=1
        if args.verbose:
            print("Data missing") #Optional data for missing data

#Update summary hash table counts for the last variant id
summary[prev_variant_id] = [correct_count, incorrect_count, nodata_count]

# Seek to beginning of file            
args.vcf.seek(0, 0)
            
for line in args.vcf:
    #Skip all lines starting with #
    if line.startswith("#"):
        continue
    
    #Initialize samples with highest GQ and the GQ value
    highest_sample = ""
    highest_gq = 0
            
    columns = line.split("\t") #Columns in vcf file
    variant_id = columns[2] #Current variant id
    gq_index = columns[8].split(":").index("GQ") #Index of current GQ value
    
    #Skip variants that aren't in vcf or four
    if variant_id not in id2gt or variant_id not in mismatch: 
        continue
    
    for s in range(9, len(columns)):
        sample_index = s - 9 #Indices of samples in vcf
        gt = columns[s].split(":")[0] #Current genoype 
        gq = columns[s].split(":")[gq_index] #Current GQ value

        #For all mismatched values: find highest one
        if mismatch[variant_id][sample_index] == True:
            if int(gq) >= highest_gq:
                highest_sample = samples[sample_index] #Update new high sample name
                highest_gq = int(gq) #Update new highest GQ value
                    
    #For all variants in summary hash table: print summary counts as well as sample with highest GQ value
    if variant_id in summary:
        print(variant_id, " ".join([str(x) for x in summary[variant_id]]), highest_sample, highest_gq)
                      
#Print table that shows counts of zeros, ones and twos in vcf (columns) and four (lines) file
print("  ", "0", "1", "2", sep = "  ")
print("0:", zeros)
print("1:", ones)
print("2:", twos)