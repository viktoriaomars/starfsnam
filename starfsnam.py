# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:08:35 2018

@author: Viktoria
"""

from __future__ import print_function

import argparse

id2gt = {} # Dictionary/hash table
samples2index = {}
samples = []
four_samples = []

parser = argparse.ArgumentParser()
parser.add_argument('vcf', type=argparse.FileType('r'), help="VCF file with variants")
parser.add_argument('four', type=argparse.FileType('r'), help="Four file with imputed genotypes")
parser.add_argument('--min_gq', type=int, default=0, help="Minimum GQ value")
parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
#args = parser.parse_args("test.vcf test.four --min_gq 30".split(" "))
args = parser.parse_args()

for line in args.vcf:
    if line.startswith("##"):
        continue
    elif line.startswith("#"):
        header = line.rstrip().split("\t")
        for s in range(9, len(header)):
            samples2index[header[s]] = s-9
            samples.append(header[s])
        continue
        
        
    columns = line.split("\t")
    variant_id = columns[2]
    id2gt[variant_id] = [-1]*(len(columns) - 9)
    gq_index = columns[8].split(":").index("GQ")
        
    for s in range(9, len(columns)):
        gt = columns[s].split(":")[0]
        gq = columns[s].split(":")[gq_index]
        if args.min_gq > int(gq):
            continue
        else:
            id2gt[variant_id][s-9] = gt.count("1")


zeros = [0, 0, 0]
ones = [0, 0, 0]
twos = [0, 0, 0]

correct_count = 0
incorrect_count = 0
nodata_count = 0
prev_variant = ""
print("Variant", "Correct", "Incorrect", "No Data", "Incorrect sample with highest GQ")
mismatch = {}
summary = {}
    
for line in args.four:
    four_columns = line.rstrip().split("\t")
    sample_no = four_columns[0]
    variant = four_columns[1]
    gt1 = four_columns[2]
    gt2 = four_columns[3]
    
    
    if variant != prev_variant:
        if prev_variant == "":
            prev_variant = variant
            mismatch[variant] = [False] * len(samples2index)
        else:
            summary[prev_variant] = [correct_count, incorrect_count, nodata_count]
            prev_variant = variant
            correct_count = 0
            incorrect_count = 0
            nodata_count = 0
            mismatch[variant] = [False] * len(samples2index)
        
        
    if sample_no not in samples2index:
        continue
        
    if variant not in id2gt:
        continue
        
    sample_index = samples2index[sample_no]
    variant_gt = id2gt[variant]
    
        
    if variant_gt[sample_index] != -1 and ((gt1 == '0') or (gt1 == '1')) and ((gt2 == '0') or (gt2 == '1')):
        gt1 = int(gt1)
        gt2 = int(gt2)
            
        if (gt1 + gt2) == 0:
            zeros[variant_gt[sample_index]] +=1
        elif (gt1 + gt2) == 1:
            ones[variant_gt[sample_index]] +=1
        else:
            twos[variant_gt[sample_index]] +=1
        
            
        if variant_gt[sample_index] == (gt1 + gt2):
            correct_count +=1
            if args.verbose:
                print(sample_no, "is the same")
        else:
            incorrect_count +=1
            mismatch[variant][sample_index] = True
            
            if args.verbose:
                print("Mismatch on", variant, "in", sample_no, sample_index)
        
    else:
        nodata_count +=1
        if args.verbose:
            print("Data missing")

summary[prev_variant] = [correct_count, incorrect_count, nodata_count]

# Seek to beginning of file            
args.vcf.seek(0, 0)
            
for line in args.vcf:
    if line.startswith("#"):
        continue
    
    
    highest_sample = ""
    highest_gq = 0
            
    columns = line.split("\t")
    variant_id = columns[2]
    gq_index = columns[8].split(":").index("GQ")
    
    if variant_id not in id2gt or variant_id not in mismatch: 
        continue
    
    for s in range(9, len(columns)):
        sample_index = s - 9
        gt = columns[s].split(":")[0]
        gq = columns[s].split(":")[gq_index]

    
        if mismatch[variant_id][sample_index] == True:
            if int(gq) >= highest_gq:
                highest_sample = samples[sample_index]
                highest_gq = int(gq)
                    
    # Gripa drasl ur usmmary hakkatöflunni
    if variant_id in summary:
        print(variant_id, " ".join([str(x) for x in summary[variant_id]]), highest_sample, highest_gq)
                      
                    
                       
            
print("  ", "0", "1", "2", sep = "  ")
print("0:", zeros)
print("1:", ones)
print("2:", twos)