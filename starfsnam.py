# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:08:35 2018

@author: Viktoria
"""

id2gt = {} # Dictionary/hash table
samples2index = {}
four_samples = []

with open("test.vcf", "r") as f:
    for line in f:
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            header = line.rstrip().split("\t")
            for s in range(9, len(header)):
                samples2index[header[s]] = s-9
                
            #print(samples2index)
            
        #print (line.rstrip().split("\t"))
        
        
        columns = line.split("\t")
        variant_id = columns[2]
        id2gt[variant_id] = [0]*(len(columns) - 9)
        
        for s in range(9, len(columns)):
            gt=columns[s].split(":")[0]
            id2gt[variant_id][s-9] = gt.count("1")
            
        #print (id2gt)
        
with open("test.four", "r") as g:
    zeros = [0, 0, 0]
    ones = [0, 0, 0]
    twos = [0, 0, 0]
    
    for line in g:
        four_columns = line.rstrip().split("\t")
        sample_no = four_columns[0]
        variant = four_columns[1]
        gt1 = four_columns[2]
        gt2 = four_columns[3]
        
        
        if sample_no not in samples2index:
            continue
        
        if variant not in id2gt:
            continue
        
        sample_index = samples2index[sample_no]
        variant_gt = id2gt[variant]
        

        if ((gt1 == '0') or (gt1 == '1')) and ((gt2 == '0') or (gt2 == '1')):
            gt1 = int(gt1)
            gt2 = int(gt2)
            
            if (gt1 + gt2) == 0:
                zeros[variant_gt[sample_index]] +=1
            elif (gt1 + gt2) == 1:
                ones[variant_gt[sample_index]] +=1
            else:
                twos[variant_gt[sample_index]] +=1
            
            
            if variant_gt[sample_index] == (gt1 + gt2):
                print(sample_no, "is the same")
            else:
                print("Mismatch on", variant, "in", sample_no)
        else:
            print("Data missing")
                     
        
        
print("  ", "0", "1", "2", sep = "  ")
print("0:", zeros)
print("1:", ones)
print("2:", twos)