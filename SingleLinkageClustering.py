#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: alexander lucaci
@Temple University
"""
# =============================================================================
# Imports
# =============================================================================
from os import listdir
from os.path import isfile, join
from datetime import datetime
import sys, multiprocessing
from multiprocessing import Pool

# =============================================================================
# Declares
# =============================================================================
startTime = datetime.now()
count = 0
numlines = 1000
families = []
datalist, thresholded_datalist = [], []
#ID_Threshold = 95.0
#GAP_Threshold = 2.0
ID_Threshold = float(sys.argv[1])
GAP_Threshold = float(sys.argv[2])
#MASTER_FAMILY, MASTER_FAMILYHolder = [], []
tag = "SLC" + str(int(ID_Threshold)) + "_" + str(int(GAP_Threshold))

# =============================================================================
# Helper functions
# =============================================================================
def dup_check(count, prot_id, id_list, mode):
    global families
    #print("checking:", id1)
    occurences = 0
    protid_list = [] 
    
    for x in range(len(families)):
        if prot_id in families[x]:
            #print(id1, x)
            occurences += 1
            #id1_dict[id1] =  id1_dict[id1] + x
            protid_list.append(x)
            
    if occurences > 1:
        print(count, mode, prot_id, id_list, occurences, protid_list)
        #dont just print, combine the families. remove duplicates
        families[protid_list[0]] = list(set(families[protid_list[0]] + families[protid_list[1]]))
        families.pop(protid_list[1])        

def Load_Datalist(fname):
    print("# Loading datalist")
    global count, datalist, thresholded_datalist
    with open(fname, "r") as f:
        while True:
            a = f.readline()        
            if a == "":
                break
            # --- DEBUG --- #
            #if count == numlines:
            #    break
            if a[0] == "#":
                #print(a[0])
                continue
            if a[:2] == "ID1":
                #print(a[2])
                print("--- error ---", a[:2])
                continue
            if count > 5: #skip the header information.    
                #print("working")
                b = a.split(",")
                if b[0] == "ID1": #skip header information
                    continue
                datalist.append(b)
                
                #threshold the values here. x% id, x% gap
                try:
                    if float(b[2]) > ID_Threshold and float(b[3]) < GAP_Threshold:
                    #print("here")
                        thresholded_datalist.append(b)
                except:
                    print(b, "--- ERROR ---")
            count += 1
	
def UsedWithin(fname): #Homage to Dr. Russ Hermansan
    global thresholded_datalist, families
    print("#  Starting to build families:", fname)

    #How many pairs do we have?
    num_datapairs = len(thresholded_datalist)
    
    #Lets set the initial family
    if len(families) == 0:
        families.append([thresholded_datalist[0][0], thresholded_datalist[0][1]])
 
    #Lets start real work. Loop through each protein datapair (ID1, ID2, PERCENT ID, PERCENT GAP)
    #Check for its occurence in all existing families
    for a in range(1, num_datapairs): #skips the first family
        id1 = thresholded_datalist[a][0]
        id2 = thresholded_datalist[a][1]
        id1_list, id2_list = [], []

        if id1 == id2: #if they are the same protein, go to the next pair.
            continue

        for b in range(len(families)):
            if id1 in families[b]:
                #print(id1, "{is in family}", b, "so append", id2, "to family")
                id1_list.append(b)
                
            if id2 in families[b]:
                #print(id2, "{is in family}", b, "so append", id1, "to family")
                id2_list.append(b)

        #print(id1, id2, id1_list, id2_list)
        job = ""
        
        #all families have been checked	
        if len(id1_list) == 0 and len(id2_list) == 0:
             #no matches, create new family
             families.append([id1,id2])
             continue
         
        if len(id1_list) == 1:
            #prevent duplicates
            if id2 in families[id1_list[0]]:
                continue
            families[id1_list[0]].append(id2)
            job += "in id1_list=1,"
        
        if len(id2_list) == 1:
            #prevent duplicates
            if id1 in families[id2_list[0]]:
                continue
            families[id2_list[0]].append(id1)
            job += "in id2_list=1,"
            
        #dup_check(id1, id1_list, "id1")
        #dup_check(id2, id2_list, "id2")

        #if len(id1_list) == 1 and len(id2_list) == 1:
            #UsedWithin(fname)
            #dup_check(id1, id1_list)
            #dup_check(id2, id2_list)

        if len(id1_list) > 1:
            #print(a, id1_list, id2_list)
            if id1_list == id2_list:
                continue
            #if the family in id1 matches the family in id2
            
            for m in id1_list[1:]:
                families[id1_list[0]].append(families[m])
           
            p_count = 0
            for p in id1_list[1:]:
                families.pop(p-p_count)
                p_count += 1
            job += "in id1_list>1,"
        
        if len(id2_list) > 1:
            #print(a, id1_list, id2_list)
            if id1_list == id2_list:
                continue
            
            for m in id2_list[1:]:
                families[id2_list[0]].append(families[m])
            
            p_count = 0
            for p in id2_list[1:]:
                families.pop(p-p_count)
                p_count += 1
          
            job += "in id2_list>1,"

        dup_check(a, id1, id1_list, "id1")
        dup_check(a, id2, id2_list, "id2")                
        #print(id1, id2, id1_list, id2_list, job)            
         
# =============================================================================
# Main program subroutine
# =============================================================================
def main(fname):
    global datalist, thresholded_datalist, families, tag
    print()
    print("# Starting SLC:", fname)
    Load_Datalist(fname) #meaning, the master datafile with all of the mafft psa results. 36 of 36.
    print("# Datalist has been loaded.", datetime.now()-startTime)   
    
    #what does the data look like?
    print("# Number of data pairs (Original):", len(datalist))
    print("# Number of data pairs (Thresholded):", len(thresholded_datalist), "of", len(datalist))
    
    #lets remove unthresholded datalist from memory
    datalist = []
           
    #The build families function
    UsedWithin(fname)
       
    print("# Building complete")
    print("# Printing families to file")
    
    import csv
    with open(tag + ".csv", "w+") as f:
        wr = csv.writer(f)
        for item in families:
            wr.writerow(item) 

           #f.write(str(families[i]) + "\n")
    
    f.close()
    
    print("## Summary Stats ##")
    print("## Number of familes:", len(families))
    print("## Total Runtime ##", datetime.now() - startTime)
    print()   

              
# =============================================================================
# Main program starts here.
# =============================================================================
mypath = "/usr/home/l/125/tuk13147/SLC/MASTER/"
pool = Pool(multiprocessing.cpu_count())
pool_args = [mypath+"MAFFT_PSA_36_of_36.csv"]

pool.starmap(main, [pool_args])

#main("MAFFT_PSA_36_of_36.csv")


# =============================================================================
# End of file
# =============================================================================
#Udiv__TRINITY_DN65382_c4_g2_i1.1 Lh3tf001413g0p [8705] [8705, 10912]
#what happens in above?
#should be to add id1 to 10912.


"""
<Sample input data>
    
#Blast file:,/home/alexander/Projects/Spiders/Blastp/Blastp_results/Aarg_TAED_pro.fa_all-vs-allothers_REDOAarg_TAED_pro_blast_db.csv,
#Species 1 FASTA:,/home/alexander/Projects/Spiders/Blastp/Fasta_Originals/Aarg_TAED_pro.fa,
#Species 2 FASTA:,/home/alexander/Projects/Spiders/Blastp/Fasta_Originals/Aarg_TAED_pro.fa,
# of pairs,641063,641063,
## Starting alignments,
ID1,ID2,PERCENT ID,PERCENT GAP,
Aarg__TRINITY_DN100_c0_g1_i1.1,Aarg__TRINITY_DN100_c0_g1_i1.1,100.0,0.0,
Aarg__TRINITY_DN100012_c0_g1_i1.1,Aarg__TRINITY_DN100012_c0_g1_i1.1,100.0,0.0,
Aarg__TRINITY_DN100032_c0_g1_i1.1,Aarg__TRINITY_DN100032_c0_g1_i1.1,100.0,0.0,

"""
