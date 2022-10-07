#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 13:44:17 2022

@author: nkleinbo
"""

from argparse import ArgumentParser
import os

#common prefix for two strings:
def common_prefix(arr1, arr2):
    prefix = [];
    i = 0
    #don't use "other entries"
    if(arr2[0] == "other entries"):
        return arr1
    while i <= len(arr1) - 1 and i <= len(arr2) - 1:    
        if (arr1[i] != arr2[i]):
            break
        prefix.append(arr1[i])
        i += 1
    return (prefix)

#lca for a list of taxids:
def get_lca_from_list(list_of_taxids, taxid_lineage_mapping):
    lca_string = taxid_lineage_mapping[list_of_taxids[0]]
    lca_array = lca_string.split("; ")
    
    #don't use "other entries", replace with second taxid, if available
    if(len(list_of_taxids) > 1 and lca_array[0] == "other entries"):
        lca_string = taxid_lineage_mapping[list_of_taxids[1]]
        lca_array = lca_string.split("; ")
    
    #print("Searching LCA")
    print (lca_string)
    
    for i in range (1, len(list_of_taxids)):
        print(taxid_lineage_mapping[list_of_taxids[i]])
        lca_array = common_prefix(lca_array, taxid_lineage_mapping[list_of_taxids[i]].split("; "))
        
    #print("Result:")
    print("; ".join(lca_array))
    return (lca_array)

#read the blast result file, compute lca and write results to new lca_file
def read_blast_results_and_store_lca(blast_file, lca_file, taxid_lineage_mapping, percent_of_best_score):
    with open (blast_file) as blast:
        with open(lca_file, "w") as lca_f:
            best_score = 0
            read_id = ""
            taxid_array = []
            best_hit = []
            for hit_line in blast:
                hit_line = hit_line.rstrip()
                hit = hit_line.split("\t")
                #skip when taxid not in mapping
                if(hit[12] not in taxid_lineage_mapping):
                    continue
                if(hit[0] != read_id):
                    #store results of previous read:
                    if(len(taxid_array) > 0):
                        print(taxid_array)
                        lca_array = get_lca_from_list(taxid_array, taxid_lineage_mapping)
                        print(lca_array)
                        lca = "; ".join(lca_array)
                        lca_name = lca_array[-1]
                        best_hit.append(lca_name)
                        best_hit.append(lca)
                        lca_f.write("\t".join(best_hit)+"\n")
                    #new read:
                    best_score = float(hit[11])
                    best_hit = hit
                    read_id = hit[0]
                    taxid_array = []
                    taxid_array.append(hit[12])
                else:
                    #further hits for read, push to array, when within score range:
                    if(float(hit[11]) >= (percent_of_best_score / 100.0) * best_score and hit[12] not in taxid_array):
                        taxid_array.append(hit[12])
                
            #store last LCA
            lca_array = get_lca_from_list(taxid_array, taxid_lineage_mapping)
            lca = "; ".join(lca_array)
            lca_name = lca_array[-1]
            best_hit.append(lca_name)
            best_hit.append(lca)
            lca_f.write("\t".join(best_hit)+"\n")
                
    
#read the taxid to lineage mapping file
def read_lineage_mapping(path_to_mapping_file):
    taxid_lineage_mapping = {}
    taxid_name_mapping = {}
    count = 0
    print ("Reading lineage mapping file ... ")
    with open(path_to_mapping_file) as f:
        for line in f:
            count  += 1
            line = line.strip()
            spl = line.split("|")
            taxid_name_mapping[spl[0].strip()] = spl[1].strip()
            taxid_lineage_mapping[spl[0].strip()] = spl[2].strip()+" "+spl[1].strip()
    print ("Done.")
    return (taxid_lineage_mapping, taxid_name_mapping)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        metavar="BLAST output", help="Path to the BLAST or DIAMOND input file, needs to be in tabular format without headers. Query id is expected in column 1, Accession number in column 2, taxid in column 13")
    parser.add_argument("-o", "--output", dest="output", help="Path to output file", required=True,
                        metavar="output file")
    parser.add_argument("-l", "--lineage", dest="lineage", help="Path to lineage file (fullnamelineage.dmp if extracted from NCBI new_taxdump.tar file)", required=True,
                        metavar="lineage file")
    parser.add_argument("-p", "--percent", dest="percent", help="Percent of best score to be considered in LCA computation)", required=False, default=100.0, type=float,
                        metavar="percent of best score")
    parser.add_argument("-f", "--force", dest="force", help="Forces overwriting the output file if it already exists, otherwise the program quits when the file exists.", required=False,
                        action="store_true")
#    parser.add_argument("-q", "--query-col", dest="query_col", default=1, type=int, metavar="column of query id",
#                        help="The column number in the input file containing the query id (beginning with 1)")
#    parser.add_argument("-t", "--taxa-col", dest="taxa_col", default=2, type=int, metavar="column of taxid",
#                        help="")
    
    args = parser.parse_args()
    
    #check if output file exists and remove if --force is set
    if(os.path.exists(args.output)):
        if(args.force):
            os.system("rm "+args.output)
        else:
            exit("The specified output file "+args.output+" already exists.")
    
    #filter the BLAST  results for bad hits:
    input_filtered = ""
    if(os.path.exists(args.input)):
        input_filtered = args.input + ".filtered"
        os.system("grep -v 'environmental' "+args.input+" | grep -v 'uncultured' | grep -v 'unidentified' > "+input_filtered)
    else:
        exit("The specified input file "+args.input+" does not exist.")

    taxid_lineage_mapping = {}
    taxid_name_mapping = {}
    if(os.path.exists(args.lineage)):
        (taxid_lineage_mapping, taxid_name_mapping) = read_lineage_mapping(args.lineage)
        #print ("Skipped reading lineage")
    else:
        exit("The specified lineage file "+args.lineage+" does not exist.")
    
    
    #compute LCA and create a new file
    print("Computing LCA...")
    read_blast_results_and_store_lca(input_filtered, args.output, taxid_lineage_mapping, args.percent)
    print("Done.")
    

 