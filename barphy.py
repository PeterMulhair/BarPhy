#!/usr/bin/env python3
import pandas as pd
import os
import sys
import glob
import toytree
import toyplot
import toyplot.pdf
import numpy as np
import shutil
import argparse
import subprocess
from datetime import datetime
from subprocess import call as unix
from collections import defaultdict, Counter
from ete3 import Tree
from Bio import SeqIO
from joblib import Parallel, delayed

# Author: Peter Mulhair
# Date: 14/01/2020 

'''
Pipeline to build phylogenetic trees to address
query barcode species from DToL data that have 
been assigned to the incorrect species from the
BOLD database
'''

date = datetime.today().strftime('%Y_%m_%d')

def create_directories(argv=None):
    #Create directories for output
    try:
        os.mkdir('barphy_output_' + date)
    except:
        print('WARNING: barphy_output exists, overwriting previous run now.\n')
        shutil.rmtree('barphy_output_' + date)
        os.mkdir('barphy_output_' + date)

    os.mkdir('barphy_output_' + date + '/blast_output')
    os.mkdir('barphy_output_' + date + '/phylo_output')

    
def read_input(fasta_files: str, indet_file):
    #Read fasta and identifier input files from user args
    plateID_info = defaultdict(list)
    if indet_file.endswith('.xlsx'):
        info_file = indet_file.split('.xlsx')[0]
        data_xls = pd.read_excel(indet_file, dtype=str, index_col=None)
        data_xls.to_csv(info_file + '.csv', encoding='utf-8', index=False)
        info_file = info_file + '.csv'
    elif indet_file.endswith('.csv'):
        info_file = indet_file
    else:
        print('ERROR: Excel or csv file required as input')
        return None
    
    with open(info_file) as f:
        for line in f:
            lines = line.split(',')
            plate = lines[0]
            sample = lines[1]
            ID = lines[2]
            plateID = plate + '_' + ID
            order = lines[3]
            family = lines[4]
            genus = lines[5]
            species = lines[6]
            info_list = order, family, genus, species, sample
            plateID_info[plateID].append(info_list)

    fasta_file_list = [i.strip('\n') for i in open(fasta_files)]

    return fasta_file_list, plateID_info

def blast_run(fasta_files: str, indet_info):
    '''
    Execute BLAST search and parsing
    '''
    ##Check if blast is installed localy
    if shutil.which('blastn') is None:
        sys.exit('BLAST not installed locally')
    else:
        input_files = read_input(fasta_files, indet_info)
        plateID_info = input_files[1]
        #For each fasta file run blast search against BOLD database
        print('Running BLAST search for barcode ID...\n')
        for fasta in input_files[0]:
            fasta_name = fasta.split('.')[0]
            unix('blastn -query ' + fasta + ' -db BOLD_db -evalue 1e-5 -num_threads 4 -max_target_seqs 20 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out barphy_output_' + date + '/blast_output/' + fasta_name + '_blastout.tsv', shell=True)

        ##Parse blast output
        plateID_hit = defaultdict(list)
        for blastout in glob.glob('barphy_output_' + date + '/blast_output/*tsv'):
            with open(blastout) as f:
                for line in f:
                    lines = line.split('\t')
                    plateID = lines[0]
                    BOLD_hit = lines[1]
                    percentID = lines[3]
                    species_hit = BOLD_hit.split('|')[1]
                    try:
                        species_hit = species_hit.split('_')[0] + '_' + species_hit.split('_')[1]#Correct for abberant BOLD sp ID names
                    except:
                        continue
                    if '_sp' not in species_hit:
                        hit_info = species_hit, percentID
                        plateID_hit[plateID].append(hit_info)
        
        ##Find queries with failed barcode runs or whos genus or species are absent from BOLD db
        issue_sp_list = []#List of species that pipeline cannot ID 
        barcode_ID_list = []
        for fasta in input_files[0]:
            with open(fasta) as f:
                for record in SeqIO.parse(f, 'fasta'):
                    ID = record.description
                    barcode_ID_list.append(ID)
        BOLD_sp_list = []
        BOLD_genus_list = []
        with open('../../src/barphy/BOLD_db.fasta') as f:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.description
                sp = ID.split('|')[1]
                BOLD_sp_list.append(sp)
                genus = sp.split('_')[0]
                BOLD_genus_list.append(genus)

        for plateIDs in plateID_info.keys():
            spIDs = plateID_info[plateIDs][0][3]
            sampleIDs_issue = plateID_info[plateIDs][0][4]
            if plateIDs not in barcode_ID_list:
                issue_sp_list.append(plateIDs)
            else:
                spIDs = spIDs.replace(' ', '_')
                if spIDs not in BOLD_sp_list:
                    if spIDs.split('_')[0] not in BOLD_genus_list:#Check is genus in present in BOLD db
                        issue_sp_list.append(plateIDs)
                        print(plateIDs, sampleIDs_issue, spIDs, 'Genus absent from BOLD db')
                    else:#Check if species is present in BOLD db
                        print(plateIDs, sampleIDs_issue, spIDs, 'Species absent from BOLD db')
                        issue_sp_list.append(plateIDs)
                
                
        ##Write unambiguous IDs from blast to file, return ambiguous IDs for phylo identification
        sp_mismatch_dict = defaultdict(list)
        with open('barphy_output_' + date + '/species_identification_' + date + '.csv', 'a+') as outF:
            for plateID, hits in plateID_hit.items():
                species_info = plateID_info[plateID]
                speciesID = species_info[0][3]
                speciesID = speciesID.replace(' ', '_')
                
                species_hits = []
                pass_hits = defaultdict(list)
                all_hits = {}
                for hit in hits:
                    species_hit = hit[0]
                    PID = hit[1]
                    species_hits.append(species_hit)
                    all_hits[species_hit] = PID
                    if float(PID) >= 95:
                        sp_PID_hit = species_hit, PID
                        pass_hits[species_hit].append(PID)
                        
                species_count = Counter(species_hits)

                #Write straightforward blast IDs to file
                if (speciesID == hits[0][0]) and (len(set(species_hits)) == 1) and (float(hits[0][1]) >= 95):
                    PID_hit = hits[0][1]
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')
                    
                elif (speciesID == hits[0][0]) and (float(hits[0][1]) >= 95):
                    PID_hit = hits[0][1]
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')
                    
                elif (len(hits) > 1) and (speciesID == hits[1][0]) and (float(hits[1][1]) >= 95) and (len(set(species_hits[1:])) == 1):
                    PID_hit = hits[0][1]
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')

                elif (speciesID in pass_hits.keys()) and (len(set(pass_hits.keys())) == 1):
                    PID_hit = max(pass_hits.values())
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')
                    
                elif (speciesID in pass_hits.keys()) and (len(set(pass_hits.keys())) == 2) and ((species_info[0][0] in set(pass_hits.keys())) or (species_info[0][1] in set(pass_hits.keys()))):
                    PID_hit = max(pass_hits.values())
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')
                    
                elif (speciesID in pass_hits.keys()) and (len(set(pass_hits.keys())) == 3) and (('Insecta' in set(pass_hits.keys()))) and ((species_info[0][0] in set(pass_hits.keys())) or (species_info[0][1] in set(pass_hits.keys()))):
                    PID_hit = max(pass_hits.values())
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')

                elif (speciesID in pass_hits.keys()) and (len(pass_hits[speciesID])/sum(([len(x) for x in pass_hits.values()])) >= 0.9):
                    #print(hits, 'res7')
                    PID_hit = max(pass_hits.values())
                    outF.write(plateID + ',' + species_info[0][-1] + ',' + species_info[0][0] + ',' + species_info[0][1] + ',' + species_info[0][2] + ',' + speciesID + '\n')
                    
                    
                else:#For mismatches, get level of clade to search against in a phylo ID 
                    if plateID not in issue_sp_list:

                        if speciesID.split('_')[0] in str(hits):
                            sp_query_list = []
                            if len(pass_hits) > 0:
                                for query_hit in pass_hits.keys():
                                    if '_' in query_hit:
                                        sp_query_list.append(query_hit)
                            else:
                                for query_hit in species_hits:
                                    if '_' in query_hit:
                                        sp_query_list.append(query_hit)
                                        
                                    
                            sp_query_list = set(sp_query_list)
                            sp_mismatch_dict[plateID] = list(sp_query_list)

                        elif (len(set(all_hits.keys())) == 2) and ((species_info[0][0] in set(all_hits.keys())) or (species_info[0][1] in set(all_hits.keys()))):
                            sp_query_search = speciesID.split('_')[0]
                            sp_mismatch_dict[plateID].append(sp_query_search)
                            
                        else:
                            print(plateID, speciesID, species_hits)
                            
    ##Write specimens being sent for tree ID to tsv files in queries outdir
    if len(sp_mismatch_dict) > 0:#If there are ambigbuous spIDs, send for phylo ID
        os.makedirs('barphy_output_' + date + '/phylo_output/queries', exist_ok=True)
        with open('barphy_output_' + date + '/phylo_output/queries/BOLD_queries.tsv', 'w') as outF:
            for k, v in sp_mismatch_dict.items():
                v = set(v)
                outF.write(k + '\t')
                for query in v:
                    outF.write(query + ',')
                outF.write('\n')


def data_download(fasta_files: str, indet_info):
    '''
    Function called when tree identification
    required. Downloads and filters barcode
    data from boldsystems.org for queries species.
    '''

    pwd = os.getcwd()
    
    input_files = read_input(fasta_files, indet_info)

    query_dict = defaultdict(list)
    with open('barphy_output_' + date + '/phylo_output/queries/BOLD_queries.tsv') as f:
        for line in f:
            q_list = []
            lines = line.split('\t')
            plate = lines[0]
            query_list = lines[1].strip()
            queries = query_list.split(',')
            for q in set(queries[:-1]):
                q_list.append(q)
            query_dict[plate] = q_list
    
    fasta_ID_dict = {}
    for fasta in input_files[0]:
        with open(fasta) as f:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id
                seq = str(record.seq)
                if ID in query_dict.keys():
                    fasta_ID_dict[ID] = seq

    plateID_info = input_files[1]
    
    os.chdir('barphy_output_' + date + '/phylo_output/queries')

    
    #Download data for query species from BOLD database
    for plateID, query_search in query_dict.items():
        species_info = plateID_info[plateID]
        speciesID = species_info[0][3]
        speciesID = speciesID.replace(' ', '_').replace('.', '')
        specimen_ID = species_info[0][-1]

        plate_seq = fasta_ID_dict[plateID]

        try:
            os.mkdir(specimen_ID)
        except:
            shutil.rmtree(specimen_ID)
            os.mkdir(specimen_ID)
            
        os.chdir(specimen_ID)
        
        with open(specimen_ID + '_' + speciesID + '.fas', 'w') as outF:
            outF.write('>' + plateID + '_' + specimen_ID + '_query' + '\n' + plate_seq + '\n')

        #Download barcode sequences from BOLD
        print('\nDownloading barcode data from BOLD (boldsystems.org)...')
        specID = speciesID.split('_')[0].title() + '%20' + speciesID.split('_')[1]
        unix('wget -q --show-progress  http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + specID, shell=True)
        for query_searchID in query_search:
            query_searchID = query_searchID.split('_')[0].title() + '%20' + query_searchID.split('_')[1]
            unix('wget -q --show-progress  http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + query_searchID, shell=True)

        for bold in glob.glob('sequence*'):
            query_fasta = bold.split('=')[1]
            query_fasta = query_fasta.replace(' ', '_')
            bold = bold.replace(' ', '\ ')
            unix('mv ' + bold + ' ' + query_fasta + '.fas', shell=True)

        #Combine BOLD barcode seqs with DToL query barcode seq
        unix('cat *fas >> ' + specimen_ID + '_query.fasta', shell=True)

        #Remove duplicate sequences from the fasta file
        with open(specimen_ID + '_query.fasta') as f, open(specimen_ID + '_query.fa', 'w') as outF:
            bold_id_dict = {}
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.description
                seq = str(record.seq)
                if ID not in bold_id_dict.keys():
                    bold_id_dict[ID] = seq
                    
            for k, v in bold_id_dict.items():
                if ('COI' in k) or ('query' in k):
                    outF.write('>' + k + '\n' + v + '\n')
                    
        os.chdir('../')
        
    os.chdir(pwd)
                    

def midpoint_root(tree):
    '''
    Function to root trees produced from
    tree_build, using ete3 for midpoint
    rooting
    '''
    t = Tree(tree, format=1)
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)
    #Write rooted tree to file
    t.write(format=1, outfile = tree + ".rooted")
      
    
def tree_build(specimenIDs):
    '''
    Created multiple sequence alignments from
    barcode fasta files using mafft. Create
    phylogenetic tree using IQTree.
    '''
    pwd = os.getcwd()
            
    os.chdir('barphy_output_' + date + '/phylo_output/queries/' + specimenIDs)
            
    #Align using mafft and infer gene tree with IQTree
    for fa in glob.glob("*fa"):
        unix('sed -i "s/ /_/g" ' + fa, shell=True)
        unix('mafft --quiet ' + fa + ' > ' + fa.split('.')[0] + '.mft', shell=True)  
        #unix('mafft --maxiterate 1000 --localpair ' + fa + ' > ' + fa.split('.')[0] + '.mft', shell=True)
        #unix('iqtree -quiet -s ' + fa.split('.')[0] + '.mft', shell=True)
        unix('iqtree -quiet -m GTR+G -s ' + fa.split('.')[0] + '.mft', shell=True)

    #Midpoint rooting to root the gene tree
    for tree in glob.glob("*.treefile"):
        midpoint_root(tree)

    #Create pdf file with tree image
    for rooted_tree in glob.glob("*.rooted"):
        
        rtre = toytree.tree(rooted_tree)
        Nnodes = rtre.nnodes
        
        colorlist = ["#de2d26" if ("query" in tip) or ("DToL" in tip) else "#000000" for tip in rtre.get_tip_labels()]
        if Nnodes < 80:
            canvas, axes, mark  = rtre.draw(tip_labels_align=True, tip_labels_colors=colorlist, width=1000, height=600, tip_labels_style={"font-size": "15px"});
        elif Nnodes < 200:
            canvas, axes, mark  = rtre.draw(tip_labels_align=True, tip_labels_colors=colorlist, width=1000, height=700, tip_labels_style={"font-size": "15px"});
        elif Nnodes < 600:
            #canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, edge_widths=0.1, layout='c', edge_type='p', width=800, height=800, tip_labels_style={"font-size": "2px"});
            canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, tip_labels_align=True, width=1000, height=5000, tip_labels_style={"font-size": "15px"});
        else:
            canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, tip_labels_align=True, width=1000, height=8000, tip_labels_style={"font-size": "15px"});
            #canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, edge_widths=0.1, layout='c', edge_type='p', width=600, height=600, tip_labels_style={"font-size": "1px"});
            
        toyplot.pdf.render(canvas, specimenIDs + "_tree.pdf")
                            
    os.chdir(pwd)


create_directories()    
#blast_run('fastas.txt', 'WW_12_Barcoding_Data.csv')
blast_run('fastas.txt', 'DToL_Sample_Sheet_C9X0BBWH.xlsx')

#Check if there are queries for phylo ID
if os.path.isfile('barphy_output_' + date + '/phylo_output/queries/BOLD_queries.tsv'):
    data_download('fastas.txt', 'DToL_Sample_Sheet_C9X0BBWH.xlsx')

    plate_specimen = {}
    with open('DToL_Sample_Sheet_C9X0BBWH.csv') as f:
        for line in f:
            lines = line.split(',')
            plate = lines[0]
            sample = lines[1]
            ID = lines[2]
            plateID = plate + '_' + ID
            plate_specimen[plateID] = sample

    specimen_list = []
    with open('barphy_output_' + date + '/phylo_output/queries/BOLD_queries.tsv') as f:
        for line in f:
            lines = line.split('\t')
            plateId = lines[0]
            specimenID = plate_specimen[plateId]
            specimen_list.append(specimenID)

    print('\n')
    print('Constructing MSAs and phylogenetic trees...')
    print('(If this step fails, see log files in barphy_output_' + date + '/phylo_output/queries/ for what went wrong)')
    print('\n')
    Parallel(n_jobs=10)(delayed(tree_build)(specID) for specID in specimen_list)












###########################
'''
if args.barcode:#If parsing barcode excel file, run pipeline
    data_download(args.barcode)
    
    Parallel(n_jobs=10)(delayed(tree_build)('output/' + query) for query in query_list)
    print('\n\n\nComplete; output directories for batch job  can be found in output/')
    
else:#If using specific queries, run pipeline
    query_sp = args.query[0][0]
    query_sp2 = args.query[0][1]
    
    query_sp = query_sp.split('/')[-1]
    
    data_download_user(query_sp, query_sp2)
    
    if query_sp2.split('_')[1] == '':
        tree_build('output/' + query_sp.split('.fa')[0] + '_genus_query')
        print('\nComplete; output files can be found in output/' + query_sp.split('.fa')[0] + '_genus_query')
    else:
        tree_build('output/' + query_sp.split('.fa')[0] + '_query')
        print('\nComplete; output files can be found in output/' + query_sp.split('.fa')[0] + '_query')
'''        
