#!/usr/bin/env python3
import pandas as pd
import os
import glob
import toytree
import toyplot
import toyplot.pdf
import numpy as np
import shutil
import argparse
from subprocess import call as unix
from ete3 import Tree
from Bio import SeqIO
from joblib import Parallel, delayed

# Author: Peter Mulhair
# Date: 14/01/2020 
# Usage: python3 barcode_queries.py --barcode <files_to_consider>
# Alternative usage: python3 barcode_queries.py --query <species.fasta> <species>

'''
Pipeline to build phylogenetic trees to address
query barcode species from DToL data that have 
been assigned to the incorrect species from the
BOLD database
'''


parse = argparse.ArgumentParser()

parse.add_argument("--barcode", type=str, help="name barcode results file to parse")
parse.add_argument("--query", action="append", nargs="+")

args = parse.parse_args()

query_list = []
#func_pass = 0


def data_download(barcode_results):
    '''
    Function called when using --barcode flag. Parse DToL 
    barcode excel sheet to find cases where species 
    identifation do not match barcode sequence.Download 
    barcode data from boldsystems.org for queries species.
    '''
    
    #Convert excel format to csv for parsing
    if args.barcode.lower().endswith('.xlsx'):
        barcode_results = args.barcode.split('.xlsx')[0]
        data_xls = pd.read_excel(args.barcode, dtype=str, index_col=None)
        data_xls.to_csv(barcode_results + '.csv', encoding='utf-8', index=False)
    elif args.barcode.lower().endswith('.csv'):
        barcode_results = args.barcode.split('.csv')[0]
    else:
        print('Error:')
        print('Excel or csv file required as input for --barcode')
        return None

    os.makedirs('output', exist_ok=True)
    os.makedirs('queries', exist_ok=True)
    #Download data for query species from BOLD database
    with open(barcode_results + '.csv') as f:
        next(f)
        for line in f:
            lines = line.split(',')
            specimen_ID = lines[2]
            expected_sp = lines[6].title() + '_' + lines[7]
            result_sp = lines[10].title() + '_' + lines[11]
            DNA = lines[15]
           
            flag = lines[14]
            if flag == 'C':#Currently just parsing cases where barcode is flagged as yellow and has at least genus level hit
                if result_sp.split('_') != ['', '']:
                    with open('barphy_results.csv', 'a+') as outBar:
                        outBar.write(specimen_ID + ',' + expected_sp + ',' + result_sp + '\n')

                    if expected_sp.split('_')[1] == 'sp':
                        with open('queries/' + specimen_ID + '_' + expected_sp + '.fasta', 'w') as outF_B:
                            outF_B.write('>' + expected_sp + '|' + specimen_ID + '\n' + DNA + '\n')

                    else:
                        #Save specimen IDs to list
                        query_list.append(specimen_ID)

                        try:
                            os.mkdir('output/' + specimen_ID)
                        except:
                            shutil.rmtree('output/' + specimen_ID)
                            os.mkdir('output/' + specimen_ID)

                        os.chdir('output/' + specimen_ID)

                        with open(specimen_ID + '_query.fasta', 'w') as outF:
                            outF.write('>' + expected_sp + '_DToL' + '\n' + DNA + '\n')

                        expected_sp_BOLD = lines[6].title() + '%20' + lines[7]
                        result_sp_BOLD = lines[10].title() + '%20' + lines[11]

                        #Download barcode sequences from BOLD
                        print('Downloading barcode data from BOLD (boldsystems.org)...')
                        unix('wget -q --show-progress  http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + expected_sp_BOLD, shell=True)
                        unix('wget -q --show-progress http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + result_sp_BOLD, shell=True)

                        for bold in glob.glob('sequence*'):
                            sp_fasta = bold.split('=')[1]
                            sp_fasta = sp_fasta.replace(' ', '_')
                            bold = bold.replace(' ', '\ ')
                            unix('mv ' + bold + ' ' + sp_fasta + '.fas', shell=True)

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
                                if ('COI' in k) or ('DToL' in k):
                                    outF.write('>' + k + '\n' + v + '\n')

                        os.chdir('../../')
                    
            elif flag == 'B':
                with open('barphy_results.csv', 'a+') as outBar:
                    outBar.write(specimen_ID + ',' + expected_sp + ',' + result_sp + '\n')
                with open('queries/' + specimen_ID + '_' + expected_sp + '.fasta', 'w') as outF_B:
                    outF_B.write('>' + expected_sp + '|' + specimen_ID + '\n' + DNA + '\n')
                             


def data_download_user(query, species):
    '''
    Function called when using --query flag. Requires fasta
    file of query species barcode, along with name of species
    or genus to search against from BOLD database. Gets barcode
    sequence for the query species and related species specified
    '''
    
    os.makedirs('output', exist_ok=True)
    sp_query = query.split('.fa')[0]
    if len(species) > 1:
        try:
            os.mkdir('output/' + sp_query + '_genus_query')
        except:
            shutil.rmtree('output/' + sp_query + '_genus_query')
            os.mkdir('output/' + sp_query + '_genus_query')
        os.chdir('output/' + sp_query + '_genus_query')
    else:
        if species.split('_')[1] == '':
            try:
                os.mkdir('output/' + sp_query + '_genus_query')
            except:
                shutil.rmtree('output/' + sp_query + '_genus_query')
                os.mkdir('output/' + sp_query + '_genus_query')
            os.chdir('output/' + sp_query + '_genus_query')
        else:
            try:
                os.mkdir('output/' + sp_query + '_query')
            except:
                shutil.rmtree('output/' + sp_query + '_query')
                os.mkdir('output/' + sp_query + '_query')
            os.chdir('output/' + sp_query + '_query')

                                                                            

    with open('../../queries/' + query) as f, open(sp_query + '_query.fasta', 'w') as outF:
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.description
            seq = str(record.seq)
            
            outF.write('>' + ID + '_query' + '\n' + seq + '\n')
            
    expected_sp_BOLD = sp_query.split('_')[0].title() + '%20' + sp_query.split('_')[1]
    if len(species) > 1:
        for sp_name in species:
            sister_sp_BOLD = sp_name.split('_')[0].title() + '%20' + sp_name.split('_')[1]
            unix('wget http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + sister_sp_BOLD, shell=True)
    else:
        sister_sp_BOLD = species.split('_')[0].title() + '%20' + species.split('_')[1]
        #Download barcode sequences from BOLD
        #unix('wget http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + expected_sp_BOLD, shell=True)
        unix('wget http://www.boldsystems.org/index.php/API_Public/sequence?taxon=' + sister_sp_BOLD, shell=True)
    
    for bold in glob.glob('sequence*'):
        sp_fasta = bold.split('=')[1]
        sp_fasta = sp_fasta.replace(' ', '_')
        bold = bold.replace(' ', '\ ')
        unix('mv ' + bold + ' ' + sp_fasta + '.fas', shell=True)

    #Combine BOLD barcode seqs with DToL query barcode seq 
    unix('cat *fas >> ' + sp_query + '_query.fasta', shell=True)

    #Remove duplicate sequences from the fasta file
    with open(sp_query + '_query.fasta') as f, open(sp_query + '_query.fa', 'w') as outF:
        bold_id_dict = {}
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.description
            seq = str(record.seq)
            if ID not in bold_id_dict.keys():
                bold_id_dict[ID] = seq
                
        for k, v in bold_id_dict.items():
            if ('COI' in k) or ('query' in k):
                outF.write('>' + k + '\n' + v + '\n')
                                                                                                                                                                                                                                            
    os.chdir('../../')


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
      
    
def tree_build(query):
    '''
    Created multiple sequence alignments from
    barcode fasta files using mafft. Create
    phylogenetic tree using IQTree
    '''
    pwd = os.getcwd()
    os.chdir(query)
    
    #Align using mafft and infer gene tree with IQTree
    for fa in glob.glob("*fa"):
        unix('sed -i "s/ /_/g" ' + fa, shell=True)
        print('\n')
        print('Constructing MSA and phylogenetic tree...')
        print('(If this step fails, see log files for what went wrong)')
        print('\n')
        unix('mafft --quiet ' + fa + ' > ' + fa.split('.')[0] + '.mft', shell=True)
        unix('trimal -in ' + fa.split('.')[0] + '.mft -out ' + fa.split('.')[0] + '.trim.mft -gappyout', shell=True)
        #unix('mafft --maxiterate 1000 --localpair ' + fa + ' > ' + fa.split('.')[0] + '.mft', shell=True)
        #unix('iqtree -s ' + fa.split('.')[0] + '.trim.mft', shell=True)
        unix('iqtree -quiet -m GTR+G -s ' + fa.split('.')[0] + '.trim.mft', shell=True)

    #Midpoint rooting to root the gene tree
    for tree in glob.glob("*.treefile"):
        midpoint_root(tree)
    
    #Create pdf file with tree image
    for rooted_tree in glob.glob("*.rooted"):
        queryID = query.split('/')[-1]

        rtre = toytree.tree(rooted_tree)
        Nnodes = rtre.nnodes

        colorlist = ["#de2d26" if ("query" in tip) or ("DToL" in tip) else "#000000" for tip in rtre.get_tip_labels()]
        #colorlist = ["#de2d26" if "DToL" in tip else "#000000" for tip in rtre.get_tip_labels()]
        if Nnodes < 80:
            canvas, axes, mark  = rtre.draw(tip_labels_align=True, tip_labels_colors=colorlist, width=1000, height=1000, tip_labels_style={"font-size": "15px"});
        elif Nnodes < 600:
            #canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, edge_widths=0.1, layout='c', edge_type='p', width=800, height=800, tip_labels_style={"font-size": "2px"});
            canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, tip_labels_align=True, width=1000, height=5000, tip_labels_style={"font-size": "10px"});
        else:
            canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, tip_labels_align=True, width=1000, height=8000, tip_labels_style={"font-size": "5px"});
            #canvas, axes, mark  = rtre.draw(tip_labels_colors=colorlist, edge_widths=0.1, layout='c', edge_type='p', width=600, height=600, tip_labels_style={"font-size": "1px"});
            
        toyplot.pdf.render(canvas, queryID + "_tree.pdf")

        #unix('Rscript ' + pwd + '/plot_tree.R -t ' + rooted_tree + ' -o ' + queryID + '_tree.pdf > /dev/null 2>&1', shell=True)
        #unix('Rscript ' + pwd + '/plot_tree.R -t ' + rooted_tree + ' -o ' + queryID + '_tree.pdf', shell=True)
    os.chdir(pwd)

#tree_build('output/gbroad_wasp_bracon_genus_query')

if args.barcode:#If parsing barcode excel file, run pipeline
    data_download(args.barcode)
    
    Parallel(n_jobs=10)(delayed(tree_build)('output/' + query) for query in query_list)
    print('\n\n\nComplete; output directories for batch job  can be found in output/')
    
else:#If using specific queries, run pipeline
    query_sp = args.query[0][0]
    query_sp2 = args.query[0][1:]
    
    query_sp = query_sp.split('/')[-1]

    data_download_user(query_sp, query_sp2)

    if len(query_sp2) > 1:
        tree_build('output/' + query_sp.split('.fa')[0] + '_genus_query')
        print('\nComplete; output files can be found in output/' + query_sp.split('.fa')[0] + '_genus_query')
    else:
        if query_sp2.split('_')[1] == '':
            tree_build('output/' + query_sp.split('.fa')[0] + '_genus_query')
            print('\nComplete; output files can be found in output/' + query_sp.split('.fa')[0] + '_genus_query')
        else:
            tree_build('output/' + query_sp.split('.fa')[0] + '_query')
            print('\nComplete; output files can be found in output/' + query_sp.split('.fa')[0] + '_query')
        
