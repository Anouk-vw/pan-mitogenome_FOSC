import re
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from ete3 import NCBITaxa
import pandas as pd
import numpy as np
import seaborn as sns
import io
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def reorder_lists(*lists):
    """
    Make sure the list of genes starts with NAD2
    """
    # Find the index of 'NAD2' in each list
    nad2_indices = []
    for lst in lists:
        try:
            nad2_indices.append(lst.index('nad2'))
        except ValueError:
            nad2_indices.append(0)
    # Determine the maximum index of 'NAD2'
    #max_index = max(nad2_indices)

    # Reorder each list, starting from NAD2
    reordered_list = []
    for lst, nad2_index in zip(lists, nad2_indices):
        if nad2_index != 0:
            reordered_list = lst[nad2_index:] + lst[:nad2_index]
        else:
            reordered_list = lst
        # Add additional elements to reach the length of the longest list
        #reordered_list += [None] * (max_index - len(reordered_list))
        #reordered_lists.append(reordered_list)

    return reordered_list

def orientation_similarity(list1, reversed_list1, target_list):
    """
    check forward/reverse strand
    Make sure the list orientation is most similar to the others
    for list 1 and the reverse of list 1, find most simialr to target list
    """
    # Calculate similarity scores for forward orientation
    forward_score = sum(1 for x, y in zip(list1, target_list) if x == y) / len(target_list)

    # Calculate similarity scores for reverse orientation of list2
    #reversed_list1 = list1[::-1]
    reverse_score = sum(1 for x, y in zip(reversed_list1, target_list) if x == y) / len(target_list)

    if forward_score > reverse_score:
        return False
    elif forward_score < reverse_score:
        return True
    elif forward_score == reverse_score:
        return True

def add_empy_genes(list1, target_list):
    temp_list = list1.copy()
    discrepancy = len(target_list) - len(list1)
    if discrepancy > 0:
        #unequal list size
        for x_add in range(0, discrepancy):
                for i, disc in enumerate(range(0, discrepancy+1)):
                    temp_list.insert(i, 'x')
                    if orientation_similarity(list1, temp_list, target_list):
                        list1 = temp_list.copy()
                        discrepancy = discrepancy -1
                    else:
                        #keep previous list
                        temp_list = list1.copy()
    print(len(list1), discrepancy)
    return list1

def fetch_genome_size(faidx_in):
    with open(faidx_in) as faidx:
        for line in faidx:
            return int(line.split('\t')[1])

def get_taxids(acc_file):
    
    ncbi = NCBITaxa()
    acc_tax_dict = {}
    with open(acc_file) as acc_tax:
        for line in acc_tax:
            acc, tax = line.strip().split('\t')
            acc_tax_dict[acc] = tax
    return acc_tax_dict

def fetch_genus(ncbi_taxid, level):
    # Create an instance of NCBITaxa
    ncbi = NCBITaxa()

    # Fetch taxonomic lineage
    lineage = ncbi.get_lineage(ncbi_taxid)

    # Get taxonomic ranks for each taxid in the lineage
    lineage_ranks = ncbi.get_rank(lineage)

    # Find the taxid for genus
    genus_taxid = None
    for taxid in lineage:
        rank = lineage_ranks[taxid]
        if rank == level:
            genus_taxid = taxid
            break

    # If genus taxid is found, return its name
    if genus_taxid:
        genus_name = ncbi.get_taxid_translator([genus_taxid])[genus_taxid]
        return genus_name
    else:
        return None

def get_leave_order(acc_tax_dict):
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(acc_tax_dict.values())
    leaf_order = []
    for leaf in tree:
        leaf_order.append(str(leaf).replace('-', '').strip())
    return leaf_order
    
def get_elements(gff_files, Fusarium_order):
    elements_dict = {} 
    for gff_file in gff_files:
        with open(gff_file) as gff:
            genome = re.match('(.*\.[0-9])_[A-Z]', gff_file.split('/')[-1]).group(1)
            elements = []
            for line in gff:
                if 'mRNA' in line:
                    try:
                        m = re.match('(.*Name=)(.*);tr(.*)', line)
                        if m.group(2) in Fusarium_order and m.group(2) not in elements:
                            elements.append(m.group(2))
                        else:
                            continue
                    except AttributeError:
                        continue

            #should genes be reversed?
            if orientation_similarity(reorder_lists(elements), reorder_lists(elements[::-1]), Fusarium_order):
                elements_dict[genome] = reorder_lists(elements[::-1])
            else:
                elements_dict[genome] = reorder_lists(elements)
    return elements_dict

def make_heatmap(elements_dict, acc_tax_dict, Fusarium_order):
    leaf_order = get_leave_order(acc_tax_dict)
    heatmap_dict = {}
    for genome, genes in elements_dict.items():
        color_gene= []
        genes_to_plot = []
        #sizes = []
        heatmap_dict[genome] = []
        for g in genes:
            try:
                color_gene.append(color_dict[g])
                genes_to_plot.append(Fusarium_order.index(g))
                #sizes.append(10)
            except KeyError:
                continue

        if len(genes_to_plot) != len(Fusarium_order):
            absent_genes = len(Fusarium_order)- len(genes_to_plot)
            genes_to_plot.extend([np.nan for x in range(0,absent_genes)])
            
        heatmap_dict[genome] = genes_to_plot

    hm_plot = pd.DataFrame(heatmap_dict).T
    ordered = {}
    for genome in hm_plot.index:
        if genome.split('.')[0] in acc_tax_dict.keys():
            taxid_genome = acc_tax_dict[genome.split('.')[0]]
            try:
                order = leaf_order.index(taxid_genome)
                ordered[order] = genome
            except ValueError:
                continue
    reindex_dict = dict(sorted(ordered.items()))
    hm_plot_re = hm_plot.reindex(reindex_dict.values())
    genus_index = []
    for x in reindex_dict.values():
        genus_index.append(fetch_genus(acc_tax_dict[x.split('.')[0]], 'class'))
    hm_plot_re.index = genus_index
    return hm_plot_re

def plot_heatmap(heatmap_to_plot):
    remove_scer = [False if 'cerevisiae' in x else True for x in heatmap_to_plot.index]
    remove_scer[remove_scer.index(False)] = True
    
    genus_index = []
    for x in reindex_dict.values():
        genus_index.append(fetch_genus(acc_tax_dict[x.split('.')[0]], 'class'))

    g = sns.clustermap(hm_plot_re.iloc[2:,0:14][remove_scer[2:]], cmap = 'rainbow',
                   row_cluster = False, col_cluster = False, square = True, linecolor = 'white')#, figsize = (40,400))

    g.ax_heatmap.yaxis.tick_left()
    plt.show()
    return
    
    
if __name__ == '__main__':
    gff_files = [f'/home/anouk/anouk2/mt_foc_2/graph_analysis/fonesca_genomes/annotation/Supplementary_FileS1/{x}' for x in os.listdir('/home/anouk/anouk2/mt_foc_2/graph_analysis/fonesca_genomes/annotation/Supplementary_FileS1/') if x.endswith('gff3')]
    Fusarium_order = ['nad2',
                     'nad3',
                     'atp9',
                     'cox2',
                     'nad4L',
                     'nad5',
                     'cob',
                     'cox1',
                     'nad1',
                     'nad4',
                     'atp8',
                     'atp6',
                     'cox3',
                     'nad6']
 
    color_dict = {'nad2': (0.5, 0.0, 1.0, 1.0),
                 'nad3': (0.3588235294117647, 0.2199463578396686, 0.9938591368952737, 1.0),
                 'atp9': (0.21764705882352942, 0.42912060877260894, 0.9755119679804366, 1.0),
                 'cox2': (0.07647058823529412, 0.6172782212897929, 0.9451838281608196, 1.0),
                 'nad4L': (0.0725490196078431, 0.7829276104921027, 0.9005867023006374, 1.0),
                 'nad5': (0.21372549019607845, 0.9005867023006374, 0.8469582108244671, 1.0),
                 'cob': (0.3549019607843137, 0.9741386021045101, 0.7829276104921028, 1.0),
                 'cox1': (0.503921568627451, 0.9999810273487268, 0.7049255469061472, 1.0),
                 'nad1': (0.6450980392156862, 0.9741386021045102, 0.622112816721474, 1.0),
                 'nad4': (0.7862745098039214, 0.9005867023006376, 0.5316594672504363, 1.0),
                 'atp8': (0.9274509803921569, 0.7829276104921029, 0.43467642176596505, 1.0),
                 'atp6': (1.0, 0.6172782212897929, 0.3265387128400833, 1.0),
                 'cox3': (1.0, 0.4291206087726091, 0.2199463578396687, 1.0),
                 'nad6': (1.0, 0.21994635783966857, 0.11065268189150082, 1.0), 'x': (0,0,0,0)}

    acc_tax_dict = get_taxids('/home/anouk/anouk2/mt_foc_2/graph_analysis/fonesca_genomes/annotation/accession_w_taxid.txt')
    elements_dict = get_elements(gff_files, Fusarium_order)   
    #print([len(v) for v in elements_dict.values()])
    plot_heatmap(make_heatmap(elements_dict, acc_tax_dict, Fusarium_order))
