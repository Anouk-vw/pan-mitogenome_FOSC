"""
Input : GFA
output: 
- core/unique/acc nodes
- Length of different categories
- Number of different variants
- plots of different variants

(mpi) anouk@tertiair:~/anouk2/mt_foc_2:python Scripts/graph_stats.py graph_v3/output/mito_final/pggb_out/all_genomes.fasta.fc357dd.7bdde5a.9a43c27.smooth.1.gfa graph_v3/output/mito_final/pggb_out/depth.csv graph_v3/output/mito_final/pggb_out/deconstructed.vcf 

"""

from sys import argv 
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

def get_node_lengths(graph, num_genomes):
    node_length_list = []
    core_length = 0
    core_length_nodup = 0
    with open(graph) as gfa_file:
        for line in gfa_file:
            if line.startswith('S'):
                s, id, seq, dp, rc = line.split('\t')
                s, i, depth = dp.split(':')
                node_length_list.append(len(seq))
    return node_length_list


def coreness(node_depth, num_genomes, graph):
    """
    node_depth = .csv file from odgi depth -d
    
    """
    #calculate node depth with odgi
    node_depth_df = pd.read_csv(node_depth, sep = '\t', index_col = 0)
    node_depth_df['node_length'] = get_node_lengths(graph, num_genomes)
    
    #[1] == uniq
    unique_length = list(node_depth_df.groupby('depth.uniq').sum()['node_length'])[1]
    unique_count = list(node_depth_df.groupby('depth.uniq').count()['node_length'])[1]
    #[-1] == core
    core_length = list(node_depth_df.groupby('depth.uniq').sum()['node_length'])[-1]
    core_count = list(node_depth_df.groupby('depth.uniq').count()['node_length'])[-1]
    #[1:-1] == accessory
    accessory_length = node_depth_df.groupby('depth.uniq').sum()[1:-1].node_length.sum()
    accessory_count = node_depth_df.groupby('depth.uniq').count()[1:-1].node_length.sum()
    
    print(f'core_len {core_length} core_nodes {core_count} \n\
    accessory_len {accessory_length} accessory_nodes {accessory_count} \n\
    unique_len {unique_length} unique_nodes {unique_count}')
    
    return 

"""
Parse variants
"""

def parse_vcf(vcf_path):
    """
    vcf file obtained with
    vg autoindex ... 
    vg deconstruct ... 
    """
    lines = []
    with open(vcf_path) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            else:
                line_list = line.strip().split('\t')
                genome = line_list[0]
                pos = line_list[1]
                nodes = line_list[2]
                ref= line_list[3]
                alt = line_list[4]
                x = line_list[5]
                strand = line_list[6]
                info = line_list[7]
                gt = line_list[7:]
                lines.append(line.strip().split('\t'))
    return lines

def extract_variant(line):
    #genome, pos, nodes, ref, alt, x, strand, info, gt  = line
    nodes = line[2]
    ref = line[3]
    alt = line[4]
    nodes_l = re.split('>|<', nodes)[1:]
    if len(ref) == 1 and len(alt) == 1:
        variant = 'SNP'
        len_var = 1
    elif len(ref) > 10 or len(alt) > 10 or (len(ref) > 10 and len(alt) > 10):
        variant = 'indel'
        len_var = len(alt)
    else:
        variant = 'small_indel'
        len_var = len(alt)
    return (nodes_l, variant, len_var)

#is variant LVR of Core? -> make lists
def list_vars(vcf_lines, LVR_node):
    core_lens = []
    LVR_lens = []
    core = []
    LVR  = []
    variants = []
    for line in vcf_lines:
        nodes, variant, len_var = extract_variant(line)
        variants.append(variant)
        if int(nodes[0]) > LVR_node or int(nodes[1]) > LVR_node:
            LVR_lens.append(len_var)
            LVR.append(variant)
        else:
            core_lens.append(len_var)
            core.append(variant)
    print(np.unique(variants, return_counts = True), len(variants), np.unique(core,return_counts = True), np.unique(LVR,return_counts = True))
    return LVR, core, LVR_lens, core_lens

def plot_variants(vcf_path, LVR_node):
    lines_vcf = parse_vcf(vcf_path)
    num_LVR, num_core, len_LVR, len_core = list_vars(lines_vcf, LVR_node)
    
    plt.subplots(figsize = (4,3))
    plt.bar([x*0.18 for x in range(0,6,2)], [x / 31.918 for x in np.unique(num_LVR, return_counts = True)[1]], width = 0.18, color = 'grey')
    plt.bar([x*0.18 for x in range(1,7,2)], [x / 50.788 for x in np.unique(num_core, return_counts = True)[1]], width = 0.18, color = 'lightgrey')
    plt.xticks([x*0.2 for x in range(0,6,2)], ['SNP', 'small indel\n(<10 bp)', 'indel'])
    plt.box(False)
    plt.grid()
    plt.savefig('variant_per_kbp_LVR_core.pdf')
    
    plt.subplots(figsize = (4,3))
    plt.bar([x*0.2 for x in range(0,3)], [x for x in np.unique(num_LVR + num_core, return_counts = True)[1]], width = 0.15, color = 'lightgrey')
    plt.xticks([x*0.2 for x in range(0,3)], ['SNP', 'small indel\n(<10 bp)', 'indel'])
    plt.box(False)
    plt.grid()
    plt.savefig('variant_distribution.pdf')
    
if __name__ == '__main__':
    graph_file = argv[1]
    depth_file = argv[2]
    vcf_file = argv[3]
    
    coreness(depth_file, 472, graph_file)
    plot_variants(vcf_file,6750)
