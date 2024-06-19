## Create pangenome graph
Graph_v3 contains the scripts to create the pangenome graph

Config.yaml: set parameters and input files

Use run.sh to run construct the graph


Returns a graph in .gfa format, dgi statistics and a node presence absence matrix.

## analyse graph
scripts contains the scripts used to analyse the pangenome graph.

- PCA_graph.py: create a PCA from the presence absence matrix

- graph_stats.py: list number of core, unique and accessory nodes and variants in the graph.
input: matrix.csv, depth CSV from odgi depth -d and VCF from vg deconstruct -a
    
- gene_synteny_plot.py: plot gene synteny of multiple mitogenome GFF files
