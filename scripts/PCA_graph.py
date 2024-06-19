"""
Script to make the PCA from the node-matrix
"""
from sys import argv
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.colors as pltcol
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def read_matrix(matrix_file):
    """
    Initiate matrix
    """
    matrix = pd.read_csv(matrix_file)
    matrix.index = [x.split('.fas')[0] for x in matrix.index]
    return matrix

def Run_PCA(matrix):
    """
    Run PCA
    """
    #PCA on full matrix
    pca_all = PCA(n_components=20)
    projected = pca_all.fit_transform(abs(matrix.iloc[:,:]))
    return projected, pca_all

def group_clusters(matrix, PCA_projection, betweenxmin, betweenxmax, betweenymin, betweenymax, category, value):
    """
    Group clusters in PCA in a category based on boundaries
    """
    matrix.loc[(PCA_projection[:,0] > betweenxmin) & (PCA_projection[:,0] < betweenxmax) &\
    (PCA_projection[:,1] < betweenymax) & (PCA_projection[:,1] > betweenymin), category] = value
    return matrix

def color_PCA(matrix):
    """
    Get list of colors based on categories in matrix
        Categories can be defined by group_clusters
        
    TODO: Categories and colors are given in the function
    """
    dict_col_coreClust = {'core 1':'#ECD2DE','core 2':'#CBE5BE','core 3':'#A7CAEB', 'core 4':'#95AB91', 'core 5':'#F2B2D0'}
    dict_col_LVRClust = {'LVR 1':'#8A98CD','LVR 2':'#F04D6D','LVR 3': '#FCCF4C'}

    color_core = []
    for x in list(matrix['core']):
        try:
            color_core.append(dict_col_coreClust[x])
        except KeyError:
            color_core.append('black')
        
    color_LVR = []
    for x in list(matrix['LVR']):
        try:
            color_LVR.append(dict_col_LVRClust[x])
        except KeyError:
            color_LVR.append('black')
    return color_LVR, color_core

def plot_PCA(PCA_projections, pcas, colors):
    """
    Plot PCAs in the list, in a row
    """
    fig, ax = plt.subplots(1,len(PCA_projections), figsize = (15,5))
    for i, pp in enumerate(PCA_projections):
        ax[i].scatter(pp[:, 0], pp[:, 1], edgecolor='none', alpha=0.5, c= color_LVR,
                cmap=plt.cm.get_cmap('Paired', 6))
        ax[i].set_xlabel(f'PC 1 ({round(pcas[i].explained_variance_ratio_[0],2) * 100}%)')
        ax[i].set_ylabel(f'PC 3 ({round(pcas[i].explained_variance_ratio_[1],2) * 100}%)')
        ax[i].set_title('matrix all')
        ax[i].set(frame_on=False)
        ax[i].grid(True, color='grey', linestyle=(0,(5,5)), linewidth=0.5, alpha = 0.5)
    fig.tight_layout()
    plt.show()
    return

def make_figure(matrix_raw):
    #specify splitting node
    LV_matrix = matrix.loc[:,'6784':]
    core_matrix = matrix.loc[:,:'6784']
    
    #run all PCAs
    core_project, pca_core = Run_PCA(core_matrix)
    all_project, pca_all = Run_PCA(matrix)
    LVR_project, pca_LVR = Run_PCA(LV_matrix)
    
    #add all core categories
    matrix = group_clusters(matrix, core_project, -10, -2.5, -5, 5, 'core', 'core 1')
    matrix = group_clusters(matrix, core_project, 7, 12, 5, 20, 'core', 'core 2')
    matrix = group_clusters(matrix, core_project, 10, 20, -10, 0, 'core', 'core 3')
    matrix = group_clusters(matrix, core_project, 0, 5, 2.5, 7.5, 'core', 'core 4')
    
    #add all LVR categories
    matrix = group_clusters(matrix, LVR_project, -15, 0, -10, 10, 'LVR', 'LVR 1')
    matrix = group_clusters(matrix, LVR_project, 50, 70, -15, 10, 'LVR', 'LVR 2')
    matrix = group_clusters(matrix, LVR_project, 20, 40, 40, 55, 'LVR', 'LVR 3')
    
    #specify color for categories
    color_LVR, color_core = color_PCA(matrix)
    
    #plot PCA
    #to plot PCA with core-colors, add color_core instead of color_LVR
    plot_PCA([LVR_project, core_project, all_project], [pca_LVR, pca_core, pca_all], color_LVR)
    
    return
    
if __name__ == '__main__':
    matrix_in = argv[1]
    matrix = read_matrix(matrix_in)
    
    make_figure(matrix)

