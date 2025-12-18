import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
import re
from Bio import Phylo
import matplotlib.patches as mpatches
import networkx as nx

def generate_plots(df, output_prefix):
    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    distances = df.where(mask).stack().values
    
    plt.figure(figsize=(10, 6))
    sns.histplot(distances, bins=20, kde=True, color="skyblue")
    plt.title("Distribution of Allelic Distances")
    plt.xlabel("Allelic Distance")
    plt.savefig(f"{output_prefix}_histogram.png")
    plt.close()
    
    n_samples = df.shape[0]
    fig_dim = max(12, n_samples * 0.6)
    plt.figure(figsize=(fig_dim, fig_dim))
    font_size = 10 if n_samples < 20 else 8
    
    sns.heatmap(
        df, cmap="viridis", annot=True, fmt="d", square=True,
        cbar_kws={"shrink": 0.8}, annot_kws={"size": font_size}
    )
    plt.title("Allelic Distance Heatmap")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_heatmap.png")
    plt.close()

    plt.figure(figsize=(8, 6))
    sns.violinplot(y=distances, color="lightblue", inner="quartile")
    plt.title("Plot of Allelic Distances")
    plt.ylabel("Allelic Distance")
    plt.savefig(f"{output_prefix}_violin.png")
    plt.close()

def get_group_colors(tree, meta_df):
    group_map = {}
    for _, row in meta_df.iterrows():
        if 'mlst_st' in row and pd.notna(row['mlst_st']) and str(row['mlst_st']) != "Unknown":
            grp = str(row['mlst_st'])
        elif pd.notna(row.get('patient_id')) and str(row.get('patient_id')) != "NA":
            grp = str(row.get('patient_id'))
        else:
            grp = "Unknown"
            
        group_map[row['sample_id']] = grp

    unique_groups = sorted(list(set(group_map.values())))
    
    if len(unique_groups) <= 10:
        palette = sns.color_palette("bright", len(unique_groups))
    else:
        palette = sns.color_palette("husl", len(unique_groups))
        
    color_map_mpl = {grp: color for grp, color in zip(unique_groups, palette)}
    return group_map, color_map_mpl, unique_groups

def plot_rectangular_tree(tree, group_map, color_map_mpl, unique_groups, output_file):
    def to_bio_color(mpl_color):
        return tuple(int(x * 255) for x in mpl_color)
    
    color_map_bio = {k: to_bio_color(v) for k, v in color_map_mpl.items()}
    gray_bio = (128, 128, 128)

    def color_clade(clade):
        if clade.is_terminal():
            l = group_map.get(clade.name, "Unknown")
            c = color_map_bio.get(l, gray_bio)
            clade.color = c
            return c
        else:
            child_colors = [color_clade(c) for c in clade]
            first_color = child_colors[0]
            if all(c == first_color for c in child_colors):
                clade.color = first_color
                return first_color
            else:
                clade.color = gray_bio 
                return gray_bio

    color_clade(tree.root)

    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False,
        label_func=lambda x: x.name if x.is_terminal() else "",
        branch_labels=lambda c: int(c.branch_length) if c.branch_length else ""
    )
    
    terminals = tree.get_terminals()
    for i, clade in enumerate(terminals):
        y_pos = i + 1
        x_pos = tree.distance(tree.root, clade)
        l = group_map.get(clade.name, "Unknown")
        c_mpl = color_map_mpl.get(l, "gray")
        ax.scatter(x_pos, y_pos, color=c_mpl, s=80, zorder=10, edgecolors='white', linewidth=0.5)

    handles = [mpatches.Patch(color=color_map_mpl[c], label=f"{c}") for c in unique_groups]
    plt.legend(handles=handles, title="Sequence Type", loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    
    plt.title("GrapeTree MST (Rectangular)", fontsize=14)
    plt.xlabel("Allelic Distance (Differing Loci)", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_force_directed_tree(tree, group_map, color_map_mpl, unique_groups, output_file):
    G = nx.Graph()
    
    def add_edges(clade):
        for child in clade:
            weight = child.branch_length if child.branch_length else 1.0
            G.add_edge(clade, child, weight=weight)
            add_edges(child)
            
    add_edges(tree.root)
    
    try:
        pos = nx.kamada_kawai_layout(G, weight='weight')
    except:
        pos = nx.spring_layout(G)
        
    plt.figure(figsize=(12, 12))
    
    terminals = [n for n in G.nodes() if hasattr(n, 'is_terminal') and n.is_terminal()]
    internals = [n for n in G.nodes() if n not in terminals]
    
    nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.4, width=1)
    
    edge_labels = {}
    for u, v, d in G.edges(data=True):
        w = d.get('weight', 0)
        if w > 0:
            edge_labels[(u, v)] = int(w)
            
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=7, font_color='black')
    
    nx.draw_networkx_nodes(G, pos, nodelist=internals, node_size=15, node_color='lightgray')
    
    for group in unique_groups:
        node_list = [n for n in terminals if group_map.get(n.name, "Unknown") == group]
        if node_list:
            color = color_map_mpl.get(group, "gray")
            nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_size=150, node_color=[color], label=group)
        
    labels = {n: n.name for n in terminals}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_color='black')
    
    plt.axis('off')
    plt.legend(title="Sequence Type", loc='best', frameon=False)
    plt.title("GrapeTree MST", fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_phylo_trees(tree_file, meta_df, output_prefix):
    try:
        tree = Phylo.read(tree_file, "newick")
        tree.ladderize()
    except Exception as e:
        print(f"Error reading tree file: {e}")
        return

    group_map, color_map, unique_groups = get_group_colors(tree, meta_df)
    
    plot_rectangular_tree(tree, group_map, color_map, unique_groups, f"{output_prefix}_rectangular.png")
    plot_force_directed_tree(tree, group_map, color_map, unique_groups, f"{output_prefix}_radial.png")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help="Path to distance_matrix.tsv")
    parser.add_argument('--metadata', required=True, help="Path to metadata.tsv")
    parser.add_argument('--tree', required=False, help="Path to grapetree.nwk")
    parser.add_argument('--threshold', type=int, default=12, help="Allele threshold")
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', index_col=0)
    df.index.name = "Sample"
    
    meta_df = pd.read_csv(args.metadata, sep='\t')

    generate_plots(df, "stats")
    
    if args.tree:
        generate_phylo_trees(args.tree, meta_df, "grapetree_viz")

if __name__ == "__main__":
    main()