import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

# =========================
# 1. Style and font settings
# =========================
sns.set(style="white")

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 18,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
})

plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.it'] = 'Arial:italic'

# =========================
# 2. Define labels and data
# =========================
mh_labels_ppi = ["MH1", "MH2", "MH3", "MH4", "MH5", "MH6", "MH7"]
g_labels = ["EIS1", "EIS2", "EIS3", "EIS4", "EIS5", "EIS6"]

mh_labels_pdi = ["MH1", "MH2", "MH3", "MH4", "MH5", "MH6", "MH7"]
d_labels = ["EIS1", "EIS2", "EIS3", "EIS4", "EIS5", "EIS6"]

datasets_ppi = [
    {'title': r'$\mathbf{J(A,B)}$', 'data': [("MH2","EIS2",0.071,1),("MH5","EIS5",0.250,2),("MH7","EIS3",0.125,1)], 'mh_labels': mh_labels_ppi, 'cluster_labels': g_labels, 'cmap': 'Oranges'},
    {'title': r'$\mathbf{J(A,B_{\pm1})}$', 'data': [("MH1","EIS2",0.040,1),("MH2","EIS2",0.154,4),("MH4","EIS4",0.091,2),("MH5","EIS5",0.111,2),("MH7","EIS3",0.300,3),("MH7","EIS4",0.037,1)], 'mh_labels': mh_labels_ppi, 'cluster_labels': g_labels, 'cmap': 'Oranges'},
    {'title': r'$\mathbf{J(A,B_{5\mathrm{\AA}})}$', 'data': [("MH1","EIS1",0.021,1),("MH1","EIS2",0.016,1),("MH2","EIS2",0.117,7),("MH2","EIS5",0.038,2),("MH3","EIS1",0.021,1),("MH3","EIS2",0.033,2),("MH4","EIS4",0.033,2),("MH5","EIS2",0.032,2),("MH5","EIS5",0.040,2),("MH7","EIS1",0.019,1),("MH7","EIS3",0.174,4),("MH7","EIS4",0.047,3),("MH7","EIS6",0.052,3)], 'mh_labels': mh_labels_ppi, 'cluster_labels': g_labels, 'cmap': 'Oranges'}
]

datasets_pdi = [
    {'title': r'$\mathbf{J(A,B)}$', 'data': [("MH1","EIS5",0.250,1),("MH2","EIS6",0.100,1),("MH5","EIS4",0.071,1),("MH6","EIS4",0.071,1)], 'mh_labels': mh_labels_pdi, 'cluster_labels': d_labels, 'cmap': 'Blues'},
    {'title': r'$\mathbf{J(A,B_{\pm1})}$', 'data': [("MH1","EIS5",0.125,1),("MH2","EIS6",0.133,2),("MH4","EIS1",0.273,3),("MH5","EIS4",0.045,1),("MH6","EIS4",0.045,1)], 'mh_labels': mh_labels_pdi, 'cluster_labels': d_labels, 'cmap': 'Blues'},
    {'title': r'$\mathbf{J(A,B_{5\mathrm{\AA}})}$', 'data': [("MH4","EIS1",0.150,3),("MH5","EIS4",0.040,2),("MH6","EIS4",0.061,3),("MH1","EIS5",0.059,1),("MH2","EIS6",0.138,4)], 'mh_labels': mh_labels_pdi, 'cluster_labels': d_labels, 'cmap': 'Blues'}
]

ppi_spatial_datasets = [
    ("3Å", [("MH1","EIS2",0.031,1),("MH2","EIS2",0.121,4),("MH3","EIS1",0.045,1),("MH3","EIS2",0.031,1),("MH4","EIS4",0.071,2),("MH5","EIS5",0.083,2),("MH7","EIS3",0.286,4),("MH7","EIS4",0.030,1),("MH7","EIS6",0.033,1)]),
    ("5Å", [("MH1","EIS1",0.021,1),("MH1","EIS2",0.016,1),("MH2","EIS2",0.117,7),("MH2","EIS5",0.038,2),("MH3","EIS1",0.021,1),("MH3","EIS2",0.033,2),("MH4","EIS4",0.033,2),("MH5","EIS2",0.032,2),("MH5","EIS5",0.040,2),("MH7","EIS1",0.019,1),("MH7","EIS3",0.174,4),("MH7","EIS4",0.047,3),("MH7","EIS6",0.052,3)]),
    ("7Å", [("MH1","EIS1",0.014,1),("MH2","EIS1",0.027,2),("MH3","EIS1",0.029,2),("MH7","EIS1",0.027,2),("MH1","EIS2",0.012,1),("MH2","EIS2",0.084,7),("MH3","EIS2",0.024,2),("MH5","EIS2",0.023,2),("MH7","EIS3",0.194,6),("MH4","EIS4",0.027,2),("MH7","EIS4",0.051,4),("MH2","EIS5",0.042,3),("MH5","EIS5",0.043,3),("MH6","EIS6",0.025,2),("MH7","EIS6",0.063,5)])
]

pdi_spatial_datasets = [
    ("3Å", [("MH4","EIS1",0.273,3),("MH5","EIS4",0.036,1),("MH6","EIS4",0.036,1),("MH1","EIS5",0.091,1),("MH2","EIS6",0.125,2)]),
    ("5Å", [("MH4","EIS1",0.150,3),("MH5","EIS4",0.040,2),("MH6","EIS4",0.061,3),("MH1","EIS5",0.059,1),("MH2","EIS6",0.138,4)]),
    ("7Å", [("MH4","EIS1",0.125,3),("MH7","EIS2",0.059,2),("MH4","EIS3",0.019,1),("MH5","EIS4",0.045,3),("MH6","EIS4",0.045,3),("MH1","EIS5",0.045,1),("MH2","EIS6",0.098,4),("MH3","EIS6",0.025,1)])
]

# ============================================================
# Core plotting function: left MH labels + shared colorbar on the right
# ============================================================
def plot_combined_heatmaps(datasets, filename, cmap_name, is_spatial=False):
    """
    datasets: list of dict or list of (cutoff, data) tuples
    filename: output file name
    cmap_name: 'Oranges' or 'Blues'
    is_spatial: True for spatial cutoff plots
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Get uniform labels
    if is_spatial:
        mh_labels = mh_labels_ppi if 'EIS' in filename else mh_labels_pdi
        cluster_labels = g_labels if 'EIS' in filename else d_labels
    else:
        mh_labels = datasets[0]['mh_labels']
        cluster_labels = datasets[0]['cluster_labels']

    for idx, item in enumerate(datasets):
        ax = axes[idx]

        # Parse data
        if is_spatial:
            cutoff, data = item
            title = rf'$\mathbf{{J(A,B_{{{cutoff}}})}}$'
        else:
            data = item['data']
            title = item['title']

        # Build matrix
        J = np.zeros((len(mh_labels), len(cluster_labels)))
        counts = np.zeros_like(J)
        for entry in data:
            mh, cl, val, c = entry
            i = mh_labels.index(mh)
            j = cluster_labels.index(cl)
            J[i, j] = val
            counts[i, j] = c

        # Build annotations
        annot = np.empty_like(J).astype(object)
        for i in range(J.shape[0]):
            for j in range(J.shape[1]):
                annot[i, j] = f"{int(counts[i, j])}" if J[i, j] > 0 else ""

        # Draw heatmap
        sns.heatmap(
            J, annot=annot, fmt="", cmap=cmap_name,
            vmin=0, vmax=1, linewidths=0.5, linecolor='gray',
            cbar=False, ax=ax
        )

        # Thicken subplot borders
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(2.0)
            spine.set_edgecolor('black')

        ax.set_title(title, fontsize=18, pad=10, fontweight='bold')
        ax.set_xticklabels(cluster_labels, rotation=0, ha='center')

        # Only the first subplot retains left y-axis labels; others hidden
        if idx == 0:
            ax.set_yticklabels(mh_labels, rotation=0, fontsize=16)
            ax.tick_params(axis='y', length=0)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis='y', length=0)

    # ---- Layout adjustments ----
    # left=0.06: reserve space for left MH labels
    # right=0.88: reserve space for right colorbar
    # wspace=0.10: compact subplot spacing
    plt.subplots_adjust(left=0.06, right=0.88, wspace=0.10)

    # ---- Shared colorbar (placed on the far right) ----
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_name)
    sm.set_array([])

    cbar = fig.colorbar(
        sm, ax=axes.ravel().tolist(),
        label='Jaccard similarity score',
        orientation='vertical',
        pad=0.02,
        shrink=0.8
    )
    cbar.ax.tick_params(labelsize=16)

    plt.savefig(filename, dpi=300, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    print(f"Saved: {filename}")

# ============================================================
# Execute plotting
# ============================================================
print("Starting to generate all heatmaps...")

plot_combined_heatmaps(datasets_ppi, 'PPI_combined_heatmap.tif', 'Oranges')
plot_combined_heatmaps(datasets_pdi, 'PDI_combined_heatmap.tif', 'Blues')
plot_combined_heatmaps(ppi_spatial_datasets, 'PPI_spatial_cutoffs.tif', 'Oranges', is_spatial=True)
plot_combined_heatmaps(pdi_spatial_datasets, 'PDI_spatial_cutoffs.tif', 'Blues', is_spatial=True)

print("\nAll heatmaps successfully generated!")