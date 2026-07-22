"""
PPI ΔΔG prediction correlation analysis
Plots predicted vs experimental binding free energy changes for protein-protein interactions.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error


def plot_correlation(file_path):
    """
    Generate correlation plots for multiple prediction methods against experimental ΔΔG.

    Parameters
    ----------
    file_path : str
        Path to the Excel file containing the benchmark dataset.
    """
    # Read only the required columns
    required_cols = ['experimental_DDG', 'DDMut-PPI', 'MutaBind2', 'SAAMBE-3D', 'RDE-network', 'average']
    df = pd.read_excel(file_path, sheet_name="Sheet1", usecols=required_cols)

    # Data cleaning: replace common placeholders and convert to numeric
    df_clean = (
        df
        .apply(pd.to_numeric, errors='coerce')
        .dropna(how='any')
    )

    # Set plotting style
    sns.set_theme(style="white", font="Times New Roman")
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['mathtext.fontset'] = 'cm'

    # Create 2x3 subplot grid (remove the 6th subplot)
    fig, axs = plt.subplots(2, 3, figsize=(8, 5))
    axs = axs.flatten()
    fig.delaxes(axs[5])

    # Define prediction methods
    methods = [
        {'col': 'DDMut-PPI', 'title': 'DDMut-PPI'},
        {'col': 'MutaBind2', 'title': 'MutaBind2'},
        {'col': 'SAAMBE-3D', 'title': 'SAAMBE-3D'},
        {'col': 'RDE-network', 'title': 'RDE-network'},
        {'col': 'average', 'title': 'average'}
    ]

    # Experimental values
    exp = df_clean['Experimental_DDG']

    # Color scheme (consistent orange/blue)
    scatter_color = '#fdae61'
    line_color = '#d6604d'

    # Loop over each method
    for i, method in enumerate(methods):
        ax = axs[i]
        pred = df_clean[method['col']]

        # Compute statistics
        pcc, _ = pearsonr(exp, pred)
        rmse = mean_squared_error(exp, pred, squared=False)

        # Fit regression line
        slope, intercept = np.polyfit(exp, pred, 1)

        # Scatter plot with regression line
        sns.regplot(
            x=exp, y=pred,
            ax=ax,
            line_kws={'color': line_color, 'lw': 1.5},
            scatter_kws={'s': 30, 'alpha': 0.7, 'edgecolor': 'w', 'color': scatter_color}
        )

        # Add statistics text (top-left)
        text_content = f'$y = {slope:.2f}x {intercept:+.2f}$\nRMSE = {rmse:.2f}\nPCC = {pcc:.3f}'
        ax.text(
            0.05, 0.95, text_content,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
        )

        ax.set_title(f'{method["title"]}', fontsize=13)
        ax.set_xlabel('Experimental ΔΔG (kcal/mol)', fontsize=10)
        ax.set_ylabel('Predicted ΔΔG (kcal/mol)', fontsize=10)

    plt.tight_layout()
    plt.savefig('protein_protein_correlation.tif', bbox_inches='tight', dpi=600,
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.show()

    # Print summary statistics
    print(f"Number of valid data points: {len(df_clean)}")
    print(f"Experimental ΔΔG range: [{exp.min():.2f}, {exp.max():.2f}]")

    for method in methods:
        pred = df_clean[method['col']]
        pcc, _ = pearsonr(exp, pred)
        print(f"{method['title']}: PCC = {pcc:.3f}")


# Example usage
if __name__ == "__main__":
    plot_correlation("SourceData_PPI_Benchmark_Dataset.xlsx")