import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error


def plot_correlation(file_path):
    # Read only the required five columns
    required_cols = ['experimental_DDG', 'mCSM-NA', 'PremPDI', 'SAMPDI-3Dv2', 'MutPNI']
    df = pd.read_excel(file_path, sheet_name="Nabe_DNA", usecols=required_cols)

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

    # Create 2x2 subplot grid
    fig, axs = plt.subplots(2, 2, figsize=(8, 7))
    axs = axs.flatten()

    # Define predictor methods
    methods = [
        {'col': 'mCSM-NA', 'title': 'mCSM-NA'},
        {'col': 'PremPDI', 'title': 'PremPDI'},
        {'col': 'SAMPDI-3Dv2', 'title': 'SAMPDI-3Dv2'},
        {'col': 'MutPNI', 'title': 'MutPNI'}
    ]

    # Experimental values
    exp = df_clean['DDG']

    # Color scheme
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

        # Annotate with statistics (top-left)
        text_content = f'$y = {slope:.2f}x {intercept:+.2f}$\nRMSE = {rmse:.2f}\nPCC = {pcc:.3f}'
        ax.text(
            0.05, 0.95, text_content,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none')
        )

        ax.set_title(f'{method["title"]}', fontsize=13)
        ax.set_xlabel('Experimental ΔΔG (kcal/mol)', fontsize=10)
        ax.set_ylabel('Predicted ΔΔG (kcal/mol)', fontsize=10)

    plt.tight_layout()
    plt.savefig('protein_dna_correlation.tif', bbox_inches='tight', dpi=600,
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.show()

    # Print summary statistics
    print(f"Number of valid data points: {len(df_clean)}")
    print(f"Experimental ΔΔG range: [{exp.min():.2f}, {exp.max():.2f}]")

    for method in methods:
        pred = df_clean[method['col']]
        pcc, _ = pearsonr(exp, pred)
        print(f"{method['title']}: PCC = {pcc:.3f}")


# Example usage – adjust the file path as needed
plot_correlation("SourceData_PDI_Benchmark_Dataset.xlsx")