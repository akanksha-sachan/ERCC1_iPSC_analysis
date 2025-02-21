import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scienceplots
from statannotations import Annotator

def calculate_senescence_score(data, up_genes=None, down_genes=None):
    """
    Calculate senescence signature scores for each sample based on DEGs from gene-sets
    """
    if up_genes is None and down_genes is None:
        raise ValueError("At least one of up_genes or down_genes must be provided")
    # Set index to gene names if not already done
    if "NAME" in data.columns:
        data = data.set_index("NAME")
    if "Description" in data.columns:
        data = data.drop("Description", axis=1)
    # Check which genes are present in the expression data
    available_genes = set(data.index)
    # Process up-regulated genes
    if up_genes is not None:
        up_genes_set = set(up_genes)
        up_genes = list(up_genes_set.intersection(available_genes))
        print(f"Using {len(up_genes)} upregulated genes")
    # Process down-regulated genes
    if down_genes is not None:
        down_genes_set = set(down_genes)
        down_genes = list(down_genes_set.intersection(available_genes))
        print(f"Using {len(down_genes)} downregulated genes")
    # Check if we have enough genes to proceed
    if (up_genes is None or len(up_genes) == 0) and (
        down_genes is None or len(down_genes) == 0
    ):
        raise ValueError(
            "No genes from the gene sets were found in the expression data"
        )
    sample_names = data.columns
    expr_matrix = data.select_dtypes(include=[np.number])
    # Z-standardize the expression values across samples
    z_standardized = (
        expr_matrix - expr_matrix.mean(axis=1).values.reshape(-1, 1)
    ) / expr_matrix.std(axis=1).values.reshape(-1, 1)
    # Calculate scores for each sample
    scores = pd.Series(0, index=sample_names)
    # Calculate total gene set size for normalization
    total_genes = 0
    if up_genes:
        total_genes += len(up_genes)
    if down_genes:
        total_genes += len(down_genes)
    # Calculate combined score with size normalization
    if up_genes:
        up_score = z_standardized.loc[up_genes].sum()
        scores += up_score
    if down_genes:
        down_score = z_standardized.loc[down_genes].sum()
        scores -= down_score  # Subtract because these are down-regulated
    # Normalize by square root of gene set size
    scores = scores / np.sqrt(total_genes)
    # Final z-score normalization across samples
    scores = (scores - scores.mean()) / scores.std()
    # Convert to DataFrame with meaningful column name
    scores_df = pd.DataFrame(scores, columns=["senescence_score"])
    return scores_df


def map_sample_to_category(sample_name):
    """
    Map sample names to their experimental categories
    """
    if sample_name in ["AV01", "AV02", "AV03", "AV04", "AV05"]:
        return "Ctrl Low-Fat diet"
    elif sample_name in ["AV06", "AV07", "AV08", "AV09", "AV10"]:
        return "Trf2 KO Low-Fat diet"
    elif sample_name in ["AV11", "AV12", "AV13", "AV14", "AV15"]:
        return "Ctrl Western diet"
    elif sample_name in ["AV16", "AV17", "AV18", "AV19", "AV20"]:
        return "Trf2 KO Western diet"
    return None


################### PLotting ###################

# Set global style parameters
plt.style.use(['science', 'no-latex'])  # Requires SciencePlots package
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'figure.dpi': 300,
    'figure.figsize': (6, 4)  # More compact size
})
sns.set_theme(context="paper", style="ticks", palette="deep")

def plot_scores_scatter(scores_df, category_colors, title=None, ylabel=None):
    """
    Create a scatter plot of senescence scores with colored labels by category
    """
    plt.figure(figsize=(12, 6))
    # Sort by score for clustering effect
    scores_sorted = scores_df.sort_values("senescence_score")
    # Create scatter plot
    plt.scatter(
        range(len(scores_sorted)), scores_sorted["senescence_score"], alpha=0
    )  # Make points invisible
    # Add colored labels
    for i, (idx, row) in enumerate(scores_sorted.iterrows()):
        plt.annotate(
            idx,  # sample name
            (i, row["senescence_score"]),
            color=category_colors[row["category"]],
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            fontsize=10,
        )
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], color=color, label=cat, marker="o")
        for cat, color in category_colors.items()
    ]
    plt.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1, 0.5))
    # Customize plot
    plt.grid(True, alpha=0.3)
    plt.xlabel("Samples (ordered by score)")
    plt.ylabel(ylabel if ylabel else "Senescence Score")
    plt.title(title if title else "Senescence Scores Distribution")
    plt.xticks([])
    plt.tight_layout()
    plt.show()

def permutation_test(group1, group2, n_permutations=10000):
    """
    Perform a permutation test to compute p-value for the difference between two groups.

    Parameters:
    -----------
    group1, group2 : array-like
        The scores for the two groups being compared.
    n_permutations : int, optional
        Number of permutations to generate the null distribution.

    Returns:
    --------
    observed_diff : float
        The observed mean difference between the groups.
    p_value : float
        The p-value based on the permutation test.
    """
    from scipy.stats import ttest_ind
    observed_diff = np.mean(group1) - np.mean(group2)
    combined = np.concatenate([group1, group2])
    null_diffs = []
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_group1 = combined[: len(group1)]
        perm_group2 = combined[len(group1) :]
        null_diffs.append(np.mean(perm_group1) - np.mean(perm_group2))
    null_diffs = np.array(null_diffs)
    p_value = np.mean(np.abs(null_diffs) >= np.abs(observed_diff))
    return observed_diff, p_value, null_diffs

def plot_scores_simple(scores_df, category_colors, title=None, ylabel=None):
    """
    Plot boxplot with swarmplot overlay in Prism-style formatting
    """
    plt.figure(figsize=(6, 4))  # Controlled by rcParams but explicit here
    # Create plot with updated aesthetics
    ax = sns.boxplot(
        x="category",
        y="senescence_score",
        data=scores_df,
        palette=category_colors,
        width=0.6,
        linewidth=0.7,
        fliersize=0  # Hide default outliers
    )
    # Add individual points with improved styling
    sns.swarmplot(
        x="category", 
        y="senescence_score", 
        data=scores_df,
        color=".2",  # Dark gray
        size=3.5,
        alpha=0.8,
        edgecolor="none"
    )
    # Add horizontal grid lines
    ax.yaxis.grid(True, linestyle='--', alpha=0.4)
    ax.set_axisbelow(True)
    # Customize plot elements
    plt.title(title if title else "Scores by Category", 
             fontsize=11, pad=10)
    plt.ylabel(ylabel if ylabel else "Score", 
              fontsize=10, labelpad=8)
    plt.xlabel("")  # Remove x-axis label
    # Rotate x-ticks and adjust alignment
    plt.xticks(rotation=35, ha='right', rotation_mode='anchor')
    # Clean up borders
    sns.despine(offset=5, trim=True)
    # Adjust layout with tight margins
    plt.tight_layout(pad=1.5)
    plt.show()

# Update your significance plotting function with statannotations
def plot_scores_with_significance(
    scores_df, category_colors, title=None, ylabel=None, n_permutations=1000
):
    """
    Plot scores with Prism-style significance annotations
    """
    plt.figure(figsize=(6, 4))
    # Create base plot
    ax = sns.boxplot(
        x="category",
        y="senescence_score",
        data=scores_df,
        palette=category_colors,
        width=0.6,
        linewidth=0.7
    )
    # Add swarmplot
    sns.swarmplot(
        x="category", 
        y="senescence_score", 
        data=scores_df,
        color=".2",
        size=3.5,
        alpha=0.8
    )
    # Configure annotations
    pairs = [
        ("Ctrl Low-Fat diet", "Trf2 KO Low-Fat diet"),
        ("Ctrl Western diet", "Trf2 KO Western diet"),
        ("Ctrl Low-Fat diet", "Ctrl Western diet"),
        ("Trf2 KO Low-Fat diet", "Trf2 KO Western diet"),
    ]
    # Set up statistical annotations
    annotator = Annotator(
        ax=ax,
        pairs=pairs,
        data=scores_df,
        x="category",
        y="senescence_score",
        order=scores_df['category'].unique()
    )
    # Configure annotation style
    annotator.configure(
        text_format='star', 
        loc='outside',
        line_height=0.02,
        line_offset=0.1,
        text_offset=1.5,
        fontsize=9,
        line_width=0.7
    )
    # Add annotations
    annotator.apply_and_annotate()
    # Final styling
    plt.title(title if title else "Scores by Category", fontsize=11)
    plt.ylabel(ylabel if ylabel else "Score", fontsize=10)
    plt.xlabel("")
    plt.xticks(rotation=35, ha='right')
    sns.despine(offset=5, trim=True)
    plt.tight_layout(pad=1.5)
    plt.show()