import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations import Annotator

def calculate_senescence_score(data, up_genes=None, down_genes=None):
    """
    Calculate senescence signature scores for each sample based on gene-sets with directionality
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


def map_sample_to_category_trf2(sample_name):
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

def map_sample_to_category_ercc1(sample_name):
    """
    Map sample names to their experimental categories
    """
    if sample_name in ["AV01", "AV02", "AV03", "AV04", "AV05", "AV06", "AV07"]:
        return "WT"
    elif sample_name in ["AV08", "AV09", "AV10", "AV11", "AV12", "AV13", "AV14"]:
        return "Ercc1 LKO"
    return None

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

def calculate_pairwise_significance(data, groups):
    """
    Calculate pairwise significance between all groups
    Returns a dictionary of p-values and significance levels
    """
    from scipy import stats
    results = {}
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            group1 = data[data['category'] == groups[i]]['senescence_score']
            group2 = data[data['category'] == groups[j]]['senescence_score']
            
            # Perform Mann-Whitney U test
            statistic, pvalue = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            
            # Add significance stars
            if pvalue < 0.001:
                sig = '***'
            elif pvalue < 0.01:
                sig = '**'
            elif pvalue < 0.05:
                sig = '*'
            else:
                sig = 'ns'
                
            results[(i, j)] = {'p-value': pvalue, 'significance': sig}
    
    return results

def create_overlap_data(genesets):
    names = list(genesets.keys())
    percent_matrix = pd.DataFrame(0.0, index=names, columns=names)
    count_matrix = pd.DataFrame(0, index=names, columns=names)
    
    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            # Skip diagonal and lower triangle
            if i >= j:
                continue
                
            set1 = set(genesets[name1])
            set2 = set(genesets[name2])
            
            intersection = len(set1.intersection(set2))
            min_size = min(len(set1), len(set2))
            
            # Store both percentage and count
            percent_matrix.loc[name1, name2] = (intersection / min_size * 100) if min_size > 0 else 0
            count_matrix.loc[name1, name2] = intersection
    
    return percent_matrix, count_matrix

################### Plotting ###################

# Set global style parameters
plt.style.use(['ggplot'])  # Requires SciencePlots package
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

def plot_violin_box_combo(data, x_var, y_var, title=None, x_ticks=None, palette=None, rotation=45):
    """
    Create a combined violin-box plot with consistent colors for all elements
    """
    plt.clf()
    # Reduce figure width and adjust spacing
    fig, ax = plt.subplots(figsize=(5, 6))
    
    # Adjust plot margins
    plt.subplots_adjust(left=0.15, right=0.85, bottom=0.1, top=0.9)

    # Calculate y-axis limits based on data
    y_min = data[y_var].min()
    y_max = data[y_var].max()
    
    # Add padding and round to nearest 0.5
    y_min = np.floor(y_min * 2) / 2
    y_max = np.ceil(y_max * 2) / 2
    
    # Set y-axis limits and ticks
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))  # Set tick intervals to 0.5

    # Create violin plot in the background
    violin = sns.violinplot(
        data=data, x=x_var, y=y_var,
        palette=palette, inner=None,
        linewidth=0, saturation=1.0,
        alpha=0.3, width=0.4, cut=0
    )

    # Get the unique categories in the order they appear
    categories = data[x_var].unique()

    # Create box plot with correct colors from the start
    box_plot = sns.boxplot(
        data=data, x=x_var, y=y_var,
        width=0.4, linewidth=1.2,
        flierprops={'marker': ' '},
        showmeans=False,
        boxprops={
            'facecolor': 'none',
            'edgecolor': 'none'
        },
        whiskerprops={'color': 'none'},
        medianprops={'color': 'none'},
        showcaps=False,
        ax=ax
    )

    # Count number of boxes and lines per box
    num_boxes = len(categories)
    lines_per_box = len(ax.lines) // num_boxes

    # Update box plot colors after creation
    for i, (name, box) in enumerate(zip(categories, ax.patches)):
        color = palette[name]
        
        # Create filled box with transparency
        box.set_facecolor(color)
        box.set_edgecolor('none')
        box.set_alpha(0.3)
        box.set_zorder(1)
        
        # Create box edges with full opacity
        import matplotlib.patches as mpatches
        path = box.get_path()
        edges = mpatches.PathPatch(
            path,
            facecolor='none',
            edgecolor=color,
            linewidth=1.2,
            alpha=1.0,
            zorder=2
        )
        ax.add_patch(edges)
        
        # Get and color all lines for this box
        box_lines = ax.lines[i * lines_per_box : (i + 1) * lines_per_box]
        for line in box_lines:
            line.set_color(color)
            line.set_alpha(1.0)
            line.set_linewidth(1.2)
            line.set_zorder(2)

    # Add individual points on top
    sns.stripplot(
        data=data, x=x_var, y=y_var,
        palette=palette, size=6,
        alpha=1.0, linewidth=0,
        jitter=0.2, zorder=3
    )
    
    # Calculate significance
    categories = data[x_var].unique()
    significance_info = calculate_pairwise_significance(data, categories)

    # Add significance bar
    def add_significance_bar(start, end, height, p_value, sig_symbol):
        bar_height = height
        bar_tips = 0.05
        
        # Draw the bar
        ax.plot([start, start, end, end], 
                [bar_height, bar_height + bar_tips, bar_height + bar_tips, bar_height],
                color='black', linewidth=0.8)
        
        # Add text
        text = f'p = {p_value:.4f} {sig_symbol}'
        ax.text((start + end) * 0.5, bar_height + bar_tips, 
                text, ha='center', va='bottom', fontsize=8)

    # Get current y limits
    current_ymin, current_ymax = ax.get_ylim()
    bar_height = current_ymax + 0.15

    # Add significant bars (p < 0.05 only)
    for (group1_idx, group2_idx), sig_data in significance_info.items():
        if sig_data['significance'] != 'ns':  # Only show significant comparisons
            add_significance_bar(
                group1_idx, 
                group2_idx, 
                bar_height,
                sig_data['p-value'],
                sig_data['significance']
            )
            bar_height += 0.15  # Increment height for next bar

    # Adjust y-axis limits to accommodate bars
    ax.set_ylim(current_ymin, bar_height + 0.1)

    if title:
        plt.title(title, pad=20)

    if x_ticks is None:
        ax.set_xticks([])
        ax.spines['bottom'].set_visible(False)
    else:
        ax.set_xticks(range(len(x_ticks)))
        ax.set_xticklabels(x_ticks, rotation=rotation, ha='right')
        plt.setp(ax.get_xticklabels(), rotation=rotation, ha='right')  # Add this line
        ax.spines['bottom'].set_visible(True)

    # Configure ticks and spines with thinner lines
    ax.minorticks_off()
    ax.tick_params(axis='both', which='minor', bottom=False, top=False, left=False, right=False)
    ax.tick_params(axis='x', which='major', top=False)
    ax.tick_params(axis='y', which='major', right=False, width=0.8)
    
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_tick_params(width=0.8)
    
    plt.setp(ax.get_yticklabels(), weight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.yaxis.grid(False)
    
    sns.despine(offset=5, trim=True, bottom=(x_ticks is None), right=True)
    
    # Force rotation of x-tick labels
    if x_ticks is not None:
        plt.setp(ax.get_xticklabels(), rotation=rotation, ha='right')
    
    plt.close()
    
    return fig

def plot_overlap_heatmap(percent_matrix, title, ax, is_up=True):
    # Create masks for lower triangle and diagonal
    mask_lower = np.tril(np.ones(percent_matrix.shape), k=-1).astype(bool)
    mask_diagonal = np.eye(percent_matrix.shape[0], dtype=bool)
    
    # First plot the diagonal with grey color
    sns.heatmap(percent_matrix,
                mask=~mask_diagonal,  # Only show diagonal
                cmap=['lightgrey'],   # Use grey color for diagonal
                cbar=False,
                ax=ax)
    
    # Create custom colormaps from white to green/blue
    from matplotlib.colors import LinearSegmentedColormap
    if is_up:
        colors = ['white', '#00441b']  # White to dark green
        cmap = LinearSegmentedColormap.from_list('custom_green', colors)
    else:
        colors = ['white', '#08519c']  # White to dark blue
        cmap = LinearSegmentedColormap.from_list('custom_blue', colors)
    
    # Then plot the upper triangle with data
    sns.heatmap(percent_matrix, 
                mask=mask_lower | mask_diagonal,  # Hide lower triangle and diagonal
                annot=True,  # Show values in cells
                fmt='.1f',   # Format as float with 1 decimal
                cmap=cmap,
                vmin=0,      # Set minimum value to 0
                vmax=50,     # Set maximum value to 50
                cbar_kws={'label': 'Overlap (%)', 'ticks': [0, 10, 20, 30, 40, 50]},
                ax=ax)
    
    ax.set_title(title)
    # Rotate x-axis labels for better readability
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')