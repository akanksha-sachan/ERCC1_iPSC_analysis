import pandas as pd
from collections import defaultdict
from scipy.stats import ttest_ind
from itertools import combinations
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def read_series_matrix(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    return lines

def parse_sample_metadata(lines):
    """Parses !Sample_* metadata into a dict."""
    metadata = defaultdict(list)
    for line in lines:
        if line.startswith("!Sample_"):
            parts = line.strip().split("\t")
            key = parts[0][8:]  # remove "!Sample_"
            values = parts[1:]
            metadata[key].append(values)
        elif line.startswith("!series_matrix_table_begin"):
            break
    # Merge duplicates
    merged_metadata = {
        key: list(map(lambda x: " || ".join(x), zip(*val_lists)))
        for key, val_lists in metadata.items()
    }
    return merged_metadata

def parse_characteristics(characteristics_list):
    """Parses the 'characteristics_ch1' strings into structured columns."""
    def extract_fields(entry):
        parts = entry.split(" || ")
        out = {}
        for p in parts:
            if ": " in p:
                k, v = p.split(": ", 1)
                out[k.strip().strip('"').lower()] = v.strip().strip('"')
        return out
    return [extract_fields(c) for c in characteristics_list]

def clean_metadata_df(metadata_dict):
    """Turns parsed metadata into a cleaned DataFrame with labels."""
    sample_ids = metadata_dict.get("title", [f"Sample_{i+1}" for i in range(len(next(iter(metadata_dict.values()))))])
    characteristics = metadata_dict["characteristics_ch1"]
    
    parsed = parse_characteristics(characteristics)
    df = pd.DataFrame(parsed)
    df.columns = [c.strip().strip('"').lower() for c in df.columns]
    df["sample_id"] = sample_ids
    df["sex"] = df["sex"].str.strip().str.lower().str.capitalize()
    df["disease_state"] = df["disease state"].str.strip().str.lower()
    
    # extract label from sample_id (e.g., "HA 01" → HA_1)
    def extract_label(sample):
        match = re.search(r'([A-Z]+)\s*0*(\d+)', sample)
        if match:
            return f"{match.group(1)}_{int(match.group(2))}"
        return "UNKNOWN"
    
    df["label"] = df["sample_id"].str.strip().str.strip('"').apply(extract_label)
    return df[["sample_id", "sex", "disease_state", "label"]]

def split_expression_by_sex(expression_df, metadata_df):
    """Splits expression matrix into male and female using label column."""
    expression_df = expression_df.set_index("Name")
    metadata_df["label"] = metadata_df["label"].str.strip()
    metadata_df["sex"] = metadata_df["sex"].str.strip().str.capitalize()

    female_labels = metadata_df.loc[metadata_df["sex"] == "Female", "label"]
    male_labels = metadata_df.loc[metadata_df["sex"] == "Male", "label"]

    female_labels = [l for l in female_labels if l in expression_df.columns]
    male_labels = [l for l in male_labels if l in expression_df.columns]

    female_expr = expression_df[female_labels].reset_index()
    male_expr = expression_df[male_labels].reset_index()
    
    return female_expr, male_expr

def load_gmt(filepath):
    """
    Load a GMT file into a dictionary of {pathway_name: list_of_genes}
    """
    gene_sets = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathway = parts[0]
                genes = parts[2:]  # Skip the description (parts[1])
                gene_sets[pathway] = genes
    return gene_sets

def calculate_geneset_score(data, up_genes=None, down_genes=None, mapping_df=None):
    """
    Calculate geneset signature scores for each sample based on gene sets.
    Optionally attaches disease state if mapping_df is provided.
    """
    import numpy as np
    import pandas as pd

    if up_genes is None and down_genes is None:
        raise ValueError("At least one of up_genes or down_genes must be provided")

    # Set index to gene names if not already done
    if "Name" in data.columns:
        data = data.set_index("Name")
    if "Description" in data.columns:
        data = data.drop("Description", axis=1)

    # Check which genes are present
    available_genes = set(data.index)

    if up_genes is not None:
        up_genes = list(set(up_genes).intersection(available_genes))
        print(f"Using {len(up_genes)} upregulated genes")

    if down_genes is not None:
        down_genes = list(set(down_genes).intersection(available_genes))
        print(f"Using {len(down_genes)} downregulated genes")

    if (not up_genes) and (not down_genes):
        raise ValueError("No genes from gene sets found in expression data")

    sample_names = data.columns
    expr_matrix = data.select_dtypes(include=[np.number])

    # Z-standardize expression across samples (per gene)
    z = (expr_matrix - expr_matrix.mean(axis=1).values.reshape(-1, 1)) / expr_matrix.std(axis=1).values.reshape(-1, 1)

    # Score calculation
    scores = pd.Series(0, index=sample_names, dtype=float)
    total_genes = (len(up_genes) if up_genes else 0) + (len(down_genes) if down_genes else 0)

    if up_genes:
        scores += z.loc[up_genes].sum()
    if down_genes:
        scores -= z.loc[down_genes].sum()

    scores /= np.sqrt(total_genes)  # Normalize by gene set size
    scores = (scores - scores.mean()) / scores.std()  # Z-score across samples

    scores_df = pd.DataFrame({"geneset_score": scores})

    # If mapping is provided, merge in disease state
    if mapping_df is not None:
        scores_df = scores_df.reset_index().rename(columns={"index": "label"})
        mapping_df = mapping_df[["label", "disease_state"]].drop_duplicates()
        scores_df = scores_df.merge(mapping_df, on="label", how="left").set_index("label")

    return scores_df

def plot_geneset_scores_scatter(scores_df, category_col="disease_state", palette=None, title=None, ylabel=None):
    """
    Create a scatter plot of gene set scores with colored sample labels by category.
    
    Parameters:
    - scores_df: DataFrame with index = sample label, and columns: 'geneset_score' and [category_col]
    - category_col: str, column to color by (e.g., 'disease_state')
    - palette: dict, mapping of category → color (REQUIRED)
    - title: str, optional plot title
    - ylabel: str, optional y-axis label
    """
    import matplotlib.pyplot as plt

    assert palette is not None, "You must provide a palette dictionary to ensure consistent coloring."

    plt.figure(figsize=(12, 6))
    
    # Sort by gene set score for visual ordering
    scores_sorted = scores_df.sort_values("geneset_score")

    # Plot invisible points to fix axes
    plt.scatter(
        range(len(scores_sorted)), 
        scores_sorted["geneset_score"], 
        alpha=0  # invisible markers
    )

    # Annotate each sample
    for i, (sample, row) in enumerate(scores_sorted.iterrows()):
        category = row[category_col]
        plt.annotate(
            sample,
            (i, row["geneset_score"]),
            color=palette[category],
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            fontsize=10,
        )

    # Custom legend
    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="none", label=cat, markerfacecolor=color, markersize=8)
        for cat, color in palette.items()
    ]
    plt.legend(handles=legend_elements, title=category_col, loc="center left", bbox_to_anchor=(1, 0.5))

    # Formatting
    plt.grid(True, alpha=0.3)
    plt.xlabel("Samples (ordered by score)")
    plt.ylabel(ylabel if ylabel else "Gene Set Score", weight="bold")
    plt.title(title if title else "Gene Set Scores per Sample", pad=10)
    plt.xticks([])
    plt.tight_layout()
    plt.show()

def plot_geneset_score_violin_box(scores_df, category_col="disease_state", score_col="geneset_score", palette=None, title=None):
    """
    Plot geneset scores across categories using combined violin-box-strip plot.
    
    Parameters:
    - scores_df: DataFrame with columns [category_col, score_col]
    - category_col: column with categorical grouping (e.g. disease_state)
    - score_col: column with numeric gene set score
    - palette: dict of colors for each category (REQUIRED)
    - title: optional plot title
    """

    assert palette is not None, "You must pass a color palette dictionary for consistent category coloring."

    def calculate_pairwise_significance(data, x_var, y_var):
        results = {}
        groups = data[x_var].unique()
        group_indices = {group: i for i, group in enumerate(groups)}
        for g1, g2 in combinations(groups, 2):
            x1 = data[data[x_var] == g1][y_var]
            x2 = data[data[x_var] == g2][y_var]
            stat, p = ttest_ind(x1, x2, equal_var=False)
            if p < 0.001: sig = '***'
            elif p < 0.01: sig = '**'
            elif p < 0.05: sig = '*'
            else: sig = 'ns'
            results[(group_indices[g1], group_indices[g2])] = {"p-value": p, "significance": sig}
        return results

    plt.clf()
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(left=0.15, right=0.85, bottom=0.1, top=0.9)

    # Violin plot
    sns.violinplot(
        data=scores_df, x=category_col, y=score_col,
        palette=palette, inner=None, linewidth=0, alpha=0.3, width=0.5, cut=0, ax=ax
    )

    # Box plot
    sns.boxplot(
        data=scores_df, x=category_col, y=score_col,
        palette=palette, fliersize=0, linewidth=1.2, width=0.3, ax=ax
    )

    # Strip plot
    sns.stripplot(
        data=scores_df, x=category_col, y=score_col,
        palette=palette, size=6, jitter=0.25, linewidth=0, ax=ax, zorder=3
    )

    # Significance bars
    significance = calculate_pairwise_significance(scores_df, category_col, score_col)
    y_max = scores_df[score_col].max()
    bar_height = y_max + 0.5
    bar_gap = 0.3

    for (i, j), result in significance.items():
        if result["significance"] != "ns":
            ax.plot([i, i, j, j],
                    [bar_height, bar_height + 0.05, bar_height + 0.05, bar_height],
                    color='black', linewidth=1.0)
            ax.text((i + j) / 2, bar_height + 0.07,
                    f"{result['significance']} (p={result['p-value']:.3g})",
                    ha='center', va='bottom', fontsize=9)
            bar_height += bar_gap

    ax.set_ylim(scores_df[score_col].min() - 0.5, bar_height + 0.1)
    ax.set_xlabel("")
    ax.set_ylabel("Gene Set Score", weight='bold')
    ax.set_title(title or "Gene Set Score by Category", pad=15)

    sns.despine()
    plt.tight_layout()
    plt.show()
