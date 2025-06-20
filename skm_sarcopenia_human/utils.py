import pandas as pd
from collections import defaultdict
import re

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
