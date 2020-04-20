"""
Given a file report TSV from encodeproject.org containing cloud_metadata, flatten the
cloud metadata JSON structure and return a new report with the flattened data.
"""

import json
import pandas as pd

df = pd.read_csv("file_report_2020_4_15_2h_50m.tsv", sep="\t", header=1)
# Report data is not valid JSON, need to replace single quotes. Works in this case
df["Cloud metadata"] = df["Cloud metadata"].str.replace("'", '"')
flattened_cloud = pd.json_normalize(df["Cloud metadata"].apply(json.loads))
merged = pd.concat([df, flattened_cloud], axis=1)
merged.drop("Cloud metadata", axis=1, inplace=True)
merged.to_csv("file_report_2020_4_15_2h_50m_expanded_cloud.tsv", sep="\t", index=False)
