#!/usr/bin/env python

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
if len(sys.argv) != 3:
	print("Usage: python plot_kmer_freq.py <input_file> <output_pdf>")
	sys.exit(1)
else:
	file_path = sys.argv[1]
	output_pdf = sys.argv[2]

# Read the file without header, skip any trailing text like '>250'
df = pd.read_csv(file_path, sep=r'\s+',skiprows=5, header=None, names=['Kmer Depth','Frequency'], comment='>')
base_name = os.path.basename(file_path).split('.')[0]

# Find the peak
df_tail = df.iloc[0:248] 
max_idx = df_tail['Frequency'].idxmax()
max_point = df.iloc[max_idx]

# Sort by Kmer Depth to ensure plotting order
#df = df.sort_values('Kmer Depth')

# Plot
plt.figure(figsize=(8, 6))
sns.lineplot(data=df, x='Kmer Depth', y='Frequency', marker='o', linewidth=0.1)

# peak point labeling
plt.scatter(max_point['Kmer Depth'], max_point['Frequency'], color='red', zorder=3)
plt.text(
    max_point['Kmer Depth'], max_point['Frequency'],
    f"Max: ({int(max_point['Kmer Depth'])}, {int(max_point['Frequency'])})",
    fontsize=10, color='red', ha='left', va='bottom'
)
# add labels
plt.xlabel('Kmer Depth')
plt.ylabel('Frequency')
plt.title(f'Kmer Frequency Distribution - {base_name}')
plt.tight_layout()

plt.savefig(output_pdf)
print(f"Plot saved to {output_pdf}")
