#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys

# Read data from the file
x_vals = []
y_vals = []

with open(sys.argv[1], "r") as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) == 2:
            x, y = map(int, parts)
            x_vals.append(x)
            y_vals.append(y)

# Find the first index where Y stops increasing
stop_index = None
for i in range(1, len(y_vals)):
    if y_vals[i] == y_vals[i - 1]:
        stop_index = i
        break

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(x_vals, y_vals, linestyle='-', marker='o', markersize=3)

# Highlight the point where Y stops increasing
if stop_index is not None:
    plt.plot(x_vals[stop_index], y_vals[stop_index], 'ro')  # red dot
    plt.annotate(f"({x_vals[stop_index]}, {y_vals[stop_index]})",
                 (x_vals[stop_index], y_vals[stop_index]),
                 textcoords="offset points", xytext=(10,-15), ha='left', fontsize=10)

name = sys.argv[1].split(".")[0] + ".distinct_kmer"
# Customize plot
plt.xlabel("reads number",fontsize=12)
plt.ylabel("number of distinct kmers",fontsize=12)
plt.title(sys.argv[1].split(".")[0])
plt.xticks([])  # remove x-axis ticks
plt.grid(False)

# Show the plot
plt.tight_layout()
plt.savefig(name+".pdf")


