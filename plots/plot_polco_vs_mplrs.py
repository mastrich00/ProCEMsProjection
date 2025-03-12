import csv
import matplotlib.pyplot as plt
from collections import defaultdict

def process_file(filename):
    data = defaultdict(lambda: defaultdict(float))
    iterations = set()
    
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0].startswith('Iter'):
                try:
                    iter_num = int(row[0].split()[1])
                    step = row[1]
                    runtime = float(row[2])
                    data[step][iter_num] = runtime
                    iterations.add(iter_num)
                except (IndexError, ValueError):
                    continue
    
    # Verify all iterations 0-22 exist for each step
    valid_steps = []
    all_iters = set(range(23))
    for step, iters in data.items():
        if set(iters.keys()) == all_iters:
            valid_steps.append(step)
    
    return {step: data[step] for step in valid_steps}, iterations

# Process both files
mplrs_steps, mplrs_iters = process_file('testResults/measurements/mplrs_excel_120t_5ss/excel_times.csv')
polco_steps, polco_iters = process_file('testResults/measurements/polco_excel_120t_5ss/excel_times.csv')

# Find common steps with complete data
common_steps = set(mplrs_steps.keys()) & set(polco_steps.keys())
common_steps = sorted(common_steps, key=lambda x: list(mplrs_steps.keys()).index(x))

# Create plot
n_rows = (len(common_steps) + 1) // 2
fig, axs = plt.subplots(n_rows, 2, figsize=(15, 5*n_rows))
fig.suptitle('Runtime Comparison Between mplrs and polco by Processing Step', y=1.02, fontsize=14)

# Flatten axes array if needed
if len(common_steps) > 1:
    axs = axs.ravel()

# Plot each common step
for idx, step in enumerate(common_steps):
    ax = axs[idx] if len(common_steps) > 1 else axs
    
    # Get data for both tools
    mplrs_vals = [mplrs_steps[step][i] for i in sorted(mplrs_steps[step])]
    polco_vals = [polco_steps[step][i] for i in sorted(polco_steps[step])]
    iterations = sorted(mplrs_steps[step].keys())
    
    # Plot lines
    ax.plot(iterations, mplrs_vals, marker='o', label='mplrs', color='#ff7f0e')
    ax.plot(iterations, polco_vals, marker='s', label='polco', color='#1f77b4')
    
    # Formatting
    ax.set_title(step, fontsize=12)
    ax.set_xlabel('Iteration', fontsize=10)
    ax.set_ylabel('Runtime (seconds)', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Set y-axis to log scale if values vary widely
    max_val = max(max(mplrs_vals), max(polco_vals))
    min_val = min(min(mplrs_vals), min(polco_vals))
    if max_val / min_val > 100:
        ax.set_yscale('log')

# Hide empty subplots
for idx in range(len(common_steps), len(axs)):
    axs[idx].axis('off')

plt.tight_layout()
plt.show()