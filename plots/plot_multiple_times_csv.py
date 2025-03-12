import csv
import matplotlib.pyplot as plt
from collections import defaultdict

# Define CSV files and subplot titles
csv_files = [
    'testResults/measurements/polco_excel_120t_10ss/excel_times.csv',
    'testResults/measurements/polco_excel_120t_30ss/excel_times.csv',
    'testResults/measurements/polco_excel_120t_60ss/excel_times.csv',
    'testResults/measurements/polco_excel_120t_ALLss/excel_times.csv'
]
titles = ['Step Size = 10', 'Step Size = 30', 'Step Size = 60', 'Step Size = ALL']

# Define COMPLETE color reference list (including all possible steps)
full_color_reference = [
    "Initial Setup",
    "Setup Elimination-Matrix",
    "Redund Projection-Cone",
    "Enumeration Projection-Cone",
    "Redund Rays",
    "Build projected Cone",
    "Redund projected Cone",
    "Enumerating ProCEMs",
    "Redund proCEMs",
    "Postprocessing proCEMs"
]
colorblind_palette = [
    "#4C72B0",   # Blue
    "#DD8452",   # Orange
    "#9970AB",   # Muted Purple (replacing green)
    "#8DA0CB",   # Light Blue (replacing red)
    "#8172B3",   # Dark Purple
    "#CCB974",   # Mustard
    "#64B5CD",   # Cyan-Blue
    "#DA8BC3",   # Pink
    "#7C7C7C",   # Gray
    "#F0E442"    # Yellow
]
color_map = {label: colorblind_palette[i] for i, label in enumerate(full_color_reference)}


# Define VISIBLE legend structure (excluding commented steps)
legend_groups = [
    ("1. Initial\ Setup:", ["Initial Setup"]),
    ("2. Projection:", [
        "Setup Elimination-Matrix",
        "Enumeration Projection-Cone",
        "Build projected Cone",
        "Redund projected Cone"
    ]),
    ("3. Enum.\ ProCEMs:", [
        "Enumerating ProCEMs"
    ])
]

# Create legend handles using the full color map
legend_handles = []
legend_labels = []
for group_name, items in legend_groups:
    # Add group header
    legend_handles.append(plt.Rectangle((0,0), 1, 1, fc='none', ec='none', alpha=0))
    legend_labels.append(f"$\mathbf{{{group_name}}}$")
    
    # Add items with original colors
    for item in items:
        legend_handles.append(plt.Rectangle((0,0), 1, 1, color=color_map[item]))
        legend_labels.append(f"  {item}")

# Create subplot grid with figure title
fig, axs = plt.subplots(2, 2, figsize=(18, 12))
# fig.suptitle(r'$\mathbf{polco:\ Duration\ of\ Processing\ Phases\ -\ A.\ thaliana\ plastid\ onto\ sugar\ and\ starch\ metabolism\ (120\ threads,\ different\ step\ sizes)}$', 
#              fontsize=16, y=0.98)
             
fig.suptitle(
    r'$\mathbf{Duration\ of\ Processing\ Phases\ -\ A.\ thaliana\ plastid\ onto\ sugar\ and\ starch\ metabolism}$',
    fontsize=16,
    y=0.96,  # Original y=0.98 -> 0.96
    verticalalignment='top'
)

# Add second subtitle
fig.text(
    0.5,  # x-position (centered)
    0.93,  # y-position (below main title)
    r'$\mathbf{polco\ (threads=120,\ step\ sizes=[10,\ 30,\ 60,\ All])}$',
    ha='center',
    va='top',
    fontsize=12,
    alpha=0.8
)

axs = axs.ravel()

# Process each CSV file and plot
for idx, (csv_file, title) in enumerate(zip(csv_files, titles)):
    # Data reading and processing (same as before)
    data = []
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            phase = row[0]
            if phase not in ('Start', 'End'):
                data.append(row)
    
    phases_order = []
    seen_phases = set()
    groups = defaultdict(list)
    for row in data:
        phase, action, duration_str = row[0], row[1], row[2]
        duration = float(duration_str)
        if phase not in seen_phases:
            phases_order.append(phase)
            seen_phases.add(phase)
        groups[phase].append((action, duration))
    
    # Prepare plot data
    durations_list = []
    labels_list = []
    for phase in phases_order:
        actions = groups[phase]
        durations = [d for (a, d) in actions]
        labels = [a for (a, d) in actions]
        durations_list.append(durations)
        labels_list.append(labels)
    
    # Plotting with reference styling
    ax = axs[idx]
    for i, phase in enumerate(phases_order):
        durations = durations_list[i]
        labels = labels_list[i]
        bottom = 0
        for dur, label in zip(durations, labels):
            ax.bar(i, dur, bottom=bottom, color=color_map[label], edgecolor='black')
            bottom += dur
    
    # Subplot configuration
    ax.set_xticks(range(len(phases_order)))
    ax.set_xticklabels(phases_order, rotation=45, ha='right', fontsize=8)
    ax.set_title(title, fontsize=10)
    # ax.set_yscale('log')
    ax.set_ylim(bottom=0.01, top=55)
    ax.grid(True, axis='y', linestyle='--', alpha=0.3)
    
    # Axis label handling
    if idx in [0, 2]:
        ax.set_ylabel('Duration (seconds)', fontsize=9)

# Create unified legend with reference styling
legend = fig.legend(
    legend_handles,
    legend_labels,
    title=r'$\mathbf{Processing\ Stages}$',
    title_fontsize=11,
    bbox_to_anchor=(0.85, 0.65),
    loc='center left',
    frameon=True,
    edgecolor='black',
    facecolor='white',
    framealpha=0.95,
    borderpad=1.2,
    labelspacing=1.2,
    handlelength=2,
    fontsize=9,
    handletextpad=0.6
)

# Adjust legend frame
legend.get_frame().set_linewidth(0.7)
# legend.get_frame().set_boxstyle("Round", pad=0.2, rounding_size=0.3)

# Final layout adjustments
plt.subplots_adjust(
    left=0.08,
    right=0.82,
    bottom=0.12,
    top=0.88,
    hspace=0.35,
    wspace=0.25
)

plt.show()