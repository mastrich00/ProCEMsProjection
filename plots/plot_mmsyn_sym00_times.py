import csv
from matplotlib import pyplot as plt
from collections import defaultdict

# Read and filter data
data = []
with open('testResults/measurements/polco_mmsyn_sm00_120t_2ss/mmsyn00_times.csv', 'r') as f:
# with open('testResults/measurements/mplrs_ecoli_120t_2ss/ecoli_csv.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        phase = row[0]
        if phase not in ('Start', 'End'):
            data.append(row)

# Define hierarchical legend structure
legend_groups = [
    ("1. Initial\ Setup", ["Initial Setup"]),
    ("2. Projection:", [
        "Setup Elimination-Matrix",
        "Redund Projection-Cone",
        "Enumeration Projection-Cone",
        "Redund Rays",
        "Build projected Cone",
        "Redund projected Cone"
    ]),
    # ("3. Enum.\ ProCEMs:", [
    #     "Enumerating ProCEMs",
    #     "Redund proCEMs",
    #     "Postprocessing proCEMs"
    # ])
]

# # Create color mapping
# colors = plt.cm.tab20.colors
# color_map = {}
# color_idx = 0
# for _, items in legend_groups:
#     for item in items:
#         if item not in color_map:
#             color_map[item] = colors[color_idx % len(colors)]
#             color_idx += 1
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

# Create legend handles and labels with styling
legend_handles = []
legend_labels = []
for group_name, items in legend_groups:
    # Add bold group header with larger font
    legend_handles.append(plt.Rectangle((0,0), 1, 1, fc='none', ec='none', alpha=0))
    legend_labels.append(f"$\mathbf{{{group_name}}}$")  # LaTeX bold
    
    # Add sub-items with normal font
    for item in items:
        legend_handles.append(plt.Rectangle((0,0), 1, 1, color=color_map[item]))
        legend_labels.append(f"  {item}")

# Group data by phase and preserve order
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

# Prepare data for plotting
durations_list = []
labels_list = []
for phase in phases_order:
    actions = groups[phase]
    durations = [d for (a, d) in actions]
    labels = [a for (a, d) in actions]
    durations_list.append(durations)
    labels_list.append(labels)

# Plot stacked bars with logarithmic y-axis
fig, ax = plt.subplots(figsize=(12, 7))
x = range(len(phases_order))

# Add main figure title with bold formatting
# fig.suptitle(r'$\mathbf{Duration\ of\ Processing\ Phases\ -\ e\ coli\ core\ (onto\ exchange\ reactions)}$', 
#              fontsize=14, y=.95, weight='bold')

fig.suptitle(
    r'$\mathbf{Duration\ of\ Processing\ Phases\ -\ mmsyn\_sm00\ onto\ exchange\ reactions}$',
    fontsize=16,
    y=0.96,  # Original y=0.98 -> 0.96
    verticalalignment='top'
)

# Add second subtitle
fig.text(
    0.5,  # x-position (centered)
    0.93,  # y-position (below main title)
    r'$\mathbf{polco\ (threads=120,\ step\ size=2)}$',
    ha='center',
    va='top',
    fontsize=12,
    alpha=0.8
)

for i, phase in enumerate(phases_order):
    durations = durations_list[i]
    labels = labels_list[i]
    bottom = 0
    for dur, label in zip(durations, labels):
        ax.bar(i, dur, bottom=bottom, color=color_map[label], edgecolor='black')
        bottom += dur

# Configure axes and labels
ax.set_xticks(x)
ax.set_xticklabels(phases_order, rotation=45, ha='right')
ax.set_ylabel('Duration (seconds)', fontsize=10)
ax.set_xlabel('Processing Steps', fontsize=10)

# Remove the axes title since we now have figure title
ax.set_title('')

legend = ax.legend(
    legend_handles,
    legend_labels,
    title=r'$\mathbf{Processing\ Stages}$',
    title_fontsize=12,
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    frameon=True,
    edgecolor='black',
    facecolor='white',  # Keep background white
    framealpha=0.95,
    borderpad=1.2,
    labelspacing=1.5,
    handlelength=2.5,
    fontsize=10,
    handletextpad=0.8
)

# Optional: Further refine the frame appearance
legend.get_frame().set_linewidth(0.7)  # Match value to linewidth above
# legend.get_frame().set_boxstyle("Round", pad=0.2, rounding_size=0.3)

# Add finer grid lines
ax.grid(True, which='both', axis='y', linestyle='--', linewidth=0.5, alpha=0.7)

# Adjust layout with better spacing for title and legend
plt.tight_layout(rect=[0, 0, 0.85, 0.98])  # Adjusted for title space

plt.show()
