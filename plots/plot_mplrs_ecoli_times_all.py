import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
import numpy as np
from collections import defaultdict

# Define CSV files and subplot titles
csv_files = [
    'testResults/measurements/mplrs_ecoli_120t_2ss/ecoli_csv.csv',
    'testResults/measurements/polco_ecoli_120t_2ss/ecoli_csv.csv'
]
titles = ['mplrs: Threads = 120, Step Size = 2', 'polco: Threads = 120, Step Size = 2']

# Define COMPLETE color reference list and create a color mapping using a colorblind-friendly palette
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

# Define hierarchical legend groups
legend_groups = [
    ("1. Initial\\ Setup", ["Initial Setup"]),
    ("2. Projection:", [
        "Setup Elimination-Matrix",
        "Redund Projection-Cone",
        "Enumeration Projection-Cone",
        "Redund Rays",
        "Build projected Cone",
        "Redund projected Cone"
    ]),
    ("3. Enum.\\ ProCEMs:", [
        "Enumerating ProCEMs",
        "Redund proCEMs",
        "Postprocessing proCEMs"
    ])
]
legend_handles = []
legend_labels = []
for group_name, items in legend_groups:
    # Group header (invisible patch)
    legend_handles.append(plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', alpha=0))
    legend_labels.append(f"$\\mathbf{{{group_name}}}$")
    for item in items:
        legend_handles.append(plt.Rectangle((0, 0), 1, 1, color=color_map[item]))
        legend_labels.append(f"  {item}")

# Create a figure with two rows of subplots per CSV file.
# The top row will show high values (y from 250 to 300) and the bottom row low values (y from 0 to 40).
nfiles = len(csv_files)
fig, axs = plt.subplots(2, nfiles, figsize=(12, 10), sharex=True,
                        gridspec_kw={'height_ratios': [1, 2]})
fig.subplots_adjust(hspace=0.05)  # Reduce vertical space between subplots

# Process each CSV file and plot the stacked bars on both axes
for idx, (csv_file, title) in enumerate(zip(csv_files, titles)):
    # Read CSV data (ignoring 'Start' and 'End' markers)
    data = []
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] not in ('Start', 'End'):
                data.append(row)
    
    # Group data by phase so each phase corresponds to one stacked bar
    phases_order = []
    groups = defaultdict(list)
    for row in data:
        phase, action, duration_str = row[0], row[1], row[2]
        duration = float(duration_str)
        if phase not in phases_order:
            phases_order.append(phase)
        groups[phase].append((action, duration))
    
    # Prepare durations and labels for each phase
    durations_list = []
    labels_list = []
    for phase in phases_order:
        actions = groups[phase]
        durations = [d for (a, d) in actions]
        labels = [a for (a, d) in actions]
        durations_list.append(durations)
        labels_list.append(labels)
    
    # Select the two axes for this CSV file
    ax_top = axs[0, idx]   # Upper subplot: high values (250 to 300)
    ax_bot = axs[1, idx]   # Lower subplot: low values (0 to 40)
    
    # Plot each stacked bar on both axes.
    # The entire bar is drawn on each, and the axes clip what’s outside.
    for i, phase in enumerate(phases_order):
        cumulative = 0
        for dur, label in zip(durations_list[i], labels_list[i]):
            ax_top.bar(i, dur, bottom=cumulative, color=color_map[label],
                       edgecolor='black', linewidth=0.5)
            ax_bot.bar(i, dur, bottom=cumulative, color=color_map[label],
                       edgecolor='black', linewidth=0.5)
            cumulative += dur

    # Set y-limits for broken axes
    ax_bot.set_ylim(0, 43)
    ax_top.set_ylim(275, 300)
    
    # Format x-axis and titles
    ax_bot.set_xticks(range(len(phases_order)))
    ax_bot.set_xticklabels(phases_order, rotation=45, ha='right', fontsize=8)
    ax_top.set_title(title, fontsize=10)
    ax_bot.set_ylabel('Duration (seconds)', fontsize=10)
    
    # Add grid lines
    ax_top.grid(True, axis='y', linestyle='--', alpha=0.3)
    ax_bot.grid(True, axis='y', linestyle='--', alpha=0.3)
    
    # Hide the spines between ax_top and ax_bot for a cleaner look
    ax_top.spines['bottom'].set_visible(False)
    ax_bot.spines['top'].set_visible(False)
    ax_top.tick_params(labelbottom=False)  # Remove x-axis tick labels on top
    
    # Define a fixed break marker length in pixels (the full diagonal line length)
    fixed_length_pixels = 10

    # --- Top Axis (break on the bottom border) ---
    # Get the bounding box (in display coordinates) for conversion.
    bbox_top = ax_top.get_window_extent()
    # Compute the fixed pixel lengths converted to fractions of the axis width and height.
    dx_top = (fixed_length_pixels / np.sqrt(2)) / bbox_top.width
    dy_top = (fixed_length_pixels / np.sqrt(2)) / bbox_top.height

    # Draw 45° diagonal markers at the left and right on the bottom edge.
    # Left marker: centered at (0,0)
    ax_top.plot([0 - dx_top, 0 + dx_top], [0 - dy_top, 0 + dy_top],
                transform=ax_top.transAxes, color='black', clip_on=False)
    # Right marker: centered at (1,0)
    ax_top.plot([1 - dx_top, 1 + dx_top], [0 - dy_top, 0 + dy_top],
                transform=ax_top.transAxes, color='black', clip_on=False)

    # --- Bottom Axis (break on the top border) ---
    bbox_bot = ax_bot.get_window_extent()
    dx_bot = (fixed_length_pixels / np.sqrt(2)) / bbox_bot.width
    dy_bot = (fixed_length_pixels / np.sqrt(2)) / bbox_bot.height

    # Draw 45° diagonal markers at the left and right on the top edge.
    # Left marker: centered at (0,1)
    ax_bot.plot([0 - dx_bot, 0 + dx_bot], [1 - dy_bot, 1 + dy_bot],
                transform=ax_bot.transAxes, color='black', clip_on=False)
    # Right marker: centered at (1,1)
    ax_bot.plot([1 - dx_bot, 1 + dx_bot], [1 - dy_bot, 1 + dy_bot],
                transform=ax_bot.transAxes, color='black', clip_on=False)

# Overall figure title
fig.suptitle(r'$\mathbf{Comparison\ mplrs\ vs.\ polco\ -\ e\ coli\ core\ onto\ exchange\ reactions}$',
             fontsize=16, y=0.95)

# Add the hierarchical legend on the right side
fig.legend(handles=legend_handles, labels=legend_labels, title='Processing Stages', 
           title_fontsize=11, bbox_to_anchor=(0.85, 0.75), loc='center left', 
           fontsize=9, handlelength=2, handletextpad=0.6)

plt.subplots_adjust(left=0.15, right=0.82, bottom=0.08, top=0.90)
plt.show()
