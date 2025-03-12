import matplotlib.pyplot as plt

# Data from the table
step_sizes = [5, 10, 30, 60, 111]
wall_clock_times = [279.73, 274.79, 138.61, 82.33, 39.82]

# Create figure with matching style
fig, ax = plt.subplots(figsize=(12, 7))

# Plot line with markers using tab20 color scheme
main_color = plt.cm.tab20.colors[0]
ax.plot(step_sizes, wall_clock_times, 
        marker='o', 
        color=main_color, 
        linestyle='-', 
        linewidth=2,
        markersize=8,
        markerfacecolor='white',
        markeredgewidth=1.5)

# Configure logarithmic scale as in original style
ax.set_yscale('log')
ax.set_xscale('log')

fig.suptitle(
    r'$\mathbf{Wall\ Clock\ Time\ vs.\ Step\ Size}$',
    fontsize=20,
    y=0.96,  # Original y=0.98 -> 0.96
    verticalalignment='top'
)

# Add second subtitle
fig.text(
    0.5,  # x-position (centered)
    0.93,  # y-position (below main title)
    r'$\mathbf{A.\ thaliana\ plastid,\ polco,\ 120\ threads}$',
    ha='center',
    va='top',
    fontsize=16,
    alpha=0.8
)

# Add titles and labels with LaTeX bold formatting
# fig.suptitle(r'$\mathbf{Wall\ Clock\ Time\ vs.\ Step\ Size}$', 
#              fontsize=14, y=0.95, weight='bold')
ax.set_xlabel('Step Size', fontsize=14)
ax.set_ylabel('Wall Clock Time (seconds)', fontsize=14)
ax.set_ylim(0, 300)

# Format x-axis ticks
ax.set_xticks(step_sizes)
ax.set_xticklabels(step_sizes)#, rotation=0, ha='right')

# Configure grid matching original style
ax.grid(True, which='both', 
        axis='both', 
        linestyle='--', 
        linewidth=0.5, 
        alpha=0.7)

# Add value annotations
for x, y in zip(step_sizes, wall_clock_times):
    ax.text(x+3, y+2, f'{y:.1f}', 
            ha='center', 
            va='bottom',
            fontsize=13,
            rotation=0)

# Adjust layout and spacing
plt.tight_layout()

plt.show()