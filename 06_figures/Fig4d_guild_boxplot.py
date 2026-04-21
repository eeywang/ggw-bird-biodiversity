"""
Figure 4d — Guild boxplot, 180mm NEE-compliant
Scatter + box style; NEE 7pt fonts; Senegal -> Nigeria -> Ethiopia order.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
import numpy as np
import os

plt.rcParams.update({
    'font.family': 'Helvetica', 'font.size': 7, 'font.weight': 'regular',
    'axes.linewidth': 0.5, 'axes.edgecolor': 'black', 'axes.labelcolor': 'black',
    'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.color': 'black', 'ytick.color': 'black',
    'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'text.color': 'black', 'legend.fontsize': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
try:
    import matplotlib.font_manager as fm
    if not any('helvetica' in f.lower() for f in fm.findSystemFonts()):
        plt.rcParams['font.family'] = 'Arial'
except Exception: pass

COL_GGW = '#3d89be'; COL_CONTROL = '#da6852'
guilds = ['Woodland', 'Open-habitat', 'Wetland']
guild_cols = {'Woodland': 'Woodland_change', 'Open-habitat': 'Openhabitat_change', 'Wetland': 'Wetland_change'}
countries = ['Senegal', 'Nigeria', 'Ethiopia']

base_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\results'
files = {c: f'{base_dir}\\{c}\\guild_richness\\{c}_matched_with_guild_richness.csv' for c in countries}

all_data = {}
for ctry, path in files.items():
    df = pd.read_csv(path)
    df['Group'] = df['treatment'].map({1: 'GGW', 0: 'Non-GGW'})
    all_data[ctry] = df

# Layout: 3 equal segments
total_width = 12.0; seg = total_width / 3.0
seg_centers = [seg * 0.5, seg * 1.5, seg * 2.5]
country_starts = {countries[i]: seg_centers[i] - 1 for i in range(3)}
sep_positions = [seg, seg * 2]
group_offset = [-0.20, 0.20]; box_width = 0.32
np.random.seed(42)

W_180 = 180 / 25.4; H_BOX = 80 / 25.4
fig, ax = plt.subplots(1, 1, figsize=(W_180, H_BOX), facecolor='white')
for spine in ax.spines.values():
    spine.set_visible(True); spine.set_linewidth(0.5); spine.set_color('black')
ax.set_facecolor('white')
ax.tick_params(axis='both', direction='in', length=2.5, width=0.4)

all_vals = []
for ctry in countries:
    for col in guild_cols.values():
        all_vals.extend(all_data[ctry][col].dropna().tolist())
y_pad = 1.5
ylim = (max(min(all_vals) - y_pad, -16), max(max(all_vals) + y_pad, 16))
ax.set_ylim(ylim); ax.set_xlim(0, total_width)
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
ax.axhline(0, color='black', lw=0.4, ls='--', dashes=(4, 3), zorder=1)

tick_positions = []; tick_labels_list = []
for ctry in countries:
    start = country_starts[ctry]; df = all_data[ctry]
    for j, guild in enumerate(guilds):
        col = guild_cols[guild]; x_base = start + j
        for g_idx, (group, color) in enumerate([('GGW', COL_GGW), ('Non-GGW', COL_CONTROL)]):
            subset = df[df['Group'] == group][col].dropna()
            x_center = x_base + group_offset[g_idx]
            ax.boxplot(subset, positions=[x_center], widths=box_width, patch_artist=True,
                       showfliers=False, zorder=3,
                       medianprops=dict(color='black', linewidth=0.7),
                       whiskerprops=dict(color='#808080', linewidth=0.5),
                       capprops=dict(color='#808080', linewidth=0.5),
                       boxprops=dict(facecolor=color, edgecolor=color, alpha=0.50, linewidth=0.5))
            jitter = np.random.uniform(-0.10, 0.10, size=len(subset))
            ax.scatter(x_center + jitter, subset.values, s=3, color=color,
                       alpha=0.45, edgecolors='none', zorder=4, rasterized=True)
        tick_positions.append(x_base); tick_labels_list.append(guild)

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels_list, fontsize=6, rotation=30, ha='right')
for sx in sep_positions:
    ax.axvline(sx, color='#cccccc', lw=0.5, ls='-', zorder=0)

label_y = ylim[1] + (ylim[1] - ylim[0]) * 0.04
for ctry, cx in zip(countries, seg_centers):
    ax.text(cx, label_y, ctry, ha='center', va='bottom', fontsize=7, color='black', clip_on=False)

ax.set_ylabel('\u0394 species richness', fontsize=7, labelpad=3)
leg_handles = [
    Patch(facecolor=COL_GGW, edgecolor=COL_GGW, alpha=0.6, label='GGW (treatment)'),
    Patch(facecolor=COL_CONTROL, edgecolor=COL_CONTROL, alpha=0.6, label='Non-GGW (control)'),
]
ax.legend(handles=leg_handles, loc='upper left', frameon=True, framealpha=0.9,
          edgecolor='#cccccc', fontsize=6, handletextpad=0.4, borderpad=0.4)
plt.tight_layout(rect=[0, 0, 1, 0.94], pad=0.5)
plt.show()

save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    output_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures'
    os.makedirs(output_dir, exist_ok=True)
    for fmt in ['pdf', 'png']:
        fig.savefig(os.path.join(output_dir, f'Fig4d_guild_boxplot_180mm.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white')
    print('\n  Saved: Fig4d_guild_boxplot_180mm.pdf / .png')
plt.close('all')
