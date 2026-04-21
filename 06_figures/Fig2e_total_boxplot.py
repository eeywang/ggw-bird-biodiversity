"""
Figure 2e — Scenario decomposition boxplot (180mm strip)
Nature Ecology & Evolution COMPLIANT

Only TOTAL change. Country order: Senegal -> Nigeria -> Ethiopia (west-to-east).
Each country section: GGW vs Non-GGW boxplot pair.
180mm strip to align under 88mm maps in Illustrator.
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import os

W_180 = 180 / 25.4
H_BOX =  65 / 25.4

P = {
    'font.size': 6, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 5.5,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.3, 'ytick.major.width': 0.3,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'spine_lw': 0.4, 'tick_len': 3,
    'median_lw': 1.0, 'whisker_lw': 0.5, 'cap_lw': 0.5, 'box_lw': 0.6,
    'zero_lw': 0.4, 'sep_lw': 0.5,
    'scatter_s': 3.5, 'scatter_alpha': 0.40, 'jitter': 0.08,
    'box_width': 0.28, 'group_offset': [-0.25, 0.25],
    'fs_ylabel': 7, 'fs_country': 7, 'fs_legend': 5.5,
}

COL_GGW     = '#43a959'
COL_CONTROL = '#fecb7a'
countries = ['Senegal', 'Nigeria', 'Ethiopia']

# ---- DATA ----
base_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\results'
files = {
    'Nigeria':  os.path.join(base_dir, 'Nigeria',  'Nigria_matched_with_scenarios.xlsx'),
    'Senegal':  os.path.join(base_dir, 'Senegal',  'Senegal_matched_with_scenarios.xlsx'),
    'Ethiopia': os.path.join(base_dir, 'Ethiopia', 'Ethiopia_matched_with_scenarios.xlsx'),
}

all_data = {}
for ctry, path in files.items():
    df = pd.read_excel(path)
    df['Group'] = df['treatment'].map({1: 'GGW', 0: 'Non-GGW'})
    all_data[ctry] = df

# ---- STYLE ----
plt.rcParams.update({
    'font.family': 'Helvetica', 'font.size': P['font.size'],
    'font.weight': 'regular',
    'axes.linewidth': P['axes.linewidth'], 'axes.edgecolor': 'black',
    'axes.labelcolor': 'black', 'axes.labelsize': P['axes.labelsize'],
    'xtick.labelsize': P['xtick.labelsize'], 'ytick.labelsize': P['ytick.labelsize'],
    'xtick.color': 'black', 'ytick.color': 'black',
    'xtick.major.width': P['xtick.major.width'], 'ytick.major.width': P['ytick.major.width'],
    'xtick.major.size': P['xtick.major.size'], 'ytick.major.size': P['ytick.major.size'],
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'text.color': 'black', 'legend.fontsize': P['legend.fontsize'],
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

try:
    import matplotlib.font_manager as fm
    system_fonts = [f.lower() for f in fm.findSystemFonts()]
    if not any('helvetica' in f for f in system_fonts):
        plt.rcParams['font.family'] = (
            'Arial' if any('arial' in f for f in system_fonts) else 'DejaVu Sans')
except Exception:
    pass

# ---- LAYOUT ----
total_x = 6.0
seg = total_x / 3.0
seg_centers = [seg * (i + 0.5) for i in range(3)]
sep_positions = [seg, seg * 2]
np.random.seed(42)

# ---- FIGURE ----
fig, ax = plt.subplots(1, 1, figsize=(W_180, H_BOX), facecolor='white')

for spine in ax.spines.values():
    spine.set_visible(True); spine.set_linewidth(P['spine_lw']); spine.set_color('black')
ax.set_facecolor('white')
ax.tick_params(axis='both', direction='in', length=P['tick_len'], width=P['xtick.major.width'])

all_vals = []
for ctry in countries:
    all_vals.extend(all_data[ctry]['prob_total'].dropna().tolist())
y_pad = 1.5
ylim = (min(all_vals) - y_pad, max(all_vals) + y_pad)
ax.set_ylim(ylim); ax.set_xlim(0, total_x)
ax.axhline(0, color='black', lw=P['zero_lw'], ls='--', zorder=1)

for i, ctry in enumerate(countries):
    df = all_data[ctry]
    x_center = seg_centers[i]
    for g_idx, (group, color) in enumerate([('GGW', COL_GGW), ('Non-GGW', COL_CONTROL)]):
        subset = df[df['Group'] == group]['prob_total'].dropna()
        x_pos = x_center + P['group_offset'][g_idx]
        ax.boxplot(subset, positions=[x_pos], widths=P['box_width'], patch_artist=True,
                   showfliers=False, zorder=3,
                   medianprops=dict(color='black', linewidth=P['median_lw']),
                   whiskerprops=dict(color='#808080', linewidth=P['whisker_lw']),
                   capprops=dict(color='#808080', linewidth=P['cap_lw']),
                   boxprops=dict(facecolor=color, edgecolor=color,
                                 alpha=0.50, linewidth=P['box_lw']))
        jitter = np.random.uniform(-P['jitter'], P['jitter'], size=len(subset))
        ax.scatter(x_pos + jitter, subset.values, s=P['scatter_s'], color=color,
                   alpha=P['scatter_alpha'], edgecolors='none', zorder=4)

ax.set_xticks(seg_centers)
ax.set_xticklabels(countries, fontsize=P['fs_country'])
for sx in sep_positions:
    ax.axvline(sx, color='#cccccc', lw=P['sep_lw'], ls='-', zorder=0)
ax.set_ylabel('\u0394 Expected species richness', fontsize=P['fs_ylabel'])

leg_handles = [
    Patch(facecolor=COL_GGW, edgecolor=COL_GGW, alpha=0.6, label='GGW (treatment)'),
    Patch(facecolor=COL_CONTROL, edgecolor=COL_CONTROL, alpha=0.6, label='Non-GGW (control)'),
]
ax.legend(handles=leg_handles, loc='upper left', frameon=True, framealpha=0.9,
          edgecolor='#cccccc', fontsize=P['fs_legend'], handletextpad=0.3)

plt.tight_layout()
plt.show()

save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    output_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures'
    os.makedirs(output_dir, exist_ok=True)
    for fmt, kw in [('pdf', {}), ('png', {}),
                     ('tiff', {'pil_kwargs': {'compression': 'tiff_lzw'}})]:
        fig.savefig(os.path.join(output_dir, f'Fig2e_total_boxplot_180mm.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white', **kw)
    print(f'\n\u2705 Saved: Fig2e_total_boxplot_180mm')
plt.close()
