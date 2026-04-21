"""
Extended Data Fig. 1e — Scenario pathway decomposition boxplot (180mm)
3 faceted subplots (Climate | Vegetation | Interaction), each with independent Y-axis.
Within each: 3 countries (Senegal -> Nigeria -> Ethiopia), GGW vs Non-GGW boxes.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
import numpy as np
import os

W_180 = 180 / 25.4; H_BOX = 65 / 25.4
P = {
    'font.size': 6, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 5.5,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.3, 'ytick.major.width': 0.3,
    'ytick.minor.width': 0.3, 'xtick.major.size': 3, 'ytick.major.size': 3,
    'ytick.minor.size': 1.5,
    'spine_lw': 0.4, 'tick_len': 3,
    'median_lw': 0.9, 'whisker_lw': 0.55, 'cap_lw': 0.55, 'box_lw': 0.55,
    'zero_lw': 0.5, 'sep_lw': 0.5,
    'scatter_s': 4, 'scatter_alpha': 0.25, 'jitter': 0.07,
    'box_width': 0.24, 'group_offset': [-0.18, 0.18],
    'fs_xtick': 6, 'fs_ylabel': 7, 'fs_subtitle': 7, 'fs_legend': 6,
}
COL_GGW = '#43a959'; COL_CONTROL = '#fecb7a'
components = ['Climate', 'Vegetation', 'Interaction']
component_cols = {'Climate': 'prob_climate', 'Vegetation': 'prob_vegetation', 'Interaction': 'prob_interaction'}
countries = ['Senegal', 'Nigeria', 'Ethiopia']

base_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\results'
files = {
    'Senegal':  os.path.join(base_dir, 'Senegal',  'Senegal_matched_with_scenarios.xlsx'),
    'Nigeria':  os.path.join(base_dir, 'Nigeria',  'Nigria_matched_with_scenarios.xlsx'),
    'Ethiopia': os.path.join(base_dir, 'Ethiopia', 'Ethiopia_matched_with_scenarios.xlsx'),
}
all_data = {}
for ctry in countries:
    df = pd.read_excel(files[ctry])
    df['Group'] = df['treatment'].map({1: 'GGW', 0: 'Non-GGW'})
    all_data[ctry] = df

plt.rcParams.update({
    'font.family': 'Helvetica', 'font.size': P['font.size'], 'font.weight': 'regular',
    'axes.linewidth': P['axes.linewidth'], 'axes.edgecolor': 'black', 'axes.labelcolor': 'black',
    'xtick.direction': 'in', 'ytick.direction': 'in', 'text.color': 'black',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
try:
    import matplotlib.font_manager as fm
    if not any('helvetica' in f.lower() for f in fm.findSystemFonts()):
        plt.rcParams['font.family'] = 'Arial'
except Exception: pass

np.random.seed(42)
X_POS = [1.0, 2.0, 3.0]; XLIM = (0.4, 3.6)
fig, axes = plt.subplots(1, 3, figsize=(W_180, H_BOX), sharey=False, facecolor='white',
                          gridspec_kw={'wspace': 0.28})

for ax_idx, (ax, comp) in enumerate(zip(axes, components)):
    col = component_cols[comp]
    for spine in ax.spines.values():
        spine.set_visible(True); spine.set_linewidth(P['spine_lw']); spine.set_color('black')
    ax.set_facecolor('white')
    ax.tick_params(axis='both', which='major', direction='in', length=P['tick_len'], width=P['axes.linewidth'])
    ax.tick_params(axis='x', length=0)

    comp_vals = []
    for ctry in countries: comp_vals.extend(all_data[ctry][col].dropna().tolist())
    lo, hi = min(comp_vals), max(comp_vals); rng = hi - lo
    ax.set_ylim(lo - rng*0.10, hi + rng*0.10); ax.set_xlim(*XLIM)
    for step in [0.5, 1, 2, 2.5, 5, 10, 20]:
        if rng / step <= 6: break
    ax.yaxis.set_major_locator(ticker.MultipleLocator(step))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(step / 2))
    ax.axhline(0, color='black', lw=P['zero_lw'], ls='--', dashes=(4,3), zorder=1)

    for i, ctry in enumerate(countries):
        x_base = X_POS[i]; df = all_data[ctry]
        for g_idx, (group, color) in enumerate([('GGW', COL_GGW), ('Non-GGW', COL_CONTROL)]):
            subset = df[df['Group'] == group][col].dropna()
            x_center = x_base + P['group_offset'][g_idx]
            ax.boxplot(subset, positions=[x_center], widths=P['box_width'], patch_artist=True,
                       showfliers=False, zorder=3,
                       medianprops=dict(color='black', linewidth=P['median_lw']),
                       whiskerprops=dict(color='#606060', linewidth=P['whisker_lw']),
                       capprops=dict(color='#606060', linewidth=P['cap_lw']),
                       boxprops=dict(facecolor=color, edgecolor=color, alpha=0.55, linewidth=P['box_lw']))
            jitter = np.random.uniform(-P['jitter'], P['jitter'], size=len(subset))
            ax.scatter(x_center+jitter, subset.values, s=P['scatter_s'], color=color,
                       alpha=P['scatter_alpha'], edgecolors='none', zorder=2, rasterized=True)

    ax.set_xticks(X_POS); ax.set_xticklabels(countries, fontsize=P['fs_xtick'], rotation=30, ha='right')
    ax.set_title(comp, fontsize=P['fs_subtitle'], fontweight='regular', pad=4)
    if ax_idx == 0: ax.set_ylabel('\u0394 expected species richness', fontsize=P['fs_ylabel'], labelpad=4)

leg_h = [Patch(facecolor=COL_GGW, edgecolor=COL_GGW, alpha=0.6, label='GGW (treatment)'),
         Patch(facecolor=COL_CONTROL, edgecolor=COL_CONTROL, alpha=0.6, label='Non-GGW (control)')]
fig.legend(handles=leg_h, loc='upper center', bbox_to_anchor=(0.5, 1.00), ncol=2,
           frameon=False, fontsize=P['fs_legend'], handletextpad=0.4, columnspacing=1.5)
plt.tight_layout(rect=[0, 0, 1, 0.93], pad=0.5, w_pad=1.0)
plt.show()

save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    output_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures'
    os.makedirs(output_dir, exist_ok=True)
    for fmt, kw in [('pdf',{}), ('png',{}), ('tiff',{'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(os.path.join(output_dir, f'EDFig1e_scenario_pathway_boxplot_180mm.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white', **kw)
    print('\n\u2705 Saved: EDFig1e_scenario_pathway_boxplot_180mm')
plt.close()
