"""
Extended Data Fig. 2 | Ensemble SDM performance distributions
Density ridgeline: 3 countries (Senegal / Nigeria / Ethiopia)
2 panels: TSS (a) and AUC (b). 180mm wide, NEE style.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os

FILES = {
    'Senegal':  r'D:\OneDrive - University of Leeds\biodiversity_manscript\results\NEE_SI\all_species_results\senegal_all_species_results.csv',
    'Nigeria':  r'D:\OneDrive - University of Leeds\biodiversity_manscript\results\NEE_SI\all_species_results\Nigeria_all_species_summary_291.csv',
    'Ethiopia': r'D:\OneDrive - University of Leeds\biodiversity_manscript\results\NEE_SI\all_species_results\ethiopia_all_species_results.csv',
}
OUTPUT_DIR = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures\ED'
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {'Senegal': '#E91E63', 'Nigeria': '#009688', 'Ethiopia': '#FF9800'}

data = {}
for country, path in FILES.items():
    df = pd.read_csv(path)
    if 'n_occurrence' in df.columns: df = df.rename(columns={'n_occurrence': 'n_occ'})
    if 'n_occ' in df.columns: df = df[df['n_occ'] >= 50]
    data[country] = df
    print(f"  {country}: {len(df)} species, TSS={df['mean_TSS'].mean():.3f}, AUC={df['mean_AUC'].mean():.3f}")

plt.rcParams.update({
    'font.family': 'Helvetica', 'font.size': 7, 'font.weight': 'regular',
    'text.color': 'black', 'axes.labelcolor': 'black', 'axes.linewidth': 0.5,
    'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'legend.fontsize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
try:
    import matplotlib.font_manager as fm
    if not any('helvetica' in f.lower() for f in fm.findSystemFonts()):
        plt.rcParams['font.family'] = 'Arial'
except Exception: pass

countries = ['Senegal', 'Nigeria', 'Ethiopia']
metrics = [('mean_TSS', 'True skill statistic (TSS)'), ('mean_AUC', 'Area under ROC curve (AUC)')]
fig, axes = plt.subplots(1, 2, figsize=(180/25.4, 75/25.4), facecolor='white',
                          gridspec_kw={'wspace': 0.30})

for ax_idx, (metric, xlabel) in enumerate(metrics):
    ax = axes[ax_idx]; ridge_height = 1.0; overlap = 0.5
    for i, country in enumerate(countries):
        vals = data[country][metric].dropna().values; n = len(vals); col = COLORS[country]
        x_min = 0.1 if metric == 'mean_TSS' else 0.5; x_max = 1.0
        x_grid = np.linspace(x_min, x_max, 300)
        try:
            kde = gaussian_kde(vals, bw_method='silverman'); density = kde(x_grid)
        except: density = np.zeros_like(x_grid)
        if density.max() > 0: density = density / density.max() * ridge_height * 0.8
        y_offset = (len(countries) - 1 - i) * (ridge_height - overlap)
        ax.fill_between(x_grid, y_offset, y_offset + density, color=col, alpha=0.35, linewidth=0)
        ax.plot(x_grid, y_offset + density, color=col, lw=0.8)
        ax.plot([x_min, x_max], [y_offset, y_offset], color=col, lw=0.3, alpha=0.5)
        ax.text(x_min - 0.02, y_offset + ridge_height * 0.3, f'{country}\n(n = {n})',
                ha='right', va='center', fontsize=7, color='black')
        med = np.median(vals); med_idx = np.argmin(np.abs(x_grid - med))
        ax.plot([med, med], [y_offset, y_offset + density[med_idx]], color='black', lw=0.5, ls='--', alpha=0.7)
        ax.text(med, y_offset + density[med_idx] + 0.02, f'{med:.2f}', ha='center', va='bottom', fontsize=6, color='black')

    thresh = 0.4 if metric == 'mean_TSS' else 0.7
    label = f'TSS = {thresh}' if metric == 'mean_TSS' else f'AUC = {thresh}'
    ax.axvline(thresh, color='#B71C1C', lw=0.5, ls=':', zorder=0)
    ax.text(thresh + 0.01, ax.get_ylim()[1] * 0.95, label, fontsize=6, color='black', va='top')
    ax.set_xlabel(xlabel, fontsize=8, color='black'); ax.set_yticks([])
    ax.spines['left'].set_visible(False); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5); ax.spines['bottom'].set_color('black')
    lbl = 'a' if ax_idx == 0 else 'b'
    ax.text(-0.10, 1.05, lbl, transform=ax.transAxes, fontsize=9, fontweight='bold', color='black', va='bottom')

plt.tight_layout(); plt.show()
save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    for fmt, kw in [('pdf',{}), ('png',{}), ('tiff',{'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(os.path.join(OUTPUT_DIR, f'ED_Fig2_SDM_performance.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white', **kw)
    print('\n\u2705 Saved: ED_Fig2_SDM_performance')
plt.close()
