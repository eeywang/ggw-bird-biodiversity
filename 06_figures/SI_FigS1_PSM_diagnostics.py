"""
Supplementary Fig. S1 | PSM matching diagnostics
3 rows (Senegal / Nigeria / Ethiopia) x 2 cols (Love plot | PS distribution)
180mm width. PS: dotted = before, solid = after; blue = GGW, orange = Non-GGW.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import os

COUNTRIES = {
    'Senegal': {
        'full': r'E:\2026.1.8_biomod2\SEN_covariates\Senegal_Adm4_covariates_corrected.csv',
        'matched': r'E:\2026.1.8_biomod2\senegal_results\doubly_robust_50threshold\Senegal_1to1_matched_data_50threshold.csv',
        'covars': ['temp_mean','prec_mean','elev_mean','slope_mean','ndvi_mean','log_pop','lulc_cropland','lulc_diversity'],
        'n_pairs': 80
    },
    'Nigeria': {
        'full': r'E:\2026.1.8_biomod2\NGA_covariates\Nigeria_LGA_covariates_full.csv',
        'matched': r'E:\2026.1.8_biomod2\nigeria_results\PSM_1to1\Nigeria_1to1_matched_data.csv',
        'covars': ['temp_mean','prec_mean','prec_cv','roughness_mean','ndvi_mean','log_pop','lulc_cropland','lulc_diversity'],
        'n_pairs': 270
    },
    'Ethiopia': {
        'full': r'E:\2026.1.8_biomod2\ETH_covariates\Ethiopia_Adm3_covariates_corrected.csv',
        'matched': r'E:\2026.1.8_biomod2\ethiopia_results\doubly_robust_50threshold\Ethiopia_1to1_matched_data_50threshold.csv',
        'covars': ['temp_mean','prec_mean','prec_cv','slope_mean','ndvi_mean','log_pop','lulc_cropland','lulc_diversity'],
        'n_pairs': 88
    }
}

VAR_LABELS = {
    'temp_mean': 'Mean temperature', 'prec_mean': 'Mean precipitation',
    'prec_cv': 'Precipitation CV', 'roughness_mean': 'Topographic roughness',
    'elev_mean': 'Elevation', 'slope_mean': 'Slope',
    'ndvi_mean': 'NDVI (2000\u20132006)', 'log_pop': 'Population density (log)',
    'lulc_cropland': 'Cropland proportion', 'lulc_diversity': 'Land-use diversity'
}
OUTPUT_DIR = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures\SI'
os.makedirs(OUTPUT_DIR, exist_ok=True)

COL_BEFORE = '#999999'; COL_AFTER = '#CC3311'; COL_TREAT = '#0077BB'; COL_CTRL = '#EE7733'

def calc_smd(t, c):
    p = np.sqrt((t.std()**2 + c.std()**2) / 2)
    return (t.mean() - c.mean()) / p if p > 0 else 0

def smooth_kde(data, x_grid, bw_factor=1.0):
    if len(data) < 5: return np.zeros_like(x_grid)
    try:
        kde = gaussian_kde(data, bw_method='silverman')
        kde.set_bandwidth(kde.factor * bw_factor)
        return kde(x_grid)
    except: return np.zeros_like(x_grid)

# Data
all_balance = {}; all_ps = {}
for country, cfg in COUNTRIES.items():
    full = pd.read_csv(cfg['full']).dropna(subset=['treatment'] + cfg['covars'])
    matched = pd.read_csv(cfg['matched']).dropna(subset=['treatment'] + cfg['covars'])
    balance = []
    for var in cfg['covars']:
        t_b, c_b = full[full['treatment']==1][var], full[full['treatment']==0][var]
        t_a, c_a = matched[matched['treatment']==1][var], matched[matched['treatment']==0][var]
        balance.append({'var': var, 'label': VAR_LABELS.get(var, var),
                        'smd_before': calc_smd(t_b, c_b), 'smd_after': calc_smd(t_a, c_a)})
    all_balance[country] = balance
    X = full[cfg['covars']].values; y = full['treatment'].values
    sc = StandardScaler(); lr = LogisticRegression(max_iter=2000, random_state=42)
    lr.fit(sc.fit_transform(X), y)
    ps_f = lr.predict_proba(sc.transform(X))[:, 1]
    X_m = matched[cfg['covars']].values; ps_m = lr.predict_proba(sc.transform(X_m))[:, 1]
    y_m = matched['treatment'].values
    all_ps[country] = {'before': {'treat': ps_f[y==1], 'ctrl': ps_f[y==0]},
                        'after': {'treat': ps_m[y_m==1], 'ctrl': ps_m[y_m==0]}}
    print(f"  {country}: T={sum(y==1)}, C={sum(y==0)} -> {cfg['n_pairs']} pairs")

plt.rcParams.update({
    'font.family': 'Helvetica', 'font.size': 6, 'font.weight': 'regular',
    'axes.linewidth': 0.4, 'axes.labelsize': 7, 'axes.titlesize': 7.5,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'legend.fontsize': 5.5, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
try:
    import matplotlib.font_manager as fm
    if not any('helvetica' in f.lower() for f in fm.findSystemFonts()):
        plt.rcParams['font.family'] = 'Arial'
except Exception: pass

W = 180 / 25.4; H = 200 / 25.4
fig, axes = plt.subplots(3, 2, figsize=(W, H), facecolor='white',
                          gridspec_kw={'width_ratios': [1.15, 1], 'hspace': 0.35, 'wspace': 0.30})
ctry_list = ['Senegal', 'Nigeria', 'Ethiopia']
panel_labels = [('a','b'), ('c','d'), ('e','f')]

for row, country in enumerate(ctry_list):
    bal = all_balance[country]; ps_data = all_ps[country]; n_pairs = COUNTRIES[country]['n_pairs']

    # Love plot
    ax_love = axes[row, 0]
    labels = [b['label'] for b in bal]; smd_b = [b['smd_before'] for b in bal]; smd_a = [b['smd_after'] for b in bal]
    y_pos = np.arange(len(labels))
    for th in [-0.25, -0.1, 0.1, 0.25]:
        ls = '--' if abs(th)==0.25 else ':'; cl = '#CC0000' if abs(th)==0.25 else '#FF8800'
        ax_love.axvline(th, color=cl, lw=0.4, ls=ls, zorder=0)
    ax_love.axvline(0, color='black', lw=0.3, zorder=0)
    ax_love.scatter(smd_b, y_pos, c=COL_BEFORE, s=18, marker='o', edgecolors=COL_BEFORE, lw=0.3, zorder=3, label='Before matching', alpha=0.7)
    ax_love.scatter(smd_a, y_pos, c=COL_AFTER, s=22, marker='D', edgecolors=COL_AFTER, lw=0.3, zorder=4, label='After matching')
    for i in range(len(labels)):
        ax_love.plot([smd_b[i], smd_a[i]], [y_pos[i]]*2, color='#cccccc', lw=0.4, zorder=1)
    ax_love.set_yticks(y_pos); ax_love.set_yticklabels(labels, fontsize=5.5)
    ax_love.set_xlabel('Standardised mean difference', fontsize=6.5); ax_love.invert_yaxis()
    max_abs = max(max(abs(s) for s in smd_b), max(abs(s) for s in smd_a))
    ax_love.set_xlim(-max(1.0, np.ceil(max_abs*2)/2+0.5), max(1.0, np.ceil(max_abs*2)/2+0.5))
    ax_love.set_title(f'{country} (n = {n_pairs} pairs)', fontsize=7.5, loc='left', pad=4)
    ax_love.text(-0.18, 1.05, panel_labels[row][0], transform=ax_love.transAxes, fontsize=9, fontweight='bold', va='bottom')
    if row == 0: ax_love.legend(loc='lower right', frameon=True, framealpha=0.9, edgecolor='#cccccc', fontsize=5, handletextpad=0.3, markerscale=0.8)

    # PS distribution
    ax_ps = axes[row, 1]; x_grid = np.linspace(0, 1, 500)
    for grp, col, ls_b in [('treat', COL_TREAT, ':'), ('ctrl', COL_CTRL, ':')]:
        lbl_before = f"{'GGW' if grp=='treat' else 'Non-GGW'} (before)"
        lbl_after = f"{'GGW' if grp=='treat' else 'Non-GGW'} (after)"
        kde_b = smooth_kde(ps_data['before'][grp], x_grid, 1.2)
        kde_a = smooth_kde(ps_data['after'][grp], x_grid, 1.2)
        ax_ps.fill_between(x_grid, kde_b, color=col, alpha=0.10)
        ax_ps.plot(x_grid, kde_b, color=col, lw=0.6, ls=':', label=lbl_before)
        ax_ps.plot(x_grid, kde_a, color=col, lw=1.2, label=lbl_after)
    ax_ps.set_xlabel('Propensity score', fontsize=6.5); ax_ps.set_ylabel('Density', fontsize=6.5)
    ax_ps.set_xlim(0, 1)
    ax_ps.set_title(country, fontsize=7.5, loc='left', pad=4)
    ax_ps.text(-0.15, 1.05, panel_labels[row][1], transform=ax_ps.transAxes, fontsize=9, fontweight='bold', va='bottom')
    if row == 0: ax_ps.legend(loc='upper center', frameon=True, framealpha=0.9, edgecolor='#cccccc', fontsize=5, ncol=2, columnspacing=0.8)
    ax_ps.text(0.97, 0.93, f'n = {n_pairs} pairs', transform=ax_ps.transAxes, fontsize=5, color='#666666', ha='right', va='top')

plt.tight_layout(); plt.show()
save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    for fmt, kw in [('pdf',{}), ('png',{}), ('tiff',{'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(os.path.join(OUTPUT_DIR, f'SI_FigS1_PSM_diagnostics.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white', **kw)
    print('\n\u2705 Saved: SI_FigS1_PSM_diagnostics')
plt.close()
