"""
Figure 4e — Guild-specific causal effects (Doubly Robust ATT)
Nature Ecology & Evolution COMPLIANT, 180 x 70mm
Country order: Senegal -> Nigeria -> Ethiopia (west-to-east)
Panel label not added in code — added manually in Illustrator.
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os

W_180 = 180 / 25.4; H_FULL = 70 / 25.4
P = {
    'font.size': 6, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 5.5,
    'axes.linewidth': 0.4, 'xtick.major.width': 0.3, 'ytick.major.width': 0.3,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'ms_main': 6, 'ms_edge': 0.8, 'lw_ci': 0.7, 'lw_cap': 0.5, 'cap_w': 0.10,
    'lw_zero': 0.4, 'fs_label': 5.5, 'fs_n': 5, 'fs_title': 7, 'fs_xtick': 6,
    'fs_ylabel': 7, 'fs_legend': 5.5, 'spine_lw': 0.4, 'tick_len': 3,
}

GUILD_COL = {'Woodland': '#a9d18e', 'Open-habitat': '#ffc000', 'Wetland': '#8fb0f7'}
guilds = ['Woodland', 'Open-habitat', 'Wetland']
countries = ['Senegal', 'Nigeria', 'Ethiopia']

data = {
    'Senegal': {
        'Woodland':     {'att': -0.5086, 'lo': -1.1814, 'hi':  0.1642, 'p': 0.14049, 'n': 80},
        'Open-habitat': {'att': -1.5483, 'lo': -2.5443, 'hi': -0.5523, 'p': 0.00272, 'n': 80},
        'Wetland':      {'att': -0.9640, 'lo': -1.6304, 'hi': -0.2977, 'p': 0.00519, 'n': 80},
    },
    'Nigeria': {
        'Woodland':     {'att':  1.1826, 'lo':  0.3529, 'hi':  2.0122, 'p': 0.00540, 'n': 262},
        'Open-habitat': {'att': -0.0605, 'lo': -0.6304, 'hi':  0.5093, 'p': 0.83512, 'n': 262},
        'Wetland':      {'att':  1.1297, 'lo':  0.7753, 'hi':  1.4842, 'p': 8.61e-10,'n': 262},
    },
    'Ethiopia': {
        'Woodland':     {'att': -0.2589, 'lo': -0.9340, 'hi':  0.4162, 'p': 0.45333, 'n': 88},
        'Open-habitat': {'att': -1.0290, 'lo': -2.3398, 'hi':  0.2819, 'p': 0.12580, 'n': 88},
        'Wetland':      {'att': -0.2221, 'lo': -0.5463, 'hi':  0.1021, 'p': 0.18120, 'n': 88},
    },
}

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

fig, axes = plt.subplots(1, 3, figsize=(W_180, H_FULL), facecolor='white', sharey=True)
all_lo = [data[c][g]['lo'] for c in countries for g in guilds]
all_hi = [data[c][g]['hi'] for c in countries for g in guilds]
y_pad = 0.8; ylim = (min(all_lo) - y_pad, max(all_hi) + y_pad)
x_pos = np.arange(len(guilds))

for idx, (ctry, ax) in enumerate(zip(countries, axes)):
    for spine in ax.spines.values():
        spine.set_visible(True); spine.set_linewidth(P['spine_lw']); spine.set_color('black')
    ax.set_facecolor('white')
    ax.tick_params(axis='both', direction='in', length=P['tick_len'], width=P['xtick.major.width'])
    ax.set_ylim(ylim)
    ax.axhline(0, color='black', lw=P['lw_zero'], ls='--', zorder=1)

    for j, guild in enumerate(guilds):
        d = data[ctry][guild]; col = GUILD_COL[guild]; sig = d['p'] < 0.05
        ax.plot([x_pos[j]]*2, [d['lo'], d['hi']], color=col, lw=P['lw_ci'], zorder=2, alpha=0.7)
        for y in [d['lo'], d['hi']]:
            ax.plot([x_pos[j]-P['cap_w'], x_pos[j]+P['cap_w']], [y]*2,
                    color=col, lw=P['lw_cap'], zorder=2, alpha=0.7)
        ax.plot(x_pos[j], d['att'], 'o', color=col, ms=P['ms_main'],
                mec='white', mew=P['ms_edge'], zorder=4)
        sign = '+' if d['att'] > 0 else ''
        star = '***' if d['p']<0.001 else ('**' if d['p']<0.01 else ('*' if d['p']<0.05 else ''))
        label = f"{sign}{d['att']:.2f}{star}" if sig else f"{d['att']:.2f} (ns)"
        y_txt = d['hi'] + 0.10 if d['att'] >= 0 else d['lo'] - 0.10
        va = 'bottom' if d['att'] >= 0 else 'top'
        ax.text(x_pos[j], y_txt, label, ha='center', va=va, fontsize=P['fs_label'], color='black')

    n_val = data[ctry][guilds[0]]['n']
    ax.text(1, ylim[0]+0.06, f'n = {n_val}', ha='center', va='bottom', fontsize=P['fs_n'], color='#888888')
    ax.set_xticks(x_pos); ax.set_xticklabels(guilds, fontsize=P['fs_xtick'])
    ax.set_xlim(-0.6, len(guilds)-0.4); ax.set_title(ctry, fontsize=P['fs_title'], pad=6)

axes[0].set_ylabel('Doubly robust ATT\n(\u0394 expected species)', fontsize=P['fs_ylabel'])
leg = [Line2D([0],[0], marker='o', color='w', markerfacecolor=GUILD_COL[g],
              ms=P['ms_main']*0.7, mec='white', mew=0.5, label=g) for g in guilds]
fig.legend(handles=leg, loc='lower center', bbox_to_anchor=(0.5, -0.02),
           frameon=False, fontsize=P['fs_legend'], ncol=3, handletextpad=0.3, columnspacing=1.5)
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.show()

save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    output_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures'
    os.makedirs(output_dir, exist_ok=True)
    for fmt, kw in [('pdf',{}), ('png',{}), ('tiff',{'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(os.path.join(output_dir, f'Fig4e_guild_ATT_180mm.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white', **kw)
    print('\n  Saved: Fig4e_guild_ATT_180mm')
plt.close('all')
