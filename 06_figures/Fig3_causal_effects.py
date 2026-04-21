"""
Figure 3 — Causal effects of the Great Green Wall on bird species richness
Nature Ecology & Evolution COMPLIANT

Layout (3x2 grid):
  Row 1: (a) Overall DR ATT | (b1) Naive->DR dumbbell | (b2) Descriptive->DR dumbbell
  Row 2: (c) Senegal bars   | (c) Nigeria bars         | (c) Ethiopia bars

Country order: Senegal -> Nigeria -> Ethiopia (west-to-east)
Marker convention: filled circle = DR ATT, diamond = Naive PSM, triangle = descriptive
"""

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
import numpy as np
import os

# ---- DIMENSIONS ----
FIG_W_MM, FIG_H_MM = 180, 120
FIG_W_IN, FIG_H_IN = FIG_W_MM / 25.4, FIG_H_MM / 25.4

# ---- FONT SIZES ----
FS = {'axis': 7, 'tick': 6, 'title': 7, 'panel': 8, 'label': 5.5,
      'label_sm': 5, 'sample': 5, 'annot': 5, 'legend': 5.5, 'footer': 5,
      'pathway': 5.5, 'country': 6}

# ---- LINE/MARKER ----
LW = {'spine': 0.4, 'zero': 0.4, 'ci': 0.7, 'cap': 0.5, 'connect': 0.8,
      'bar_edge': 0.4, 'err': 0.6, 'tick': 0.3}
MS = {'main': 6, 'dr': 5, 'est1': 4, 'edge': 0.8}
ERR_CAPSIZE, ERR_CAPTHICK = 2, 0.4
BAR_WIDTH, DODGE, CAP_B = 0.55, 0.15, 0.06

# ---- COLOURS ----
PAL = {
    'Nigeria':  {'dark': '#009688', 'mid': '#4DB6AC'},
    'Senegal':  {'dark': '#E91E63', 'mid': '#F06292'},
    'Ethiopia': {'dark': '#FF9800', 'mid': '#FFB74D'},
}
PATH_COL = {'Climate': '#42A5F5', 'Vegetation': '#66BB6A',
            'Interaction': '#FFA726', 'Total': '#78909C'}

# ---- STYLE ----
plt.rcParams.update({
    'font.family': 'Helvetica', 'font.size': 6, 'font.weight': 'regular',
    'axes.linewidth': LW['spine'], 'axes.edgecolor': 'black', 'axes.labelcolor': 'black',
    'axes.labelsize': FS['axis'], 'axes.titlesize': FS['title'], 'axes.titleweight': 'regular',
    'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'],
    'xtick.color': 'black', 'ytick.color': 'black',
    'xtick.major.width': LW['tick'], 'ytick.major.width': LW['tick'],
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'text.color': 'black', 'legend.fontsize': FS['legend'],
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
try:
    import matplotlib.font_manager as fm
    if not any('helvetica' in f.lower() for f in fm.findSystemFonts()):
        plt.rcParams['font.family'] = 'Arial'
except Exception: pass

# ---- DATA ----
countries = ['Senegal', 'Nigeria', 'Ethiopia']

dr_data = {
    'Senegal':  {'att': -3.30, 'lo': -5.58, 'hi': -1.03, 'p': 0.005, 'n': 80},
    'Nigeria':  {'att':  2.19, 'lo':  0.66, 'hi':  3.72, 'p': 0.005, 'n': 262},
    'Ethiopia': {'att': -1.59, 'lo': -3.92, 'hi':  0.73, 'p': 0.181, 'n': 88},
}
naive_att = {
    'Senegal':  {'att': -3.74, 'lo': -6.20, 'hi': -1.28},
    'Nigeria':  {'att':  7.45, 'lo':  4.50, 'hi': 10.40},
    'Ethiopia': {'att': -2.00, 'lo': -4.80, 'hi':  0.80},
}
desc_att = {
    'Senegal':  {'att':  0.76, 'lo': -0.58, 'hi':  2.10},
    'Nigeria':  {'att':  3.04, 'lo':  2.41, 'hi':  3.66},
    'Ethiopia': {'att': -2.26, 'lo': -5.25, 'hi':  0.73},
}
pathway_c = {
    'Senegal': {
        'Climate':     {'att': -2.16, 'lo': -4.10, 'hi': -0.22, 'p': 0.03},
        'Vegetation':  {'att': -0.92, 'lo': -2.50, 'hi':  0.66, 'p': 0.25},
        'Interaction': {'att': -0.23, 'lo': -0.44, 'hi': -0.02, 'p': 0.03},
        'Total':       {'att': -3.30, 'lo': -5.58, 'hi': -1.03, 'p': 0.005},
    },
    'Nigeria': {
        'Climate':     {'att':  0.08, 'lo': -0.85, 'hi':  1.01, 'p': 0.50},
        'Vegetation':  {'att':  2.28, 'lo':  0.78, 'hi':  3.78, 'p': 0.003},
        'Interaction': {'att': -0.17, 'lo': -0.30, 'hi': -0.04, 'p': 0.008},
        'Total':       {'att':  2.19, 'lo':  0.66, 'hi':  3.72, 'p': 0.005},
    },
    'Ethiopia': {
        'Climate':     {'att': -1.27, 'lo': -3.50, 'hi':  0.96, 'p': 0.26},
        'Vegetation':  {'att': -0.26, 'lo': -1.20, 'hi':  0.68, 'p': 0.59},
        'Interaction': {'att': -0.06, 'lo': -0.35, 'hi':  0.23, 'p': 0.68},
        'Total':       {'att': -1.59, 'lo': -3.92, 'hi':  0.73, 'p': 0.181},
    },
}

# ---- HELPERS ----
def style_box(ax):
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(LW['spine']); s.set_color('black')
    ax.set_facecolor('white')
    ax.tick_params(axis='both', direction='in', length=3, width=LW['tick'])

# ---- LAYOUT ----
fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN), facecolor='white')
gs = GridSpec(2, 3, hspace=0.50, wspace=0.42, top=0.92, bottom=0.10, left=0.08, right=0.97)
ax_a  = fig.add_subplot(gs[0, 0])
ax_b1 = fig.add_subplot(gs[0, 1])
ax_b2 = fig.add_subplot(gs[0, 2])
ax_c1 = fig.add_subplot(gs[1, 0])
ax_c2 = fig.add_subplot(gs[1, 1])
ax_c3 = fig.add_subplot(gs[1, 2])

# ---- PANEL (a): Overall DR ATT ----
style_box(ax_a)
x_pos = np.arange(len(countries))
ylim_a = (min(dr_data[c]['lo'] for c in countries) - 1.8,
          max(dr_data[c]['hi'] for c in countries) + 1.8)
ax_a.set_ylim(ylim_a)
for i, ctry in enumerate(countries):
    d = dr_data[ctry]; col = PAL[ctry]['dark']
    ax_a.plot([x_pos[i]]*2, [d['lo'], d['hi']], color=col, lw=LW['ci'], zorder=2, alpha=0.7)
    for y in [d['lo'], d['hi']]:
        ax_a.plot([x_pos[i]-0.10, x_pos[i]+0.10], [y]*2, color=col, lw=LW['cap'], zorder=2, alpha=0.7)
    ax_a.plot(x_pos[i], d['att'], 'o', color=col, ms=MS['main'], mec='white', mew=MS['edge'], zorder=4)
    sign = '+' if d['att'] > 0 else ''
    star = '**' if d['p']<0.01 else ('*' if d['p']<0.05 else '')
    lbl = f"{sign}{d['att']:.2f}{star}" if d['p']<0.05 else f"{d['att']:.2f} (ns)"
    y_t = d['hi']+0.18 if d['att']>=0 else d['lo']-0.18
    ax_a.text(x_pos[i], y_t, lbl, ha='center', va='bottom' if d['att']>=0 else 'top',
              fontsize=FS['label'], color='black')
    ax_a.text(x_pos[i], ylim_a[0]+0.12, f'n = {d["n"]}', ha='center', va='bottom',
              fontsize=FS['sample'], color='#888888')

ax_a.axhline(0, color='black', lw=LW['zero'], ls='--', zorder=1)
ax_a.set_xticks(x_pos); ax_a.set_xticklabels(countries, fontsize=FS['country'])
ax_a.set_ylabel('Doubly robust ATT\n(\u0394 expected species)', fontsize=FS['axis'])
ax_a.set_xlim(-0.6, len(countries)-0.4)
ax_a.text(-0.12, 1.06, 'a', transform=ax_a.transAxes, fontsize=FS['panel'], fontweight='bold', va='top')
ax_a.set_title('Overall causal effect', fontsize=FS['title'], pad=8)

# ---- PANEL (b): Dumbbell comparisons ----
x_pos_b = np.arange(len(countries)) * 1.4

def draw_paired_dot(ax, est1_data, est2_data, subtitle, est1_marker='D', annotation_func=None):
    style_box(ax)
    ax.axhline(0, color='black', lw=LW['zero'], ls='--', zorder=1)
    for i, ctry in enumerate(countries):
        x = x_pos_b[i]; col = PAL[ctry]['dark']
        d2 = est2_data[ctry]; v1 = est1_data[ctry]['att']; v2 = d2['att']
        x_e1, x_dr = x - DODGE, x + DODGE
        ax.plot([x_e1, x_dr], [v1, v2], color=col, lw=LW['connect'], alpha=0.18, zorder=2)
        # DR CI (solid)
        ax.plot([x_dr]*2, [d2['lo'], d2['hi']], color=col, lw=LW['ci']*0.8, alpha=0.7, zorder=3)
        for y in [d2['lo'], d2['hi']]:
            ax.plot([x_dr-CAP_B, x_dr+CAP_B], [y]*2, color=col, lw=LW['cap'], alpha=0.7, zorder=3)
        # est1 CI (dashed)
        if 'lo' in est1_data[ctry]:
            e1 = est1_data[ctry]
            ax.plot([x_e1]*2, [e1['lo'], e1['hi']], color=col, lw=LW['ci']*0.6, ls='--', alpha=0.45, zorder=3)
            for y in [e1['lo'], e1['hi']]:
                ax.plot([x_e1-CAP_B, x_e1+CAP_B], [y]*2, color=col, lw=LW['cap']*0.8, ls='--', alpha=0.45, zorder=3)
        ax.plot(x_e1, v1, marker=est1_marker, ms=MS['est1'], mfc='white', mec=col, mew=1.0, zorder=5)
        ax.plot(x_dr, v2, marker='o', ms=MS['dr'], mfc=col, mec='white', mew=MS['edge'], zorder=6)
        s1 = '+' if v1>0 else ''; s2 = '+' if v2>0 else ''
        ax.text(x_e1-CAP_B-0.06, v1, f'{s1}{v1:.2f}', ha='right', va='center', fontsize=FS['label_sm'], color=col, alpha=0.50)
        ax.text(x_dr+CAP_B+0.06, v2, f'{s2}{v2:.2f}', ha='left', va='center', fontsize=FS['label_sm'], color=col)
    ax.set_xticks(x_pos_b); ax.set_xticklabels(countries, fontsize=FS['country'])
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.7, x_pos_b[-1]+0.7)
    ax.set_ylabel('\u0394 Expected species richness', fontsize=FS['axis'])
    ax.set_title(subtitle, fontsize=FS['title'], pad=8)
    if annotation_func: annotation_func(ax)

# b1: Naive PSM -> DR
def annotate_b1(ax):
    ng_idx = countries.index('Nigeria'); x_ng = x_pos_b[ng_idx]
    ng_naive, ng_dr = naive_att['Nigeria']['att'], dr_data['Nigeria']['att']
    mid_y = (ng_naive + ng_dr) / 2
    br, bl = x_ng + 0.75, x_ng + 0.85
    ax.plot([br,bl,bl,br], [ng_naive,ng_naive,ng_dr,ng_dr], color='#888888', lw=0.6, alpha=0.70, zorder=1)
    ax.text(bl+0.08, mid_y, '71% climate\nconfounding', ha='left', va='center',
            fontsize=FS['annot'], color='#666666', style='italic', linespacing=1.2)

draw_paired_dot(ax_b1, naive_att, dr_data, 'Na\u00efve PSM \u2192 Doubly robust ATT', 'D', annotate_b1)
yp = 1.8
all_y1 = [naive_att[c]['att'] for c in countries] + [dr_data[c]['lo'] for c in countries] + \
         [dr_data[c]['hi'] for c in countries] + [naive_att[c].get('lo',naive_att[c]['att']) for c in countries] + \
         [naive_att[c].get('hi',naive_att[c]['att']) for c in countries]
ax_b1.set_ylim(min(all_y1)-yp, max(all_y1)+yp)
ax_b1.text(-0.14, 1.06, 'b', transform=ax_b1.transAxes, fontsize=FS['panel'], fontweight='bold', va='top')

# b2: Descriptive -> DR
def annotate_b2(ax):
    sn_idx = countries.index('Senegal'); x_sn = x_pos_b[sn_idx]
    sn_desc, sn_dr = desc_att['Senegal']['att'], dr_data['Senegal']['att']
    mid_y = (sn_desc + sn_dr) / 2
    br, bl = x_sn + 0.7, x_sn + 0.8
    ax.plot([br,bl,bl,br], [sn_desc,sn_desc,sn_dr,sn_dr], color='#888888', lw=0.6, alpha=0.70, zorder=1)
    ax.text(bl+0.08, mid_y, 'Sign reversal\n(+0.76 \u2192 \u22123.30)', ha='left', va='center',
            fontsize=FS['annot'], color='#666666', style='italic', linespacing=1.2)

draw_paired_dot(ax_b2, desc_att, dr_data, 'Within-zone descriptive \u2192 Doubly robust ATT', '^', annotate_b2)
all_y2 = [desc_att[c]['att'] for c in countries] + [desc_att[c]['lo'] for c in countries] + \
         [desc_att[c]['hi'] for c in countries] + [dr_data[c]['lo'] for c in countries] + [dr_data[c]['hi'] for c in countries]
ax_b2.set_ylim(min(all_y2)-yp, max(all_y2)+yp)

# Shared b legend
leg_b = [
    (Line2D([0],[0], marker='o', color='w', markerfacecolor='#555', markeredgecolor='white', ms=4, mew=0.5), 'Doubly robust ATT'),
    (Line2D([0],[0], marker='D', color='w', markerfacecolor='white', markeredgecolor='#555', ms=4, mew=0.8), 'Na\u00efve PSM estimate'),
    (Line2D([0],[0], marker='^', color='w', markerfacecolor='white', markeredgecolor='#555', ms=4, mew=0.8), 'Within-zone descriptive'),
    (Line2D([0],[0], color='#555', lw=0.8, ls='-'), '95% CI (causal)'),
    (Line2D([0],[0], color='#555', lw=0.6, ls='--', alpha=0.5), '95% CI (reference)'),
]
ax_b2.legend([h for h,_ in leg_b], [l for _,l in leg_b], loc='lower center',
             bbox_to_anchor=(-0.20, -0.30), frameon=False, fontsize=FS['legend'],
             ncol=5, handletextpad=0.3, columnspacing=0.8)

# ---- PANEL (c): Pathway bars ----
pathway_order = ['Climate', 'Vegetation', 'Interaction', 'Total']
pathway_labels = ['Climate\npathway', 'Vegetation\npathway', 'Interaction', 'Total\neffect']
bar_colors = [PATH_COL[p] for p in pathway_order]

def draw_pathway_panel(ax, ctry, panel_label=None):
    style_box(ax); data = pathway_c[ctry]; x = np.arange(len(pathway_order))
    all_lo_c = min(d['lo'] for c in pathway_c.values() for d in c.values())
    all_hi_c = max(d['hi'] for c in pathway_c.values() for d in c.values())
    ylim = (min(all_lo_c-0.8, -0.5), max(all_hi_c+0.8, 0.5))
    ax.axhline(0, color='black', lw=LW['zero'], ls='--', zorder=1)
    for j, pw in enumerate(pathway_order):
        d = data[pw]; sig = d['p']<0.05; col = bar_colors[j]
        kw = dict(width=BAR_WIDTH, color=col, edgecolor=col, linewidth=LW['bar_edge'], zorder=3)
        if sig: ax.bar(x[j], d['att'], alpha=0.85, **kw)
        else: ax.bar(x[j], d['att'], alpha=0.22, hatch='///', **kw)
        ax.errorbar(x[j], d['att'], yerr=[[d['att']-d['lo']], [d['hi']-d['att']]],
                    fmt='none', ecolor=col, elinewidth=LW['err'], capsize=ERR_CAPSIZE,
                    capthick=ERR_CAPTHICK, zorder=4, alpha=0.85)
        sign = '+' if d['att']>0 else ''
        star = '***' if d['p']<0.001 else ('**' if d['p']<0.01 else ('*' if d['p']<0.05 else ''))
        lbl = f"{sign}{d['att']:.2f}{star}" if sig else f"{sign}{d['att']:.2f} (ns)"
        y_t = d['hi']+0.12 if d['att']>=0 else d['lo']-0.12
        ax.text(x[j], y_t, lbl, ha='center', va='bottom' if d['att']>=0 else 'top',
                fontsize=FS['label'], color='black')
    ax.set_xticks(x); ax.set_xticklabels(pathway_labels, fontsize=FS['pathway'])
    ax.set_ylim(ylim); ax.set_title(ctry, fontsize=FS['title'], pad=8)
    if panel_label:
        ax.text(-0.14, 1.06, panel_label, transform=ax.transAxes,
                fontsize=FS['panel'], fontweight='bold', va='top', color='black')

draw_pathway_panel(ax_c1, 'Senegal', panel_label='c')
draw_pathway_panel(ax_c2, 'Nigeria')
draw_pathway_panel(ax_c3, 'Ethiopia')
ax_c1.set_ylabel('ATT (\u0394 expected species richness)', fontsize=FS['axis'])

# Bottom legend
leg_h = [Patch(facecolor=PATH_COL[p], alpha=0.85, edgecolor=PATH_COL[p], lw=0.4, label=f'{p} pathway' if p!='Total' else 'Total effect')
         for p in pathway_order]
leg_h.append(Patch(facecolor='white', edgecolor='#888888', hatch='///', lw=0.4, label='Not significant'))
fig.legend(handles=leg_h, loc='lower center', bbox_to_anchor=(0.5, 0.005), frameon=False,
           fontsize=FS['legend'], ncol=5, handletextpad=0.3, columnspacing=1.2, handlelength=1.2)
fig.text(0.5, -0.005, 'Error bars = 95% CI.  Solid = p < 0.05.  Hatched = ns.  *p<0.05  **p<0.01  ***p<0.001',
         ha='center', fontsize=FS['footer'], color='#777777')

plt.show()

save = input('\nSave? (y/n): ').strip().lower()
if save == 'y':
    output_dir = r'D:\OneDrive - University of Leeds\biodiversity_manscript\figures'
    os.makedirs(output_dir, exist_ok=True)
    for fmt, kw in [('pdf',{}), ('png',{}), ('tiff',{'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(os.path.join(output_dir, f'Fig3_causal_effects_revised.{fmt}'),
                    bbox_inches='tight', dpi=600, facecolor='white', **kw)
    print(f'\n\u2705 Saved: Fig3_causal_effects_revised ({FIG_W_MM}x{FIG_H_MM}mm)')
plt.close()
