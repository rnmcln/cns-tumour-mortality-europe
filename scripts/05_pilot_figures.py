"""Step 5 (pilot): figures. Low-ink style, monospace type, Nord greys with one accent colour."""
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from config import PILOT  # noqa: E402

T = PILOT / 'tables'
F = PILOT / 'figures'
F.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
INK, MID, LIGHT, RULE = '#2E3440', '#4C566A', '#AEB4BE', '#D8DEE9'
ACCENT, ACCENT2, WARM = '#5E81AC', '#88C0D0', '#BF616A'
plt.rcParams.update({
    'font.family': 'monospace', 'font.monospace': ['DejaVu Sans Mono', 'Liberation Mono'],
    'font.size': 6.5, 'axes.labelsize': 6.5, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'axes.edgecolor': MID, 'axes.labelcolor': INK, 'xtick.color': MID, 'ytick.color': MID,
    'text.color': INK, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'lines.linewidth': 0.8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 300, 'svg.hashsalt': 'p3-pilot'})


def save(fig, name):
    for ext, meta in (('pdf', {'CreationDate': None, 'ModDate': None}), ('png', {})):
        fig.savefig(F / f'{name}.{ext}', bbox_inches='tight', facecolor='white', pad_inches=0.03, metadata=meta)
    plt.close(fig)


def label(ax, s):
    ax.text(-0.12, 1.04, s, transform=ax.transAxes, fontsize=7.5, fontweight='bold', va='bottom')


# ---------------------------------------------------------------- Figure P1: trends
tr = pd.read_csv(T / 'trends_eurostat.csv')
eu = [u for u in tr.unit.unique() if u.startswith('Europe')][0]
s = tr[(tr.unit == eu) & (tr.sex == 'both')].set_index('year')
ag = pd.read_csv(T / 'age_specific_rates_eu.csv')
fig, axs = plt.subplots(1, 2, figsize=(170 * MM, 60 * MM), gridspec_kw={'wspace': 0.35})
ax = axs[0]
for col, c, lab in (('deaths', INK, 'Deaths'), ('asmr', ACCENT, 'Age-standardised rate')):
    v = s[col] / s.loc[2011, col] * 100
    ax.plot(v.index, v.values, color=c, marker='o', ms=2)
    ax.text(2023.3, v.loc[2023], f'{lab}\n{v.loc[2023]:.0f}', color=c, va='center', fontsize=6)
ax.axhline(100, color=RULE, lw=0.5, zorder=0)
ax.set_xlim(2010.5, 2029)
ax.set_xticks([2011, 2015, 2019, 2023])
ax.set_ylabel('Index, 2011 = 100')
label(ax, 'a')
ax = axs[1]
x = np.arange(len(ag))
ax.axhline(0, color=RULE, lw=0.5, zorder=0)
ax.vlines(x, 0, ag['pct_change'], color=LIGHT, lw=1.6)
ax.plot(x, ag['pct_change'], 'o', color=ACCENT, ms=2.5)
ax.set_xticks(x[::2])
ax.set_xticklabels([f'{a}' for a in ag.age[::2]])
ax.set_xlabel('Age group (lower bound, years)')
ax.set_ylabel('Change in rate, 2011-13 to 2021-23 (%)')
label(ax, 'b')
save(fig, 'figP1_trends')

# ---------------------------------------------------- Figure P2: decomposition by country
dec = pd.read_csv(T / 'decomposition_2011_13_to_2021_23.csv')
d = dec[dec.measure == 'deaths'].copy()
for c in ('population_size', 'age_structure', 'age_specific_rates', 'change'):
    d[c + '_pct'] = d[c] / d.base_2011_2013 * 100
d = d.sort_values('change_pct')
d = pd.concat([d[d.unit != eu], d[d.unit == eu]])
fig, ax = plt.subplots(figsize=(100 * MM, 120 * MM))
y = np.arange(len(d))
for col, c, lab, off in (('population_size_pct', LIGHT, 'Population size', -0.22),
                         ('age_structure_pct', ACCENT2, 'Age structure', 0.0),
                         ('age_specific_rates_pct', WARM, 'Age-specific rates', 0.22)):
    ax.barh(y + off, d[col], height=0.22, color=c, label=lab)
ax.plot(d.change_pct, y, 'D', color=INK, ms=2.2, label='Net change')
ax.axvline(0, color=MID, lw=0.5)
ax.set_yticks(y)
ax.set_yticklabels(d.unit)
for t in ax.get_yticklabels():
    if t.get_text() == eu:
        t.set_fontweight('bold')
ax.set_xlabel('Contribution to change in deaths, 2011-13 to 2021-23 (% of baseline)')
ax.legend(frameon=False, loc='lower right', fontsize=5.5)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
save(fig, 'figP2_decomposition')

# ---------------------------------------------------------- Figure P3: projections
pr = pd.read_csv(T / 'projections_eu.csv')
fig, axs = plt.subplots(1, 2, figsize=(170 * MM, 62 * MM), gridspec_kw={'wspace': 0.32})
for ax, m, ylab, scale in ((axs[0], 'deaths', 'Deaths (thousands)', 1e3), (axs[1], 'yll', 'YLL (thousands)', 1e3)):
    obs = pr[pr.scenario == 'observed']
    ax.plot(obs.year, obs[m] / scale, color=INK, marker='o', ms=1.8, label='Observed')
    for sc, c in (('A: constant rates 2017-2019', LIGHT), ('A: constant rates 2021-2023', ACCENT),
                  ('B: attenuated trend fitted 2011-2023', WARM)):
        b = pr[(pr.scenario == sc) & (pr.variant == 'BSL')]
        ax.plot(b.year, b[m] / scale, color=c, label=sc.replace('A: ', 'A, ').replace('B: ', 'B, '))
        if sc.startswith('A: constant rates 2021'):
            rng_ = pr[pr.scenario == sc].groupby('year')[m].agg(['min', 'max'])
            ax.fill_between(rng_.index, rng_['min'] / scale, rng_['max'] / scale, color=ACCENT, alpha=0.15, lw=0)
    ax.set_ylabel(ylab)
    ax.set_xticks([2011, 2020, 2030, 2040])
label(axs[0], 'a')
label(axs[1], 'b')
axs[0].legend(frameon=False, fontsize=5.3, loc='upper left')
save(fig, 'figP3_projections')

# ------------------------------------------------------ Figure P4: certification audit
r = pd.read_csv(T / 'audit_benign_uncertain_ratio_who.csv')
pv = r.pivot(index='geo', columns='period', values='ratio_D_to_C')
pv = pv.sort_values('2013-2019')
fig, ax = plt.subplots(figsize=(85 * MM, 110 * MM))
y = np.arange(len(pv))
for col, c, mk in (('2005-2012', LIGHT, 'o'), ('2013-2019', ACCENT, 'o'), ('2020-2024', WARM, 's')):
    ax.plot(pv[col], y, mk, color=c, ms=2.6, label=col)
ax.set_yticks(y)
ax.set_yticklabels(pv.index)
ax.set_xlabel('Deaths D32-D33, D42-D43 per death C70-C72')
ax.legend(frameon=False, fontsize=5.5, loc='lower right')
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
for yy in y:
    ax.axhline(yy, color=RULE, lw=0.3, zorder=0)
save(fig, 'figP4_certification')
print('figures written')
