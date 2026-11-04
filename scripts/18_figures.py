"""Step 18: figures.

Design: one style for all figures. Inter (sans-serif), 7 pt text, thin axes without boxes or gridlines, greyscale
with one muted accent colour, direct labelling instead of legends, point estimates with 95% intervals where
available. Encoding used throughout:
  combined codes (C70-C72, D32-D33, D42-D43)        black, solid
  malignant codes (C70-C72)                           grey, solid
  benign/uncertain codes, D43 in particular          accent (muted rust)
Lancet number style (midline decimal point, minus sign) in all text and tick labels.

Main figures
  Figure 1  components of the change in malignant deaths (Europe); age-standardised rates by definition
  Figure 2  country ranks on both definitions, 2020-22; change in age-standardised rates by country
  Figure 3  the behaviour-code component: countries, subcodes, registry comparisons (England, Belgium)
  Figure 4  observed and projected deaths
Supplementary figures S1-S7 (ratios, back-test errors, observed/expected 2020-23, D43 in four countries,
regions, chained components, components by country).
"""
import matplotlib
import matplotlib.patheffects
import matplotlib.ticker as mt
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from config import OUTPUT, PROCESSED, NAMES  # noqa: E402

T = OUTPUT / 'tables'
F = OUTPUT / 'figures'
F.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
BLACK, DARK, MID, LIGHT, FAINT = '#1F1F1F', '#555555', '#8F8F8F', '#C8C8C8', '#ECECEC'
ACCENT, ACCENT_LIGHT = '#B0573A', '#E8CBBF'
MINUS, DOT = '−', '·'
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Inter', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'axes.edgecolor': DARK, 'axes.labelcolor': BLACK, 'xtick.color': DARK, 'ytick.color': DARK,
    'text.color': BLACK, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 0, 'ytick.minor.size': 0,
    'lines.linewidth': 0.75, 'lines.markersize': 3, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': False, 'legend.frameon': False, 'axes.unicode_minus': True,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'paper3', 'savefig.dpi': 600})
OUT = {'COMP': ('Combined codes', BLACK), 'MAL': ('Malignant (C70–C72)', MID),
       'BU': ('Benign/uncertain (D32–D33, D42–D43)', ACCENT)}


def lancet(s):
    return str(s).replace('-', MINUS).replace('.', DOT)


def num(x, d=0, sign=False):
    s = f'{x:+,.{d}f}' if sign else f'{x:,.{d}f}'
    return s.replace('-', MINUS).replace('.', DOT)


def tick_fmt(ax, axis='y', d=None):
    def f(v, _):
        if d is None:
            s = f'{v:g}'
        else:
            s = f'{v:.{d}f}'
        return s.replace('-', MINUS).replace('.', DOT)
    (ax.yaxis if axis == 'y' else ax.xaxis).set_major_formatter(mt.FuncFormatter(f))


def years_axis(ax, ticks):
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f'{int(v)}'))


def save(fig, name):
    for ext, meta in (('pdf', {'CreationDate': None, 'ModDate': None}), ('svg', {'Date': None}), ('png', {})):
        fig.savefig(F / f'{name}.{ext}', bbox_inches='tight', facecolor='white', pad_inches=0.04, metadata=meta)
    plt.close(fig)


def panel_label(ax, s, x=-0.12, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, va='bottom', ha='left')


def spread(ys, gap):
    """Spread label positions (data units) so that neighbours are at least `gap` apart, keeping their order."""
    order = np.argsort(ys)
    v = np.array(ys, dtype=float)[order]
    for _ in range(200):
        moved = False
        for i in range(1, len(v)):
            if v[i] - v[i - 1] < gap:
                m = (v[i] + v[i - 1]) / 2
                v[i - 1], v[i] = m - gap / 2, m + gap / 2
                moved = True
        if not moved:
            break
    out = np.empty_like(v)
    out[order] = v
    return out


def end_labels(ax, x, items, gap, log=False, dx=0.4, fs=6.5):
    """Direct labels at the right end of lines: items = [(y, text, colour)] or [(y, text, colour, x_end)].
    Labels are spread vertically to avoid overlap; lines that end before x get a thin leader to their label."""
    ys = np.array([np.log10(it[0]) if log else it[0] for it in items])
    pos = spread(ys, gap)
    for it, p in zip(items, pos):
        y, text, col = it[:3]
        ly = 10 ** p if log else p
        if len(it) > 3 and it[3] < x:
            ax.plot([it[3] + 0.3, x + dx * 0.8], [y, ly], color=LIGHT, lw=0.4, clip_on=False)
        ax.text(x + dx, ly, text, color=col, fontsize=fs, va='center', ha='left')


# --------------------------------------------------------------------------------------------- figure 1
def fig1():
    """A: components of the change in malignant deaths, Europe, with bootstrap 95% intervals;
    B: age-standardised rates by definition, 2003-2022."""
    d = pd.read_csv(T / 'decomp_primary.csv').set_index(['geo', 'outcome'])
    b = pd.read_csv(T / 'decomp_primary_bootstrap.csv').set_index(['geo', 'outcome', 'component'])
    m, c, cod = d.loc[('Europe', 'MAL')], d.loc[('Europe', 'COMP')], d.loc[('Europe', 'CODING')]
    rows = [('Population size', m.population_size, b.loc[('Europe', 'MAL', 'population_size')], DARK),
            ('Age–sex structure', m.age_sex_structure, b.loc[('Europe', 'MAL', 'age_sex_structure')], DARK),
            ('Combined-code mortality', c.rates, b.loc[('Europe', 'COMP', 'rates')], DARK),
            ('Behaviour-code component', cod.rates, b.loc[('Europe', 'CODING', 'rates')], ACCENT),
            ('Total change', m.deaths_1 - m.deaths_0, b.loc[('Europe', 'MAL', 'change')], BLACK)]
    assert abs(sum(r[1] for r in rows[:4]) - rows[4][1]) < 1
    fig = plt.figure(figsize=(170 * MM, 58 * MM))
    ax = fig.add_axes([0.17, 0.16, 0.36, 0.76])
    y = np.array([4, 3, 2, 1, -0.2])
    for yy, (lab, v, ci, col) in zip(y, rows):
        ax.hlines(yy, ci.boot_lo, ci.boot_hi, color=col, lw=0.75)
        ax.plot(v, yy, 'o', color=col, ms=3.2, mfc=col if lab == 'Total change' else 'white', mew=0.75)
        ax.text(max(v, ci.boot_hi) + 180, yy, num(v, 0, True), va='center', ha='left', fontsize=6.5, color=col)
    ax.axvline(0, color=LIGHT, lw=0.5, zorder=0)
    ax.axhline(0.4, color=FAINT, lw=0.5, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(-4200, 7600)
    ax.set_xticks([-4000, -2000, 0, 2000, 4000, 6000])
    tick_fmt(ax, 'x')
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: num(v, 0)))
    ax.set_xlabel('Change in malignant deaths per year, 2003–05 to 2020–22')
    panel_label(ax, 'A', -0.42)
    s = pd.read_csv(T / 'asr_series.csv')
    eu = s[(s.series == 'Set L (WHO)') & (s.unit == 'Set L')].set_index('year')
    ax = fig.add_axes([0.66, 0.16, 0.22, 0.76])
    items = []
    for o in ['COMP', 'MAL', 'BU']:
        lab, col = OUT[o]
        ax.plot(eu.index, eu[f'{o}_asr'], color=col, lw=0.9)
        items.append((eu[f'{o}_asr'].iloc[-1], {'COMP': 'Combined codes', 'MAL': 'Malignant', 'BU': 'Benign/uncertain'}[o], col))
    end_labels(ax, eu.index.max(), items, gap=0.55)
    ax.set_ylim(0, 10)
    tick_fmt(ax, 'y')
    ax.set_ylabel('Age-standardised rate per 100,000')
    ax.set_xlim(2003, eu.index.max())
    years_axis(ax, [2003, 2008, 2013, 2018, 2022])
    panel_label(ax, 'B', -0.3)
    save(fig, 'fig1_overview')


# --------------------------------------------------------------------------------------------- figure 2
def fig2():
    """A: country ranks by malignant and combined-code age-standardised rates, 2020-22; B: change in both rates."""
    r = pd.read_csv(T / 'comparison_country_ranks.csv')
    r = r[r.period == 1]
    tr = pd.read_csv(T / 'comparison_trend_direction.csv')
    n = len(r)
    fig = plt.figure(figsize=(170 * MM, 112 * MM))
    ax = fig.add_axes([0.13, 0.09, 0.20, 0.85])
    for t in r.itertuples():
        sig = abs(t.rank_diff) >= 5 and (t.rank_diff_lo > 0 or t.rank_diff_hi < 0)
        col = (ACCENT if t.rank_diff > 0 else DARK) if sig else LIGHT
        ax.plot([0, 1], [t.rank_MAL, t.rank_COMP], color=col, lw=0.9 if sig else 0.6, zorder=3 if sig else 1)
        ax.plot([0, 1], [t.rank_MAL, t.rank_COMP], 'o', color=col, ms=2, zorder=3 if sig else 1)
        tc = BLACK if sig else MID
        ax.text(-0.06, t.rank_MAL, NAMES[t.geo], ha='right', va='center', fontsize=6.3, color=tc)
        ax.text(1.06, t.rank_COMP, NAMES[t.geo], ha='left', va='center', fontsize=6.3, color=tc)
    ax.set_ylim(n + 0.6, 0.4)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Malignant', 'Combined codes'], fontsize=6.5)
    ax.xaxis.tick_top()
    ax.set_yticks([])
    for sp in ['left', 'bottom', 'top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis='x', length=0)
    panel_label(ax, 'A', -0.62, 1.0)
    ax = fig.add_axes([0.66, 0.09, 0.31, 0.85])
    tr = tr.assign(d=tr.pct_MAL - tr.pct_COMP).sort_values('d', ascending=False).reset_index(drop=True)
    y = np.arange(len(tr))
    for yy, t in zip(y, tr.itertuples()):
        ax.hlines(yy - 0.17, t.pct_MAL_lo, t.pct_MAL_hi, color=MID, lw=0.6)
        ax.hlines(yy + 0.17, t.pct_COMP_lo, t.pct_COMP_hi, color=BLACK, lw=0.6)
    ax.plot(tr.pct_MAL, y - 0.17, 'o', color=MID, mfc='white', ms=3, mew=0.7)
    ax.plot(tr.pct_COMP, y + 0.17, 's', color=BLACK, ms=2.6)
    ax.axvline(0, color=LIGHT, lw=0.5, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([NAMES[g] for g in tr.geo], fontsize=6.3)
    ax.set_ylim(len(tr) - 0.4, -1.6)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    t0 = tr.iloc[0]
    ax.text(t0.pct_MAL, -0.95, 'Malignant', color=MID, fontsize=6.3, ha='center', va='bottom')
    ax.text(t0.pct_COMP, -0.95, 'Combined codes', color=BLACK, fontsize=6.3, ha='center', va='bottom')
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: num(v, 0)))
    ax.set_xlabel('Change in age-standardised rate, 2003–05 to 2020–22 (%)')
    panel_label(ax, 'B', -0.36, 1.0)
    save(fig, 'fig2_comparisons')


# --------------------------------------------------------------------------------------------- figure 3
HALO = [matplotlib.patheffects.withStroke(linewidth=2.2, foreground='white')]


def place_codes(fig, ax, r, sz, lines):
    """Label every country with its code, set directly against its own circle. For each country, candidate
    positions around the circle (twelve directions, then the same directions slightly further out) are tried in
    order; a position is accepted if its box overlaps no other label, no other circle, no panel annotation and
    neither reference line. Countries are placed largest first; a country that cannot be placed is moved to the
    front of the order and the placement is repeated."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.column_stack([r.bu_pct, r.mal_pct]))
    pt = fig.dpi / 72
    rad = np.sqrt(np.asarray(sz) / np.pi) * pt                      # circle radius in pixels
    circles = [(x - q, y - q, x + q, y + q) for (x, y), q in zip(pts, rad)]
    fixed = [(b.x0, b.y0, b.x1, b.y1) for b in (t.get_window_extent(rend) for t in ax.texts)]
    ab = ax.get_window_extent(rend)
    lpts = []                                                       # dense points along the reference lines
    for ln in lines:
        x, y = ln.get_data()
        tt = np.linspace(0, 1, 2000)
        lpts.append(ax.transData.transform(np.column_stack([x[0] + tt * (x[-1] - x[0]), y[0] + tt * (y[-1] - y[0])])))
    lpts = np.vstack(lpts)

    def hit(b, boxes):
        return any(b[0] < c[2] and b[2] > c[0] and b[1] < c[3] and b[3] > c[1] for c in boxes)

    def on_line(b):
        m = 0.8 * pt
        return bool(((lpts[:, 0] > b[0] - m) & (lpts[:, 0] < b[2] + m) &
                     (lpts[:, 1] > b[1] - m) & (lpts[:, 1] < b[3] + m)).any())

    angles = [0, 180, 90, 270, 30, 150, 330, 210, 60, 120, 300, 240]

    def try_place(i, taken):
        others = [c for j, c in enumerate(circles) if j != i]
        for extra in (1.0, 2.5):
            q = rad[i] + extra * pt
            for a in angles:
                ca, sa = np.cos(np.radians(a)), np.sin(np.radians(a))
                ha = 'left' if ca > 0.3 else 'right' if ca < -0.3 else 'center'
                va = 'bottom' if sa > 0.3 else 'top' if sa < -0.3 else 'center'
                t = ax.annotate(r.geo.iloc[i], (r.bu_pct.iloc[i], r.mal_pct.iloc[i]), xytext=(q * ca / pt, q * sa / pt),
                                textcoords='offset points', fontsize=6.3, color=BLACK, ha=ha, va=va, zorder=4)
                e = t.get_window_extent(rend)
                b = (e.x0 - 0.5, e.y0 - 0.5, e.x1 + 0.5, e.y1 + 0.5)
                inside = b[0] >= ab.x0 and b[2] <= ab.x1 and b[1] >= ab.y0 and b[3] <= ab.y1
                if inside and not hit(b, others) and not hit(b, taken) and not on_line(b):
                    return t, b
                t.remove()
        return None

    order = list(np.argsort(-np.asarray(r.COMP0)))
    for _ in range(len(order)):
        taken, placed, failed = list(fixed), [], None
        for i in order:
            res = try_place(i, taken)
            if res is None:
                failed = i
                break
            placed.append(res[0])
            taken.append(res[1])
        if failed is None:
            return
        for t in placed:
            t.remove()
        order.remove(failed)
        order.insert(0, failed)
    raise AssertionError(f'no free label positions for all countries in figure 3A ({r.geo.iloc[failed]})')


def fig3():
    cc = pd.read_csv(T / 'coding_evidence_country_components.csv')
    bd = pd.read_csv(T / 'coding_evidence_bounds.csv')
    sc = pd.read_csv(T / 'coding_evidence_subcodes.csv')
    r = cc[cc.reported]
    fig = plt.figure(figsize=(170 * MM, 150 * MM))
    # A: countries
    ax = fig.add_axes([0.09, 0.57, 0.38, 0.39])
    sz = 4 + 70 * r.COMP0 / r.COMP0.max()
    ax.scatter(r.bu_pct, r.mal_pct, s=sz, facecolor='white', edgecolor=DARK, lw=0.6, zorder=3)
    xs = np.array([-28.0, 12.0])
    l_unit = ax.plot(xs, -xs, color=MID, lw=0.6, ls=(0, (3, 2)), zorder=1)[0]
    ax.text(-13.0, 14.6, 'Slope ' + MINUS + '1', color=MID, fontsize=6.3, ha='left', va='center', path_effects=HALO)
    sl = bd[bd.estimate.str.startswith('regression')].iloc[0].slope
    w = r.COMP0
    xm, ym = np.average(r.bu_pct, weights=w), np.average(r.mal_pct, weights=w)
    l_reg = ax.plot(xs, ym + sl * (xs - xm), color=ACCENT, lw=0.8, zorder=2)[0]
    ax.text(-27.2, -19.0, f'Weighted regression,\nslope {num(sl, 2)}', color=ACCENT,
            fontsize=6.3, ha='left', va='center', path_effects=HALO)
    place_codes(fig, ax, r, sz, [l_unit, l_reg])
    ax.axhline(0, color=FAINT, lw=0.5, zorder=0)
    ax.axvline(0, color=FAINT, lw=0.5, zorder=0)
    ax.set_xlim(-28, 12)
    ax.set_ylim(-25, 16)
    tick_fmt(ax, 'x'); tick_fmt(ax, 'y')
    ax.set_xlabel('Benign/uncertain-coded rate component (% of 2003–05 CNS tumour deaths)')
    ax.set_ylabel('Malignant rate component (%)')
    panel_label(ax, 'A', -0.2)
    # B: subcodes
    ax = fig.add_axes([0.70, 0.57, 0.27, 0.39])
    t = sc[sc.geo == 'C4 total'].set_index('outcome').rates
    items = [('C71_9', 'C71.9, brain, unspecified', DARK), ('MAL_SPEC', 'Other C70–C72', DARK),
             ('D43_2', 'D43.2, brain, unspecified', ACCENT), ('D43_01', 'D43.0–D43.1', ACCENT),
             ('D43_OTH', 'Other D43', ACCENT), ('D32', 'D32, meninges', MID), ('D33', 'D33', MID), ('D42', 'D42', MID)]
    y = np.arange(len(items))
    for yy, (k, lab, col) in zip(y, items):
        ax.hlines(yy, 0, t[k], color=col, lw=0.6)
        ax.plot(t[k], yy, 'o', color=col, ms=3.2)
        ax.text(t[k] + (90 if t[k] >= 0 else -90), yy, num(t[k], 0, True), va='center',
                ha='left' if t[k] >= 0 else 'right', fontsize=6.3, color=col)
    ax.set_yticks(y)
    ax.set_yticklabels([x[1] for x in items])
    for lab, (_, _, col) in zip(ax.get_yticklabels(), items):
        lab.set_color(col if col != MID else DARK)
    ax.set_ylim(len(items) - 0.4, -0.6)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.axvline(0, color=LIGHT, lw=0.5, zorder=0)
    ax.set_xlim(-2700, 1700)
    ax.set_xticks([-2000, -1000, 0, 1000])
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: num(v, 0)))
    ax.set_xlabel('Rate component (deaths per year)')
    panel_label(ax, 'B', -0.62)
    # C: deaths per registration, England and Belgium
    ax = fig.add_axes([0.09, 0.07, 0.30, 0.38])
    er = pd.read_csv(T / 'registry_england_ratios.csv')
    bc = pd.read_csv(T / 'registry_belgium_counts.csv')
    items = []
    for g, col, lab in [('D43', ACCENT, 'D43'), ('MAL', MID, 'Malignant')]:
        e = er[(er.group == g) & er.year.between(2009, 2017)].sort_values('year')
        ax.plot(e.year, e.mi_ratio, '-o', color=col, ms=2.4, lw=0.75)
        items.append((e.mi_ratio.iloc[-1], f'{lab}, England', col, e.year.max()))
        bg = {'D43': 'D43', 'MAL': 'MAL'}[g]
        b = bc[bc.group == bg].sort_values('year')
        ratio = b.deaths / b.registrations
        ax.plot(b.year, ratio, '--o', color=col, ms=2.4, mfc='white', lw=0.75, dashes=(3, 1.5))
        items.append((ratio.iloc[-1], f'{lab}, Belgium', col))
    ax.set_yscale('log')
    ax.set_yticks([0.5, 1, 2, 3])
    ax.yaxis.set_minor_locator(mt.NullLocator())
    tick_fmt(ax, 'y')
    ax.set_ylim(0.45, 3.6)
    ax.axhline(1, color=FAINT, lw=0.5, zorder=0)
    ax.set_xlim(2004, 2020)
    years_axis(ax, [2004, 2008, 2012, 2016, 2020])
    end_labels(ax, 2020, items, gap=0.07, log=True, dx=0.5)
    ax.set_ylabel('Deaths per registration')
    panel_label(ax, 'C', -0.2)
    # D: site unspecified, England
    ax = fig.add_axes([0.66, 0.07, 0.20, 0.38])
    e = er[(er.group == 'C71') & er.year.between(2009, 2017)].sort_values('year')
    items = []
    for c, col, lab in [('c719', DARK, 'C71.9'), ('d432', ACCENT, 'D43.2')]:
        ax.plot(e.year, 100 * e[f'{c}_share_dth'], '-o', color=col, ms=2.4, lw=0.75)
        ax.plot(e.year, 100 * e[f'{c}_share_reg'], '--o', color=col, ms=2.4, mfc='white', lw=0.75, dashes=(3, 1.5))
        items += [(100 * e[f'{c}_share_dth'].iloc[-1], f'{lab}, deaths', col),
                  (100 * e[f'{c}_share_reg'].iloc[-1], f'{lab}, registrations', col)]
    end_labels(ax, 2017, items, gap=6, dx=0.4)
    ax.set_ylim(0, 100)
    ax.set_xlim(2009, 2017)
    years_axis(ax, [2009, 2011, 2013, 2015, 2017])
    ax.set_ylabel('Coded to unspecified site (%)')
    panel_label(ax, 'D', -0.38)
    save(fig, 'fig3_nonmalignant')


# --------------------------------------------------------------------------------------------- figure 4
def fig4():
    from lib.data import load_panel
    a = pd.read_csv(T / 'projections_aggregate.csv')
    a = a[a.variant != 'EP25'].copy()
    s5 = pd.read_csv(T / 'projections_S5_latest_rates.csv')
    unst = pd.read_csv(T / 'bapc_unstable.csv')
    p = load_panel()
    cs = pd.read_csv(OUTPUT.parent / 'protocol' / 'country_sets.csv')
    L = cs[cs.set_L].geo.tolist()
    p = p[p.geo.isin(L) & p.year.between(2003, 2023) & p['pop'].notna()]
    p['COMP'] = p.MAL + p.D32 + p.D33 + p.D42 + p.D43
    obs = p.groupby('year')[['MAL', 'COMP']].sum()
    obs85 = p[p.age >= 85].groupby('year').MAL.sum()
    fig = plt.figure(figsize=(180 * MM, 64 * MM))
    axs = [fig.add_axes([0.07 + k * 0.335, 0.13, 0.15, 0.80]) for k in range(3)]
    styles = {'M2 damped log-linear': ((0, (1, 1.2)), 'Damped trend'), 'M3 age-drift': ((0, (4, 1.5, 1, 1.5)), 'Age-drift'),
              'M5 BAPC': ((0, (5, 2)), 'Bayesian APC')}
    for ax, o in [(axs[0], 'MAL'), (axs[1], 'COMP')]:
        t = a[(a.geo == 'Set L') & (a.outcome == o)]
        bsl = t[(t.model == 'M1 constant') & (t.variant == 'BSL')].set_index('year')
        rng_ = t[t.model == 'M1 constant'].groupby('year').deaths.agg(['min', 'max'])
        ax.fill_between(rng_.index, rng_['min'] / 1000, rng_['max'] / 1000, color=FAINT, lw=0)
        ax.fill_between(bsl.index, bsl.pi_lo / 1000, bsl.pi_hi / 1000, color=LIGHT, lw=0)
        ax.plot(bsl.index, bsl.deaths / 1000, color=BLACK, lw=0.9)
        items = [(bsl.deaths.iloc[-1] / 1000, 'Constant, 2015\u201319', BLACK, bsl.index.max())]
        s = s5[s5.outcome == o].groupby('year').deaths.sum()
        ax.plot(s.index, s / 1000, color=ACCENT, lw=0.8)
        items.append((s.iloc[-1] / 1000, 'Constant, latest', ACCENT, s.index.max()))
        for mname, (ls, lab) in styles.items():
            if mname == 'M5 BAPC' and (unst.outcome == o).any():
                continue
            u = t[(t.model == mname) & (t.variant == 'BSL')].set_index('year')
            if len(u):
                ax.plot(u.index, u.deaths / 1000, color=DARK, ls=ls, lw=0.7)
                items.append((u.deaths.iloc[-1] / 1000, lab, DARK, u.index.max()))
        yy = obs.index[obs.index <= 2022]
        ax.plot(yy, obs.loc[yy, o] / 1000, 'o', color=BLACK, ms=1.8, mew=0)
        lo = min(obs.loc[yy, o].min(), bsl.pi_lo.min()) / 1000
        ax.set_ylim(np.floor(lo / 5) * 5, np.ceil(rng_['max'].max() / 1000 / 5) * 5)
        ax.set_xlim(2003, 2040)
        end_labels(ax, 2040, items, gap=(ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.07, dx=1.2, fs=6.0)
        years_axis(ax, [2005, 2020, 2035])
        ax.yaxis.set_major_locator(mt.MultipleLocator(5))
        tick_fmt(ax, 'y')
        ax.set_ylabel({'MAL': 'Malignant deaths per year (thousands)', 'COMP': 'Combined-code deaths per year (thousands)'}[o])
    ax = axs[2]
    t = a[(a.geo == 'Set L') & (a.outcome == 'MAL') & (a.variant == 'BSL')]
    items = []
    for mname, col, ls, lab in [('M1 constant', BLACK, '-', 'Constant, 2015\u201319'),
                                ('M3 age-drift', DARK, (0, (4, 1.5, 1, 1.5)), 'Age-drift')]:
        u = t[t.model == mname].set_index('year')
        ax.plot(u.index, u.deaths_85 / 1000, color=col, ls=ls, lw=0.8)
        items.append((u.deaths_85.iloc[-1] / 1000, lab, col, u.index.max()))
    yy = obs85.index[obs85.index <= 2022]
    ax.plot(yy, obs85.loc[yy] / 1000, 'o', color=BLACK, ms=1.8, mew=0)
    ax.set_ylim(0, None)
    ax.set_xlim(2003, 2040)
    end_labels(ax, 2040, items, gap=ax.get_ylim()[1] * 0.07, dx=1.2, fs=6.0)
    years_axis(ax, [2005, 2020, 2035])
    ax.yaxis.set_major_locator(mt.MultipleLocator(1))
    tick_fmt(ax, 'y')
    ax.set_ylabel('Malignant deaths at 85+ years (thousands)')
    for k, ax in enumerate(axs):
        panel_label(ax, 'ABC'[k], -0.42)
    save(fig, 'fig4_projections')


# --------------------------------------------------------------------------------------------- supplementary
def small_multiples(n, width=170, height=50, wspace=0.12):
    fig, axs = plt.subplots(1, n, figsize=(width * MM, height * MM), sharey=True, gridspec_kw={'wspace': wspace})
    return fig, axs


def facet_title(ax, s):
    ax.text(0, 1.03, s, transform=ax.transAxes, fontsize=6.8, ha='left', va='bottom')


def figS1():
    bu = pd.read_csv(T / 'bu_by_country_period.csv')
    rep = pd.read_csv(OUTPUT.parent / 'protocol' / 'country_sets.csv')
    show = set(rep[rep.country_level_reporting].geo)
    bu = bu[bu.geo.isin(show)]
    piv = bu.pivot(index='geo', columns='period', values='ratio_BU_MAL').sort_values('2020-2022', ascending=False)
    fig, ax = plt.subplots(figsize=(95 * MM, 120 * MM))
    y = np.arange(len(piv))
    ax.hlines(y, piv.min(axis=1), piv.max(axis=1), color=FAINT, lw=1.6)
    marks = [('2003-2005', 'o', 'white', DARK), ('2011-2013', 'o', LIGHT, MID), ('2020-2022', 'o', BLACK, BLACK)]
    for per, mk, fc, ec in marks:
        ax.plot(piv[per], y, mk, mfc=fc, mec=ec, ms=3, mew=0.6)
    ref = piv.loc['FR']                         # well-separated values; labels placed in a header row above the plot
    for per, _, _, ec in marks:
        ax.text(ref[per], -1.1, per[:4] + '\u2013' + per[-2:], fontsize=6, ha='center', va='bottom', color=ec)
    ax.set_yticks(y)
    ax.set_yticklabels([NAMES[g] for g in piv.index])
    ax.set_ylim(len(piv) - 0.4, -1.8)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    tick_fmt(ax, 'x')
    ax.set_xlabel('Benign/uncertain-coded deaths per malignant death')
    save(fig, 'figS1_ratio')


def figS2():
    b = pd.read_csv(T / 'backtest_detail.csv')
    b = b[~b.covid_target]
    models = ['M1 constant', 'M2 damped log-linear', 'M3 age-drift', 'M5 BAPC']
    names = ['Constant\nrates', 'Damped\nlog-linear', 'Age-drift', 'Bayesian\nAPC']
    fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 55 * MM), sharey=True, gridspec_kw={'wspace': 0.12})
    for ax, o in [(axs[0], 'MAL'), (axs[1], 'COMP')]:
        for j, h in enumerate([5, 10]):
            for i, m in enumerate(models):
                v = b[(b.outcome == o) & (b.model == m) & (b.horizon == h)].ape.to_numpy()
                if not len(v):
                    continue
                x = i + (j - 0.5) * 0.3
                q = np.percentile(v, [25, 50, 75])
                col = MID if h == 5 else BLACK
                ax.vlines(x, q[0], q[2], color=col, lw=0.8)
                ax.plot(x, q[1], 'o', mfc='white' if h == 5 else BLACK, mec=col, ms=3, mew=0.7)
                if i == 1 and o == 'MAL':
                    ax.text(x, q[2] + 0.8, f'{h} years', fontsize=6, ha='center', va='bottom', color=col)
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels(names, fontsize=6.3)
        ax.tick_params(axis='x', length=0)
        facet_title(ax, {'MAL': 'Malignant deaths', 'COMP': 'Combined-code deaths'}[o])
    hi = max(np.percentile(b[(b.outcome == o) & (b.model == m) & (b.horizon == h)].ape, 75)
             for o in ['MAL', 'COMP'] for m in models for h in [5, 10]
             if len(b[(b.outcome == o) & (b.model == m) & (b.horizon == h)]))
    axs[0].set_ylim(0, np.ceil(hi * 1.12 / 10) * 10)
    axs[0].set_ylabel('Absolute percentage error (%)')
    tick_fmt(axs[0], 'y')
    save(fig, 'figS2_backtest')


def figS3():
    oe = pd.read_csv(T / 'trend_obs_exp.csv')
    oe = oe[(oe.geo == 'Set L') & (oe.year < 3000)]
    fig, ax = plt.subplots(figsize=(95 * MM, 58 * MM))
    items = []
    for k, o in enumerate(['COMP', 'MAL', 'BU']):
        t = oe[oe.outcome == o].set_index('year')
        x = t.index + (k - 1) * 0.16
        lab, col = OUT[o]
        ax.vlines(x, t.pi_lo / t.expected, t.pi_hi / t.expected, color=LIGHT if o != 'BU' else ACCENT_LIGHT, lw=2.4)
        ax.plot(x, t.ratio, 'o', color=col, ms=3, mfc=col if o != 'MAL' else 'white', mew=0.7)
        items.append((t.ratio.iloc[-1], {'COMP': 'Combined codes', 'MAL': 'Malignant', 'BU': 'Benign/uncertain'}[o], col))
    ax.axhline(1, color=LIGHT, lw=0.5, zorder=0)
    end_labels(ax, oe.year.max() + 0.16, items, gap=0.018, dx=0.15)
    years_axis(ax, [2020, 2021, 2022, 2023])
    ax.set_xlim(2019.6, 2023.4)
    tick_fmt(ax, 'y')
    ax.set_ylabel('Observed / expected deaths')
    save(fig, 'figS3_obs_exp')


def figS4():
    """Age-standardised rates of C70-C72, D43 and D32 deaths, 2003-2022, in the four countries with the
    largest behaviour-code components (ordered by the size of the component)."""
    from lib.data import load_panel, collapse, asr
    d = pd.read_csv(T / 'decomp_primary.csv')
    cod = d[(d.outcome == 'CODING') & (d.geo.str.len() == 2)].set_index('geo').rates.sort_values(ascending=False)
    top = list(cod.index[:4])
    p = load_panel()
    grid = [0, 15] + list(range(20, 90, 5))
    q = collapse(p[p.geo.isin(top) & p.year.between(2003, 2022)], grid, by=('geo', 'year'),
                 value_cols=['MAL', 'D32', 'D43', 'pop'])
    fig, axs = small_multiples(4, height=52)
    ymax = 0
    for k, (ax, gname) in enumerate(zip(axs, top)):
        t = q[q.geo == gname]
        yrs = sorted(t.year.unique())
        items = []
        for col, c, lab in [('MAL', MID, 'C70–C72'), ('D43', ACCENT, 'D43'), ('D32', DARK, 'D32')]:
            v = [asr(t[t.year == y], col, 'pop', grid) for y in yrs]
            ymax = max(ymax, max(v))
            ax.plot(yrs, v, color=c, lw=0.8, ls='-' if col != 'D32' else (0, (3, 1.5)))
            items.append((v[-1], lab, c))
        if k == len(top) - 1:
            end_labels(ax, yrs[-1], items, gap=0.45, dx=0.6, fs=6)
        facet_title(ax, NAMES[gname])
        years_axis(ax, [2005, 2013, 2021])
    axs[0].set_ylim(0, ymax * 1.05)
    tick_fmt(axs[0], 'y')
    axs[0].set_ylabel('Age-standardised rate per 100,000')
    save(fig, 'figS4_d43_countries')


def figS5():
    s = pd.read_csv(T / 'asr_series.csv')
    s = s[s.series == 'Set L (WHO)']
    regions = ['Northern Europe', 'Western Europe', 'Southern Europe', 'Eastern Europe']
    fig, axs = small_multiples(4, height=52)
    for k, (ax, rg) in enumerate(zip(axs, regions)):
        t = s[s.unit == rg].set_index('year')
        items = []
        for o in ['COMP', 'MAL', 'BU']:
            lab, col = OUT[o]
            ax.plot(t.index, t[f'{o}_asr'], color=col, lw=0.8)
            items.append((t[f'{o}_asr'].iloc[-1], {'COMP': 'Combined codes', 'MAL': 'Malignant', 'BU': 'Benign/uncertain'}[o], col))
        if k == 3:
            end_labels(ax, t.index.max(), items, gap=0.7, dx=0.6, fs=6)
        facet_title(ax, rg)
        years_axis(ax, [2005, 2013, 2021])
    axs[0].set_ylim(0, None)
    tick_fmt(axs[0], 'y')
    axs[0].set_ylabel('Age-standardised rate per 100,000')
    save(fig, 'figS5_regions')


def figS6():
    c = pd.read_csv(T / 'decomp_chained.csv')
    ch = pd.read_csv(PROCESSED / 'coding_changes.csv')
    rep = pd.read_csv(OUTPUT.parent / 'protocol' / 'country_sets.csv')
    L = set(rep[rep.set_L].geo)
    ch = ch[ch.geo.isin(L) & ch.year.between(2004, 2022)]
    e = c[c.geo == 'Europe']
    fig = plt.figure(figsize=(170 * MM, 62 * MM))
    axs = [fig.add_axes([0.09, 0.14, 0.26, 0.8]), fig.add_axes([0.645, 0.14, 0.26, 0.8])]
    for ax, o in [(axs[0], 'MAL'), (axs[1], 'COMP')]:
        t = e[e.outcome == o].set_index('year').sort_index()
        comps = [('population_size', DARK, (0, (1, 1.2)), 'Population size'),
                 ('age_sex_structure', DARK, (0, (4, 1.5)), 'Age–sex structure'),
                 ('rates', MID, '-', 'Malignant mortality' if o == 'MAL' else 'Combined-code mortality')]
        items = []
        for k, col, ls, lab in comps:
            ax.plot(t.index, t[k], color=col, ls=ls, lw=0.8)
            items.append((t[k].iloc[-1], lab, col))
        if o == 'MAL':
            cod = e[e.outcome == 'CODING'].set_index('year').sort_index().rates
            ax.plot(cod.index, cod, color=ACCENT, lw=0.8)
            items.append((cod.iloc[-1], 'Behaviour-code component', ACCENT))
        tot = t[['population_size', 'age_sex_structure', 'rates']].sum(axis=1)
        ax.plot(t.index, tot, color=BLACK, lw=1.0)
        items.append((tot.iloc[-1], 'Total change', BLACK))
        ax.axhline(0, color=LIGHT, lw=0.5, zorder=0)
        yl = ax.get_ylim()
        for _, r in ch.iterrows():
            ax.plot([r.year, r.year], [yl[0], yl[0] + 0.035 * (yl[1] - yl[0])],
                    color=BLACK if r.source == 'documented' else MID, lw=0.5)
        ax.set_ylim(yl)
        end_labels(ax, 2022, items, gap=(yl[1] - yl[0]) * 0.06, dx=0.4, fs=6.2)
        ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: num(v, 0)))
        ax.set_ylabel('Cumulative change since 2003\n(deaths per year)')
        ax.set_xlim(2003, 2022)
        years_axis(ax, [2005, 2010, 2015, 2020])
    panel_label(axs[0], 'A', -0.36, 1.04); panel_label(axs[1], 'B', -0.36, 1.04)
    save(fig, 'figS6_chained')


def figS7():
    """Components of the change in malignant deaths by country (% of 2003-05 deaths), small multiples."""
    d = pd.read_csv(T / 'decomp_primary.csv')
    b = pd.read_csv(T / 'decomp_primary_bootstrap.csv')
    rep = pd.read_csv(OUTPUT.parent / 'protocol' / 'country_sets.csv')
    show = set(rep[rep.country_level_reporting & rep.set_L].geo)
    m = d[d.outcome == 'MAL'].set_index('geo')
    c = d[d.outcome == 'CODING'].set_index('geo').rates
    comp = d[d.outcome == 'COMP'].set_index('geo').rates
    geos = [g for g in m.index if g in show] + ['Europe']
    base = m.loc[geos, 'deaths_0']
    tab = pd.DataFrame({'size': m.loc[geos, 'population_size'], 'structure': m.loc[geos, 'age_sex_structure'],
                        'residual': comp.loc[geos], 'coding': c.loc[geos]}).div(base, axis=0) * 100
    tab['total'] = tab.sum(axis=1)
    bt = b[(b.outcome == 'MAL') & (b.component == 'change')].set_index('geo')
    tab['lo'] = bt.loc[tab.index, 'boot_lo'] / base * 100
    tab['hi'] = bt.loc[tab.index, 'boot_hi'] / base * 100
    eu = tab.loc[['Europe']]
    tab = pd.concat([eu, tab.drop('Europe').sort_values('total', ascending=False)])
    cols = [('size', 'Population\nsize', DARK), ('structure', 'Age–sex\nstructure', DARK),
            ('residual', 'Combined-code\nmortality', DARK), ('coding', 'Behaviour-code\ncomponent', ACCENT), ('total', 'Total\nchange', BLACK)]
    fig, axs = plt.subplots(1, 5, figsize=(170 * MM, 120 * MM), sharey=True, gridspec_kw={'wspace': 0.18})
    y = np.arange(len(tab)).astype(float)
    y[1:] += 0.8
    for ax, (k, lab, col) in zip(axs, cols):
        if k == 'total':
            ax.hlines(y, tab.lo, tab.hi, color=col, lw=0.6)
        ax.hlines(y, 0, tab[k], color=FAINT, lw=0.6)
        ax.plot(tab[k], y, 'o', color=col, ms=2.6)
        ax.axvline(0, color=LIGHT, lw=0.5, zorder=0)
        ax.set_title(lab, fontsize=6.5, loc='left', pad=4)
        ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: num(v, 0)))
        ax.tick_params(axis='y', length=0)
        ax.spines['left'].set_visible(False)
        ax.xaxis.set_major_locator(mt.MaxNLocator(3))
    axs[0].set_yticks(y)
    axs[0].set_yticklabels(['All 26 countries' if g == 'Europe' else NAMES[g] for g in tab.index])
    axs[0].set_ylim(y.max() + 0.6, -0.6)
    fig.text(0.5, 0.04, 'Change in malignant deaths (% of 2003–05 deaths)', ha='center', fontsize=7)
    save(fig, 'figS7_countries')


if __name__ == '__main__':
    for f in (fig1, fig2, fig3, fig4, figS1, figS2, figS3, figS4, figS5, figS6, figS7):
        f()
    print('figures written to', F)
