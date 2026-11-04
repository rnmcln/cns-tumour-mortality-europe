"""Step 6: extract coding-system and ICD-update information from the Eurostat country metadata
(hlth_cdeath_simscd_<cc>.htm, saved in data/raw/metadata). Output: data/interim/coding_metadata_blocks.md
(verbatim blocks per country, for building the table of documented coding changes)."""
import html
import re

from config import RAW, INTERIM, COUNTRIES

SRC = RAW / 'metadata'
KEYS = ['Classification and updates applied by years', '18.3.2. Automated Coding',
        '18.3.3. Underlying cause of death', '15.2. Comparability - over time']
STOP = re.compile(r'^(\d+\.\d+(\.\d+)?\.\s|Classification and updates applied by years$|Annex$|Related metadata$)')


def lines(path):
    t = path.read_text(encoding='utf-8', errors='ignore')
    t = re.sub(r'<(script|style)[^>]*>.*?</\1>', ' ', t, flags=re.S)
    s = html.unescape(re.sub('<[^>]+>', '\n', t))
    return [ln.strip() for ln in s.split('\n') if ln.strip() and 'accordion arrow' not in ln]


out = ['# Coding-system and ICD-update information from Eurostat country metadata', '',
       'Verbatim extracts (accessed: see data/raw/metadata/ACCESSED.txt). One block per country and section.', '']
for cc in COUNTRIES:
    f = SRC / f'hlth_cdeath_simscd_{cc.lower()}.htm'
    L = lines(f)
    out.append(f'## {cc}')
    for key in KEYS:
        idx = [i for i, ln in enumerate(L) if ln.startswith(key)]
        if not idx:
            continue
        i = idx[-1] if key != '15.2. Comparability - over time' else idx[-1]
        block = [L[i]]
        for ln in L[i + 1:i + 200]:
            if STOP.match(ln) and not ln.startswith(key):
                break
            block.append(ln)
        out.append(f'**{key}**: ' + ' | '.join(block[1:]))
        out.append('')
(INTERIM / 'coding_metadata_blocks.md').write_text('\n'.join(out))
print('written', len(out))


# ---------------------------------------------------------------- structured table, 2011 onwards
import pandas as pd  # noqa: E402

ANNEX = SRC / 'annexes'
YEAR = re.compile(r'^(19|20)\d{2}$')


def records_from_tokens(tokens):
    """Split a flat token list into per-year records (a 4-digit year starts a record)."""
    recs, cur = [], None
    for t in tokens:
        if YEAR.match(t):
            if cur:
                recs.append(cur)
            cur = [t]
        elif cur is not None:
            cur.append(t)
    if cur:
        recs.append(cur)
    return recs


def html_block(L, key, stop_prefix):
    idx = [i for i, ln in enumerate(L) if ln.startswith(key)]
    if not idx:
        return []
    out = []
    for ln in L[idx[-1] + 1:]:
        if ln.startswith(stop_prefix) or (re.match(r'^\d+\.\d+(\.\d+)?\.\s', ln) and not ln.startswith(key)):
            break
        out.append(ln)
    return out


def annex_table(cc, sheet_prefix):
    files = sorted(ANNEX.glob(f'hlth_cdeath_simscd_{cc.lower()}_an_*.xlsx'))
    for f in files:
        xl = pd.read_excel(f, sheet_name=None, header=None, dtype=str)
        for name, df in xl.items():
            if name.strip().startswith(sheet_prefix):
                rows = []
                for _, r in df.iterrows():
                    vals = [str(v).strip() for v in r.values if str(v).strip() not in ('', 'nan')]
                    if vals and re.match(r'^(19|20)\d{2}(\.0)?$', vals[0]):
                        vals[0] = vals[0][:4]
                        rows.append(vals)
                return rows, f.name
    return None, None


rows = []
for cc in COUNTRIES:
    L = lines(SRC / f'hlth_cdeath_simscd_{cc.lower()}.htm')
    icd, src_icd = annex_table(cc, '3.2')
    if not icd:
        icd = records_from_tokens(html_block(L, 'Classification and updates applied by years', '3.3'))
        src_icd = 'html'
    auto, src_auto = annex_table(cc, '18.3.2')
    if not auto:
        auto = records_from_tokens(html_block(L, '18.3.2. Automated Coding', '18.3.3'))
        src_auto = 'html'
    d_icd = {r[0]: ' / '.join(r[1:]) for r in icd if 2000 <= int(r[0]) <= 2025}
    d_auto = {r[0]: ' / '.join(r[1:]) for r in auto if 2000 <= int(r[0]) <= 2025}
    for y in sorted(set(d_icd) | set(d_auto)):
        rows.append({'geo': cc, 'year': int(y), 'icd_and_update': d_icd.get(y, ''),
                     'automated_coding': d_auto.get(y, ''), 'source_icd': src_icd, 'source_auto': src_auto})
tab = pd.DataFrame(rows).drop_duplicates(['geo', 'year'], keep='last').sort_values(['geo', 'year'])
tab.to_csv(INTERIM / 'coding_metadata_by_year.csv', index=False)

# candidate change points: years in which the automated-coding status or the coding-system family
# (Iris, ACME, MICAR, STYX, MIKADO, MUSE, own, manual) differs from the previous documented year;
# software version updates alone are not counted
FAMILIES = ['iris', 'acme', 'micar', 'styx', 'mikado', 'muse', 'own application', 'manual']


def family(txt):
    t = txt.lower()
    if not t.strip():
        return None
    if t.startswith('no') and not any(f in t for f in FAMILIES[:-1]):
        return 'no automated coding'
    if 'stopped using iris' in t:
        return 'manual'
    return '+'.join(f for f in FAMILIES if f in t) or 'automated, unspecified'


chg = []


def html_rows(cc):
    """Rows of table 18.3.2 in the metadata page, as lists of cell texts."""
    t = (SRC / f'hlth_cdeath_simscd_{cc.lower()}.htm').read_text(encoding='utf-8', errors='ignore')
    i = t.find('18.3.2. Automated Coding')
    seg = t[i:t.find('18.3.3', i)]
    out = []
    for r in re.findall(r'<tr[^>]*>(.*?)</tr>', seg, flags=re.S):
        cells = [html.unescape(re.sub('<[^>]+>', '', c)).strip() for c in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', r, flags=re.S)]
        out.append(cells)
    return out


# Belgium reports two regions in separate columns (Flanders and Brussels; Wallonia): change points per region
for region, (use_col, sys_col) in (('Flanders and Brussels', (1, 3)), ('Wallonia', (2, 4))):
    prev = None
    for r in html_rows('BE'):
        if len(r) == 5 and YEAR.match(r[0]):
            use, system = r[use_col].lower(), r[sys_col].lower()
            cur = 'no automated coding' if use.startswith('no') else (family(system) if system else None)
            if cur and prev and cur != prev:
                chg.append({'geo': 'BE', 'year': int(r[0]), 'from': f'{region}: {prev}', 'to': f'{region}: {cur}'})
            if cur:
                prev = cur
for g, s in tab.groupby('geo'):
    if g == 'BE':
        continue
    prev = None
    for _, r in s.iterrows():
        cur = family(r.automated_coding)
        if cur and prev and cur != prev:
            chg.append({'geo': g, 'year': r.year, 'from': prev, 'to': cur})
        if cur:
            prev = cur
pd.DataFrame(chg).to_csv(INTERIM / 'coding_change_points.csv', index=False)
print(tab.groupby('geo').year.agg(['min', 'max', 'count']).to_string())
print(pd.DataFrame(chg).to_string())
