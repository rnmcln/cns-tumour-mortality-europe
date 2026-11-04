"""Step 2c (needs internet): cancer registration data for the registry comparisons (steps 17e and 17f).

England: Office for National Statistics, Cancer registration statistics, England, annual reference tables for
diagnosis years 2009-2017, and ONS mid-year population estimates, detailed time series 2001-2020 (MYEB1).
Belgium: Belgian Cancer Registry, Primary brain and other central nervous system tumours in adults in Belgium
2004-2020 (Brussels, 2023), full report (PDF).

Files are written to data/raw/registry_england/ and data/raw/registry_belgium/ under the names used by the
pipeline and checked against the SHA-256 checksums of the files used for the published analysis
(data/SHA256SUMS_raw_inputs.txt). A mismatch is reported, not overwritten silently: a provider may have
replaced a file since the analysis was run.
"""
import hashlib
import urllib.request

from config import RAW, ROOT

ONS = ('https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/'
       'datasets/cancerregistrationstatisticscancerregistrationstatisticsengland/')
UA = 'cns-tumour-mortality-europe/1.0 (research data download)'
FILES = {
    'registry_england/ons_2009.xls': ONS + '2009/cancerarv200_tcm77-230877.xls',
    'registry_england/ons_2010.xls': ONS + '2010/cancerarv2010_tcm77-267683.xls',
    'registry_england/ons_2011.xls': ONS + '2011/cancerstatisticsregistrations2011metadataandreferencetables_tcm77-315737.xls',
    'registry_england/ons_2012.xls': ONS + '2012/copyofcancerstatisticsregistrations2012finalxlsversion1_tcm77-367713.xls',
    'registry_england/ons_2013.xls': ONS + '2013/referencetablearv2013c44_tcm77-410534.xls',
    'registry_england/ons_2014.xls': ONS + '2014/cancerregistrationstatistics20162014datareferencetablesfinalinc.table13.xls',
    'registry_england/ons_2015.xls': ONS + '2015/cancerregistrations2015final22.05.2017.xls',
    'registry_england/ons_2016.xls': ONS + '2016/cancerregistrations2016final.xls',
    'registry_england/ons_2017.xlsx': ONS + '2017/finaltables14.xlsx',
    'registry_england/ons_pop_detailed_2001_2020.zip': (
        'https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/'
        'datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland/mid2001tomid2020detailedtimeseries/'
        'ukdetailedtimeseries2001to2020.zip'),
    'registry_belgium/bcr_cns_adults_2004_2020.pdf': 'https://kankerregister.org/sites/default/files/2024/bcrpub_neuroonco_en_0.pdf',
}


def expected():
    out = {}
    for line in (ROOT / 'data' / 'SHA256SUMS_raw_inputs.txt').read_text().splitlines():
        h, p = line.split(maxsplit=1)
        out[p.strip().removeprefix('data/raw/')] = h
    return out


def main():
    exp = expected()
    for rel, url in FILES.items():
        dest = RAW / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={'User-Agent': UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            dest.write_bytes(r.read())
        h = hashlib.sha256(dest.read_bytes()).hexdigest()
        status = 'OK' if exp.get(rel) == h else 'CHECKSUM DIFFERS FROM THE FILE USED IN THE ANALYSIS'
        print(f'{rel}: {status}')


if __name__ == '__main__':
    main()
