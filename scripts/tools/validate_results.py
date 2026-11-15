"""Portable checks of archived manuscript results, independent of local paths."""
from pathlib import Path
from collections import Counter
import csv,json,sys
ROOT=Path(__file__).resolve().parents[2]
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
    checks=[]
    def check(name,actual,expected):
        checks.append((name,actual==expected));print(('PASS' if actual==expected else 'FAIL'),name,actual)
    for dataset,total in [('16S',3652),('ITS',2310)]:check('metadata '+dataset,len(rows(ROOT/'data/metadata'/f'metadata_{dataset}.txt')),total)
    for dataset,want in [('16s',{'NRH':12190,'generalist':1095,'intermediate':3813,'specialist':5747}),('its',{'NRH':2860,'generalist':234,'intermediate':522,'specialist':1302})]:
        check('host breadth '+dataset,dict(Counter(x['category'] for x in rows(ROOT/'results/tables/gos'/f'gos_classified_{dataset}.txt'))),want)
    for name,count in [('sensitivity_robust17.tsv',17),('sensitivity_discordant9.tsv',9)]:check(name,len(rows(ROOT/'server_result/function_prediction/result_tables'/name)),count)
    check('output map',len(rows(ROOT/'manifests/figure_map.tsv')),84)
    for n in range(1,8):check('Table S'+str(n),(ROOT/'results/tables'/f'Table S{n}.xlsx').is_file(),True)
    return 0 if all(ok for _,ok in checks) else 1
if __name__=='__main__':sys.exit(main())
