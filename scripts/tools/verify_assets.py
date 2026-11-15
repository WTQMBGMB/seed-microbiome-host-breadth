"""Verify the released processed files without needing the original workstation."""
from pathlib import Path
import argparse,csv,hashlib,sys
ROOT=Path(__file__).resolve().parents[2]
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=ROOT);a=p.parse_args()
    rows=list(csv.DictReader((ROOT/'manifests/release_files.tsv').open(encoding='utf-8'),delimiter='\t'))
    errors=[]
    for row in rows:
        f=a.root/row['path']
        if not f.is_file():errors.append('missing: '+row['path']);continue
        if f.stat().st_size!=int(row['bytes']) or digest(f)!=row['sha256']:errors.append('mismatch: '+row['path'])
    if errors:
        print('\n'.join(errors));return 1
    print(f'ASSET CHECK PASS: {len(rows)} processed input/result files');return 0
if __name__=='__main__':sys.exit(main())
