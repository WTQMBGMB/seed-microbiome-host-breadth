"""Download published Zenodo data, verify SHA-256, and extract only data paths."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,shutil,stat,sys,urllib.request,urllib.parse,zipfile,zlib
ROOT=Path(__file__).resolve().parents[2]
ALLOWED={'data','server_result','results','manuscript_tables'}
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def extract(archive,root):
    root=root.resolve()
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            rel=PurePosixPath(info.filename)
            mode=info.external_attr>>16
            if rel.is_absolute() or '..' in rel.parts or not rel.parts or rel.parts[0] not in ALLOWED or '\\' in info.filename or stat.S_ISLNK(mode):
                raise ValueError('Unsafe archive member: '+info.filename)
            p=(root/Path(*rel.parts)).resolve();p.relative_to(root)
            if info.is_dir():p.mkdir(parents=True,exist_ok=True);continue
            if p.exists():
                crc=0
                with p.open('rb') as f:
                    for b in iter(lambda:f.read(1024*1024),b''):crc=zlib.crc32(b,crc)
                if p.stat().st_size==info.file_size and crc==info.CRC:continue
                raise FileExistsError(f'{p} differs from this release; use a clean project copy')
            p.parent.mkdir(parents=True,exist_ok=True)
            with z.open(info) as src,p.open('wb') as dst:shutil.copyfileobj(src,dst,1024*1024)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record');parser.add_argument('--archives',type=Path,help='Use local archives instead of network');parser.add_argument('--root',type=Path,default=ROOT);a=parser.parse_args()
    cfg=json.loads((ROOT/'release-assets.json').read_text());record=a.record or cfg.get('zenodo_record');download=a.archives or ROOT/'.local_state/downloads';download.mkdir(parents=True,exist_ok=True)
    remote={}
    if not a.archives:
        if not record:parser.error('A published Zenodo record is required: pass --record https://zenodo.org/records/ID')
        u=urllib.parse.urlparse(str(record))
        if u.scheme!='https' or u.netloc!='zenodo.org':parser.error('Use a public https://zenodo.org/records/ID URL')
        rid=u.path.rstrip('/').split('/')[-1]
        if not rid.isdigit():parser.error('Record ID must be numeric')
        with urllib.request.urlopen('https://zenodo.org/api/records/'+rid,timeout=90) as response:metadata=json.load(response)
        remote={x['key']:x for x in metadata['files']}
    for item in cfg['archives']:
        path=download/item['name']
        if not path.exists():
            if a.archives:raise FileNotFoundError(path)
            file=remote.get(item['name'])
            if not file:raise FileNotFoundError(item['name']+' not present in record')
            url=file['links'].get('self') or file['links'].get('content')
            parsed=urllib.parse.urlparse(url)
            if parsed.scheme!='https' or parsed.hostname!='zenodo.org':raise ValueError('Unexpected download origin')
            temp=path.with_suffix('.partial')
            with urllib.request.urlopen(url,timeout=180) as src,temp.open('wb') as dst:shutil.copyfileobj(src,dst,1024*1024)
            temp.replace(path)
        if path.stat().st_size!=item['bytes'] or digest(path)!=item['sha256']:raise ValueError('Archive checksum mismatch: '+item['name'])
        print('Verified:',item['name'],flush=True);extract(path,a.root)
    print('Download and extraction complete. Run python scripts/tools/verify_assets.py.')
if __name__=='__main__':main()
