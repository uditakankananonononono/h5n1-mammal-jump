#!/usr/bin/env python3
"""PB2 mammal-vs-avian enrichment analysis per locked pb2/GATES.md."""
import json, re, hashlib, math
from collections import defaultdict, Counter
from scipy.stats import fisher_exact

START, END = '2021-01-01', '2026-09-23'
MINLEN = math.ceil(759*0.95)  # 722
MAMMALIA, AVES = '40674', '8782'

def norm_date(d):
    mons={m:i+1 for i,m in enumerate(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])}
    m=re.match(r'(\d{1,2})-(\w{3})-(\d{4})', d)
    if m: return f"{m.group(3)}-{mons[m.group(2)]:02d}-{int(m.group(1)):02d}"
    m=re.match(r'(\w{3})-(\d{4})', d)
    if m: return f"{m.group(2)}-{mons[m.group(1)]:02d}-15"
    m=re.match(r'(\d{4})$', d)
    if m: return f"{m.group(1)}-06-30"
    return ''

seqs={}
hdr=None; seq=[]
def flush():
    if hdr is None: return
    m=re.match(r'>lcl\|(\S+?)_prot_(\S+?)_\d+ \[gene=PB2\].*\[location=([^]]*)\]', hdr)
    if not m: return
    seqs[m.group(1)]=''.join(seq).replace('*','')
for line in open('data/pb2_cds_aa.fna', errors='replace'):
    if line.startswith('>'):
        flush(); hdr=line.rstrip(); seq=[]
    else: seq.append(line.strip())
flush()
print('PB2 protein records:', len(seqs))

meta=json.load(open('results/gb_meta.json'))
linmap=json.load(open('results/taxlineage.json'))

def classify(tid):
    l=linmap.get(tid,{})
    lins=set(l.get('lineage_ids',[]))
    if AVES in lins: return 'avian','avian'
    if MAMMALIA in lins:
        if tid=='9606': return 'mammal','human'
        if tid=='9913': return 'mammal','cattle'
        if '9681' in lins: return 'mammal','felids'
        if '9655' in lins: return 'mammal','mustelids'
        if lins & {'9706','34893','9708'}: return 'mammal','pinnipeds'
        return 'mammal','other_mammals'
    return None,None

kept=[]; drops=Counter()
for acc,prot in seqs.items():
    m=meta.get(acc)
    if not m: drops['no_meta']+=1; continue
    d=norm_date(m.get('collection_date',''))
    if not d or d<START or d>END: drops['date']+=1; continue
    tid=m.get('taxon','')
    cls,grp=classify(tid)
    if not cls: drops['host_unknown_or_excluded']+=1; continue
    if tid=='10090': drops['lab_mouse']+=1; continue
    if len(prot)<MINLEN: drops['short']+=1; continue
    kept.append({'acc':acc,'cls':cls,'grp':grp,'date':d,'host':m.get('host',''),'taxon':tid,'prot':prot})
print('kept:', len(kept), dict(drops))
print('by class:', Counter(k['cls'] for k in kept))
print('mammal groups:', Counter(k['grp'] for k in kept if k['cls']=='mammal'))

byseq=defaultdict(lambda:{'n_m':0,'n_a':0,'grps':Counter(),'accs':[]})
for k in kept:
    e=byseq[k['prot']]
    if k['cls']=='mammal': e['n_m']+=1
    else: e['n_a']+=1
    e['grps'][k['grp']]+=1
    if len(e['accs'])<10: e['accs'].append(k['acc'])
uniq=list(byseq.items())
Nm=sum(1 for k in kept if k['cls']=='mammal'); Na=len(kept)-Nm
print('unique proteins:', len(uniq), 'isolates: mammal',Nm,'avian',Na)

L=max(len(p) for p,_ in uniq)
ref=[]
for i in range(L):
    c=Counter()
    for p,e in uniq:
        if i<len(p): c[p[i]]+=e['n_a']+e['n_m']
    ref.append(c.most_common(1)[0][0])
refseq=''.join(ref)

def run_tests(mask_groups=None):
    nm=na=0
    per=defaultdict(lambda:[Counter(),Counter()])
    for p,e in uniq:
        if mask_groups and e['grps']:
            dom=e['grps'].most_common(1)[0][0]
            if dom in mask_groups:
                e2m=0; e2a=e['n_a']
            else:
                e2m=e['n_m']; e2a=e['n_a']
        else:
            e2m=e['n_m']; e2a=e['n_a']
        nm+=e2m; na+=e2a
        for i,a in enumerate(p[:L]):
            if e2m: per[i][0][a]+=e2m
            if e2a: per[i][1][a]+=e2a
    tests=[]
    for i in range(L):
        mam,avi=per[i]
        r=ref[i]
        for a in (set(mam)|set(avi))-{r}:
            m1=mam[a]; m0=nm-m1; a1=avi[a]; a0=na-a1
            if m1==0: continue
            _,pval=fisher_exact([[m1,m0],[a1,a0]])
            orh=((m1+0.5)*(a0+0.5))/((m0+0.5)*(a1+0.5))
            tests.append({'pos':i+1,'ref':r,'alt':a,'m1':m1,'a1':a1,'p':pval,'or':orh,
                          'mf':m1/nm if nm else 0,'af':a1/na if na else 0})
    tests.sort(key=lambda t:t['p'])
    M=len(tests)
    for rank,t in enumerate(tests,1): t['q']=t['p']*M/rank
    prev=1
    for t in reversed(tests):
        t['q']=min(t['q'],prev); prev=t['q']
    return tests,nm,na

tests,nm,na=run_tests()
sig=[t for t in tests if t['q']<0.05 and t['mf']>t['af']]
print('tests:',len(tests),'sig(mammal-enriched):',len(sig))
for t in sig[:40]:
    print(t['ref']+str(t['pos'])+t['alt'],'m',t['m1'],'a',t['a1'],'q',f"{t['q']:.1e}",'OR',f"{t['or']:.1f}",'mf',f"{t['mf']:.3f}",'af',f"{t['af']:.4f}")

groups=['human','cattle','felids','mustelids','pinnipeds','other_mammals']
jack={}
for g in groups:
    jt,jnm,jna=run_tests(mask_groups={g})
    jack[g]={(t['pos'],t['alt']):t['q'] for t in jt}
    print('jackknife without',g,'nm',jnm,'tests',len(jt))

json.dump({'Nm':Nm,'Na':Na,'kept':len(kept),'unique':len(uniq),'ref':refseq,
           'drops':dict(drops),'tests':tests,
           'jack':{g:{f'{p}|{a}':q for (p,a),q in d.items()} for g,d in jack.items()}},
          open('results/analysis_raw.json','w'))
print('saved results/analysis_raw.json')
