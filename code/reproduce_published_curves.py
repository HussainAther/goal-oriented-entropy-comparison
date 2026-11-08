"""Compare independently integrated posterior quantiles to authors' archived CSVs.
First run export_posterior_draws.R on macOS from study root, then:
 python reproduce_published_curves.py --archive data.zip --exports ENT03N_exports --out reproduction.csv
This is a descriptive reproduction, NOT held-out validation.
"""
import argparse, zipfile, io, pathlib
import pandas as pd
import numpy as np
from model import integrate

p=argparse.ArgumentParser();p.add_argument('--archive',default='/mnt/data/data.zip');p.add_argument('--exports',required=True);p.add_argument('--out',default='reproduction.csv');p.add_argument('--draw-limit',type=int,default=1000)
a=p.parse_args()

def group_params(row,group,cs):
    if group=='monocultures':
        s={'Bc':1,'Ce':2,'Le':3}[cs[:2]]
        rep=int(cs[2:]); key=lambda v:f'{v}[{s},{rep}]'
        q=lambda v:f'{v}[{s}]'
        return ({'b':row['b'],'K':np.exp(row['Klog']),'eA':np.exp(row['eAlog']),
                 'a':np.exp(row[q('alog')]),'h':np.exp(row[q('hlog')]),'eR':row[q('eR')],'m':row[q('m')]},key)
    i=int(cs);key=lambda v:f'{v}[{i}]'
    return ({'b':row['b'],'K':np.exp(row['Klog']),'eA':np.exp(row['eAlog']),
             'a':np.exp(row['alog']),'h':np.exp(row['hlog']),'eR':row['eR'],'m':row['m']},key)

recs=[]
with zipfile.ZipFile(a.archive) as z:
 for group in ('monocultures','polycultures'):
    fitfile=pathlib.Path(a.exports)/f'{group}_1000_draws.csv'
    post=pd.read_csv(fitfile).head(a.draw_limit)
    pred=pd.read_csv(io.BytesIO(z.read(f'data/bayesian_predictions_{group}.csv')))
    for cs,dd in pred.groupby('Chemostat',sort=False):
      cs=str(cs)
      dd=dd.iloc[:200].reset_index(drop=True)  # source archive repeats identical 200-point blocks five times
      for seg in ('pre','post'):
        # pre integrates from first reported point through exactly zero; CSV repeats zero at each segment.
        # ensure aligned with original 100-point grid: -17..0 / -13..0 (pre) and 0..19 / 0..12 (post)
        t0=(-17 if group=='monocultures' else -13) if seg=='pre' else 0
        t1=(19 if group=='monocultures' else 12) if seg=='post' else 0
        ts=np.linspace(t0,t1,100)
        samples=[];failed=0
        for _,row in post.iterrows():
          pars,key=group_params(row,group,cs)
          pref='B0' if seg=='pre' else 'BT'
          try:
            y0=np.exp([row[key(pref+v+'log')] for v in ['N','A','R']])
            samples.append(integrate(y0,pars,ts))
          except (KeyError,ValueError,RuntimeError,OverflowError,FloatingPointError): failed+=1
        if failed: print('WARNING failed draws',group,cs,seg,failed)
        if not samples: continue
        arr=np.asarray(samples)
        for j,t in enumerate(ts):
          # source CSV contains pre and post with repeated t=0. Match by segment and position,
          # not by time alone. At pre t=0, archived rows include a second t=0 post row.
          for idx,lab in ((0,'N'),(1,'A'),(2,'R')):
            calc=np.quantile(arr[:,j,idx],[.05,.5,.95])
            source_col={'N':'N','A':'Biovolume','R':'Biovolume'}[lab]
            suffix='' if lab=='N' else ('_A' if lab=='A' else '_R')
            cols=[f'{source_col}.Q{k}{suffix}' for k in [1,2,3]]
            # Original pre/post blocks are length 100 each; identify via array order.
            rows=dd.iloc[0:100] if seg=='pre' else dd.iloc[100:200]
            if j>=len(rows): continue
            ref=rows.iloc[j][cols].to_numpy(dtype=float)
            for q,(x,y) in enumerate(zip(calc,ref)):
              recs.append({'group':group,'chemostat':cs,'segment':seg,'time':t,'component':lab,
                           'quantile':[.05,.5,.95][q],'calculated':x,'published':y,
                           'relative_discrepancy':abs(x-y)/max(1.0,abs(y))})
result=pd.DataFrame(recs);result.to_csv(a.out,index=False)
if len(result):
 print(result.groupby(['group','component']).relative_discrepancy.agg(['count','median','max']).to_string())
 print('overall median relative discrepancy',result.relative_discrepancy.median())
else:print('No comparisons were produced')
