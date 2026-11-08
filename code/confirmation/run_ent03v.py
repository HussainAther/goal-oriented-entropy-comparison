"""ENT-03V: fixed 20-micromol N, 1% entropy-warning confirmation on unused posterior rows.
New draws excluded from prior ENT-03P/Q/U sample. Published Wojcik three-state ODE.
Caution: all posterior draws conditioned on all chemostats; not independent biological validation.
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score,brier_score_loss,log_loss
from model import rhs,jacobian
from derivatives import hessian_direction
R=Path(__file__).resolve().parent
old=pd.read_csv(R/'discovery_cases.csv');old=old[old.delta_N==20].copy()
assert len(old)==192
old['key']=old.group+'|'+old.chemostat.astype(str)
rng=np.random.default_rng(2026110803)
T=np.linspace(0,5,41);dn=20.
# fixed: six new draws per group, each applied to twelve chemostats
chosen={}
for grp in ['monocultures','polycultures']:
 dd=pd.read_csv(R/'posterior_inputs'/f'{grp}_1000_draws.csv')
 used=set(old.loc[old.group==grp,'draw'].astype(int))
 pool=np.array(sorted(set(range(len(dd)))-used),dtype=int)
 chosen[grp]=list(np.sort(rng.choice(pool,size=6,replace=False)))
print('Preselected untouched zero-based draw indices:',chosen,flush=True)
res=[]; fail=[]
for grp in chosen:
 draws=pd.read_csv(R/'posterior_inputs'/f'{grp}_1000_draws.csv')
 for idx in chosen[grp]:
  d=draws.iloc[idx]
  for j in range(12):
   if grp=='monocultures':
    si=j//4+1;ri=j%4+1;cs=['Bc','Ce','Le'][si-1]+str(ri)
    get=lambda s:d[f'{s}[{si}]'];ix=lambda s:d[f'{s}[{si},{ri}]']
   else:
    cs=str(j+1).zfill(3);get=lambda s:d[s];ix=lambda s:d[f'{s}[{j+1}]']
   try:
    p=dict(b=float(d.b),K=float(np.exp(d.Klog)),eA=float(np.exp(d.eAlog)),a=float(np.exp(get('alog'))),h=float(np.exp(get('hlog'))),eR=float(get('eR')),m=float(get('m')))
    y0=np.exp([ix('BTNlog'),ix('BTAlog'),ix('BTRlog')]);S=np.diag(y0);Sinv=np.diag(1/y0)
    z0=np.r_[y0,1.,0.,0.,0.,0.,0.,np.eye(3).ravel()]
    def f(t,z):
     X=z[:3];u=z[3:6];w=z[6:9];phi=z[9:].reshape(3,3);J=jacobian(X,p)
     return np.r_[rhs(t,X,p),J@u,J@w+.5*hessian_direction(X,u,p),(J@phi).ravel()]
    sol=solve_ivp(f,(0,5),z0,t_eval=T,rtol=1e-8,atol=1e-8,method='DOP853')
    if not sol.success:raise RuntimeError('reference: '+sol.message)
    X=sol.y.T[:,:3];u=sol.y.T[:,3:6];w=sol.y.T[:,6:9];phi=sol.y.T[:,9:].reshape(-1,3,3)
    P1=X+dn*u;P2=P1+dn**2*w
    if np.any(P1<=0) or np.any(P2<=0):raise ValueError('negative Taylor state')
    ypert=y0.copy();ypert[0]+=dn
    pert=solve_ivp(lambda t,y:rhs(t,y,p),(0,5),ypert,t_eval=T,rtol=1e-8,atol=1e-8,method='DOP853')
    if not pert.success:raise RuntimeError('perturbed: '+pert.message)
    actual=pert.y.T
    def D(v):return np.sum((v*np.log(v/X)-v+X)/y0,axis=1)
    truth=D(actual);d1=D(P1);d2=D(P2);valid=truth>1e-10;tv=T>0
    Derr=float(np.mean(np.abs(d1[valid]-truth[valid])/truth[valid]))
    scorr=float(np.mean(np.linalg.norm(dn**2*w[tv]/y0,axis=1)/np.maximum(np.linalg.norm(dn*u[tv]/y0,axis=1),1e-14)))
    ecorr=float(np.mean(np.abs(d2[tv]-d1[tv])/np.maximum(d1[tv],1e-12)))
    rx=[]
    for k in range(len(T)):
     Jscaled=Sinv@jacobian(X[k],p)@S
     rx.append(np.linalg.eigvalsh((Jscaled+Jscaled.T)/2)[-1])
    res.append(dict(group=grp,chemostat=cs,draw=int(idx),delta_N=dn,D_error1=Derr,label=int(Derr>.01),entropy_correction_ratio=ecorr,state_correction_ratio=scorr,reactivity_max=float(max(rx)),D_error2=float(np.mean(np.abs(d2[valid]-truth[valid])/truth[valid])),state_error1=float(np.mean(np.sqrt(np.mean(((P1-actual)/np.maximum(X,y0*.01))**2,axis=1))))))
   except Exception as ex:fail.append(dict(group=grp,chemostat=cs,draw=int(idx),error=str(ex)))
  print('Completed',grp,'draw',idx,'cases',len(res),'failures',len(fail),flush=True)
new=pd.DataFrame(res);new.to_csv(R/'confirmation_cases.csv',index=False)
pd.DataFrame(fail,columns=['group','chemostat','draw','error']).to_csv(R/'failures.csv',index=False)
pd.DataFrame([(g,x) for g,ids in chosen.items() for x in ids],columns=['group','zero_based_draw']).to_csv(R/'heldout_draw_indices.csv',index=False)
assert len(new)==144 and len(fail)==0, 'Incomplete confirmation; check failure file.'
# DISCOVERY calibrated ONCE using all original 192 configurations; no LOCO on the new data
spec={'baseline':[], 'entropy':['entropy_correction_ratio'],'state':['state_correction_ratio'],'both':['entropy_correction_ratio','state_correction_ratio'],'state_plus_reactivity':['state_correction_ratio','reactivity_max']}
metrics=[];pred=[]
for name,cols in spec.items():
 if name=='baseline':prob=np.full(len(new),float((old.D_error1>.01).mean()))
 else:
  def matrix(df):
   return np.column_stack([np.log10(np.maximum(df[c].to_numpy(float),1e-15)) if c.endswith('_ratio') else df[c].to_numpy(float) for c in cols])
  clf=make_pipeline(StandardScaler(),LogisticRegression(C=1.,max_iter=2000))
  clf.fit(matrix(old),(old.D_error1>.01).astype(int));prob=clf.predict_proba(matrix(new))[:,1]
 y=new.label.to_numpy();metrics.append(dict(method=name,train_n=len(old),test_n=len(new),test_positives=int(y.sum()),auc=float(roc_auc_score(y,prob)),brier=float(brier_score_loss(y,prob)),logloss=float(log_loss(y,np.clip(prob,1e-10,1-1e-10)))))
 for i,row in new.iterrows():pred.append(dict(group=row.group,chemostat=row.chemostat,draw=int(row.draw),label=int(row.label),method=name,probability=float(prob[i])))
pd.DataFrame(metrics).to_csv(R/'confirmation_metrics.csv',index=False)
pr=pd.DataFrame(pred);pr.to_csv(R/'confirmation_predictions.csv',index=False)
wide=pr.pivot(index=['group','chemostat','draw','label'],columns='method',values='probability').reset_index();wide['key']=wide.group+'|'+wide.chemostat
keys=wide.key.unique();blocks={k:wide[wide.key==k] for k in keys};rng2=np.random.default_rng(2026110804);diff=[]
for i in range(2000):
 boot=pd.concat([blocks[k] for k in rng2.choice(keys,size=len(keys),replace=True)],ignore_index=True)
 yy=boot.label.to_numpy()
 out={'iteration':i,'brier_entropy_minus_state':brier_score_loss(yy,boot.entropy)-brier_score_loss(yy,boot.state)}
 if len(np.unique(yy))==2:out['auc_entropy_minus_state']=roc_auc_score(yy,boot.entropy)-roc_auc_score(yy,boot.state)
 else:out['auc_entropy_minus_state']=np.nan
 diff.append(out)
pd.DataFrame(diff).to_csv(R/'cluster_bootstrap.csv',index=False)
print(pd.DataFrame(metrics).to_string(index=False))
print('paired diffs',pd.DataFrame(diff)[['brier_entropy_minus_state','auc_entropy_minus_state']].quantile([.025,.5,.975]).to_string())
