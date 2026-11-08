"""ENT-04B fixed-threshold posterior-draw confirmation, matched-information diagnostics.
Source: ENT-03V discovery and confirmation indices. All predictors from X,u,w only.
"""
from pathlib import Path
import sys, time
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,brier_score_loss,log_loss
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'dependencies'))
from model import rhs,jacobian
from derivatives import hessian_direction
BASE=P/'dependencies'
T=np.linspace(0.,5.,41);eps=20.
old=pd.read_csv(BASE/'discovery_cases.csv');old=old[old.delta_N==20].copy()
new=pd.read_csv(BASE/'confirmation_cases.csv').copy()
assert len(old)==192 and len(new)==144

def compute_one(group,cs,idx,draws):
 d=draws.iloc[int(idx)]
 if group=='monocultures':
  si=['Bc','Ce','Le'].index(cs[:2])+1;ri=int(cs[-1]);get=lambda s:d[f'{s}[{si}]'];ix=lambda s:d[f'{s}[{si},{ri}]']
 else:
  ri=int(cs);get=lambda s:d[s];ix=lambda s:d[f'{s}[{ri}]']
 pars=dict(b=float(d.b),K=float(np.exp(d.Klog)),eA=float(np.exp(d.eAlog)),a=float(np.exp(get('alog'))),h=float(np.exp(get('hlog'))),eR=float(get('eR')),m=float(get('m')))
 y0=np.exp([ix('BTNlog'),ix('BTAlog'),ix('BTRlog')]);z0=np.r_[y0,1.,0.,0.,0.,0.,0.]
 def f(t,z):
  X=z[:3];u=z[3:6];w=z[6:9];J=jacobian(X,pars)
  return np.r_[rhs(t,X,pars),J@u,J@w+.5*hessian_direction(X,u,pars)]
 sol=solve_ivp(f,(0,5),z0,t_eval=T,rtol=1e-8,atol=1e-8,method='DOP853')
 if not sol.success:raise RuntimeError(sol.message)
 X=sol.y[:3].T;u=sol.y[3:6].T;w=sol.y[6:9].T
 p1=X+eps*u;h=eps**2*w;p2=p1+h
 if (X<=0).any() or (p1<=0).any() or (p2<=0).any():raise ValueError('non-positive Taylor state')
 # original reference-normalized generalised I divergence; scalar normalizer cancels feature ratios
 D=lambda y,ref:np.sum((y*np.log(y/ref)-y+ref)/y0,axis=1)
 d1=D(p1,X);d2=D(p2,X)
 S=np.sum(h*np.log(p1/X)/y0,axis=1) # first variation of goal at p1
 C=D(p2,p1) # nonnegative exact Bregman remainder
 assert np.max(np.abs((S+C)-(d2-d1))) < 2e-4, (np.max(np.abs((S+C)-(d2-d1))),np.max(np.abs(d2-d1)))
 pos=T>0; denom=np.maximum(d1[pos],1e-12)
 fisher=np.sqrt(np.sum(h[pos]**2/(p1[pos]*y0),axis=1))/np.maximum(np.sqrt(np.sum((p1[pos]-X[pos])**2/(p1[pos]*y0),axis=1)),1e-12)
 # signed GO adjoint equals signed goal gradient action of integrated second variation
 # not a separately solved continuous adjoint; mathematically equivalent endpoint sensitivity
 return dict(entropy_exact=float(np.mean(np.abs(S[pos]+C[pos])/denom)),
             goal_linear=float(np.mean(np.abs(S[pos])/denom)),
             goal_plus_curvature=float(np.mean(np.abs(S[pos]+C[pos])/denom)),
             curvature_only=float(np.mean(C[pos]/denom)),
             fisher_norm=float(np.mean(fisher)),
             fisher_alignment=float(np.mean(np.abs(np.sum((p1[pos]-X[pos])*h[pos]/(p1[pos]*y0),axis=1))/np.maximum(np.sum((p1[pos]-X[pos])**2/(p1[pos]*y0),axis=1),1e-12))),
             state_norm=float(np.mean(np.linalg.norm(h[pos]/y0,axis=1)/np.maximum(np.linalg.norm(eps*u[pos]/y0,axis=1),1e-14))),
             max_identity_error=float(np.max(np.abs((S+C)-(d2-d1)))))

def features(df,label):
 out=[];drawsets={g:pd.read_csv(BASE/f'{g}_1000_draws.csv') for g in ['monocultures','polycultures']}
 for i,r in df.iterrows():
  feat=compute_one(r.group,str(r.chemostat),int(r.draw),drawsets[r.group]);out.append(dict(group=r.group,chemostat=str(r.chemostat),draw=int(r.draw),D_error1=float(r.D_error1),label=int(float(r.D_error1)>.01),**feat))
  if (len(out)%24)==0:print(label,len(out),'/',len(df),flush=True)
 result=pd.DataFrame(out);result.to_csv(P/f'{label}_features.csv',index=False);return result
A=features(old,'discovery');B=features(new,'confirmation')
# Validate preserved originally published correction feature (denominator numerical handling may differ)
for tag,base,now in [('discovery',old,A),('confirmation',new,B)]:
 delta=np.abs(now.entropy_exact.to_numpy()-base.entropy_correction_ratio.to_numpy())
 print(tag,'entropy feature absolute discrepancy median',np.median(delta),'max',np.max(delta),flush=True)
 # Note exact matches could differ slightly in integration tolerances; log discrepancies
pred=[];summary=[]
choices={'entropy_exact':['entropy_exact'],'goal_linear':['goal_linear'],'goal_plus_curvature':['goal_plus_curvature'],
         'fisher_norm':['fisher_norm'],'fisher_alignment':['fisher_alignment'],'state_norm':['state_norm'],
         'goal_linear_plus_curvature':['goal_linear','curvature_only']}
for name,cols in choices.items():
 def design(df):return np.log10(np.maximum(df[cols].to_numpy(float),1e-15))
 clf=make_pipeline(StandardScaler(),LogisticRegression(C=1.,max_iter=2000))
 clf.fit(design(A),A.label)
 prob=clf.predict_proba(design(B))[:,1];y=B.label.to_numpy()
 summary.append(dict(method=name,auc=roc_auc_score(y,prob),brier=brier_score_loss(y,prob),log_loss=log_loss(y,np.clip(prob,1e-12,1-1e-12))))
 for k,v in enumerate(prob):pred.append(dict(group=B.iloc[k].group,chemostat=B.iloc[k].chemostat,draw=int(B.iloc[k].draw),label=int(y[k]),method=name,probability=float(v)))
S=pd.DataFrame(summary).sort_values('brier');S.to_csv(P/'matched_metrics.csv',index=False)
pr=pd.DataFrame(pred);pr.to_csv(P/'matched_predictions.csv',index=False)
wide=pr.pivot(index=['group','chemostat','draw','label'],columns='method',values='probability').reset_index();wide['key']=wide.group+'|'+wide.chemostat
keys=wide.key.unique();blocks={key:wide[wide.key==key] for key in keys};rng=np.random.default_rng(2026110805);boot=[]
for j in range(1500):
 x=pd.concat([blocks[k] for k in rng.choice(keys,len(keys),replace=True)],ignore_index=True);y=x.label.to_numpy();r=dict(iteration=j)
 for comp in ['goal_linear','goal_plus_curvature','fisher_norm','fisher_alignment','state_norm']:
  r['brier_entropy_minus_'+comp]=brier_score_loss(y,x.entropy_exact)-brier_score_loss(y,x[comp])
  r['auc_entropy_minus_'+comp]=roc_auc_score(y,x.entropy_exact)-roc_auc_score(y,x[comp]) if len(np.unique(y))==2 else np.nan
 boot.append(r)
pd.DataFrame(boot).to_csv(P/'cluster_bootstrap.csv',index=False)
print(S.to_string(index=False));print('DONE',flush=True)
