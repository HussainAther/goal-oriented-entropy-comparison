"""ENT-03Q: Reference-trajectory-only diagnostics for first-order approximation failure.
Inputs: ENT03P posterior draws and cases. No nonlinear perturbed outcomes are used as features.
"""
import os,sys
import numpy as np,pandas as pd
from scipy.integrate import solve_ivp
from sklearn.metrics import roc_auc_score,average_precision_score
ROOT=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,ROOT)
from model import rhs,jacobian
from run_ent03p_helpers import hessian_direction
base=pd.read_csv(os.path.join(ROOT,'inputs','cases.csv'))
T=np.linspace(0,5,41)
features=[]; failures=[]
for (group,cs,draw),rows in base.groupby(['group','chemostat','draw'],sort=False):
 try:
  d=pd.read_csv(os.path.join(ROOT,'inputs',group+'_1000_draws.csv')).iloc[int(draw)]
  if group=='monocultures':
   si=['Bc','Ce','Le'].index(cs[:2])+1;ri=int(cs[2:]);get=lambda s:d[f'{s}[{si}]']; ix=lambda s:d[f'{s}[{si},{ri}]']
  else:
   ci=int(cs);get=lambda s:d[s];ix=lambda s:d[f'{s}[{ci}]']
  p=dict(b=float(d.b),K=float(np.exp(d.Klog)),eA=float(np.exp(d.eAlog)),a=float(np.exp(get('alog'))),h=float(np.exp(get('hlog'))),eR=float(get('eR')),m=float(get('m')))
  y0=np.exp([ix('BTNlog'),ix('BTAlog'),ix('BTRlog')]);S=np.diag(y0);Sinv=np.diag(1/y0)
  z0=np.r_[y0,1.,0.,0.,0.,0.,0.,np.eye(3).ravel()]
  def aug(t,z):
   X=z[:3];u=z[3:6];w=z[6:9];Phi=z[9:].reshape(3,3);J=jacobian(X,p)
   return np.r_[rhs(t,X,p),J@u,J@w+.5*hessian_direction(X,u,p), (J@Phi).ravel()]
  sol=solve_ivp(aug,(0,5),z0,t_eval=T,rtol=1e-8,atol=1e-8,method='DOP853')
  if not sol.success:raise ValueError(sol.message)
  X=sol.y.T[:,:3];u=sol.y.T[:,3:6];w=sol.y.T[:,6:9];Phi=sol.y.T[:,9:].reshape(-1,3,3)
  gains=[]; reacts=[]
  for t in range(len(T)):
   M=Sinv@Phi[t]@S
   gains.append(np.linalg.svd(M,compute_uv=False)[0])
   Jsc=Sinv@jacobian(X[t],p)@S
   reacts.append(np.linalg.eigvalsh((Jsc+Jsc.T)/2)[-1])
  valid= T>0
  for delta in rows.delta_N:
   p1=X+delta*u;p2=p1+delta**2*w
   def div(Y):
    z=Y/X
    return np.sum((Y*np.log(z)-Y+X)/y0,axis=1)
   d1=div(p1);d2=div(p2)
   f=dict(group=group,chemostat=cs,draw=draw,delta_N=delta,
       entropy_correction_ratio=float(np.mean(np.abs(d2[valid]-d1[valid])/np.maximum(d1[valid],1e-12))),
       entropy_correction_abs=float(np.mean(np.abs(d2[valid]-d1[valid]))),
       state_correction_ratio=float(np.mean(np.linalg.norm((delta**2*w[valid])/y0,axis=1)/np.maximum(np.linalg.norm((delta*u[valid])/y0,axis=1),1e-14))),
       state_second_correction=float(np.mean(np.linalg.norm(delta**2*w[valid]/y0,axis=1))),
       first_state_size=float(np.mean(np.linalg.norm(delta*u[valid]/y0,axis=1))),
       propagator_gain=float(max(gains)),
       reactivity_max=float(max(reacts)),
       reactivity_time_avg=float(np.trapezoid(reacts,T)/(T[-1]-T[0])))
   features.append(f)
 except Exception as ex:failures.append(dict(group=group,chemostat=cs,draw=draw,error=repr(ex)))
 if len(features)%48==0: print('completed',len(features),flush=True)
f=pd.DataFrame(features)
f.to_csv(os.path.join(ROOT,'reference_only_features.csv'),index=False)
pd.DataFrame(failures,columns=['group','chemostat','draw','error']).to_csv(os.path.join(ROOT,'failures.csv'),index=False)
m=base.merge(f,on=['group','chemostat','draw','delta_N'],validate='one_to_one')
m['chemostat_key']=m.group+'|'+m.chemostat
m.to_csv(os.path.join(ROOT,'joined_cases.csv'),index=False)
feat=['entropy_correction_ratio','entropy_correction_abs','state_correction_ratio','state_second_correction','first_state_size','propagator_gain','reactivity_max','reactivity_time_avg']
rows=[]
for delta,grp in m.groupby('delta_N'):
 for target,threshold in [('D_error1',.01),('D_error1',.05),('state_error1',.01),('state_error1',.05)]:
  label=(grp[target]>threshold).astype(int).to_numpy()
  for key in feat:
   score=grp[key].to_numpy()
   rows.append(dict(delta_N=delta,target=target,threshold=threshold,positive=int(label.sum()),n=len(label),predictor=key,auc=roc_auc_score(label,score) if 0<label.sum()<len(label) else np.nan,average_precision=average_precision_score(label,score) if 0<label.sum()<len(label) else np.nan,median_indicator=float(np.median(score))))
A=pd.DataFrame(rows);A.to_csv(os.path.join(ROOT,'discrimination.csv'),index=False)
print('cases',len(m),'failures',len(failures))
print(A[(A.target=='D_error1') & (A.threshold==.05)].pivot(index='predictor',columns='delta_N',values='auc').round(3).to_string())
print('label frequencies',A.groupby(['delta_N','target','threshold']).positive.first().to_string())
