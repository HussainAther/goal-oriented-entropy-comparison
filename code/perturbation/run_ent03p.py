import sys,os
import numpy as np,pandas as pd
from scipy.integrate import solve_ivp
ROOT=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,ROOT)
from model import rhs,jacobian
from runpy import run_path
# Define second directional term independently of side-effecting prior driver.
def hessian_direction(y,u,p):
 n,A,R=y; b,K,eA,a,h,eR,m=[p[k] for k in ('b','K','eA','a','h','eR','m')]
 gp=b*K/(K+n)**2;gpp=-2*b*K/(K+n)**3
 qp=a/(1+a*h*A)**2;qpp=-2*a*a*h/(1+a*h*A)**3
 FA=A*gpp*u[0]**2+2*gp*u[0]*u[1]
 FR=R*qpp*u[1]**2+2*qp*u[1]*u[2]
 return np.array([-FA/eA,FA-FR,eR*FR])
T=np.linspace(0,5,41); scales=[5.,20.,80.]
rng=np.random.default_rng(20261108)
res=[]; draws_used=[]; failed=[]
for group in ['monocultures','polycultures']:
 data=pd.read_csv(os.path.join(ROOT,'inputs',group+'_1000_draws.csv'))
 selected=np.sort(rng.choice(len(data),size=8,replace=False))
 for draw_index in selected:
  d=data.iloc[draw_index]
  for j in range(12):
   if group=='monocultures':
    si=j//4+1;ri=j%4+1;cs=['Bc','Ce','Le'][si-1]+str(ri)
    get=lambda s:d[f'{s}[{si}]']
    ix=lambda s:d[f'{s}[{si},{ri}]']
   else:
    cs=str(j+1).zfill(3)
    get=lambda s:d[s]
    ix=lambda s:d[f'{s}[{j+1}]']
   p=dict(b=float(d.b),K=float(np.exp(d.Klog)),eA=float(np.exp(d.eAlog)),a=float(np.exp(get('alog'))),h=float(np.exp(get('hlog'))),eR=float(get('eR')),m=float(get('m')))
   y0=np.exp([ix('BTNlog'),ix('BTAlog'),ix('BTRlog')]);draws_used.append((group,cs,int(draw_index)))
   def aug(t,z):
    X=z[:3];u=z[3:6];w=z[6:]
    J=jacobian(X,p)
    return np.r_[rhs(t,X,p),J@u,J@w+0.5*hessian_direction(X,u,p)]
   z0=np.r_[y0,1.,0.,0.,0.,0.,0.]
   try:
    sol=solve_ivp(aug,(0,5),z0,t_eval=T,rtol=1e-8,atol=1e-8,method='DOP853')
    if not sol.success:raise RuntimeError(sol.message)
    X,u,w=sol.y.T[:,:3],sol.y.T[:,3:6],sol.y.T[:,6:]
    for dn in scales:
     ypert=y0.copy();ypert[0]+=dn
     so=solve_ivp(lambda t,y:rhs(t,y,p),(0,5),ypert,t_eval=T,rtol=1e-8,atol=1e-8,method='DOP853')
     if not so.success:raise RuntimeError(so.message)
     actual=so.y.T; p1=X+dn*u;p2=p1+dn**2*w
     pos1=bool(np.all(p1>0));pos2=bool(np.all(p2>0))
     def D(z):
      if np.any(z<=0):return None
      return np.sum((z*np.log(z/X)-z+X)/y0,axis=1)
     dtrue=D(actual);d1=D(p1);d2=D(p2)
     valid=dtrue>1e-10
     e1=np.mean(np.sqrt(np.mean(((p1-actual)/np.maximum(X,y0*.01))**2,axis=1)))
     e2=np.mean(np.sqrt(np.mean(((p2-actual)/np.maximum(X,y0*.01))**2,axis=1)))
     err1=float(np.mean(np.abs(d1[valid]-dtrue[valid])/dtrue[valid])) if pos1 and valid.any() else np.nan
     err2=float(np.mean(np.abs(d2[valid]-dtrue[valid])/dtrue[valid])) if pos2 and valid.any() else np.nan
     res.append(dict(group=group,chemostat=cs,draw=int(draw_index),delta_N=dn,initial_N=y0[0],p1_positive=pos1,p2_positive=pos2,state_error1=e1,state_error2=e2,D_error1=err1,D_error2=err2,D2_wins=bool(err2<err1) if np.isfinite(err1) and np.isfinite(err2) else None))
   except Exception as exc:
    failed.append((group,cs,draw_index,str(exc)))
  print(group,'draw',draw_index,'done',flush=True)
metrics=pd.DataFrame(res);metrics.to_csv(os.path.join(ROOT,'cases.csv'),index=False)
pd.DataFrame(draws_used,columns=['group','chemostat','draw']).to_csv(os.path.join(ROOT,'sampled_draws.csv'),index=False)
pd.DataFrame(failed,columns=['group','chemostat','draw','error']).to_csv(os.path.join(ROOT,'failures.csv'),index=False)
rows=[]
for dn,d in metrics.groupby('delta_N'):
 eligible=d.dropna(subset=['D_error1','D_error2']).copy()
 dif=eligible.D_error1.to_numpy()-eligible.D_error2.to_numpy()
 # draw and chemostat are crossed: bootstrap chemostat as cluster, then within-chemostat draw samples
 chemkeys=eligible.group+'|'+eligible.chemostat
 keys=chemkeys.unique();sub={k:eligible.loc[chemkeys==k] for k in keys}
 boot=[];win=[]
 for _ in range(2000):
  selected_keys=rng.choice(keys,len(keys),replace=True)
  arr=np.concatenate([rng.choice((sub[k].D_error1-sub[k].D_error2).to_numpy(),len(sub[k]),replace=True) for k in selected_keys])
  boot.append(float(np.median(arr)));win.append(float(np.mean(arr>0)))
 rows.append(dict(delta_N=dn,total_cases=len(d),eligible=len(eligible),p1_invalid=int((~d.p1_positive).sum()),p2_invalid=int((~d.p2_positive).sum()),D2_wins=int(np.sum(dif>0)),D2_win_share=float(np.mean(dif>0)),median_D1=float(np.median(eligible.D_error1)),median_D2=float(np.median(eligible.D_error2)),median_difference_D1_minus_D2=float(np.median(dif)),median_diff_ci_low=float(np.quantile(boot,.025)),median_diff_ci_high=float(np.quantile(boot,.975)),win_share_ci_low=float(np.quantile(win,.025)),win_share_ci_high=float(np.quantile(win,.975)),median_state1=float(np.median(d.state_error1)),median_state2=float(np.median(d.state_error2)),state2_wins=int((d.state_error2<d.state_error1).sum())))
summary=pd.DataFrame(rows);summary.to_csv(os.path.join(ROOT,'summary.csv'),index=False)
print('\n',summary.to_string(index=False));print('failed',len(failed))
