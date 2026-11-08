"""ENT-03N: Published Wojcik et al. nutrient-algae-rotifer ODE (post- and pre-pulse).
Parameters b,K,eA,a,h,eR,m are on NATURAL scales. State is (nitrogen, algae, rotifers).
Python port of model_mono/model_poly in get_bayesian_predictions_*.R.
"""
import numpy as np
from scipy.integrate import solve_ivp

DILUTION=0.2
FEED_N=80.0

def rhs(t, state, pars):
    n, algae, rot = state
    b,K,eA,a,h,eR,m = [float(pars[k]) for k in ('b','K','eA','a','h','eR','m')]
    if K <= 0 or eA <= 0 or a < 0 or h < 0:
        raise ValueError('Invalid nonpositive saturation/conversion parameter')
    FA = b*n*algae/(K+n)
    FR = a*algae*rot/(1+a*h*algae)
    return np.array([DILUTION*(FEED_N-n)-FA/eA,
                     FA-FR-DILUTION*algae,
                     eR*FR-(m+DILUTION)*rot])

def jacobian(y,pars):
    n,A,R = np.asarray(y,dtype=float)
    b,K,eA,a,h,eR,m=[float(pars[k]) for k in ('b','K','eA','a','h','eR','m')]
    g=b*n/(K+n); gp=b*K/(K+n)**2
    q=a*A/(1+a*h*A); qp=a/(1+a*h*A)**2
    return np.array([
     [-DILUTION-A*gp/eA,-g/eA,0.0],
     [A*gp,g-qp*R-DILUTION,-q],
     [0.0,eR*qp*R,eR*q-m-DILUTION]
    ])

def integrate(y0,pars,times,rtol=1e-9,atol=1e-9):
    times=np.asarray(times,float)
    if np.any(np.diff(times)<0): raise ValueError('times must be nondecreasing')
    if np.any(np.asarray(y0)<=0): raise ValueError('initial states must be positive')
    res=solve_ivp(lambda t,y:rhs(t,y,pars),(float(times[0]),float(times[-1])),np.asarray(y0,float),
                  t_eval=times,rtol=rtol,atol=atol, method='DOP853')
    if not res.success: raise RuntimeError(res.message)
    if np.any(res.y <= 0): raise RuntimeError('non-positive trajectory')
    return res.y.T
