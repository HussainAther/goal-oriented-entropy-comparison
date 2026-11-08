import numpy as np
def hessian_direction(y,u,p):
 n,A,R=y; b,K,eA,a,h,eR,m=[p[k] for k in ('b','K','eA','a','h','eR','m')]
 gp=b*K/(K+n)**2;gpp=-2*b*K/(K+n)**3
 qp=a/(1+a*h*A)**2;qpp=-2*a*a*h/(1+a*h*A)**3
 FA=A*gpp*u[0]**2+2*gp*u[0]*u[1]
 FR=R*qpp*u[1]**2+2*qp*u[1]*u[2]
 return np.array([-FA/eA,FA-FR,eR*FR])
