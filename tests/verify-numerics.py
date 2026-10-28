"""Check a Python reproduction of model.R's RK4 scheme, not execution of R itself."""
import numpy as np
from scipy.integrate import solve_ivp
import json

def simulate(model, beta=.004, gamma=.2, m=0, v=0, a=.25, mu=.05, days=60, max_step=.1):
    states={'SIR':'SIR','SIS':'SI','SEIR':'SEIR','SIRD':'SIRD'}[model]
    mu=mu if model=='SIRD' else 0
    if model=='SIS': v=0
    init=np.array([{'S':99.,'E':0.,'I':1.,'R':0.,'D':0.}[s] for s in states])
    def f(t,y):
        q=dict(zip(states,y)); n=sum(z for k,z in q.items() if k!='D'); s=q['S']; i=q['I']; inf=beta*s*i
        z={'S':m*n*(1-v)-m*s-inf,'I':inf-(gamma+m+mu)*i}
        if model=='SIS': z['S']+=gamma*i
        else: z['R']=gamma*i-m*q['R']+m*v*n
        if model=='SEIR': z.update(E=inf-(a+m)*q['E'],I=a*q['E']-(gamma+m)*i)
        if model=='SIRD': z['D']=mu*i
        return np.array([z[s] for s in states])
    times=np.arange(0,days+.001,.25); out=[init.copy()]; y=init.copy()
    rate=max(beta*sum(init)+m,a+m,gamma+m+mu,1)
    hmax=min(max_step,.2/rate)
    for t0,t1 in zip(times[:-1],times[1:]):
        steps=int(np.ceil((t1-t0)/hmax)); h=(t1-t0)/steps
        for _ in range(steps):
            k1=f(0,y); k2=f(0,y+h*k1/2); k3=f(0,y+h*k2/2); k4=f(0,y+h*k3)
            y=y+h*(k1+2*k2+2*k3+k4)/6
        out.append(y.copy())
    out=np.array(out)
    ref=solve_ivp(f,(0,days),init,t_eval=times,method='DOP853',rtol=1e-12,atol=1e-13).y.T
    err=float(abs(out-ref).max()); mass=float(abs(out.sum(axis=1)-init.sum()).max())
    assert err<1e-5 and mass<1e-7 and out.min()>=-1e-10
    return {'model':model,'m':m,'v':v,'maximum_reference_error':err,'maximum_mass_error':mass,'minimum_state':float(out.min()),'initial_I_change':float(f(0,init)[states.index('I')])}
results=[]
for model in ['SIR','SIS','SEIR','SIRD']:
    results.append(simulate(model)); results.append(simulate(model,m=.03,v=.2))
results.append(simulate('SIR',beta=.0015)); results.append(simulate('SIR',gamma=.4))
print(json.dumps({'scope':'Python numerical reproduction only; R execution not verified','checks':results},indent=2))
