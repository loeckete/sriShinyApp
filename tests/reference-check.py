"""Independent equations and high-accuracy reference. Requires scipy, only for development."""
import csv
from scipy.integrate import solve_ivp

def f(model):
    def rhs(t,y):
        if model=='SIS':
            s,i=y; return [-.004*s*i+.2*i,.004*s*i-.2*i]
        if model=='SEIR':
            s,e,i,r=y; return [-.004*s*i,.004*s*i-.25*e,.25*e-.2*i,.2*i]
        if model=='SIRD':
            s,i,r,d=y; return [-.004*s*i,.004*s*i-.25*i,.2*i,.05*i]
        s,i,r=y; return [-.004*s*i,.004*s*i-.2*i,.2*i]
    return rhs
with open('tests/reference.csv','w') as file:
    w=csv.writer(file); w.writerow(['model','time','state','value'])
    for model,states,initial in [('SIR','SIR',[99,1,0]),('SIS','SI',[99,1]),('SEIR','SEIR',[99,0,1,0]),('SIRD','SIRD',[99,1,0,0])]:
        sol=solve_ivp(f(model),(0,60),initial,t_eval=[0,5,10,20,40,60],method='DOP853',rtol=1e-12,atol=1e-13)
        for j,t in enumerate(sol.t):
            for k,state in enumerate(states): w.writerow([model,t,state,sol.y[k,j]])
        print(model, 'reference complete')
