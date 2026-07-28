import numpy as np
from scipy.integrate import odeint, solve_ivp


def gLV(y, t, I_simul, g_simul, k_simul):
    dydt = np.zeros_like(y)
    for i in range(len(y)):
        dydt[i] = g_simul[i]* y[i] * (1- np.sum(I_simul[i,:] * y)/k_simul[i])
    return dydt

def run_lotka_volterra(y0, t, s_idx, I, g, k):
    
    s_idx=np.where(s_idx)[0].tolist()
    N=len(y0)
    y0_simul=y0[s_idx]
    I_simul=I[s_idx,:]
    I_simul=I_simul[:,s_idx]
    g_simul=g[s_idx]
    k_simul=k[s_idx]
    
    def f(t,y) : return gLV(y, t, I_simul, g_simul, k_simul)
    t_span = (t[0], t[-1])
    y = solve_ivp(f, t_span, y0_simul, method='RK23', t_eval=t)
    y=y.y[:,-1]
    y_out=np.zeros(N)
    for i in range(y.shape[0]):
        y_out[s_idx[i]]=y[i] 
    return y_out

def run_lotka_volterra_dynamics(y0, t, s_idx, I, g, k):
    
    s_idx=np.where(s_idx)[0].tolist()
    N=len(y0)
    y0_simul=y0[s_idx]
    I_simul=I[s_idx,:]
    I_simul=I_simul[:,s_idx]
    g_simul=g[s_idx]
    k_simul=k[s_idx]
    

    def f(t,y) : return gLV(y, t, I_simul, g_simul, k_simul)
    t_span = (t[0], t[-1])
    y = solve_ivp(f, t_span, y0_simul, t_eval=t)
    #time = y.t
    #y=y.y[:,-1]
    #y_out=np.zeros(N)
    #for i in range(y.shape[1]):
    #    y_out[s_idx[i]]=y[-1,i] 
    return y
