def ftauq_sub(x,ot,l,m):
    fig,ax=plt.subplots(1,3, figsize=(18,5),sharex=False, sharey=False)
    for overtones in ot: # loop through various overtones to get 3 values (f*tau,tau/tau)
        ftau = []
        tau_tau = []
        q_factor = []
        for chi in x:
            f = ringdown.qnms.get_ftau(1,chi,overtones,l,m)
            f0 = f[0]*T_MSUN
            ftau.append(f0)
            t = ringdown.qnms.get_ftau(1,chi,overtones,l,m)
            t0 = t[1]/T_MSUN
            tau_tau.append(t0)
            q = ringdown.qnms.get_ftau(1,chi,overtones,l,m)
            q0 = q[0]*q[1]*math.pi
            q_factor.append(q0)
        ax[0].plot(x,ftau,label=f'n = ${{{overtones}}}$')
        ax[0].tick_params(labelsize=20)
        ax[1].plot(x,tau_tau,label=f'n = ${{{overtones}}}$')
        ax[1].tick_params(labelsize=20)
        ax[2].plot(x,q_factor,label=f'n = ${{{overtones}}}$')
        ax[2].tick_params(labelsize=20)
        fig.subplots_adjust(wspace=.3)
    temp_string = f'{l},{m}' 
    ax[0].set_ylabel(r'$f_{'+ temp_string+',n}$$t_M$',fontsize=25)
    ax[1].set_ylabel(r'$\tau_{'+ temp_string+',n}$/$t_M$',fontsize=25)
    ax[2].set_ylabel(r'$Q_{'+ temp_string+',n}$',fontsize=25)
    for ax in ax[0],ax[1],ax[2]:
        ax.set_xlabel('$\\chi$', fontsize=25)
        ax.legend(fontsize=15)