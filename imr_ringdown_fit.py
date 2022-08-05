# creating python function with specific signature
# see https://smarie.github.io/python-makefun/ for reference

def rd_fun(modes=[(2,2,0),(2,2,1)]):
    # create function signature
    func_signature="rngdwn_multi_modes(time,M,chi"
    # loop through defined qn modes
    for mode in modes:
        l,m,n = mode
        # add an amplitude and phase parameter for each mode
        func_signature+= f",a{l}{m}{n},phi{l}{m}{n}"
    func_signature+=")"
    # kwargs = a220, phi220, a221, phi221, ... , (whatever modes you add a_lmn/phi_lmn)
    def new_fun(time,M,chi,**kwargs):
        # creating empty time array to hold rngdwn() return
        ex = np.zeros(time.shape)
        # loop through defined modes
        for mode in modes:
            l,m,n=mode
            amplitude=kwargs[f"a{l}{m}{n}"]
            phi=kwargs[f"phi{l}{m}{n}"]
            # add back to empty time array  
            ex += rngdwn(time,M,chi,l,m,n,amplitude,phi)
            #print(time,M,chi,kwargs)
        # return ex to have combined signal for all defined modes
        return ex
    # this fx now looks like what scipy.curve_fit expects...
    return makefun.create_function(func_signature,new_fun)

##################################################################################################################################

# creating ringdown waveform 

def rngdwn(time,M,chi,l,m,n,amplitude,phi):
    # using get_ftau to retrieve frequency and damping time for given values (ftau = [freq, tau])
    ftau = ringdown.qnms.get_ftau(M,chi,n,l,m)
    gamma = (ftau[1])**-1
    t0=0
        
    wf_kws = dict(
    A = amplitude,
    phi = phi,
    f = ftau[0],                                                            
    gamma = gamma,                                                          
    )
        
    def get_signal(time, A, phi, f, gamma):
        # generate sinusoid
        s = A*np.cos(2*np.pi*f*(time-t0) + phi)*np.exp(-gamma*abs(time-t0))
        return ringdown.Data(s, index=time)

    signal = get_signal(time, **wf_kws)
    
    return signal

##################################################################################################################################

def plot_IMR_curvefit_ringdown(x0,xdata,ydata,p0,bounds,approx):  
    
    # use curve_fit
    popt,pcov = curve_fit(x0,xdata,ydata,p0,bounds=bounds)
    # changing ydata & x0 into appropriate object
    ydata_array = ydata.data
    x0_exp = x0(xdata, *popt)._values
    #residual
    r = ydata_array - x0_exp
    ss_res = np.sum(r**2)
    ss_tot = np.sum((ydata_array-np.mean(ydata_array))**2)
    # r**2 value subtracting the residual from 1
    r_squared = 1 - (ss_res / ss_tot)
    
    #plot
    fig1, ax1 = plt.subplots(figsize=(10,4))
    ax2 = ax1.twiny()
    ax1.set_xlim(0,.03)
    ax2.set_xlim(0,.03)
    ax1.set_xlabel('Time (s)',fontsize=20)
    ax1.set_ylabel('Strain',fontsize=20)
    ax2.set_frame_on(False)             
    ax2.get_xaxis().tick_bottom()           
    ax2.axes.get_xaxis().set_visible(False)
    ax1.tick_params(axis='both',labelsize=15)
    ax2.tick_params(axis='both',labelsize=15)
    #plot fit
    ax1.plot(xdata, ydata, label=f'{approx} Plus',c='blue')
    ax1.plot(xdata, test(xdata, *popt),label='scipy.curve_fit',c='orange')
    ax1.legend(bbox_to_anchor=(1.008, 1.02),loc='upper right',fontsize=15)
    ax1.legend(bbox_to_anchor=(1.008, 1.02),loc='upper right',fontsize=15)
    
    return popt, r_squared

##################################################################################################################################























