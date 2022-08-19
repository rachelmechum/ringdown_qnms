# define what modes to use
def rd_fun(modes=[(2,2,0),(2,2,1)]):
    # create function signature
    func_signature="rngdwn_multi_modes(time,M,chi"
    # loop through defined modes
    for mode in modes:
        l,m,n = mode
        func_signature+= f",a{l}{m}{n},phi{l}{m}{n}"
    func_signature+=")"
    # where kwargs = a220, phi220, a221, phi221, ...
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
        # return ex to have combined signal for all defined modes
        return ex
    # this fx now looks like what scipy.curve_fit expects...
    return makefun.create_function(func_signature,new_fun)

# let "test = rd_fun" // test now takes inputs of (time,M,chi,a220,phi220,a221,phi221) --> curve_fit can play with
# this is only being done for XP approx. with 22 mode and 012 overtones
test = rd_fun(modes=[(2,2,0),(2,2,1),(2,2,2)])

#---------------------------------------------------------------------------------------------------------------------------------

def rngdwn(time,M,chi,l,m,n,amplitude,phi):
    
    # using get_ftau to retrieve frequency and damping time for given values (ftau = [freq, tau])
    ftau = ringdown.qnms.get_ftau(M,chi,n,l,m)
    gamma = (ftau[1])**-1
    t0=0
    
    # dict of values for sinusoid
    wf_kws = dict(
    A = amplitude,
    phi = phi,
    f = ftau[0],                                                            
    gamma = gamma,                                                          
    )
    
    # generate sinusoid and store data as ringdown.Data object
    s = amplitude*np.cos(2*np.pi*ftau[0]*(time-t0) + phi)*np.exp(-gamma*abs(time-t0))
    return ringdown.Data(s, index=time)

#---------------------------------------------------------------------------------------------------------------------------------

def plot_ringdown(x0,xdata,ydata,p0,bounds,approx): 
    
'''

    x0     = 'test' with what is defined in rd_fun
    xdata  = defined in XPringdown_plots function 
    ydata  = defined in XPringdown_plots function
    p0     = array of guess values for M,chi,amps,phi
    bounds = array of upper and lower bounds for M,chi,amps,phi
    approx = waveform approximant to use (here IMRPhenomXP)

'''
    
    
    # implementing curve_fit
    popt,pcov = curve_fit(x0,xdata,ydata,p0,bounds=bounds)
    
    # changing timeseries to array
    ydata_array = ydata.data
    x0_exp = x0(xdata, *popt)._values
    # calculating r**2 value
    r = ydata_array - x0_exp 
    ss_res = np.sum(r**2)
    ss_tot = np.sum((ydata_array-np.mean(ydata_array))**2)
    r_squared = 1 - (ss_res / ss_tot)
    #plot
    plt.figure(figsize=(10,4))
    ax1.set_xlim(0,.03)
    ax2.set_xlim(0,.03)
    ax1.set_xlabel('Time (s)',fontsize=20) 
    ax1.set_ylabel('Strain',fontsize=20) 
    ax2.set_frame_on(False)             
    ax2.get_xaxis().tick_bottom()           
    ax2.axes.get_xaxis().set_visible(False)
    ax1.tick_params(axis='both',labelsize=15) 
    ax2.tick_params(axis='both',labelsize=15)
    #labels
    ax1.plot(xdata, ydata, label=f'{approx} Plus')
    ax1.plot(xdata, test(xdata, *popt),label='QNM fit')
    ax1.legend(bbox_to_anchor=(1.008, 1.02),loc='upper right',fontsize=15)
    ax1.legend(bbox_to_anchor=(1.008, 1.02),loc='upper right',fontsize=15)
    
    return popt, r_squared

#---------------------------------------------------------------------------------------------------------------------------------

def XPringdown_plots(mass1,mass2,spin1x,spin1y,spin1z,spin2x,spin2y,spin2z,first_t0,last_t0,time_step,x0,p0,bounds,approx): 
    
'''

    mass1     = mass of BH1
    mass2     = mass of BH2
    spin1x    = x-component of spin1
    spin1y    = y-component of spin1
    spin1z    = z-component of spin1
    spin2x    = x-component of spin2
    spin2y    = y-component of spin2
    spin2z    = z-component of spin2
    first_t0  = earliest start time of ringdown analysis
    last_t0   = latest start time of ringdown anaylsis
    time_step = time step between first_t0 and last_t0
    x0        = 'test'
    p0        = array of guess bounds for M,chi,amp,phi
    bounds    = array of upper and lower bounds for M,chi,amp,phi
    approx    = waveform approximant to use (here IMRPhenomXP)

'''
    inclination = [0, np.pi/6, np.pi/3, np.pi/2]
    df_list=[]
    # generating waveforms for different inclinations
    for i in inclination:        
        hp, hc = get_td_waveform(approximant="IMRPhenomXP",
                                 mass1=mass1,
                                 mass2=mass2,
                                 delta_t=1.0/16384,
                                 f_lower=30,
                                 distance=1e-14,
                                 inclination=i,
                                 spin1x=spin1x,
                                 spin1y=spin1y,
                                 spin1z=spin1z,
                                 spin2x=spin2x,
                                 spin2y=spin2y,
                                 spin2z=spin2z)
        # cropping ringdown time for desired start and end times
        cropped = np.arange((-hp.start_time-(first_t0)), ((abs(hp.start_time))+(last_t0)) , time_step)
        # loopig through each timestep for cropped times
        for crop_start in cropped:
            # sample times ---> timeseries
            # any time 'cond' is used, I am cropping the data down to focus on ringdown
            hp_smp_ts = TimeSeries(hp.sample_times, delta_t=1/16384, epoch=0)
            cond = hp_smp_ts.crop(crop_start,0)                                      

            # h+/hx data ---> timeseries
            hp_data = TimeSeries(hp, delta_t=1/16384, epoch=0)  
            cond_hp = hp_data.crop(crop_start,0)
            hp_signal = cond_hp.numpy()                  

            # repeat
            hc_data = TimeSeries(hc, delta_t=1/16384, epoch=0)              
            cond_hc = hc_data.crop(crop_start,0)
            hc_signal = cond_hc.numpy()                                         

            # use for hp+hc ringdown analysis
            time = cond.numpy()            
            
            x0 = x0         # test defined in rd_fun
            xdata = time    # defined above with timeseries/array
            ydata = cond_hp # defined above with timeseries/array
            p0 = p0         # guess values for mass,chi,amplitudes and phases of each mode configuration
            bounds = bounds # create upper and lower bounds for M,chi,amp,phi
            a, r_squared = plot_ringdown(x0,xdata,ydata,p0,bounds,approx)

            my_dict = dict()
            for ind, v in enumerate(a):
                my_dict[ind] = [v]

            # append necessary values to dataframe
            my_dict[len(a)] = [r_squared]
            my_dict[len(a)+1] = [crop_start]
            my_dict[len(a)+2] = [i]
            temp_df = pd.DataFrame.from_dict(my_dict)
            temp_df.columns = ["mass", "chi",
                               "a220","phi220",
                               "a221","phi221",
                               "a222","phi222",
                               "r_squared","t0","inclination"]
            # instead of times changing for each input case, this adjusts the times
            # to subtract the actual 'start at 0 merger time'
            temp_df['t0'] = temp_df['t0'] - (abs(hp.start_time))
            df_list += [temp_df]
        
    df = pd.concat(df_list,ignore_index=True)
        
    return df