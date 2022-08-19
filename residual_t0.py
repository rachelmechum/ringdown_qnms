# XP & XPHM

def r2_t0_XP(name_of_csv,your_title):
    
    read = pd.read_csv(f'{name_of_csv}')
    t0 = read[["t0"]]
    inc = read[["inclination"]]
    r2 = read[["r_squared"]]
    mass = read[["mass"]]
    chi = read[["chi"]]
    
    #slicing times
    # these values are only useful for t0 of .006 to .011 with timestep .001
    t0_0 = t0.iloc[:17, :]
    t0_pi6 = t0.iloc[17:34, :]
    t0_pi3 = t0.iloc[34:51, :]
    t0_pi2 = t0.iloc[51:68, :]
    
    #slicing inclination
    # these values are only useful for t0 of .006 to .011 with timestep .001
    inc_0 = inc.iloc[:17, :]
    inc_pi6 = inc.iloc[17:34, :]
    inc_pi3 = inc.iloc[34:51, :]
    inc_pi2 = inc.iloc[51:68, :]

    #slicing r squared
    # these values are only useful for t0 of .006 to .011 with timestep .001
    r2_0 = r2.iloc[:17, :]
    r2_pi6 = r2.iloc[17:34, :]
    r2_pi3 = r2.iloc[34:51, :]
    r2_pi2 = r2.iloc[51:68, :]
    
    plt.figure(figsize=(10,5))
    plt.grid()
    plt.scatter(t0_0,(1-(r2_0)), marker='^',s=100,label='0')
    plt.scatter(t0_pi6,(1-(r2_pi6)), marker='*',s=150,label='$\pi/6$')
    plt.scatter(t0_pi3,(1-(r2_pi3)), marker='.',s=150,label='$\pi/3$')
    plt.scatter(t0_pi2,(1-(r2_pi2)), marker='x',s=100,label='$\pi/2$')
    plt.legend(title='Inclination',fontsize=15)
    plt.ylabel('$1-R^2$',fontsize=25)
    plt.yscale('log')
    plt.xlabel('$t_0$ (sec)',fontsize=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(f'{your_title}',fontsize=25)
    plt.show()