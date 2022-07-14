import numpy as np
import gwsurrogate
from pycbc.types import TimeSeries
import makefun

# You only need to pull the catalog once
# gwsurrogate.catalog.pull('NRHybSur3dq8')
# An example of loading the surrogate: 
# sur = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

def wrap_surrogate_to_pycbc(
    sur,
    mass_1,
    mass_2,
    spin_1x=0,
    spin_1y=0,
    spin_1z=0,
    spin_2x=0,
    spin_2y=0,
    spin_2z=0,
    luminosity_distance=1,
    incl=None,
    phi=None,
    f_low=20,
    dt=1 / 4096.,
    mode_list=None):
    """
    sur ::gwsurrogate.surrogate object
        The surrogate to use, as loaded from gwsurrogate.LoadSurrogate()
    mass_1 ::float
        The mass of the primary in MSun
    mass_2 ::float
        The mass of the secondary in MSun
    spin_1x ::float
        The dimensionless spin component in the x direction of the primary
    spin_1y ::float
        The dimensionless spin component in the y direction of the primary
    spin_1z ::float
        The dimensionless spin component in the z direction of the primary
    spin_2x ::float
        The dimensionless spin component in the x direction of the secondary
    spin_2y ::float
        The dimensionless spin component in the y direction of the secondary
    spin_2z ::float
        The dimensionless spin component in the z direction of the secondary
    luminosity_distance ::float
        The luminosity distance to the system in Mpc, defaults to 1
    incl ::float
        The inclination of the system - if not passed then modes will remain uncombined
    phi ::float
        The reference phase of the system - must be passed with incl to combine modes
    f_low ::float
        The lower frequency in Hz at which to generate the waveform
    dt ::float 
        The sampling rate for the waveform, should generally be 1 / 2^n for some n (e.g. 1 / 4096)
    mode_list ::list
        A list of modes, each a tuple (l,m) to generate the waveform with
    """
    # NR Convention
    q = mass_1 / mass_2
    chi1 = [spin_1x, spin_1y, spin_1z]
    chi2 = [spin_2x, spin_2y, spin_2z]
    total_mass = mass_1 + mass_2

    if (incl is None) != (phi is None):
        raise ValueError(
            "Either both of incl and phi should be specified, or neither should be"
        )
    
    # Get waveform from surrogate
    _, h, _ = sur(
        q,
        chi1,
        chi2,
        M=total_mass,
        f_low=f_low,
        dt=dt,
        mode_list = mode_list,
        inclination=incl,
        phi_ref=phi,
        dist_mpc=luminosity_distance,
        units='mks'
    )
    
    # If combined modes then h is just a complex valued np array
    if incl is not None:
        hp = TimeSeries(np.real(h), delta_t=dt)
        hc = TimeSeries(-np.imag(h), delta_t=dt)
    # Not combined modes --> h is a dict of complex valued arrays
    else:
        hp = dict()
        hc = dict()
        for mode, hlm in h.items():
            hp[mode] = TimeSeries(np.real(h[mode]), delta_t=dt)
            hc[mode] = TimeSeries(-np.imag(h[mode]), delta_t=dt)

    return hp, hc

# A demonstration, here we are constructing a polynomial to fit with
def polynomial_term(x, a, n):
    # make a single term in a polynomial, given the coefficient and the power
    return a * x ** n

# This function will return a combined function
def make_polynomial_function(nmax=2):
    # Setup the function's desired call signature
    func_sig = "polynomial_function(x"
    for n in range(nmax):
        func_sig += f",a{n}"
    func_sig += ")"
    print(func_sig)
    
    # make the actual function itself, using *args
    def unwrapped_polynomial(x, **kwargs):
        polynomial = 0
        for key, val in kwargs.items():
            order = float(key[1:])
            polynomial += polynomial_term(x, val, order)
        return polynomial
            
    return makefun.create_function(func_sig, unwrapped_polynomial)

# Example fitting with out artificial function
# test_x_values = np.arange(100)
# polynomial_to_fit = 4 + 5 * test_x_values + 3 * test_x_values ** 2 + 2 * test_x_values ** 3
# fitting_func = make_polynomial_function(nmax=4)
# from scipy.optimize import curve_fit
# attempt_fit = curve_fit(fitting_func, test_x_values, polynomial_to_fit, [0,0,0,0])