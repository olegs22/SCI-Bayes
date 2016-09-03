import numpy as np
from numpy.polynomial.polynomial import polyval


#class with the 21cm signal parametrizations
class model_21cm(object):
    def __init__(self, freq, model_type='gaussian'):
        self.freqs = freq
        self.model = model_type

    def __call__(self, *pars):
        if self.model == 'gaussian':
            T, nu, sigma = pars
            T_b = -T * np.exp(-(nu - self.freqs)**2 / 2. / sigma**2)
            return T_b

        elif self.model == 'tanh':
            x0, xz, xdz, T0, Tz, Tdz, J0, Jz, Jdz = pars

            v_0=1420.4057
            z = v_0 / self.freqs - 1.
            T_cmb = 2.725 * (1. + z)
            Tg = T_cmb * ((1. + z) / (1. + 150.))**2
    
    
            x_par = 0.5 * x0 * (np.tanh((xz - z) / xdz) + 1.)
            T_par = 0.5 * T0 * (np.tanh((Tz - z) / Tdz) + 1.) + Tg
            J_par = 0.5 * J0 * (np.tanh((Jz - z) / Jdz) + 1.)
    
            J_par *= 1e-12
    
            from coupling_coefficients import *
    
            cc=CouplingCoefficients()
            x_c = cc.CollisionalCouplingCoefficient(z,T_par)
            x_a = cc.RadiativeCouplingCoefficient(z,J_par)
    
            T_s = (1. + x_c + x_a) / (T_cmb**-1. + x_c * T_par**-1. + x_a * T_par**-1.)
    
            dTb = 27.0 * (1. - x_par) * np.sqrt(( 1. + z) / 10.) * (1. - T_cmb / T_s)
            return dTb


#class with the foreground model (log not base 10 log)
class foreground(object):
    def __init__(self, freq):
        self.freqs = freq

    def __call__(self, *pars):
        c0, c1, c2 = pars
        params1 = np.array([c0, c1, c2])
        T_gx = np.exp(polyval(np.log(self.freqs / 80.), params1))
        return T_gx

def radiometer(Tsys, tint, channel):
    hz_per_mhz = 1e6
    s_per_hr = 3600.

    return Tsys / np.sqrt(tint * s_per_hr * channel * hz_per_mhz)

# the functions needed for emcee taking in to account the priors

#first we define the likelihood function for each parametrization


def likelihood(pars, model, T_sky, freqs, err):
    if model == 'gaussian':
        signal_pars = pars[:3]
        fore_pars = pars[3:]

        signal = model_21cm(freqs, model)
        fore = foreground(freqs)

        Tb = signal(signal_pars[0], signal_pars[1], signal_pars[2]) * 1e-3
        Tgx = fore(fore_pars[0], fore_pars[1], fore_pars[2])
        T_model = Tb + Tgx

        p = (((T_sky - T_model) / err)**2) + np.log(2. * np.pi * err**2)
        return -0.5 * np.sum(p)

    elif model == 'tanh':
        signal_pars = pars[:9]
        fore_pars = pars[9:]

        signal = model_21cm(freqs, model)
        fore = foreground(freqs)

        Tb = signal(signal_pars) * 1e-3
        Tgx = fore(fore_pars)
        T_model = Tb + Tgx
        
        p = (((T_sky - T_model) / err)**2) + np.log(2. * np.pi * err**2)
        return -0.5 * np.sum(p)

#priors
def priors(pars, list):
    if list[0]<pars[0]<list[1] and list[2]<pars[1]<list[3] and list[4]<pars[2]<list[5]\
       and list[6]<pars[3]<list[7] and list[8]<pars[4]<list[9] and list[10]<pars[5]<list[11]\
       and list[12]<pars[6]<list[13] and list[14]<pars[7]<list[15] and list[16]<pars[8]<list[17]\
       and list[18]<pars[9]<list[19] and list[20]<pars[10]<list[21] and list[22]<pars[11]<list[23]\
       and list[24]<pars[12]<list[25]:
        return 0.0
    return -np.inf

                
def result(pars, model, T_sky, freqs, err, list):
    p = priors(pars, list)
    if not np.isfinite(p):
        return -np.inf
    return p + likelihood(pars, model, T_sky, freqs, err)


#elif model == 'tanh':
    """
if list[0]<pars[0]<list[1] and list[2]<pars[1]<list[3] and list[4]<pars[2]<list[5]\
    and list[6]<pars[3]<list[7] and list[8]<pars[4]<list[9] and list[10]<pars[5]<list[11]:
        return 0.0
    return -np.inf
"""
