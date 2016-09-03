"""

coupling_coefficients.py

Author: Jordan Mirocha
Affiliation: UCLA
Created on: Thu Apr 14 10:30:09 PDT 2016

Description: Adapted from Hydrogen class of the ARES code.

References:
Zygelman, B. 2005, ApJ, 622, 1356

"""

import numpy as np

A10 = 2.85e-15       # HI 21cm spontaneous emission coefficient 
T_star = 0.068       # Temperature difference between HI hyperfine states 

# Rate coefficients for spin de-excitation - from Zygelman (2005)

# H-H collisions.
T_HH = \
    [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, \
     25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, \
     90.0, 100.0, 200.0, 300.0, 500.0, 700.0, 1000.0, \
     2000.0, 3000.0, 5000.0, 7000.0, 10000.0]

kappa_HH = \
    [1.38e-13, 1.43e-13, 2.71e-13, 6.60e-13, 1.47e-12, 2.88e-12, \
     9.10e-12, 1.78e-11, 2.73e-11, 3.67e-11, 5.38e-11, 6.86e-11, \
     8.14e-11, 9.25e-11, 1.02e-10, 1.11e-10, 1.19e-10, 1.75e-10, \
     2.09e-10, 2.56e-10, 2.91e-10, 3.31e-10, 4.27e-10, 4.97e-10, \
     6.03e-10, 6.87e-10, 7.87e-10]
            
# H-e collisions.            
T_He = \
    [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 
     1000.0, 2000.0, 3000.0, 5000.0, 7000.0, 
     10000.0, 15000.0, 20000.0]
        
kappa_He = \
    [2.39e-10, 3.37e-10, 5.30e-10, 7.46e-10, 1.05e-9, 1.63e-9, 
     2.26e-9, 3.11e-9, 4.59e-9, 5.92e-9, 7.15e-9, 7.71e-9, 
     8.17e-9, 8.32e-9, 8.37e-9, 8.29e-9, 8.11e-9]

T_HH = np.array(T_HH)
T_He = np.array(T_He)

class CouplingCoefficients(object):
    def __init__(self):
        """
        Parameters
        ----------
        interp_method: str
            Either 'linear' or 'cubic'. The latter will only work if you
            have scipy installed.
            
        """
        
        self.tabulated_coeff = {'kappa_H': kappa_HH, 'kappa_e': kappa_He, 
            'T_H': T_HH, 'T_e': T_He}

        nH0 = 1.889e-07
        self.nH = lambda z: nH0 * (1. + z)**3

    def kappa_H(self, Tk):
        """
        Rate coefficient for spin-exchange via H-H collsions.
        """
        
        return np.interp(Tk, T_HH, kappa_HH)

    def kappa_e(self, Tk):       
        """
        Rate coefficient for spin-exchange via H-electron collsions.
        """       
        
        return np.interp(Tk, T_He, kappa_He)

    def CollisionalCouplingCoefficient(self, z, Tk):
        """
        Parameters
        ----------
        z : float
            Redshift
        Tk : float
            Kinetic temperature of the gas [K]
                
        """
                
        return self.nH(z) * self.kappa_H(Tk) * T_star / A10 \
            / (2.715 * (1. + z))
    
    def RadiativeCouplingCoefficient(self, z, Ja):
        """
        Return radiative coupling coefficient (i.e., Wouthuysen-Field effect).
        
        .. note :: Assumes S_alpha = 1.
        
        """
                
        return 1.81e11 * Ja / (1. + z)


if __name__ == '__main__':
    cc = CouplingCoefficients()
    
    import matplotlib.pyplot as pl
    
    pl.scatter(cc.tabulated_coeff['T_H'], cc.tabulated_coeff['kappa_H'], 
        color='k')
    pl.scatter(cc.tabulated_coeff['T_e'], cc.tabulated_coeff['kappa_e'], 
        color='k')
        
    # Interpolated values
    T = np.logspace(0, 4.5)
    pl.loglog(T, cc.kappa_H(T), color = 'k', ls = '-', 
        label=r'$\kappa_{10}^{\mathrm{HH}}$')
    pl.loglog(T, cc.kappa_e(T), color = 'k', ls = '--', 
        label=r'$\kappa_{10}^{\mathrm{eH}}$')
        
    pl.xlabel(r'$T / K$')
    pl.ylabel(r'$\kappa_{10}$')
    pl.ylim(1e-13, 1e-7)
    pl.legend(loc='lower right')
    
    
        