'''
A generalized precession parameter $\chip$ to interpret gravitational-wave data
Gerosa et al (2020). arXiv:2011.11948
'''

import numpy as np
import scipy
import precession # Try: pip install precession

# Code units: G=c=M=1

def ftor_PN(f, M_msun, q, chi1, chi2, theta1, theta2, deltaphi):
    '''Convert GW frequency to PN orbital separation conversion'''

    c_cgs = 2.99e10
    G_cgs = 6.67e-8
    om = np.pi * f
    M_sec = M_msun * 2e33 * G_cgs / c_cgs**3
    mom = M_sec * om
    m1 = 1 / (1+q)
    m2 = q / (1+q)
    eta = m1*m2
    ct1 = np.cos(theta1)
    ct2 = np.cos(theta2)
    ct12 = np.sin(theta1) * np.sin(theta2) * np.cos(deltaphi) + ct1 * ct2
    # Eq. 4.13, Kidder 1995. gr-qc/9506022
    r = (mom)**(-2./3.)*(1. \
                    - (1./3.)*(3.-eta)*mom**(2./3.)  \
                    - (1./3.)* ( chi1*ct1*(2.*m1**2.+3.*eta) + chi2*ct2*(2.*m2**2.+3.*eta))*mom \
                    + ( eta*(19./4. + eta/9.) -eta*chi1*chi2/2. * (ct12 - 3.*ct1*ct2 ))*mom**(4./3.)\
                    )
    return r


def omegatilde(q):
    '''Ratio between the spin frequency, leading order term. Eq (13)'''
    return q*(4*q+3)/(4+3*q)

def chip_terms(q,chi1,chi2,theta1,theta2):
    '''Two chip terms'''
    term1 = chi1*np.sin(theta1)
    term2 = omegatilde(q) * chi2*np.sin(theta2)
    return term1,term2

def chip_heuristic(q,chi1,chi2,theta1,theta2):
    '''Heuristic definition of chip. Eq (14)'''
    term1, term2 = chip_terms(q,chi1,chi2,theta1,theta2)
    return np.maximum(term1,term2)

def chip_generalized(q,chi1,chi2,theta1,theta2,deltaphi):
    '''Generalized definition of chip. Eq (15)'''
    term1, term2 = chip_terms(q,chi1,chi2,theta1,theta2)
    return (term1**2 + term2**2 + 2*term1*term2*np.cos(deltaphi))**0.5

def chip_asymptotic(q,chi1,chi2,theta1,theta2):
    '''Asymptotic definition of chip. Eq (19)'''
    term1, term2 = chip_terms(q,chi1,chi2,theta1,theta2)
    return (np.abs(term1-term2) * scipy.special.ellipe(-4*term1*term2/(term1-term2)**2) + np.abs(term1+term2) * scipy.special.ellipe(4*term1*term2/(term1+term2)**2))/np.pi



@np.vectorize
def chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,r=None,fref=None, M_msun=None):
    '''Averaged definition of chip. Eq (15) and Appendix A'''

    # Convert frequency to separation, if necessary
    if r is None and fref is None: raise ValueError
    elif r is not None and fref is not None: raise ValueError
    if r is None:
        if M_msun is None: raise ValueError
        # Eq A1
        r = ftor_PN(fref, M_msun, q, chi1, chi2, theta1, theta2, deltaphi)
    #Compute constants of motion
    L = (r**0.5)*q/(1+q)**2
    S1 = chi1/(1.+q)**2
    S2 = (q**2)*chi2/(1.+q)**2
    # Eq 5
    chieff = (chi1*np.cos(theta1)+q*chi2*np.cos(theta2))/(1+q)
    # Eq A2
    J = (L**2 + S1**2 + S2**2 + 2*L*(S1*np.cos(theta1) + S2*np.cos(theta2)) + 2*S1*S2*(np.sin(theta1)*np.sin(theta2)*np.cos(deltaphi) + np.cos(theta1)*np.cos(theta2)))**0.5
    # Solve dSdt=0. Details in arXiv:1506.03492
    Sminus,Splus=precession.Sb_limits(chieff,J,q,S1,S2,r)

    def integrand_numerator(S):
        '''chip(S)/dSdt(S)'''
        #Eq A3
        theta1ofS = np.arccos( (1/(2*(1-q)*S1)) * ( (J**2-L**2-S**2)/L - 2*q*chieff/(1+q) ) )
        #Eq A4
        theta2ofS = np.arccos( (q/(2*(1-q)*S2)) * ( -(J**2-L**2-S**2)/L + 2*chieff/(1+q) ) )
        # Eq A5
        deltaphiofS = np.arccos( ( S**2-S1**2-S2**2 - 2*S1*S2*np.cos(theta1ofS)*np.cos(theta2ofS) ) / (2*S1*S2*np.sin(theta1ofS)*np.sin(theta2ofS)) )
        # Eq A6 (prefactor cancels out)
        dSdtofS = np.sin(theta1ofS)*np.sin(theta2ofS)*np.sin(deltaphiofS)/S
        # Eq (15)
        chipofS = chip_generalized(q,chi1,chi2,theta1ofS,theta2ofS,deltaphiofS)
        return chipofS/dSdtofS

    numerator = scipy.integrate.quad(integrand_numerator, Sminus, Splus)[0]

    def integrand_denominator(S):
        '''1/dSdt(S)'''
        #Eq A3
        theta1ofS = np.arccos( (1/(2*(1-q)*S1)) * ( (J**2-L**2-S**2)/L - 2*q*chieff/(1+q) ) )
        #Eq A4
        theta2ofS = np.arccos( (q/(2*(1-q)*S2)) * ( -(J**2-L**2-S**2)/L + 2*chieff/(1+q) ) )
        # Eq A5
        deltaphiofS = np.arccos( ( S**2-S1**2-S2**2 - 2*S1*S2*np.cos(theta1ofS)*np.cos(theta2ofS) ) / (2*S1*S2*np.sin(theta1ofS)*np.sin(theta2ofS)) )
        # Eq A6 (prefactor cancels out)
        dSdtofS = np.sin(theta1ofS)*np.sin(theta2ofS)*np.sin(deltaphiofS)/S

        return 1/dSdtofS

    denominator = scipy.integrate.quad(integrand_denominator, Sminus, Splus)[0]

    return numerator/denominator



if __name__ == "__main__":

    if True:
        q=0.7
        chi1=0.3
        chi2=1
        theta1=np.pi/3
        theta2=np.pi/4
        deltaphi=np.pi/5
        fref= 20 # Hz
        M_msun=60 # Msun
        print("Heuristic chip:", chip_heuristic(q,chi1,chi2,theta1,theta2))
        print("Asymptotic chip:", chip_asymptotic(q,chi1,chi2,theta1,theta2))
        print("Generalized chip:", chip_generalized(q,chi1,chi2,theta1,theta2,deltaphi))
        print("Averaged chip:", chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,fref=fref,M_msun=M_msun))

        #Alternatively:
        #r=10
        #print("Averaged chip:", chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,r=r))

    #Also works with numpy array:
    if False:
        q=np.array([0.7,0.7])
        chi1=np.array([0.3,0.3])
        chi2=np.array([1,1])
        theta1=np.array([np.pi/3,np.pi/3])
        theta2=np.array([np.pi/4,np.pi/4])
        deltaphi=np.array([np.pi/5,np.pi/5])
        fref= np.array([20,20]) # Hz
        M_msun= np.array([60,60]) # Msun
        print("Heuristic chip:", chip_heuristic(q,chi1,chi2,theta1,theta2))
        print("Asymptotic chip:", chip_asymptotic(q,chi1,chi2,theta1,theta2))
        print("Generalized chip:", chip_generalized(q,chi1,chi2,theta1,theta2,deltaphi))
        print("Averaged chip:", chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,fref=fref,M_msun=M_msun))
