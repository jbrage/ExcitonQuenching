import numpy as np
from math import pi
from scipy.interpolate import interp1d

'''
Implementation of the Bethe-Bloch stopping power formula for ions (in water)
'''

m_proton = 1.672621e-27         # [kg]
m_electron = 9.10938e-31        # [kg]
q_electron = 1.60217e-19        # [C]
c = 299792458                   # [m/s]
I_Joule = 75*q_electron         # Average ionization energy of water [Joule]
N_A = 6.0221409e+23             # [/mol]
epsilon_0 = 8.85418e-12         # [F/m]
Z_A_water = 0.555               # ratio of Z to A in water


def beta_to_gamma(beta):
    '''
    beta to gamma
    '''
    return 1.0/np.sqrt(1 - beta**2)


def beta_to_Tmax_Joule(beta, A):
    '''
    INPUT:
    - beta = v/c
    - nucleon number of the particle

    OUTPUT:
    - Max kinematic energy [joule] transferred from particle to electron
    '''
    gamma = beta_to_gamma(beta)
    m_particle = A*m_proton
    R = m_electron/m_particle
    Tmax_Joule = 2*m_electron * c**2 * beta**2 * gamma**2/(1+2*gamma*R + R**2)
    return Tmax_Joule


def E_MeV_to_beta(E_MeV_per_A, A_ion):
    '''
    INPUT:
    - kinetic energy per nucleon [MeV/A]
    - particle's nucleon number

    OUTPUT:
    - beta = v/c
    '''
    E_MeV = A_ion * E_MeV_per_A
    E_Joule = E_MeV*1e6*q_electron          # total kinetic energy [Joule]
    m_particle = A_ion*m_proton        # ion mass [kg]
    beta = np.sqrt(1. - 1./(1. + E_Joule/(m_particle*c**2))**2)
    return beta


def load_ShellCorrection_data(E_MeV_per_A, A_ion):
    '''
    Shell correction (SC): 0.1 to 10 MeV protons in water for I = 75 eV
    Data extracted from DOI 10.1088/0031-9155/54/11/012

    The shell correction is the same for a given velocity in a material;
    The proton energy is converted to a velocity and the SC is interpolated as
    a function of velocity (it drops quickly to 0 for increasing v)

    INPUT:
    - kinetic energy per nucleon [MeV/A]
    - nucleon number

    OUTPUT:
    - shell correction in liquid water for that particle energy
    '''

    # load the shell correction data and interpolate as a function of velocity
    fdir = "EQ_model_functions/"
    if __name__ == '__main__':
        fdir = ""

    fname = "%sshell_correction_water.dat" % fdir
    data = np.genfromtxt(fname, delimiter=",", dtype=float)
    E_keV, SC = data[:, 0], data[:, 1]
    E_MeV = E_keV*1e-3
    v_data = E_MeV_to_beta(E_MeV, A_ion=1)*c  # [m/s]
    interpol = interp1d(v_data, SC, fill_value="extrapolate")

    # velocity of particle in question
    v_particle = E_MeV_to_beta(E_MeV_per_A, A_ion)*c

    SC_interpol = interpol(v_particle)
    if SC_interpol < 0.0:
        SC_interpol = 0.0
    return SC_interpol


def check_BetheBloch_limit(
        E_MeV_per_A, A_ion, PRINT_WARNING, Z_A_material=7.5):

    '''
    Check if the Bethe-Bloch equation is 'applicable',
    Attix chapter 8: Stopping power, A Soft-collision term
    (although) is in good agreement with experimental data for lower energies

    A less strict (not implemented) requirement is beta > 0.02.

    INPUT:
    - kinetic energy per nucleon [MeV/A]
    - nucleon number
    - boolean PRINT_WARNING (True prints the warning)
    - Z/A of the material, defaults to liquid water

    OUTPUT:
    - Boolean validity
    '''
    beta = E_MeV_to_beta(E_MeV_per_A, A_ion)
    ratio = (Z_A_material / 137. / beta)**2
    valid_energy = True

    if ratio > 1 or beta > 0.99:
        valid_energy = False
        if PRINT_WARNING:
            msg = "\n# %0.2f MeV/A \n# => (Z/137 beta)**2 = %0.2f" % (
                                                            E_MeV_per_A, ratio)
            msg += "\n# (Z/137beta)**2 << 1 isn't satisfied\n"
            print(msg)
    return valid_energy


def BetheBlochEq(
        E_MeV_per_A, z_ion, A_ion, Z_A_material=Z_A_water,
        shell_correction=False, PRINT_WARNING=False):
    '''
    Calculates the Bethe-Bloch stopping power

    INPUT:
    - kinetic energy per nucleon [MeV/A]
    - z_ion, particle charge, multipla of the elementary charge
    - A_ion, nucleon number
    - Z_A_material, Z/A of the given material
    - shell_correction (boolean), True loads the shell corrections (liq. H2O)
    - PRINT_WARNING (boolean),True prints a warning of the particle is too slow

    OUTPUT:
    - stopping power in MeV cm^2/g
    - True/False depending on whether the particle energy is too low to a
    reliable calculation of stopping power use Bethe-Bloch
    '''
    valid_energy = check_BetheBloch_limit(E_MeV_per_A, A_ion, PRINT_WARNING)

    # OBS, The factor 1/(4*pi*epsilon_0) = 1 in (many) textbooks
    K1 = 4*pi/(m_electron*c**2)                 # 1/Joule
    K2 = (q_electron**2/(4*pi*epsilon_0))**2    # Joule^2 *m^2
    K = K1*K2                                   # Joule * m^2
    K *= N_A    # Joule * m^2 /mol => A_water has units of g/mol => Joule*m^2/g
    K /= q_electron                             # to eV *m^2/g
    K *= 1e4                                    # to eV *cm^2/g
    K *= 1e-6                                   # to MeV cm^2 /g

    beta = E_MeV_to_beta(E_MeV_per_A, A_ion)
    gamma = beta_to_gamma(beta)
    Tmax_Joule = beta_to_Tmax_Joule(beta, A_ion)

    # full equation
    front_factor = K * z_ion**2 * Z_A_material / beta**2
    term1 = 0.5*np.log(2*m_electron*c**2*beta**2*gamma**2*Tmax_Joule/I_Joule**2)
    term2 = beta**2

    if shell_correction:
        C_Z = load_ShellCorrection_data(E_MeV_per_A, A_ion)
        dEdx_MeV_cm2_g = front_factor*(term1 - term2 - C_Z)
    else:
        dEdx_MeV_cm2_g = front_factor*(term1 - term2)

    if dEdx_MeV_cm2_g < 0.0:
        dEdx_MeV_cm2_g = np.nan
    return dEdx_MeV_cm2_g, valid_energy


if __name__ == '__main__':

    '''
    Run an example
    '''
    Energies_MeV = np.logspace(0, 3, 10)

    z_ion = 4
    A_ion = 8

    print("# z = %i, A = %i" % (z_ion, A_ion))

    txt = "# E [MeV/A], \tdEdx [MeV cm^2/g] \tEq. applicable?"
    print(txt)
    for E_MeV_A in Energies_MeV:
        dEdx_MeV_cm2_g, valid_range = BetheBlochEq(
            E_MeV_A, z_ion, A_ion, Z_A_water)

        print("%0.2f,\t\t%0.2f\t\t\t%s" % (E_MeV_A,dEdx_MeV_cm2_g,valid_range))
