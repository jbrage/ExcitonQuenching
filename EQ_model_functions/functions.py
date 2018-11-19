import numpy as np
from scipy.integrate import trapz
from EQ_model_functions.Bethe_Bloch_stopping_power import BetheBlochEq, E_MeV_to_beta


def stopping_power(E_MeV_per_A, z_projectile, A_projectile, density_g_cm3, Z_A_scintillator):
    '''
    Calculate the stopping power use the Bethe-Bloch formula

    INPUT:
    - kinetic energy per nucleon [MeV/A]
    - charge of the projectile (multipla of the elementary charge)
    - nucleon number of the projectile
    - density [g/cm^3] of the scintillator
    - material specific Z/A

    too_slow (boolean), is the Bethe-Bloch eq applicable for that energy?

    OUTPUT:
    - kinetic energy per nucleon [MeV/A]
    - LET [MeV/cm]
    '''
    use_shell_correction = True
    dEdx_MeV_cm2_g = np.zeros(len(E_MeV_per_A))
    for idx, E_MeV in enumerate(E_MeV_per_A):

        params = [E_MeV, z_projectile, A_projectile, Z_A_scintillator, use_shell_correction]
        dEdx_MeV_cm2_g[idx], too_slow = BetheBlochEq(*params)

    LET_MeV_cm = dEdx_MeV_cm2_g*density_g_cm3
    return E_MeV_per_A, LET_MeV_cm


def scintillator_parameters(scint_name):
    '''
    INPUT:
    - Scintillator name (see dictionary)

    OUTPUT:
    Returns an array with the scintillator's:
    - Decay time [s] (primary component)
    - Light yield [photons/MeV] (dictionary in percent of anthracene)
    - Density [g/cm^3]
    of the scintillator.

    See e.g. https://www.crystals.saint-gobain.com/products/scintillating-fiber
    '''

    Z_A_scintillator = 0.5555 # water
    anthracene = 17400 # photons per MeV

    global scint_dic
    scint_dic = {
        'BCF-12' : [3.2e-9, 0.46*anthracene, 1.05, Z_A_scintillator],
        'BCF-60' : [7e-9,   0.41*anthracene, 1.05, Z_A_scintillator],
        'BC-400' : [2.4e-9, 0.65*anthracene, 1.03, Z_A_scintillator],
        'Pilot-U': [1.4e-9, 0.76*anthracene, 1.023, Z_A_scintillator],
    }

    error_text = "\n\n\tScintillator '%s' not implemented" % scint_name
    error_text += "\n\tAvaliable scintillators: %s\n" % [i for i in scint_dic.keys()]
    assert scint_name in scint_dic.keys(), error_text

    return scint_dic[scint_name]



# def core_radius_cm(beta, density_g_cm3):
#     '''
#     INPUT:
#     - beta = v/c
#     - scintillator density [g/cm^3]
#
#     OUTPUT:
#     - core radius [cm]
#     '''
#     r_core_cm = 11.6e-7
#     rMin_cm = r_core_cm *beta
#     return rMin_cm/density_g_cm3


def track_structure_parameters(track_structure_model, E_MeV_per_A, z_projectile, A_projectile, density_g_cm3, Z_A_scintillator = 0.555):
    '''
    Calculates the track structure radii for the ion track

    INPUT:
    - track structure model
    - kinetic energy per nucleon [MeV/A]
    - charge of the projectile (multipla of the elementary charge)
    - nucleon number of the projectile
    - density [g/cm^3] of the scintillator
    - Z/A ratio of the scintillator

    OUTPUT:
    - core radius [cm]
    - penumbra radius [cm]
    - LET [MeV/cm]
    '''

    E_MeV_per_A, LET_MeV_cm = stopping_power(E_MeV_per_A, z_projectile, A_projectile, density_g_cm3, Z_A_scintillator)
    beta = E_MeV_to_beta(E_MeV_per_A, A_projectile)

    # rMin_cm = core_radius_cm(beta, density_g_cm3)

    r_core_cm = 11.6e-7
    if track_structure_model == "Gaussian":
        r0 = 0.5e-6 #cm
        b_cm = 2*r0/np.sqrt(np.pi) # from Birks (1964)
        rMin_cm = b_cm*np.ones(len(LET_MeV_cm))
        rMax_cm = -1*np.ones(len(LET_MeV_cm))
    elif track_structure_model == "Scholz_Kraft":
        gamma = 0.05e-4
        delta = 1.7
        rMax_cm = gamma*E_MeV_per_A**delta
        rMin_cm = r_core_cm*np.ones(len(E_MeV_per_A))
    elif track_structure_model == 'Chatterjee_Schaefer':
        rMax_um = 0.768*E_MeV_per_A - 1.925*np.sqrt(E_MeV_per_A) + 1.257
        rMax_cm = rMax_um*1e-4 # to cm
        rMin_cm = r_core_cm *beta
    else:
        print("Track structure model unknown")

    rMax_cm /= density_g_cm3
    rMin_cm /= density_g_cm3

    return rMin_cm, rMax_cm, LET_MeV_cm


def integrate_signal(emissionResults):
    '''
    Integrate the fluorescence emission as a function of time

    INPUT:
    - results from solving the PDE

    OUTPUT:
    - integrated fluorescence
    '''
    _, time_s, fluorescence_emission, dt, n_tries = emissionResults
    return trapz(fluorescence_emission, time_s), dt, n_tries


def Blanc_model_parameters(track_structure):
    '''
    Returns an array with
    1) diffusion constant [cm^2/s]
    2) quenching parameter alpha [cm^3/s]
    3) quenching parameter beta [cm^6/s]

    beta = 0. in this current (lacks experimental data for high-LET)

    See
        Christensen JB and Andersen CE (2018), Phys. Med. Biol. (63) 195010
    for details
    '''

    diff_cm2_s = 5.0e-4
    track_structure_params = {
        "Gaussian" :            [diff_cm2_s, 4.5e-9, 0.0e-25],
        "Scholz_Kraft" :        [diff_cm2_s, 9.0e-8, 0.0e-25],
        "Chatterjee_Schaefer" : [diff_cm2_s, 1.5e-9, 0.0e-31]
    }

    TSnames = track_structure_params.keys()
    error_text = "\n\n\tTrack structure '%s' not implemented" % track_structure
    error_text += "\n\tAvaliable models: %s\n" % [i for i in TSnames]
    assert track_structure in TSnames, error_text

    return track_structure_params[track_structure]
