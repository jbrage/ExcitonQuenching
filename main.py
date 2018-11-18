import numpy as np
import sys
sys.path.append('./EQ_cythonized_PDE')
from evolveDensitiesCython import PDEsolver

from EQ_model_functions.functions import scintillator_parameters, track_structure_parameters, \
    Blanc_model_parameters, integrate_signal


def getQCF(scintillator_name, track_structure_name, E_MeV_per_A, z_projectile, A_projectile):
    '''
    Calculate the quenching correction factor

    INPUT:
    - scintillator name
    - track structure model
    - kinetic energy per nucleon [MeV/A]
    - projectile charge (multipla of the elementary charge)
    - nucleon number of the projectile

    OUTPUT:
    - array with the quenching correction factors for each specified energy
    '''

    # scintillator parameters
    scint_decaytime_s, light_yield, density_g_cm3, Z_A_scintillator = scintillator_parameters(scintillator_name)

    # Blanc model parameters
    diff_cm2_s, alpha_cm3_s, beta_cm6_s = Blanc_model_parameters(track_structure_name)

    # projectile track parameters
    rMin_cm, rMax_cm, LET_MeV_cm = track_structure_parameters(
                                                            track_structure_name,
                                                            E_MeV_per_A,
                                                            z_projectile,
                                                            A_projectile,
                                                            density_g_cm3,
                                                            Z_A_scintillator
                                                        )

    N_0 = light_yield * LET_MeV_cm # linear exciton density
    QCF_array = np.empty(len(E_MeV_per_A))

    for idx, (N0, rmin, rmax) in enumerate(zip(N_0, rMin_cm, rMax_cm)):
        n_tries = 1 # number of times the functions recursively calls itself to find a time step dt

        emissionResults = PDEsolver(track_structure_name,
                                        N0,
                                        rmin,
                                        rmax,
                                        scint_decaytime_s,
                                        n_tries,
                                        diff_cm2_s,
                                        alpha_cm3_s,
                                        beta_cm6_s
                                    )
        scintillator_response, dt, n_tries = integrate_signal(emissionResults)


        # reference calculation using the same time step as above
        reference_results = PDEsolver(track_structure_name,
                                        N0,
                                        rmin,
                                        rmax,
                                        scint_decaytime_s,
                                        n_tries,
                                        diff_cm2_s,
                                        dt = dt
                                    )

        reference_signal, dt, n_tries = integrate_signal(reference_results)

        QCF = reference_signal / scintillator_response
        QCF_array[idx] = QCF

    return LET_MeV_cm, QCF_array



if __name__ == "__main__":

    scintillator = "BCF-12"
    track_structure_name = "Scholz_Kraft"

    z_projectile, A_projectile = 2, 4
    E_MeV_per_A = np.linspace(10, 200, 3)

    getQCF(scintillator, track_structure_name, E_MeV_per_A, z_projectile, A_projectile)
