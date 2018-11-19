import numpy as np
from main import getQCF


if __name__ == '__main__':
    '''
    Running an exmaple with the BCF-12 scintillator exposed to
    helium ions at different energies per nucleon using the Scholz-Kraft
    track structure model

    Implement
        a new scintillator in
        - scintillator_parameters() in functions.py

        a new track structure model in
        - track_structure_parameters() in functions.py
    '''

    scintillator = "BCF-12"
    track_structure_name = "Scholz_Kraft"

    z_projectile, A_projectile = 2, 4 # helium
    E_MeV_per_A = np.linspace(2, 300, 10)


    LET_MeV_cm, QCFs = getQCF(scintillator, track_structure_name, E_MeV_per_A, z_projectile, A_projectile)

    print("# Particle: z = {:d}, A = {:d}. \n# Scintillator: {} ".format(z_projectile, A_projectile, scintillator))
    print("# E [MeV/A], LET [MeV/cm], QCFs")

    for (E, LET, QCF) in zip(E_MeV_per_A, LET_MeV_cm, QCFs):
        print("  {:9.3g},{:13.3g},{:6.4g}".format(E, LET, QCF))
