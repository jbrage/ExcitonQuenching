import matplotlib.pyplot as plt
import numpy as np
import os, sys
sys.path.append('./cython')
from evolveDensitiesCython import PDEsolver
from parameters import scintillator_parameters, track_structure_parameters, \
    Blanc_parameters, integrate_signal

def initialise_cython(scintillator, track_structure, LET_MeV_cm, print_info = True):

    # load scintillator specific parameters
    tau_s, A_MeV, density_g_cm3 = scintillator_parameters(scintillator, print_info)
    N_0 = A_MeV * LET_MeV_cm # linear exciton density

    # load track structure parameters
    TS_parameters = track_structure_parameters(track_structure, LET_MeV_cm, density_g_cm3, print_info)
    E_MeV, rMin_cm, rMax_cm, voxelSize, gridSize = TS_parameters

    # load the parameters for the Blanc model
    D_diff, alpha_bimol, beta_trimol = Blanc_parameters(track_structure, print_info)

    if print_info:
        print("\n# ===============================")
        print("# Energy = %0.3g MeV" % E_MeV)
        print("# LET = %0.3g MeV/cm" % LET_MeV_cm)
    # print("# biQuenchFactor factor = %0.3g" % biQuenchFactor)
    # print("# Ion diff factor = %0.3g*e-4" % ion_diff_factor)


    # find an appropriate time step, upper limit defined by the von Neumann criterion
    SCALE_TIME = 1.0
    emissionResults = [-1]
    while emissionResults[0] < 0:

        SCALE_TIME *= 1.05
        # stop if the calculation takes too long:
        if SCALE_TIME > 1000:
            # print ("# fucked up")
            return -1, SCALE_TIME

        emissionResults = PDEsolver(track_structure,
                                    N_0,
                                    D_diff,
                                    gridSize,
                                    voxelSize,
                                    alpha_bimol,
                                    tau_s,
                                    rMax_cm,
                                    rMin_cm,
                                    print_info,
                                    SCALE_TIME
                                    )

    scint_response = integrate_signal(emissionResults)

    # ==========================================================================
    # Calculate the fluorescence of max(LET) with \biQuenchFactor = 0
    # for a linear respone
    reference_results = PDEsolver(track_structure,
                                N_0,
                                D_diff,
                                gridSize,
                                voxelSize,
                                0., # alpha_bimol
                                tau_s,
                                rMax_cm,
                                rMin_cm,
                                print_info,
                                SCALE_TIME
                                )


    ref_signal = integrate_signal(reference_results)
    QCF = ref_signal / scint_response
    return QCF




if __name__ == "__main__":

    print_info = False
    scintillator = "BCF-12"
    track_structure = "Gaussian"
    track_structure = "Scholz_Kraft"
    track_structure = "Chatterjee_Schaefer"
    LET_MeV_cm_list = [30, 40, 50]
    QCF_list = []

    for LET_MeV_cm in LET_MeV_cm_list:
        QCF = initialise_cython(scintillator, track_structure, LET_MeV_cm, print_info)
        QCF_list.append(QCF)

    ptext ="""# Scintillator: %s
# Track structure model: %s
# LET [MeV/cm], QCF
""" % (scintillator, track_structure)

    os.system("mkdir -p results")
    fname = "results/%s_%s.txt" % (scintillator, track_structure)
    with open(fname, "w") as outfile:
        print("# LET,\tQCF: \n# ============")
        outfile.write(ptext)
        # outfile.write("# Scintillator: %s\n" % scintillator)
        # outfile.write("# Track structure model: %s\n" % track_structure)
        # outfile.write("# LET [MeV/cm], QCF\n")
        for idx, LET_MeV_cm in enumerate(LET_MeV_cm_list):
            QCF = QCF_list[idx]
            print_text = "%0.3g,\t%0.3g" % (LET_MeV_cm, QCF)
            save_text = "%0.3g,%0.3g\n" % (LET_MeV_cm, QCF)
            print(print_text)
            outfile.write(save_text)


    # BCF-12 response in proton beams from Wang et al (2012)
    # fpath = "../../Experimental_data/BCF12_all_datapoints.dat"
    # exp_data = np.genfromtxt(fpath, delimiter = ",")
    # exp_LET, exp_signal, exp_sigmas = exp_data[:,0], exp_data[:,1], exp_data[:,2]
    #
    # PRINT = True
    # PRINT = False
    #
    # factor = 1
    # # LET_list = np.asarray(exp_LET[::factor
    # LET_list = np.linspace(5, 70, 10)
    # scintillator = "BC-400"
    # scintillator = "BCF-531"
    #
    # TS = ["Gaussian"]
    # TS = ["Scholz_Kraft"]
    # TS = ["Chatterjee_Schaefer"]
    # TS = ["Gaussian", "Scholz_Kraft", "Chatterjee_Schaefer"]
    #
    # dirname = "../vary_diff_alpha_cythonized/optimised_results"
    # for track_structure in TS:
    #     print(track_structure)
    #     #fname = "optimised_results/%s_parameters.dat" % (track_structure)
    #     #data = np.genfromtxt(fname, delimiter =",")
    #     #alpha, ion_diff_factor, chi2, best_scale, best_scale_sigma = data[0], data[1], data[2], data[3], data[4]
    #     fname = "%s/%s_parameters.dat" % (dirname, track_structure)
    #     data = np.genfromtxt(fname, delimiter =",")
    #     alpha, ion_diff_factor =  data[0], data[1]
    #
    #     print("D = %g" % ion_diff_factor)
    #     print("alpha = %g" % alpha)
    #
    #     w_dir = "."
    #     time_stamp = ""
    #     # LET_list = np.linspace(min(exp_LET)*0.9, max(exp_LET)*1.1, 10, endpoint = True)
    #     if PRINT:
    #         print("\n%s track structure model." % track_structure)
    #
    #     with open("response_%s_%s.txt" % (scintillator, track_structure), 'w') as outfile:
    #
    #         for idx, LET_MeV_cm in enumerate(LET_list):
    #             response = initialise_cython(scintillator, alpha, ion_diff_factor, track_structure, LET_MeV_cm, PRINT, w_dir, time_stamp)[0]
    #             print(LET_MeV_cm, response)
    #             outfile.write("%g,%g\n" % (LET_MeV_cm, response))
