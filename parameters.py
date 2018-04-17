import numpy as np
from math import exp, pi, sqrt
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import os

def E_MeV_beta_from_LET(LET_list):
    # libamtrack calculated
    path_data = "input_data/E_LET_beta_libamtrack.dat"
    data = np.genfromtxt(path_data, delimiter = ",")
    d_LET = data[:,0]
    d_E_MeV = data[:,1]
    d_beta = data[:,2]

    # check the LET is withtin the range before interpolation
    if LET_list < min(d_LET) or LET_list > max(d_LET):
        error_text = "\nThe LET value is outside the range of [%0.3g, %0.3g] MeV/cm" \
        % (min(d_LET), max(d_LET))
        raise ValueError(error_text)

    interpolated_energy = interp1d(d_LET, d_E_MeV, kind = 'cubic')
    energy_list = interpolated_energy(LET_list)

    interpolated_beta = interp1d(d_LET, d_beta , kind = 'cubic')
    beta_list = interpolated_beta(LET_list)
    return energy_list, beta_list


def scintillator_parameters(scint_name, print_info):
    # decay time [s], light yield relative to anthracene, and density [g/cm^3]
    # see e.g. https://www.crystals.saint-gobain.com/products/scintillating-fiber

    scint_dic = {'BCF-12' : [3.2e-9, 0.46, 1.05],
                 'BC-400' : [2.4e-9, 0.65, 1.03]
                }

    error_text = "\n\nScintillator %s is not (yet) included! \
        \nCheck out: %s" % (scint_name, scint_dic.keys())
    assert scint_name in scint_dic.keys(), error_text
    tau_s, A_percentage, density_g_cm3 = scint_dic[scint_name]

    # convert the light yield to photons per MeV
    anthracene = 17400 # photons per MeV
    A_MeV = A_percentage * anthracene

    if print_info:
        print("# Using the %s scintillator" % scint_name)

    return tau_s, A_MeV, density_g_cm3


def track_structure_parameters(name, LET_MeV_cm, density_g_cm3, print_info):
    TS_names = ["Gaussian", "Scholz_Kraft", "Chatterjee_Schaefer"]
    error_text = "\n\nThe track structure %s is not supported. \
        \nCheck out: %s" % (name, TS_names)
    assert name in TS_names, error_text

    E_MeV, beta = E_MeV_beta_from_LET(LET_MeV_cm)

    if name == "Gaussian":
        r0 = 0.5e-6 #cm
        b_cm = 2*r0/sqrt(pi) # from Birks (1964)
        rMin_cm = b_cm
        rMax_cm = -1
    elif name == "Scholz_Kraft":
        gamma = 0.05e-4
        delta = 1.7
        rMin_cm = 10.0e-7
        rMax_cm = gamma*E_MeV**delta
    else:
        r_core_cm = 11.6e-7 # r_core = 11.6 nm
        rMin_cm = beta*r_core_cm
        # Chatterjee-Schaefer r max:
        rMax_cm = 0.768*E_MeV - 1.925*sqrt(E_MeV) + 1.257
        rMax_cm *= 1e-4 # to cm

    rMin_cm /= density_g_cm3
    rMax_cm /= density_g_cm3
    print("r max = ", rMax_cm)
    print("r min = ", rMin_cm)

    if name == "Gaussian":
        voxelSize = b_cm/40. # [cm] distance between two neighbouring voxels
        gridSize = int(b_cm*20./voxelSize)
    else:
        if LET_MeV_cm < 30:
            voxelSize = rMin_cm/5. # [cm] distance between two neighbouring voxels
        else:
            voxelSize = rMin_cm/10. # [cm] distance between two neighbouring voxels
        gridSize = int(2*rMax_cm*0.02/voxelSize)

    if print_info:
        print("# with the %s track structure model" % name)

    return E_MeV, rMin_cm, rMax_cm, voxelSize, gridSize


def integrate_signal(emissionResults):
    initialised_cc, time, signal = emissionResults
    integration = trapz(signal, time)
    return integration


def Blanc_parameters(track_structure, print_info):

    path = "../../fff.py"
    if os.path.isfile(path):
        Blanc_parameters = [0, 0, 0]
        if print_info:
            print("\n# Loading the Blanc parameters from: \n# %s" % path)
    else:
        # [diffusion, bimolecular quenching, trimolecular quenching]
        quench_dic = {
            "Gaussian" : [4.28532,1.3413,0.0910405],
            "Scholz_Kraft" : [5,4.5068,0.266882],
            "Chatterjee_Schaefer" : [5.97,0.606107,0.00419259]
            }
        Blanc_params = np.asarray(quench_dic[track_structure])
        if print_info:
            print("# Loading the tabulated Blanc parameters")

    # diffusion in units of [1e-5 cm^2/s]
    Blanc_params[0] *= 1e-5
    # bimol quench in units of [1e-9 cm^3/s]
    Blanc_params[1] *= 1e-9
    # trimol quench in units of [1e-25 cm^6/s]
    Blanc_params[2] *= 1e-25
    return Blanc_params
