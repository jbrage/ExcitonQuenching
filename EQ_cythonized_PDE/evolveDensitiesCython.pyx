from __future__ import division
import numpy as np
import sys
cimport numpy as np
# from math import exp
from libc.math cimport exp, sqrt, M_PI as pi, log, abs

DTYPE = np.double
ctypedef np.double_t DTYPE_t

cimport cython
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(False) # turn off checks for zero division


def find_time_step_dt(double diff_cm2_s, double voxelSize_cm):
    '''
    Von Neumann scheme as defined in:
    Dehghan M (2004), Appl. Math. Comput. 150(1), 5-19

    INPUT:
    - Exciton diffusion parameter [cm^2/s]
    - distance [cm] between to voxels
    - parameter to reduce the time step if previous step was too large

    OUTPUT:
    - time step dt [s]
    - Lax-Wendroff parameter s_x (see Daghan-paper referenced above)
    '''

    cdef double dt = 1., sx = 0., cx = 0.
    cdef bint von_neumann_expression = False

    while not von_neumann_expression:
        dt /= 1.1
        sx = diff_cm2_s*dt/(voxelSize_cm*voxelSize_cm)
        cx = 0. # not used in curent version
        # check von Neumann's criterion
        von_neumann_expression = (2*sx + cx*cx <= 1 and cx*cx <= 2*sx)
    return dt, sx



def calculate_time_steps(double decay_time_tau_s, double dt):
    '''
    Calculates the number of time steps

    INPUT:
    - decay time of the scintillator [s]
    - length of a time step [s]

    OUTPUT:
    - number of iterations in time
    '''
    cdef double scale_factor = 5
    cdef size_t computation_time_steps, min_number_steps = 100
    computation_time_steps = int(scale_factor*decay_time_tau_s/dt)

    if computation_time_steps < min_number_steps:
        computation_time_steps = min_number_steps
    return computation_time_steps


def Gaussian_distribution(int nVoxelsArray, double voxelSize_cm, double N_0):
    '''
    Gaussian distribution
    INPUT:
    - grid size (number of voxels in the 1D array)
    - distance [cm] between to voxels
    - linear exciton density [exictons per unit length]

    OUTPUT:
    - 1D array with a centered Gaussian distribution
    - number of exitons
    '''
    cdef np.ndarray[DTYPE_t, ndim=1] exctionArray = np.zeros(nVoxelsArray)
    cdef double nExcitons = 0.0, distance_cm, ExcitonDensity
    cdef double r0_cm = 0.5e-6, b_cm = 2*r0_cm/sqrt(pi)
    cdef double Gaussian_factor = N_0/(pi*b_cm*b_cm)
    cdef size_t voxel_i, mid_array = int(nVoxelsArray/2.)

    for voxel_i in range(nVoxelsArray):
        distance_cm = abs(voxel_i - mid_array)*voxelSize_cm
        ExcitonDensity = Gaussian_factor * exp( - distance_cm*distance_cm/(b_cm*b_cm))

        exctionArray[voxel_i] = ExcitonDensity
        nExcitons += ExcitonDensity
    return exctionArray, nExcitons


def amorphous_track_structure_model_distribution(str track_structure_model,
                                            int nVoxelsArray,
                                            double voxelSize_cm,
                                            double N_0,
                                            double core_radius_cm,
                                            double penumbra_radius_cm
                                           ):
    '''
    Distributes the excitons according to the selected track structure model
    and evovles the system according to the Lax-Wendroff scheme.

    If the time step is too large, the function calls itself with a smaller
    time step

    INPUT:
    - track structure name
    - grid size (number of voxels in the 1D array)
    - distance [cm] between to voxels
    - linear exciton density [exictons per unit length]
    - core and penumra radii [cm]

    OUTPUT:
    - 1D array with a centered Gaussian distribution
    - number of exitons
    '''
    cdef np.ndarray[DTYPE_t, ndim=1] exctionArray = np.zeros(nVoxelsArray)
    cdef double nExcitons_penumbra = 0., nExcitons_core = 0.
    cdef double nExcitons, distance_cm, SK_const, CS_const
    cdef size_t voxel_i, mid_array = int(nVoxelsArray/2.)

    # calculate core and penumbral exciton densities
    # (or simply return the Gaussian)
    if track_structure_model == "Scholz_Kraft":
        SK_const = N_0/(pi*(1. + 2.*log(penumbra_radius_cm/core_radius_cm)))
        core = SK_const
        penumbra = SK_const
    elif track_structure_model == "Chatterjee_Schaefer":
        CS_const = log(sqrt(exp(1))*penumbra_radius_cm/core_radius_cm)
        core = N_0/(2*pi) + N_0/(4*pi) / CS_const
        penumbra = N_0/(4*pi) / CS_const
    else:
        return Gaussian_distribution(nVoxelsArray, voxelSize_cm, N_0)

    # insert the ion track exciton densities in the array
    for voxel_i in range(nVoxelsArray):
        distance_cm = abs(voxel_i - mid_array)*voxelSize_cm
        if distance_cm <= core_radius_cm:
            # core region
            exctionArray[voxel_i] = core/(core_radius_cm*core_radius_cm)
            nExcitons_core += core/(core_radius_cm*core_radius_cm)
        elif (distance_cm > core_radius_cm and distance_cm <= penumbra_radius_cm):
            # penumbral region
            exctionArray[voxel_i] = penumbra/(distance_cm*distance_cm)
            nExcitons_penumbra += penumbra/(distance_cm*distance_cm)
        else:
            # larger than track radius
            exctionArray[voxel_i] = 0.

    nExcitons = nExcitons_core + nExcitons_penumbra
    return exctionArray, nExcitons


def PDEsolver( str track_structure_model,
               double N_0,
               double core_radius_cm,
               double penumbra_radius_cm,
               double decay_time_tau_s,
               int n_tries,
               double diff_cm2_s,
               double alpha = 0.,
               double beta = 0.,
               double dt = -1
               ):

    '''
    Solve the Blanc model equation subject to the initial condition
    given by the ion and track structure model

    INPUT:
    - track structure model
    - linear exciton density [exictons per unit length]
    - ion penumbral radius [cm]
    - ion core radius [cm]
    - decay time of the scintillator [s]
    - grid size (number of voxels in the 1D array)
    - distance [cm] between to voxels
    - exciton diffusion constant [cm^2/s]
    - bimolecular quenching parameter [cm^3/s]
    - trimolecular quenching parameter [cm^6/s]
    - parameter to scale the time step dt if necessary


    OUTPUT:
    - number of excitons initialized in the array
    - array with the time [s] for each step
    - array with the relative fluorescence emission for each step
    '''

    # grid values from stability tests
    cdef double voxelSize_cm = 3e-8
    cdef int nVoxelsArray = int(1.0e-5/voxelSize_cm)
    
    if track_structure_model == 'Chatterjee_Schaefer':
        voxelSize_cm = 0.8e-8
        nVoxelsArray = int(1.0e-5/voxelSize_cm)

    # cdef double voxelSize_cm = 10e-9
    # cdef int nVoxelsArray = int(5e-5/voxelSize_cm)

    # preallocate variables and arrays
    cdef np.ndarray[DTYPE_t, ndim=1] exctionArray = np.empty(nVoxelsArray)
    cdef np.ndarray[DTYPE_t, ndim=1] exctionArray_temp = np.empty(nVoxelsArray)
    cdef double nExcitons, events_per_s, flourescenceRate, quenchingRate
    cdef size_t i

    cdef double countUniMolQuenched, countBiMolQuenched, countTriMolQuenched, fluorescenceCounter
    cdef double n_Excitons, fluorescenceRatio = 0.5
    cdef double exctionArray_temp_voxel = 0., fluorescence = 0.,
    cdef double unimolQuenching, bimolecularQuench, trimolecularQuench, ExcitonLoss
    cdef int time_step

    # Calculate the time step dt
    cdef double sx
    cdef size_t computation_time_steps
    if dt > 0:
        sx = diff_cm2_s*dt/(voxelSize_cm*voxelSize_cm)
    else:
        dt, sx = find_time_step_dt(diff_cm2_s, voxelSize_cm)
    computation_time_steps = calculate_time_steps(decay_time_tau_s, dt)

    # preallocate the arrays for fluorescence emission as a function of time
    cdef np.ndarray[DTYPE_t, ndim=1] fluorescence_emission = np.empty(computation_time_steps)
    cdef np.ndarray[DTYPE_t, ndim=1] time_s = np.empty(computation_time_steps)

    # fluoresence and quenching rates
    events_per_s = 1./decay_time_tau_s
    flourescenceRate = events_per_s*fluorescenceRatio
    quenchingRate = events_per_s*(1-fluorescenceRatio)

    # get the Exciton densities and number of initialised excitons
    exctionArray, nExcitons = amorphous_track_structure_model_distribution(
                                            track_structure_model,
                                            nVoxelsArray,
                                            voxelSize_cm,
                                            N_0,
                                            core_radius_cm,
                                            penumbra_radius_cm
                                            )

    # iterate through the array
    for time_step in range(computation_time_steps):

        fluorescenceCounter = 0.
        for i in range(1, nVoxelsArray-1):
            # evolve according to the Lax-Wendroff scheme
            n_Excitons = exctionArray[i]

            exctionArray_temp_voxel = sx*(exctionArray[i-1] + exctionArray[i+1])
            exctionArray_temp_voxel += (1.- 2.*sx) * n_Excitons

            fluorescence = n_Excitons*flourescenceRate*dt
            unimolQuenching = n_Excitons*quenchingRate*dt
            bimolecularQuench = alpha*n_Excitons*n_Excitons*dt
            trimolecularQuench = beta*n_Excitons*n_Excitons*n_Excitons*dt

            fluorescenceCounter += fluorescence
            ExcitonLoss = fluorescence + unimolQuenching + bimolecularQuench + trimolecularQuench

            # removing more exictons than avaliable?
            # return -1 => restarts function with a different time step dt
            if ExcitonLoss > exctionArray_temp_voxel:
                n_tries += 1

                if n_tries > 20:
                    error_msg = "# Voxel is too large for the LET. \n# ... Exiting."
                    sys.exit(error_msg)

                # try recursively again with a smaller time step
                PDEsolver(track_structure_model,
                                N_0,
                                core_radius_cm,
                                penumbra_radius_cm,
                                decay_time_tau_s,
                                n_tries,
                                diff_cm2_s,
                                alpha,
                                beta,
                                dt = dt/2.
                                )

            else:
                exctionArray_temp[i] = exctionArray_temp_voxel - ExcitonLoss

        time_s[time_step] = time_step*dt
        fluorescence_emission[time_step] = fluorescenceCounter

        # update the array before next iteration
        for i in range(1, nVoxelsArray-1):
            exctionArray[i] = exctionArray_temp[i]

    return [nExcitons, time_s, fluorescence_emission, dt, n_tries]
