from __future__ import division
import numpy as np
cimport numpy as np
# from math import exp
from libc.math cimport exp, sqrt, M_PI as pi, log

DTYPE = np.double
ctypedef np.double_t DTYPE_t

cimport cython
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(True) # turn off checks for zero division

def PDEsolver( str track_structure,
               double N_0,
               double ion_diff,
               long gridSize,
               double voxelSize,
               double alpha,
               double decayConst_tau,
               double rMax_cm,
               double rMin_cm,
               bint PRINT,
               double SCALE_TIME
               ):


    computation_time_steps = -1 # calculate atomatically
    multiply_tau = 3 # multipla of the decay time (i.e. simulation time)
    scale_decay_time = multiply_tau*2
    # supported track structure model?
    assert(track_structure in ["Gaussian", "Scholz_Kraft", "Chatterjee_Schaefer"])

    # preallocate variables and arrays
    cdef np.ndarray[DTYPE_t, ndim=1] cArray = np.zeros(gridSize)
    cdef np.ndarray[DTYPE_t, ndim=1] cArray_temp = np.zeros(gridSize)
    cdef double distance_cm, ion_density, initialised_cc = 0.
    cdef int i, mid_array
    mid_array = int(gridSize/2.)

    cdef double r0, b_cm, initialised_core = 0., initialised_penumbra = 0.
    cdef double Gaussian_factor, decayTime, flourescenceRate, quenchingRate
    cdef double SK_const, CS_const, core, penumbra, dist_squared

    decayTime = 1./decayConst_tau
    fluorescenceRatio = 0.5
    flourescenceRate = decayTime*fluorescenceRatio
    quenchingRate = decayTime*(1-fluorescenceRatio)

    # calculate the densities in the 1-D array based on the chosen track structure model
    if track_structure == "Gaussian":
        r0 = 0.5e-6 #cm
        b_cm = 2*r0/sqrt(pi)
        Gaussian_factor = N_0/(pi*b_cm*b_cm)
        for i in range(gridSize):
            distance_cm = abs(i - mid_array)*voxelSize
            ion_density = Gaussian_factor * exp( - distance_cm*distance_cm/(b_cm*b_cm))

            cArray[i] = ion_density
            initialised_cc += ion_density
        initialised_cc *= voxelSize

    else:
        if track_structure == "Scholz_Kraft":
            SK_const = N_0/(pi*(1. + 2.*log(rMax_cm/rMin_cm)))
            core = SK_const
            penumbra = SK_const
        else:
            # Chatterjee-Schaefer track structure
            CS_const = log(sqrt(exp(1))*rMax_cm/rMin_cm)
            core = N_0/(2*pi) + N_0/(4*pi) / CS_const
            penumbra = N_0/(4*pi) / CS_const

        initialised_penumbra = 0.
        initialised_core = 0.

        for i in range(gridSize):
            distance_cm = abs(i - mid_array)*voxelSize
            if distance_cm <= rMin_cm:
                dist_squared = rMin_cm*rMin_cm
                cArray[i] = core/dist_squared
                initialised_core += core/dist_squared
            elif (distance_cm > rMin_cm and distance_cm <= rMax_cm):
                dist_squared = distance_cm*distance_cm
                cArray[i] = penumbra/dist_squared
                initialised_penumbra += penumbra/dist_squared
            else:
                cArray[i] = 0.

        initialised_core *= voxelSize
        initialised_penumbra *= voxelSize
        initialised_cc = initialised_core + initialised_penumbra

    # calculate appropriate time step based on the diffusion and grid spacing
    cdef double dt = 1., sx = 0., cx = 0.
    cdef bint von_neumann_expression = False

    while not von_neumann_expression:
        dt /= 1.01
        # as defined in the Deghan (2004) paper
        sx = ion_diff*dt/(voxelSize*voxelSize)
        cx = 0.
        # check von Neumann's criterion
        von_neumann_expression = (2*sx + cx*cx <= 1 and cx*cx <= 2*sx)

    # dt /= 100.

    dt /= SCALE_TIME
    # dt /= 2.

    if computation_time_steps < 0:
        computation_time_steps = int(scale_decay_time*decayConst_tau/dt)
    if computation_time_steps == 0:
        computation_time_steps = 100

    if PRINT and SCALE_TIME < 1.06:
        if alpha > 0:
            print("# Time step, dt = %0.3E s" % (dt))
            print("# Total iterations = %0.3E" % (float(computation_time_steps)*gridSize))
        else:
            print("# Reference computation ...")

    # preallocate arrays and variables
    cdef np.ndarray[DTYPE_t, ndim=1] time_list = np.zeros(computation_time_steps)
    cdef np.ndarray[DTYPE_t, ndim=1] fluorescence_list = np.zeros(computation_time_steps)
    # cdef np.ndarray[DTYPE_t, ndim=1] uniMolQuench_list = np.zeros(computation_time_steps)
    # cdef np.ndarray[DTYPE_t, ndim=1] biMolQuench_list = np.zeros(computation_time_steps)
    # cdef np.ndarray[DTYPE_t, ndim=1] triMolQuench_list = np.zeros(computation_time_steps)

    cdef double countUniMolQuenched, countBiMolQuenched, countTriMolQuenched
    cdef double countEmittedFluorescence
    cdef double cArray_temp_voxel = 0., fluorescence = 0., quenching = 0.
    cdef double bimolecularQuench = 0., trimolecularQuench = 0., check_NaN, all_losses
    cdef int time_step, quarter_length
    quarter_length = int(gridSize/4.)

    for time_step in range(computation_time_steps):
        countUniMolQuenched = 0.
        countBiMolQuenched = 0.
        countTriMolQuenched = 0.
        countEmittedFluorescence = 0.

        for i in range(1, gridSize-1):

            cArray_temp_voxel = sx*(cArray[i-1] + cArray[i+1])
            cArray_temp_voxel += (1.- 2.*sx) * cArray[i]

            fluorescence = cArray[i]*flourescenceRate*dt
            quenching = cArray[i]*quenchingRate*dt
            bimolecularQuench = alpha*cArray[i]*cArray[i]*dt
            # trimolecularQuench = beta*cArray[i]*cArray[i]*cArray[i]*dt

            # losses to fluorescence and quenching
            countEmittedFluorescence += fluorescence
            countUniMolQuenched += quenching
            countBiMolQuenched += bimolecularQuench
            # countTriMolQuenched += trimolecularQuench

            all_losses = fluorescence + quenching + bimolecularQuench
            # all_losses = fluorescence + quenching + bimolecularQuench + trimolecularQuench

            if all_losses > cArray_temp_voxel:
                return [-1]
            else:
                cArray_temp[i] = cArray_temp_voxel - all_losses

        time_list[time_step] = time_step*dt
        fluorescence_list[time_step] = countEmittedFluorescence
        # uniMolQuench_list[time_step] = countUniMolQuenched
        # biMolQuench_list[time_step] = countBiMolQuenched
        # triMolQuench_list[time_step] = countTriMolQuenched

        # update the array before next iteration
        for i in range(1, gridSize-1):
            cArray[i] = cArray_temp[i]

    return [initialised_cc, time_list, fluorescence_list]
