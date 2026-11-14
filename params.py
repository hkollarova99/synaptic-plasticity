"""
Regularly spiking (RS) neuron according to Izhikevich
Nearest-neighbour, presynaptically centered, additive metaSTDP
copyright Lubica Benuskova

Constants and model parameters.

Cross references:
  C:     reference/plastic-rs.c  (original C code this is a port of)
  rs.c:  later version of the same C code (STIM_INTENSITY, 0.05 weights)
  Paper: Benuskova & Abraham (2007) J Comput Neurosci 22:129-133
"""

# --- simulation control -------------------------------------------------------
# C:25 used 100000 ms (= 1.7 min). With SAMPLE_PERIOD = 1 min this gave only two
# plot points (t = 0 and t = 1 min), and theta_M (TAU = 30 s) had barely settled.
# 30 min gives 31 plot points. rs.c used 25200000 (= 7 h, the paper's time scale).
SIMULATION_PERIOD = 1_800_000   # total # of time steps of the simulation in ms

TIME_STEP = 1               # C:27  time step in ms (model equations assume 1 ms, see simulation.cycle_cell)
SAMPLE_PERIOD = 60000       # C:28  time step between sample points for plots
MINUTES = 60000.0           # C:29  to get time in minutes in the format of real number

AP_PEAK = 30.0              # C:31  peak of action potential for the given neuron
P = 0.007                   # C:32  probability P=0.00X corresponds to X Hz average frequency of homogenous Poisson series
MAX = 10000                 # C:33  max number of presynaptic spikes between 2 postsynaptic spikes = 10000

TAU = 30000                 # C:35  Integration period for theta_M in ms (Paper sec. 2.3 quotes 60 s)
BOLTZMAN = 0                # C:36  should a Boltzman average be used
ALPHA = 2000.0              # C:37  theta_M is scaled by this parameter (Paper: alpha = 2000)

# Input type. C:223-252 uncorrelated inputs are active; C:254-277 correlated
# (synchronous) inputs are commented out. rs.c uses the correlated ones.
CORRELATED_INPUT = False

# --- weight initialization (Initial bias) ---------------------------------------
# FIX: C:47-48 used weights of 5 (mV per input spike). With these the cell never
# fires (output rate 0 Hz), so STDP/BCM never run. The plotting code (C:396-397)
# and the commented-out hard limits (C:331-332) assume w(0) = 0.05. Following
# rs.c and Paper sec. 2.1 (Eq. 4: I = h*w*N, w(0) = 0.05, N = 150), each input
# spike now adds weight * STIM_INTENSITY = 7.5 mV.
WEIGHT_1_INIT = 0.05
WEIGHT_2_INIT = 0.05
STIM_INTENSITY = 150        # rs.c: number of engaged input fibers (Paper: N_mpp = N_lpp = 150)

# Optional hard weight limits (C:331-332 / C:355-356, commented out in C).
# Paper sec. 2.3 used 70-160 % of the initial value. Off by default to match C.
USE_WEIGHT_LIMITS = False
WEIGHT_MIN_FACTOR = 0.70
WEIGHT_MAX_FACTOR = 1.60

# --- STDP rule parameters (C:51-54, Paper sec. 2.3) ----------------------------
TAU_LTP = 20.0      # LTP tau in ms
TAU_LTD = 100.0     # LTD tau in ms
A0_LTP = 0.02       # initial LTP peak at t=0
A0_LTD = 0.01       # initial LTD peak at t=0

# --- Izhikevich's neuron parameter values for RS cell (C:57-65) ----------------
# c = -69 and d = 2 differ from Izhikevich's standard RS values (-65, 8) on
# purpose: Paper sec. 2.1 uses these.
COEFF1 = 0.04
COEFF2 = 5
COEFF3 = 140
COEFF4 = 1

A = 0.02
B = 0.2
C = -69.0
D = 2

# --- output data files used for plots (C:161-166) ------------------------------
FILE_INPUT_1 = "input1_pulse_times.dat"
FILE_INPUT_2 = "input2_pulse_times.dat"
FILE_CELL_VOLTAGE = "cell_voltage.dat"
FILE_WEIGHT_1 = "weight_1.dat"
FILE_WEIGHT_2 = "weight_2.dat"
FILE_THETA_M = "theta_M.dat"
