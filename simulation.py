"""Network stuff: inputs, neuron update, STDP training, BCM threshold and sampling.

Port of the "Network Stuff" section of reference/plastic-rs.c (C:179-401).
"C:<line>" comments refer to line numbers in that file. "Paper" refers to
Benuskova & Abraham (2007) J Comput Neurosci 22:129-133.
"""

import math
import random

from neuron import Neuron, make_neuron
from params import (
    A, A0_LTD, A0_LTP, ALPHA, AP_PEAK, B, C, COEFF1, COEFF2, COEFF3, COEFF4, CORRELATED_INPUT, D,
    MAX, MINUTES, P, SAMPLE_PERIOD, STIM_INTENSITY, TAU_LTD, TAU_LTP, TIME_STEP,
    USE_WEIGHT_LIMITS, WEIGHT_1_INIT, WEIGHT_2_INIT, WEIGHT_MAX_FACTOR, WEIGHT_MIN_FACTOR,
)


class SimulationState:
    """Holds the C global variables (C:67-78) and the open data files (C:40)."""

    def __init__(self, files):
        self.input_1 = 0        # input = 0 if there is no presynaptic spike
        self.input_2 = 0        # input = 1 if there is a presynaptic spike
        self.t1_cell = 0        # time of the previous postsynaptic spike
        self.t2_cell = 0        # time of the following postsynaptic spike
        self.index_1 = 0        # cycle index for input 1
        self.index_2 = 0        # cycle index for input 2
        self.t_input_1 = [0] * MAX          # arrays for times of presynaptic spikes
        self.t_input_2 = [0] * MAX
        self.input_1_thetaM = [0.0] * MAX   # arrays for current values of theta_M
        self.input_2_thetaM = [0.0] * MAX
        self.last_sample = 0    # sampling of variables for plotting
        self.mean_out_frequency = 0.0
        # rs.c: number of stimulated input fibers (Paper Eq. 4: N_mpp, N_lpp)
        self.n_input_1 = STIM_INTENSITY
        self.n_input_2 = STIM_INTENSITY

        # file variables for data used for plots
        self.f_input_1 = files["input_1"]
        self.f_input_2 = files["input_2"]
        self.f_cell_voltage = files["cell_voltage"]
        self.f_weight_1 = files["weight_1"]
        self.f_weight_2 = files["weight_2"]
        self.f_theta_M = files["theta_M"]

    def reset_arrays(self):
        """Array (re)initialization (C:191-192 and C:364-365)."""
        self.t_input_1 = [0] * MAX
        self.t_input_2 = [0] * MAX
        self.input_1_thetaM = [0.0] * MAX
        self.input_2_thetaM = [0.0] * MAX


def run_network(t: int, s: SimulationState) -> None:
    """Main simulation loop (C:183-212)."""
    # The update equations in cycle_cell() and the per-step spike probability P
    # are written for 1 ms steps, as in the C code (Paper: "We update the neuron
    # every 1 ms"). They are not scaled by TIME_STEP.
    if TIME_STEP != 1:
        raise ValueError("TIME_STEP must be 1 ms: the model equations are not scaled by dt")

    generate_input = spontaneous_input_correlated if CORRELATED_INPUT else spontaneous_input

    single_cell = make_neuron()

    # initialization of arrays (C:191-192)
    s.reset_arrays()

    time = 0
    while time <= t:
        generate_input(time, single_cell, s)      # C:197  generation of inputs at time t
        cycle_cell(time, single_cell, s)          # C:198  update the cell according to the current inputs

        if time == s.t2_cell and s.t1_cell > 0:
            stdp_train_cell(single_cell, s)       # C:200  train the cell according to the current inputs and output

        bcm_threshold(time, single_cell)          # C:202  update theta_M for the next cycle
        test_sample(time, single_cell, s)         # C:203  monitor variables

        time = time + TIME_STEP                   # C:205  increase time step by 1

    # C:208 divided by SIMULATION_PERIOD. This uses the period actually simulated
    # (identical when called from main.py).
    s.mean_out_frequency = s.mean_out_frequency / (float(t) / 1000.0)
    print("Average output frequency\t", end="")   # C:209-210
    print(f"{s.mean_out_frequency:f}")


def _record_input_1(time: int, n: Neuron, s: SimulationState) -> None:
    """Stores a presynaptic spike at input 1 for the next STDP update (C:230-236)."""
    if s.t1_cell > 0:
        if s.index_1 < MAX and time > s.t2_cell:
            s.t_input_1[s.index_1] = time               # times of presynaptic spikes at input_1
            s.input_1_thetaM[s.index_1] = n.theta_M     # values of theta_M when there is spike at input_1
            s.index_1 = s.index_1 + 1


def _record_input_2(time: int, n: Neuron, s: SimulationState) -> None:
    """Stores a presynaptic spike at input 2 for the next STDP update (C:245-251)."""
    if s.t1_cell > 0:
        if s.index_2 < MAX and time > s.t2_cell:
            s.t_input_2[s.index_2] = time               # times of presynaptic spikes at input_2
            s.input_2_thetaM[s.index_2] = n.theta_M     # values of theta_M when there is spike at input_2
            s.index_2 = s.index_2 + 1


def spontaneous_input(time: int, n: Neuron, s: SimulationState) -> None:
    """Generates spontaneous, totally uncorrelated Poisson inputs (C:217-252)."""
    s.input_1 = 0
    s.input_2 = 0  # C:221  initialization of inputs

    # C:224  rnd = random()/RAND_MAX, a random number from interval [0, 1)
    rnd = random.random()
    if P > rnd:
        s.input_1 = 1
        s.f_input_1.write(f"{float(time):f}\t{s.input_1:d}\n")  # C:228  times of presynaptic spikes at input_1
        _record_input_1(time, n, s)

    rnd = random.random()  # C:239
    if P > rnd:
        s.input_2 = 1
        s.f_input_2.write(f"{float(time):f}\t{s.input_2:d}\n")  # C:243  times of presynaptic spikes at input_2
        _record_input_2(time, n, s)


def spontaneous_input_correlated(time: int, n: Neuron, s: SimulationState) -> None:
    """Correlated, i.e. synchronous Poisson inputs (C:254-277, commented out in C).

    Selected with params.CORRELATED_INPUT = True (this is the default in rs.c).
    """
    s.input_1 = 0
    s.input_2 = 0

    rnd = random.random()
    if P > rnd:
        s.input_1 = 1
        s.input_2 = 1
        s.f_input_1.write(f"{float(time):f}\t{s.input_1:d}\n")
        s.f_input_2.write(f"{float(time):f}\t{s.input_2:d}\n")

        # C:264-275  both inputs are stored together
        if s.t1_cell > 0:
            if s.index_1 < MAX and s.index_2 < MAX and time > s.t2_cell:
                s.t_input_1[s.index_1] = time
                s.t_input_2[s.index_2] = time

                s.input_1_thetaM[s.index_1] = n.theta_M
                s.input_2_thetaM[s.index_2] = n.theta_M

                s.index_1 = s.index_1 + 1
                s.index_2 = s.index_2 + 1


def cycle_cell(time: int, n: Neuron, s: SimulationState) -> None:
    """Updates the neuron output according to its inputs (C:285-303).

    Izhikevich model, Paper Eq. 1-3, forward Euler with dt = 1 ms.
    """
    n.output = 0

    if n.v >= AP_PEAK:  # C:289  spike was emitted in the previous step -> reset (Paper Eq. 3)
        if s.t1_cell == 0:
            s.t1_cell = time
            s.t2_cell = time
        else:
            s.t1_cell = s.t2_cell
            s.t2_cell = time
        n.v = C
        n.u = n.u + D
        s.mean_out_frequency = s.mean_out_frequency + 1
    else:
        # C:296  dv/dt = 0.04 v^2 + 5 v + 140 - u (Paper Eq. 1)
        n.v = n.v + COEFF1 * (n.v * n.v) + COEFF2 * n.v + COEFF3 - COEFF4 * n.u
        # C:297-298  synaptic input I(t) = h(t) w N (Paper Eq. 4).
        # FIX: C adds only the weight. The weight times the number of stimulated
        # fibers is taken from rs.c, which is the convention used for w(0) = 0.05.
        if s.input_1 > 0:
            n.v = n.v + n.weight_1 * s.n_input_1
        if s.input_2 > 0:
            n.v = n.v + n.weight_2 * s.n_input_2
        if n.v >= AP_PEAK:  # C:299  capping of output spike peak and setting the variable output to 1
            n.v = AP_PEAK
            n.output = 1
        n.u = n.u + A * (B * n.v - n.u)  # C:300  du/dt = a (b v - u) (Paper Eq. 2)


def _stdp_update(weight: float, weight_init: float, t_input: list, input_thetaM: list,
                 s: SimulationState) -> float:
    """Nearest-neighbour, presynaptically centred STDP with BCM-scaled
    amplitudes for one input (C:315-337 for input 1, C:339-361 for input 2).

    Every presynaptic spike t_pre stored between the previous (t1_cell) and the
    following (t2_cell) postsynaptic spike is paired with both of them:
      LTP: dt = t2_cell - t_pre > 0   (pre before post, Paper Eq. 5)
      LTD: dt = t1_cell - t_pre < 0   (post before pre, Paper Eq. 6)
    and the weight is updated as w <- w (1 + ltp - ltd) (Paper Eq. 7).
    """
    index = 0
    while index < MAX:
        t_pre = t_input[index]
        if t_pre == 0:  # C:316 / C:340  no more presynaptic spikes at this input
            break
        theta = input_thetaM[index]

        # FIX: C:310 declares ltp, ltd, A_ltp, A_ltd without initializing them,
        # so a value from the previous presynaptic spike (or from the other input)
        # could be reused when an "if" below is false. In C this is undefined
        # behaviour, so here every presynaptic spike starts from 0.
        ltp = 0.0
        ltd = 0.0
        A_ltp = 0.0

        # C:320-321  Paper Eq. 8-9: A+(t) = A+(0) / theta_M,  A-(t) = A-(0) * theta_M
        if theta != 0:
            A_ltp = A0_LTP * (1.0 / theta)  # update of LTP peak according to the current theta_M
        A_ltd = A0_LTD * theta              # update of LTD peak according to the current theta_M

        if s.t2_cell > t_pre:
            ltp = A_ltp * math.exp(-(s.t2_cell - t_pre) / TAU_LTP)  # C:324  pre-post pairing
        if s.t1_cell < t_pre:
            ltd = A_ltd * math.exp((s.t1_cell - t_pre) / TAU_LTD)   # C:326  post-pre pairing

        weight = weight * (1 + ltp - ltd)  # C:328  overall update of the weight
        weight = _limit_weight(weight, weight_init)  # C:331-332  optional hard limits

        index = index + 1
    return weight


def _limit_weight(weight: float, weight_init: float) -> float:
    """Optional hard limits on weights (C:331-332 / C:355-356, commented out in C).

    Paper sec. 2.3 used 70-160 % of the initial value. Enabled with
    params.USE_WEIGHT_LIMITS.
    """
    if not USE_WEIGHT_LIMITS:
        return weight
    return min(max(weight, WEIGHT_MIN_FACTOR * weight_init), WEIGHT_MAX_FACTOR * weight_init)


def stdp_train_cell(n: Neuron, s: SimulationState) -> None:
    """STDP training (C:308-369)."""
    s.index_1 = 0
    s.index_2 = 0  # C:313
    n.weight_1 = _stdp_update(n.weight_1, WEIGHT_1_INIT, s.t_input_1, s.input_1_thetaM, s)  # C:315-337
    n.weight_2 = _stdp_update(n.weight_2, WEIGHT_2_INIT, s.t_input_2, s.input_2_thetaM, s)  # C:339-361

    # resetting (C:364-367)
    s.reset_arrays()
    s.index_1 = s.index_2 = 0
    s.t1_cell = s.t2_cell  # important: time of the following postsynaptic spike becomes the time of the preceding one


def bcm_threshold(time: int, n: Neuron) -> None:
    """BCM theta_M calculation (C:375-382, Paper Eq. 10).

    mem is a running average of the squared (binary) output over roughly the
    last TAU ms, and theta_M = alpha * <c>.
    """
    n.mem = (n.mem * (n.pf - 1.0) + (n.output * n.output)) / n.pf
    n.pf = n.pf * n.boltzfactor + 1.0
    n.theta_M = ALPHA * n.mem


def test_sample(time: int, n: Neuron, s: SimulationState) -> None:
    """Tests the cell every so often (C:387-401)."""
    if n.output > 0:
        s.f_cell_voltage.write(f"{time / MINUTES:f}\t{1:1d}\n")  # C:390

    # saving data point for plots every minute, i.e. 60000 ms (C:393-399)
    if time == 0 or time - s.last_sample >= SAMPLE_PERIOD:
        s.last_sample = time
        # C:396-397  % change with respect to the initial weight. C hard-codes
        # 0.05, which equals WEIGHT_*_INIT now (see params.py).
        s.f_weight_1.write(f"{time / MINUTES:f}\t{100 * (n.weight_1 - WEIGHT_1_INIT) / WEIGHT_1_INIT:6.3f}\n")
        s.f_weight_2.write(f"{time / MINUTES:f}\t{100 * (n.weight_2 - WEIGHT_2_INIT) / WEIGHT_2_INIT:6.3f}\n")
        s.f_theta_M.write(f"{time / MINUTES:f}\t{n.theta_M:f}\n")  # C:398
