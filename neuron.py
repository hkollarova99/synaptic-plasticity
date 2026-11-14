"""Neuron structure and its creation / initialization.

Port of the "Neuron Stuff" section of reference/plastic-rs.c (C:84-153).
"""

import math
from dataclasses import dataclass

from params import B, BOLTZMAN, C, TAU, TIME_STEP, WEIGHT_1_INIT, WEIGHT_2_INIT


@dataclass
class Neuron:
    """struct _neuron (C:88-102)."""
    # variables describing the neuron
    v: float            # membrane voltage in mV
    u: float            # membrane recovery variable
    output: int         # binary output
    # weights of neuron
    weight_1: float     # input 1 weight
    weight_2: float     # input 2 weight
    # BCM averaging stuff
    theta_M: float      # dynamic modification LTD-LTP threshold
    mem: float          # average of postsynaptic activity
    pf: float           # used to calculate mem
    boltzfactor: float  # used to calculate mem


def make_neuron_struct() -> Neuron:
    """Creates a new neuron and initializes it (C:110-129)."""
    return Neuron(
        v=C,
        u=B * C,
        output=0,
        weight_1=WEIGHT_1_INIT,
        weight_2=WEIGHT_2_INIT,
        theta_M=1.0,
        mem=0.0,
        pf=1.0,
        boltzfactor=0.0,
    )


def init_neuron(n: Neuron) -> None:
    """Sets up the BCM running average (C:131-138).

    With BOLTZMAN == 0, pf = 1/(1 - boltzfactor) is a fixed point of
    pf <- pf*boltzfactor + 1 (see simulation.bcm_threshold). mem is then an
    exponential moving average of the output with time constant TAU.
    """
    n.boltzfactor = math.exp(-(float(TIME_STEP) / float(TAU)))
    if BOLTZMAN == 0:
        n.pf = 1.0 / (1 - n.boltzfactor)
    else:
        n.pf = 1.0


def make_neuron() -> Neuron:
    """Creates and initializes a new neuron (C:141-149)."""
    n = make_neuron_struct()
    init_neuron(n)
    return n
