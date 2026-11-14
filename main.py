"""
Regularly spiking (RS) neuron according to Izhikevich
Nearest-neighbour, presynaptically centered, additive metaSTDP
copyright Lubica Benuskova

Port of reference/plastic-rs.c. "C:<line>" comments refer to that file.
"""

import random
import time

from params import (
    FILE_CELL_VOLTAGE, FILE_INPUT_1, FILE_INPUT_2, FILE_THETA_M, FILE_WEIGHT_1, FILE_WEIGHT_2,
    SIMULATION_PERIOD,
)
from simulation import SimulationState, run_network


def get_time() -> int:
    """Returns the microseconds part of the current real clock time,
    used for setting up the new seed of the random number generator (C:406-413, tv_usec)."""
    return time.time_ns() // 1000 % 1_000_000


def main() -> None:
    # each time the code is run, it starts with a different seed for the
    # random number generator, set according to real clock
    random.seed(get_time())  # C:157-159  srandom(z)

    # C:161-175  open the data files; they are closed automatically at the end
    with open(FILE_INPUT_1, "w") as f_input_1, \
            open(FILE_INPUT_2, "w") as f_input_2, \
            open(FILE_CELL_VOLTAGE, "w") as f_cell_voltage, \
            open(FILE_WEIGHT_1, "w") as f_weight_1, \
            open(FILE_WEIGHT_2, "w") as f_weight_2, \
            open(FILE_THETA_M, "w") as f_theta_M:
        state = SimulationState({
            "input_1": f_input_1,
            "input_2": f_input_2,
            "cell_voltage": f_cell_voltage,
            "weight_1": f_weight_1,
            "weight_2": f_weight_2,
            "theta_M": f_theta_M,
        })
        run_network(SIMULATION_PERIOD, state)  # C:168  all computations happen here


if __name__ == "__main__":
    main()
