# synaptic-plasticity

Python port of C simulation of a **regularly spiking (RS)
Izhikevich neuron** with two inputs. Its synapses learn through
**nearest-neighbour, presynaptically centred STDP**, and the LTP/LTD amplitudes
are scaled by a **BCM sliding modification threshold θ_M** (metaplasticity).

The model follows:

> Benuskova L., Abraham W.C. (2007). *STDP rule endowed with the BCM sliding
> threshold accounts for hippocampal heterosynaptic plasticity.*
> J Comput Neurosci 22:129–133. doi:10.1007/s10827-006-0002-x

The original C code is in [`reference/plastic-rs.c`](reference/plastic-rs.c).
Comments of the form `C:<line>` in the Python source point to the matching
line there.

## Requirements

- Python 3.10+ (tested with 3.13)
- Standard library only, no extra packages

## Running

```bash
python main.py
```

The program prints the mean postsynaptic firing rate, for example:

```
Average output frequency	0.870556
```

and writes six tab-separated data files to the current directory:

| File | Columns | Content |
|---|---|---|
| `input1_pulse_times.dat` | time (ms), 1 | presynaptic spike times at input 1 |
| `input2_pulse_times.dat` | time (ms), 1 | presynaptic spike times at input 2 |
| `cell_voltage.dat` | time (min), 1 | postsynaptic spike times |
| `weight_1.dat` | time (min), % change | weight 1 relative to its initial value, sampled every minute |
| `weight_2.dat` | time (min), % change | weight 2 relative to its initial value, sampled every minute |
| `theta_M.dat` | time (min), θ_M | BCM modification threshold, sampled every minute |

The random seed comes from the system clock, so every run is different.

## Project layout

| File | Purpose | Original C |
|---|---|---|
| `main.py` | seeds the RNG, opens the output files, starts the simulation | C:155-177, C:406-413 |
| `params.py` | all constants and model parameters | C:25-65 |
| `neuron.py` | `Neuron` data class and its initialization | C:84-153 |
| `simulation.py` | main loop, Poisson inputs, Izhikevich update, STDP, BCM, sampling | C:179-401 |

## Model summary

On each 1 ms step (`run_network`):

1. **Inputs** (`spontaneous_input`): each input fires a Poisson spike with
   probability `P` = 0.007, i.e. 7 Hz.
2. **Neuron** (`cycle_cell`): Izhikevich's model, paper Eq. 1–4, with
   a = 0.02, b = 0.2, c = −69, d = 2. Each input spike adds
   `weight × STIM_INTENSITY` mV.
3. **STDP** (`stdp_train_cell`): runs on every postsynaptic spike. Each
   presynaptic spike since the previous postsynaptic spike is paired with
   both its nearest neighbours:
   - LTP: `A+ · exp(−Δt/20 ms)`
   - LTD: `A− · exp(Δt/100 ms)`
   - update: `w ← w · (1 + ltp − ltd)` (paper Eq. 5–7)
4. **BCM** (`bcm_threshold`): θ_M = α · ⟨c⟩. ⟨c⟩ is an exponential running
   average of the postsynaptic spikes over `TAU` = 30 s, and α = 2000.
   The LTP and LTD amplitudes depend on it as A+ = A+(0)/θ_M and
   A− = A−(0)·θ_M (paper Eq. 8–10).

## Configuration (`params.py`)

| Parameter | Default | Meaning |
|---|---|---|
| `SIMULATION_PERIOD` | 1 800 000 | simulated time in ms (30 min). `rs.c` used 25 200 000 (7 h) |
| `WEIGHT_1_INIT`, `WEIGHT_2_INIT` | 0.05 | initial weights (paper sec. 2.1) |
| `STIM_INTENSITY` | 150 | number of input fibres engaged |
| `P` | 0.007 | per-ms spike probability of each input |
| `CORRELATED_INPUT` | `False` | `True` makes both inputs fire synchronously |
| `USE_WEIGHT_LIMITS` | `False` | `True` keeps weights within 70–160 % of their initial value (paper sec. 2.3) |
| `TAU`, `ALPHA` | 30000, 2000 | BCM averaging period (ms) and scaling |

`TIME_STEP` must stay at 1 ms because the equations are not scaled by dt.

## Differences from `plastic-rs.c`

The port keeps the original structure and order of operations. These changes
were made to fix bugs or inconsistencies:

1. **Synaptic weights.** The C code started with weights of 5 mV per spike.
   At that strength the neuron never fired (0 Hz), so no plasticity happened.
   The C plotting code (`100*(w-0.05)/0.05`) and the commented-out weight
   limits both assume w(0) = 0.05. As in the later `rs.c` and in the paper,
   the weights now start at 0.05 and are multiplied by 150 input fibres,
   i.e. 7.5 mV per input spike.
2. **Uninitialized STDP variables.** In C, `ltp`, `ltd` and `A_ltp` are
   declared without a value (C:310), so an old value could be reused when the
   pairing conditions were false. In C this is undefined behaviour. They are
   now reset to 0 for every presynaptic spike.
3. **Simulation length.** 30 min instead of 100 s. With one sample per
   minute, the C length produced only two plot points.
4. **Small changes.**
   - The mean frequency is divided by the period actually simulated.
   - Plot percentages use `WEIGHT_*_INIT` instead of the literal 0.05.
   - The correlated input generator and the hard weight limits, which are
     commented out in C, can be switched on in `params.py`.

These values still differ from the paper and were left as in the C code:

| Parameter | C code | Paper / `rs.c` |
|---|---|---|
| `TAU` | 30 s | 60 s |
| spike peak / firing threshold | 30 mV | 55 mV peak, 24 mV threshold |
| spontaneous input | uncorrelated 7 Hz | mostly correlated 8 Hz |

## Example results

Default 30 min run:

| t (min) | weight 1 (%) | weight 2 (%) | θ_M |
|---|---|---|---|
| 0 | 0.0 | 0.0 | 0.00 |
| 5 | −93.5 | +96.0 | 1.38 |
| 10 | −99.5 | +98.5 | 1.53 |
| 20 | −100.0 | +101.3 | 2.87 |
| 30 | −100.0 | +100.9 | 1.53 |

Mean output frequency: 0.87 Hz. The paper's target spontaneous rate is about 1 Hz.

With two independent inputs, the synapses compete: one weight drops to zero
and the other roughly doubles. Which one wins varies from run to run. The
paper notes that the sliding threshold "is not able to prevent synaptic
decrease to zero", which is why it used the 70–160 % limits.

Results from other settings (30 min runs):

| Setting | Output | Final weights |
|---|---|---|
| `USE_WEIGHT_LIMITS = True` | about 0.6 Hz | both near +60 % (the upper limit) |
| `CORRELATED_INPUT = True` | about 0.88 Hz | both stay within 0–3 % |
