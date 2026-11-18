# Excited Ultralight Halos around Black Holes

Code accompanying the paper

> S. Bottaro, V. Firsov, G. Perez, H. Steuer, O. Telem,
> *Capturing Ultralight Dark Matter Around Black Holes*, arXiv:XXXX.XXXXX

## Citation

If you use this code, please cite the paper above.

## License

See [`LICENSE`](LICENSE).


The pipeline has three steps:

1. **`ME_Numerics/`** (C++): Monte Carlo evaluation of the second-order matrix elements (diagrams) for transitions between bound and unbound states of the gravitational atom. Produces log files with the value and error of every matrix element for a set of focusing parameters `xi`.
2. **`ODE_Numerics/`** (C++): reads those log files, builds the coefficients of the rate equations for the occupation numbers `N_1s`, `N_2s`, `N_2p`, and solves them on a grid in `(xi, f_a)`.
3. **Notebooks** (Python): make the figures from the output of step 2.

The outputs of our own runs of steps 1 and 2 are included, so the figures can be reproduced without rerunning the expensive steps.

```
.
├── ME_Numerics/              matrix elements (step 1)
│   ├── src/main.cpp          diagram templates, xi values, MC settings
│   ├── src/coulomb.hpp       Coulomb (unbound) wavefunctions
│   ├── submit.sh             submission to an LSF cluster as a job array
│   └── CMakeLists.txt
├── ODE_Numerics/             rate equations (step 2)
│   ├── src/main.cpp
│   ├── logs/logs_<date>/     our matrix-element results (output of step 1)
│   ├── *_solve_ODE_logs_*.txt  our ODE results on the (xi, f_a) grid
│   ├── *_solve_ODE_plot_*.txt  our ODE time evolution at selected points
│   └── CMakeLists.txt
├── logs_ODE_solutions.ipynb  parameter-space figures (step 3)
└── logs_ODE_plots.ipynb      time-evolution figures (step 3)
```

## Requirements

- **ME_Numerics:** C++17 compiler, CMake >= 3.16, Boost (header-only special functions).
- **ODE_Numerics:** C++20 compiler, CMake >= 3.16, [Matplot++](https://github.com/alandefreitas/matplotplusplus) (needs gnuplot to display plots) and [indicators](https://github.com/p-ranav/indicators).
- **Notebooks:** Python 3 with `numpy`, `matplotlib`, `tqdm` and `matplotlib-label-lines`.

## 1. Matrix elements (`ME_Numerics`)

Build:

```bash
cd ME_Numerics
cmake -S . -B build
cmake --build build -j
```

Each run of the executable computes **one diagram template** for all values of `xi` listed in `vector_xi` in `src/main.cpp`, using 8 threads. The template is selected by the environment variable `LSB_JOBINDEX` (1-based index into `vector_diagram_template`, which currently has 130 entries). Results are appended to a file `simulation_<timestamp>.log_<index>` in the directory `log_dir` set in `src/main.cpp` (the directory must exist).

**On an LSF cluster**, `submit.sh` submits all templates as a job array (set `LOG_DIR` in `submit.sh` to the same directory as `log_dir`):

```bash
./submit.sh            # all templates
./submit.sh 17         # only template 17
./submit.sh 40-82      # a range
EASY_MODE=1 ./submit.sh   # quick test with 1000x less statistics
```

**Locally**, for a quick test of a single template:

```bash
LSB_JOBINDEX=1 EASY_MODE=1 ./build/Viktor_MC_integral
```

Use `EASY_MODE=0` for full statistics (`submit.sh` does this by default).

**Diagram templates** are strings of the form `a+b-c=d=a+b-c`. Bound states are written as `nlm` (with `/` standing for `m = -1`, e.g. `21/` is the 2p state with `m = -1`), unbound states as `k1, k2, ...`, and the middle entry `d` is the intermediate state connecting the two vertices. For every template and `xi`, the log file contains the template, the values of `alpha` and `xi`, the matrix element as `(re,im)`, and its MC error.

### Our results

Our matrix-element log files are in `ODE_Numerics/logs/logs_<date>/`, where the ODE code expects them. `logs_2026_09_17` is the latest complete run of all 130 templates; the other folders contain earlier runs and reruns. When the same template and `xi` appear several times, the ODE code keeps the result with the smallest error. Folders whose names do not start with `logs` (e.g. `_logs_2026_09_15`) are ignored.

## 2. Rate equations (`ODE_Numerics`)

Build in the same way:

```bash
cd ODE_Numerics
cmake -S . -B build
cmake --build build -j
```

Run from inside `ODE_Numerics`, since the code reads the matrix elements from `./logs/logs*/`:

```bash
./build/solve_ODE
```

The settings are in `main()` in `src/main.cpp`:

- `astro_regime`: the central object, `"BH"` (10 solar masses at 1 kpc from the Galactic Centre), `"IBH"` (10^3 solar masses), `"SBH"` (10^5 solar masses), `"Sun"` or `"star"` (no absorption);
- `vector_xi`, `vector_f`: the grid in `xi` and `f_a` (in GeV);
- `plotting_regime`: `"false"` (no plot), `"true"` (show a plot) or `"log"` (also write the time evolution to a `*_solve_ODE_plot_*.txt` file);
- `ASSUMED_DECAY_PHASE` (near the top of the file): `1.0` for the phases giving maximal decay, as used in the paper.

The results at the end of the evolution are written to `<astro_regime>_solve_ODE_logs_<timestamp>.txt` in the working directory.

### Our results

- `BH_solve_ODE_logs_2026-10-04_10-38-11.txt`: 10 solar-mass BH, used for the halo-density figure in the main text and the halo-mass figure in the appendix.
- `IBH_solve_ODE_logs_2026-10-04_19-26-27.txt`: 10^3 solar-mass BH, used for the intermediate-mass BH figure in the appendix.
- `SBH_*`, `Sun_*`, `star_*`: additional runs not shown in the paper.
- `BH_solve_ODE_plot_*.txt`: time evolution at selected `(xi, f_a)` points, used for the halo-density-over-time figure.

## 3. Figures

- `logs_ODE_solutions.ipynb` reads the `*_solve_ODE_logs_*.txt` files and makes the parameter-space figures (`capturing parameter space.pdf`, `BH_density_parameter_space.pdf`, `BH_mass_parameter_space.pdf`, `IBH_mass_parameter_space.pdf`, and the corresponding plots for the other regimes).
- `logs_ODE_plots.ipynb` reads the `*_solve_ODE_plot_*.txt` files and makes `density_over_time.pdf`.

The notebooks were written for Google Colab. To run them locally, remove the `drive.mount(...)` cell and set `file_path = 'ODE_Numerics/' + filename`.
