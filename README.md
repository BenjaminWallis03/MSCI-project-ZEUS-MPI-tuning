# MSCI-project-ZEUS-MPI-tuning
Monte Carlo tuning of photoproduction MPI parameters in Pythia8 using ZEUS 2007 and ZEUS 2001 data, Rivet analyses, and Professor optimisation.

---

## Environment

All simulations, event generation, and tuning procedures in this project were performed within a ROOT-enabled Linux terminal environment.  
This ensured compatibility between Pythia8, Rivet, and Professor, and provided a consistent framework for compiling and executing the analysis code.

---

## Software Versions Used

- Pythia8: 8.315
- Rivet: v4.1.2
- Professor: 2.5.6

---

## Methodology

This project consists of a full Monte Carlo tuning workflow using Pythia8, Rivet, and Professor. The procedure can be divided into several key stages:

---

### 1. Construction of Run Cards and Event Generation Framework

The first step was to define the event generation setup using Pythia8.

Custom run cards (`.run.cmnd`) were created for each dataset:
- ZEUS 2007 dijet
- ZEUS 2007 multijet
- ZEUS 2001
- Combined dataset configuration

These were based on standard Pythia8 templates, with modifications including:
- Beam configuration for photoproduction (γp collisions)
- Phase space cuts (e.g. `PhaseSpace:pTHatMin`)
- MPI parameter settings (for later tuning)

In parallel, a custom event generation executable (`nano_generate.cc`) was implemented. This program:
- Reads the `.run.cmnd` configuration file
- Initialises Pythia8 with the specified settings
- Converts generated events into HepMC3 format
- Passes events to Rivet for analysis
- Produces output in YODA format (`rivet.yoda`)

The generator was compiled using:

```
g++ scripts/nano_generate.cc -o nano_generate \
  $(pythia8-config --cxxflags --ldflags) \
  $(rivet-config --cxxflags --ldflags)
```

---

### 2. Validation of Event Generation

Before performing any tuning, the correctness of the simulation setup was verified.

This was done by studying the effect of the minimum hard scattering scale (`pTHatMin`) on the ZEUS 2007 dijet dataset.

Using the dijet run card, simulations were performed with:
- `pTHatMin = 8, 10, 12, 15 GeV`
- followed by a finer scan: `2, 4, 6, 8, 10 GeV`

The resulting Rivet distributions were compared to ZEUS data.

This step served two purposes:
- Ensuring that the generator + Rivet pipeline was functioning correctly
- Understanding how phase space cuts affect agreement with data

It was observed that higher `pTHatMin` values suppress low transverse energy regions, leading to mismodelling in MPI-sensitive observables.

---

### 3. Single-Parameter Sensitivity Studies

To understand the role of individual MPI parameters, dedicated scans were performed using the ZEUS 2007 multijet configuration.

Each parameter was varied independently while keeping others fixed:

- `MultipartonInteractions:pT0Ref`
- `MultipartonInteractions:ecmPow`
- `MultipartonInteractions:expPow`
- `BeamRemnants:primordialKT`
- `StringPT:sigma`

These studies aimed to identify:
- Which observables are sensitive to each parameter
- The qualitative effect of parameter variation on distributions

However, it was later identified that the chosen `pTHatMin` values were too high for this dataset. As a result:
- The low-energy region (most sensitive to MPI) was underpopulated
- The apparent parameter sensitivity in this region was artificially suppressed

This limitation is important when interpreting the results of these scans.

---

### 4. Professor Tuning Procedure

A full multi-parameter tuning was performed using Professor (Prof2).

Three key MPI parameters were selected:
- `pT0Ref`
- `ecmPow`
- `expPow`

These parameters control the regularisation scale and energy scaling of Multiple Parton Interactions, and are expected to have the strongest impact on MPI-sensitive observables.

---

#### Sampling the Parameter Space

For each dataset (dijet, multijet, combined, and ZEUS 2001), a parameter space was defined using a `ranges.dat` file.

Professor was then used to generate a set of parameter points within this space using:

```
prof2-sample ranges.dat -n 64 -o scan -t run.cmnd.in
```

where:
- `ranges.dat` defines the allowed parameter ranges
- `-n 64` specifies the number of sampled points
- `scan/` is the output directory containing generated run configurations
- `run.cmnd.in` is the template Pythia8 run card

This step produced a set of run cards sampling the parameter space using a Sobol sequence (default Professor sampling method).

---

#### Monte Carlo Generation

For each sampled parameter point:
- A corresponding run card was generated
- 200,000 events were simulated using Pythia8
- Rivet was used to analyse the events and produce YODA histograms

This resulted in a training dataset mapping:

```
(parameter values) → (observable distributions)
```

---

#### Interpolation

The sampled Monte Carlo results were interpolated using:

```
prof2-ipol scan/ ipol.dat
```

This step constructs polynomial response surfaces that approximate the dependence of each observable bin on the input parameters.

For specific datasets (e.g. combined tuning), separate interpolation files were created, e.g.:

```
prof2-ipol scan_combined/ ipol_combined.dat
```

---

#### Parameter Tuning

The interpolated response surfaces were then used to perform a χ² minimisation against experimental data:

```
prof2-tune ipol.dat -d /work/tune_zeus2007/refdata
```

where:
- `/work/tune_zeus2007/refdata` contains the ZEUS experimental reference data

For weighted fits (used in the combined tuning), the procedure included bin weighting:

```
prof2-tune ipol_2007_combined.dat \
  -d /work/tune_zeus2007/refdata \
  -w weights_2007_clean.dat \
  -o tune_weighted.dat
```

This step determines the optimal parameter values by minimising the χ² between the interpolated Monte Carlo predictions and the experimental measurements.

---

#### Output

The output of the tuning procedure is a set of optimised MPI parameters (the “tune”)

---

### 5. Validation of Tuned Parameters

The tuned parameter sets obtained from Professor were then validated.

For each tune:
- A new Pythia8 run was performed using the optimised parameters
- Rivet was used to generate updated distributions
- These were directly compared to ZEUS data

This was done separately for:
- ZEUS dijet observables
- ZEUS multijet observables
- Combined dataset results
- Cross-validation between different tunes

---

### 6. Plotting and Comparison

Final results were visualised using Rivet’s plotting tools.

This included:
- Overlaying Monte Carlo predictions with ZEUS data
- Comparing different parameter tunes on the same distributions
- Identifying regions of agreement and disagreement

These plots form the basis of the analysis presented in the report.

---

### Summary of Workflow

1. Define run cards and generator setup  
2. Validate simulation behaviour (pTHatMin studies)  
3. Perform single-parameter sensitivity scans  
4. Run Professor tuning (multi-parameter optimisation)  
5. Validate tuned parameter sets  
6. Produce final comparison plots

---

## Reproducibility Notes

- All run configurations are stored in `cards/`
- Generator code is provided in `scripts/`
- Analysis is performed using Rivet

---

## Future Work

- Sensitivity-weighted combined tuning  
- Inclusion of NLO matrix elements  
- Improved photon matter profile modelling  

---

## Author

Benjamin Wallis  
MSci Physics Project
