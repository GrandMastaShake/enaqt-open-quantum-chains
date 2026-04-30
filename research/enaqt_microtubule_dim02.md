# Dimension 2: QuTiP HEOM Simulation Framework for ENAQT

## Executive Summary

QuTiP 5 (released March 2024) provides the most comprehensive open-source framework for simulating Environment-Assisted Quantum Transport (ENAQT) via its Hierarchical Equations of Motion (HEOM) solver. The HEOM solver supports multiple bath types (bosonic and fermionic), mixed baths, arbitrary spectral densities via fitting, and provides access to auxiliary density operators (ADOs) for extracting environmental information. The QD3SET-1 database -- generated using QuTiP's HEOM implementation -- serves as a key benchmark resource. GPU acceleration via QuTiP-JAX shows crossover advantages for large systems, though HEOM-specific GPU support is still under development. The Julia-based HierarchicalEOM.jl offers significant speedups over QuTiP's HEOM for large systems. PIQS enables scalability to hundreds of qubits for permutationally symmetric systems.

---

## 1. QuTiP 5's HEOM Solver: Documentation and API

### 1.1 Bath Types Supported

Claim: QuTiP 5's HEOMSolver supports multiple built-in bath types including Drude-Lorentz (Matsubara and Padé expansions), UnderDamped Brownian motion, Lorentzian fermionic baths, and generic BosonicBath/FermionicBath from user-defined correlation function coefficients. Bosonic and fermionic baths can be mixed in a single problem.
Source: QuTiP 5.2 Documentation -- HEOM API
URL: https://qutip.readthedocs.io/en/latest/apidoc/heom.html
Date: 2026-01-26 (live documentation)
Excerpt: "Each bath can be specified as either an object of type Bath, BosonicBath, FermionicBath, or their subtypes, or as a tuple (env, Q), where env is an ExponentialBosonicEnvironment or an ExponentialFermionicEnvironment and Q the system coupling operator... Each bath must be either bosonic or fermionic, but bosonic and fermionic baths can be mixed."
Context: Full API documentation listing DrudeLorentzBath, DrudeLorentzPadeBath, UnderDampedBath, LorentzianBath, LorentzianPadeBath, FermionicBath, and BosonicBath classes.
Confidence: high

Claim: The built-in bath classes are mirrored by a newer "environment" API that separates bath properties from coupling operators, allowing reuse of bath definitions across different systems.
Source: QuTiP 5 Documentation -- Bosonic Environments Guide
URL: https://qutip.readthedocs.io/en/latest/guide/heom/bosonic.html
Date: 2026-01-26
Excerpt: "Using the environment API, we first create an abstract DrudeLorentzEnvironment describing the bath, and then use its functions to create exponential expansions such as the Matsubara and Pade ones... Note that the coupling operator Q is not part of the environment objects."
Context: The environment API was introduced in QuTiP 5.1.0 to provide more powerful and flexible bath definitions.
Confidence: high

### 1.2 Key API Features

Claim: The HEOMSolver API includes: dynamics via `run()`, steady-state via `steady_state()`, terminator support for improved accuracy, ADO inspection via `HierarchyADOsState`, and current calculations from first-level ADOs for fermionic baths.
Source: QuTiP 5 Documentation
URL: https://qutip.readthedocs.io/en/latest/apidoc/heom.html
Date: 2026-01-26
Excerpt: "The HEOMSolver class constructs the right-hand-side of Eq.(11), taking as input one or more baths... The HEOMSolver object can then be used to solve the coupled set of ordinary differential equations using standard libraries, or calculate the steady-state."
Context: From the QuTiP-BoFiN paper description of the API.
Confidence: high

### 1.3 Supported Spectral Densities

Claim: QuTiP's HEOM supports Drude-Lorentz, UnderDamped Brownian motion, and arbitrary spectral densities via fitting routines (Matsubara, Padé, or correlation function fitting).
Source: QuTiP 4.7 Documentation -- Introduction to HEOM
URL: https://qutip.org/docs/4.7/guide/heom/intro.html
Date: N/A (documentation)
Excerpt: "In the HEOM, for bosonic baths, one typically chooses a Drude-Lorentz spectral density: J_D = 2*lambda*gamma*omega / (gamma^2 + omega^2), or an under-damped Brownian motion spectral density: J_U = alpha^2 * Gamma * omega / [(omega_c^2 - omega^2)^2 + Gamma^2 * omega^2]."
Context: These are the standard spectral densities used in photosynthetic ENAQT simulations.
Confidence: high

---

## 2. QuTiP Examples: ENAQT and Dephasing-Assisted Transport

### 2.1 BoFiN FMO Example -- Direct ENAQT Simulation

Claim: The BoFiN example notebook 2 provides a complete HEOM simulation of the Fenna-Matthews-Olson (FMO) photosynthetic complex, demonstrating how non-Markovian environments protect against pure dephasing -- the core ENAQT mechanism.
Source: GitHub -- tehruhn/bofin
URL: https://github.com/tehruhn/bofin/blob/main/examples/example-2-FMO-example.ipynb
Date: 2020-10-18
Excerpt: "In this example notebook we outline how to employ the HEOM to solve the FMO photosynthetic complex dynamics. We aim to replicate the results in reference..."
Context: Also referenced at https://github.com/qutip/qutip-notebooks/blob/master/examples/heom/heom-2-fmo-example.ipynb as the official QuTiP FMO tutorial.
Confidence: high

Claim: The QuTiP-BoFiN paper explicitly demonstrates that a suitable non-Markovian environment can protect against pure dephasing in the FMO complex, and that system-environment entanglement oscillations persist on longer timescales than electronic coherence.
Source: Physical Review Research 5, 013181 (2023)
URL: https://link.aps.org/doi/10.1103/PhysRevResearch.5.013181
Date: 2023-03-15
Excerpt: "For the latter, we both clarify how a suitable non-Markovian environment can protect against pure dephasing, and model recent experimental results demonstrating the suppression of electronic coherence. Importantly, we show that by combining the HEOM method with the reaction coordinate method we can observe nontrivial system-environment entanglement oscillations on timescales substantially longer than electronic coherence alone."
Context: This is the definitive QuTiP-BoFiN paper by Lambert et al. (2023).
Confidence: high

### 2.2 Multiple Baths -- Photosynthesis Pattern

Claim: QuTiP's HEOM documentation includes a multiple-baths example that models energy transfer where each basis state interacts with a separate bath, directly applicable to photosynthetic complexes like FMO.
Source: QuTiP 5 Documentation -- Bosonic Environments
URL: https://qutip.readthedocs.io/en/latest/guide/heom/bosonic.html
Date: 2026-01-26
Excerpt: "In the example below we calculate the evolution of a small system where each basis state of the system interacts with a separate bath. Such an arrangement can model, for example, the Fenna-Matthews-Olson (FMO) pigment-protein complex which plays an important role in photosynthesis."
Context: Code example shows how to construct separate DrudeLorentzBath objects for each site and pass a list of baths to HEOMSolver.
Confidence: high

---

## 3. Spin-Boson Models with Pure Dephasing in QuTiP

### 3.1 Built-in Spin-Boson Example

Claim: QuTiP provides a complete spin-boson example using HEOMSolver with underdamped Brownian motion spectral density, comparing HEOM results against mesolve (Lindblad) and brmesolve (Bloch-Redfield), showing substantial differences in the non-Markovian regime.
Source: QuTiP 5 Paper (arXiv:2412.04705)
URL: https://arxiv.org/html/2412.04705v1
Date: 2024-12-06
Excerpt: "For the example of a standard spin-boson problem, we compare the output of the different master equation solvers available in QuTiP, namely, the Lindblad master equation solver (mesolve()), the Bloch-Redfield equation solver (brmesolve()) and the hierarchical equations of motion solver (HEOMSolver) for a spin coupled to a bosonic bath described by an underdamped Brownian motion spectral density... Notice that both brmesolve and mesolve coincide, but they differ substantially from the HEOM result in this deeply non-Markovian regime."
Context: Parameters used: lambda=0.5*Delta, Gamma=0.1*Delta, T=0.5*Delta, Nk=5, Nc=6, omega0=3*Delta/2. This is the canonical pure dephasing regime relevant to ENAQT.
Confidence: high

### 3.2 Spin-Boson Localization-Delocalization Transition

Claim: QuTiP 5 includes a complete example of the spin-boson localization-delocalization phase transition at T=0, demonstrating how to use correlation function fitting with arbitrary spectral densities.
Source: QuTiP 5 Paper Examples
URL: https://github.com/qutip/qutip-paper-v5-examples
Date: 2024-12-01
Excerpt: "heom_transition.py: A more involved example of how to use the fitting features to capture the spin-boson phase transition with the HEOM solver. Used to generate Fig. 12"
Context: The example uses an Ohmic spectral density with polynomial cutoff and demonstrates fitting the correlation function at T=0 where standard spectral density approaches fail.
Confidence: high

### 3.3 Pure Dephasing Model Validation

Claim: The pure dephasing model (Q=sigma_z) is the standard test case for ENAQT simulations, and the cumulant equation matches the exact HEOM solution for this case, as demonstrated using QuTiP-BoFiN.
Source: arXiv:2403.04488v2 -- Dynamics of the Non-equilibrium Spin-Boson Model
URL: https://arxiv.org/html/2403.04488v2
Date: 2024-09-14
Excerpt: "For the pure dephasing model we show that the cumulant equation matches the exact solution, as does time convolutionless second order master equation but not Bloch-Redfield (BR). Code to reproduce to simulate the cumulant equation in general has been made available by the authors."
Context: Uses QuTiP BoFiN HEOM implementation as the exact reference. Parameters: f1=f2=0, f3=1 (pure dephasing), omega0/T=4, lambda/gamma=1/4, gamma=5*omega0.
Confidence: high

---

## 4. QuTiP-JAX for GPU Acceleration

### 4.1 GPU Speedup for Quantum Dynamics

Claim: QuTiP-JAX provides GPU acceleration via the JAX data layer and Diffrax integrator. For an Ising spin chain benchmark on an NVIDIA A100 GPU, a crossover point exists where GPU outperforms CPU by up to two orders of magnitude. However, mesolve (open systems) reaches memory limits at ~11 spins, while sesolve (closed systems) scales to ~22 spins.
Source: QuTiP 5 Paper (arXiv:2412.04705)
URL: https://arxiv.org/html/2412.04705v1
Date: 2024-12-06
Excerpt: "By comparing the performance of simulating the dynamics of the Ising spin chain with the standard QuTiP CPU-bound method (adams) and the equivalent diffrax method on a GPU, we observe a threshold in system size where the GPU outperforms the CPU calculation by up to two orders of magnitude... However, even with the jaxdia format and using a state-of-the-art graphics card with 80 gigabytes of RAM, the memory limit is reached already at 11 spins [for mesolve]."
Context: The GPU memory limit is the primary bottleneck. PIQS and HEOMSolver do not yet support the JAX data layer.
Confidence: high

### 4.2 HEOM GPU Support Status

Claim: GPU implementation of HEOM is listed as a future development project for QuTiP, requiring CUDA/OpenCL knowledge. No native GPU HEOM support exists yet.
Source: QuTiP Development Ideas -- HEOM GPU
URL: https://qutip.readthedocs.io/en/v5.2.3/development/ideas/heom-gpu.html
Date: 2026-01-26
Excerpt: "The Hierarchical Equations of Motion (HEOM) method is a non-perturbative approach to simulate the evolution of the density matrix of dissipative quantum systems. The underlying equations are a system of coupled ODEs which can be run on a GPU. This will allow the study of larger systems... Expected outcomes: A version of HEOM which runs on a GPU; Performance comparison with the CPU version; Implement dynamic scaling."
Context: Listed as a "Hard" difficulty project. Mentors: Neill Lambert, Alex Pitchford, Shahnawaz Ahmed, Simon Cross.
Confidence: high

### 4.3 Dynamiqs -- Alternative GPU-Accelerated Solver

Claim: Dynamiqs is a JAX-based quantum simulation library that outperforms QuTiP on GPU by up to 60x for specific benchmarks. For the dissipative cat CNOT benchmark, Dynamiqs v0.3.0 achieves 60x speedup over QuTiP 5.0.4 via GPU acceleration.
Source: Alice & Bob Blog -- Meet Dynamiqs
URL: https://alice-bob.com/blog/dynamiqs-gpu-opensource-quantum-simulation-library/
Date: 2024-11-20
Excerpt: "For all systems that we consider here, the latest release of Dynamiqs (v0.2.2) is at least on par with QuTiP and features up to a 30x speed increase on the dissipative cat CNOT benchmark by leveraging GPU acceleration. On Dynamiqs v0.3.0 which is under development, the addition of the sparse representation of matrices provides a significant additional speedup, up to 60x on the same benchmark."
Context: Dynamiqs uses JAX for GPU acceleration and is designed for batched simulations. However, it does not currently support HEOM.
Confidence: medium

---

## 5. Transport Efficiency Calculations Using QuTiP

### 5.1 Lindblad Master Equation for Transport

Claim: QuTiP's mesolve can simulate quantum transport in open chains using Lindblad dissipators. The QuTiP 5 paper includes a complete example of transport through a quantum system with left/right leads modeled as Lindblad jump operators, solved on GPU using JAX.
Source: QuTiP 5 Paper (arXiv:2412.04705)
URL: https://arxiv.org/html/2412.04705v1
Date: 2024-12-06
Excerpt: "Below is the implementation in Python using QuTiP and JAX: [code showing] L0 = (qt.liouvillian(H) + qt.lindblad_dissipator(c_op_L) - 0.5*qt.spre(c_op_R.dag()*c_op_R) - 0.5*qt.spost(c_op_R.dag()*c_op_R))"
Context: This example models transport through a two-level system with left and right leads using Lindblad operators, solved with automatic differentiation enabled via JAX.
Confidence: high

### 5.2 HEOM for Non-Markovian Transport

Claim: QuTiP's HEOM can calculate non-equilibrium currents between system and fermionic baths from first-level ADOs, enabling exact transport efficiency calculations in the non-Markovian regime.
Source: QuTiP 4.7 Documentation -- Fermionic Environments
URL: https://qutip.org/docs/4.7/guide/heom/fermionic.html
Date: N/A
Excerpt: "The currents between the system and a fermionic bath may be calculated from the first level auxiliary density operators (ADOs) associated with the exponents of that bath. The contribution to the current into a given bath from each exponent in that bath is: Contribution from Exponent = +/- i Tr(Q^+/- . A)"
Context: Code examples show how to calculate left and right currents using `heom_current()` functions operating on ADO states.
Confidence: high

### 5.3 ENAQT Transport Efficiency -- Experimental Connection

Claim: ENAQT has been experimentally demonstrated in perovskite nanocrystal superlattices, where maximum steady-state diffusion occurs at intermediate temperatures where static disorder and thermal fluctuations balance. The Haken-Strobl-Reineker (HSR) model was used for simulation.
Source: Nature Communications (2025)
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025-02-02
Excerpt: "Most remarkably, at intermediate temperatures, the interplay between coherence and dephasing collaboratively enhances transport. We observed a maximum for the steady-state diffusion constant at temperatures where static disorder and thermal fluctuations balance each other, providing experimental evidence of ENAQT."
Context: Simulations used the Anderson Hamiltonian and Haken-Strobl-Reineker model. QuTiP could replicate these simulations using HEOM or Lindblad approaches.
Confidence: high

---

## 6. Benchmarks: QuTiP HEOM vs. Other Packages

### 6.1 HierarchicalEOM.jl (Julia) -- Significant Speedup

Claim: HierarchicalEOM.jl, a Julia framework built on QuantumToolbox.jl, achieves "a significant speedup with respect to the corresponding method in the Quantum Toolbox in Python (QuTiP)." It features memory optimization with lazy operators for large-scale systems.
Source: arXiv:2306.07522 / Communications Physics 6, 313 (2023)
URL: https://arxiv.org/abs/2306.07522
Date: 2023-06-13
Excerpt: "HierarchicalEOM.jl achieves a significant speedup with respect to the corresponding method in the Quantum Toolbox in Python (QuTiP), upon which this package is founded... The lazy operator feature enables memory-efficient representation of HEOMLS matrices, with memory savings proportional to the number of ADOs, making it possible to simulate larger systems with higher hierarchy tiers."
Context: QuantumToolbox.jl is the Julia equivalent of QuTiP. HierarchicalEOM.jl is maintained by the QuTiP organization (https://github.com/qutip/HierarchicalEOM.jl).
Confidence: high

### 6.2 DM-HEOM -- Distributed Memory Framework

Claim: DM-HEOM is a portable, scalable solver-framework for HEOM using distributed memory (MPI), developed at Zuse Institute Berlin for photoactive system simulations on HPC clusters.
Source: GitHub -- noma/dm-heom
URL: https://github.com/noma/dm-heom
Date: N/A
Excerpt: "DM-HEOM comprises a software suite to compute Open Quantum System Dynamics with the Hierarchical Equations of Motion approach... DM-HEOM: A Portable and Scalable Solver-Framework for the Hierarchical Equations of Motion, 2018 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW)"
Context: DM-HEOM is designed for large-scale distributed memory systems and has been used for FMO complex simulations. Reference: M. Noack et al., J. Comput. Chem. 39, 1779 (2018).
Confidence: high

### 6.3 GPU-HEOM (Tanimura Group) -- 400x+ Speedups

Claim: CUDA-based GPU implementations of HEOM achieve speedups of 400-458x over CPU implementations for FMO-like systems, enabling calculation of 7-site FMO models at high hierarchy depths.
Source: Tsuchimoto & Tanimura, J. Chem. Theory Comput. 11, 3859 (2015)
URL: http://theochem.kuchem.kyoto-u.ac.jp/public/TT15JCTC.pdf
Date: 2015
Excerpt: "We developed a compact code for the reduced hierarchy equations of motion (HEOM) for a graphics processor unit (GPU) that can treat the system as large as 4096 energy states."
Context: The code uses Padé spectrum decomposition (PSD) and exponential integrators. The CUDA source code is provided as supporting information.
Confidence: high

### 6.4 BoFiN-fast -- C++ Acceleration

Claim: BoFiN-fast is a hybrid C++/Python implementation of the BoFiN HEOM solver that performs RHS construction in C++, achieving similar speed to the current QuTiP 5 pure-Python implementation.
Source: QuTiP Documentation -- Previous Implementations
URL: https://qutip.readthedocs.io/en/latest/guide/heom/history.html
Date: 2026-01-26
Excerpt: "The construction of the right-hand side matrix for BoFiN was slow, so BoFiN-fast, a hybrid C++ and Python implementation, was written that performed the right-hand side construction in C++. It was otherwise identical to the pure Python version... The current implementation is a rewrite of BoFiN in pure Python. Its right-hand side construction has similar speed to BoFiN-fast, but is written in pure Python."
Context: QuTiP 5's HEOM solver matches BoFiN-fast performance while remaining in pure Python.
Confidence: high

---

## 7. QuTiP Extensions for Photosynthetic Complexes

### 7.1 QuTiP-BoFiN -- The Definitive HEOM Library

Claim: QuTiP-BoFiN (Phys. Rev. Research 5, 013181, 2023) is the definitive numerical HEOM library integrated with QuTiP, with applications in light-harvesting, quantum control, and single-molecule electronics.
Source: Physical Review Research
URL: https://link.aps.org/doi/10.1103/PhysRevResearch.5.013181
Date: 2023-03-15
Excerpt: "Here we present a numerical library in Python, integrated with the powerful QuTiP platform, which implements the HEOM for both bosonic and fermionic environments. We demonstrate its utility with a series of examples... a study of the dynamics of energy transfer in the Fenna-Matthews-Olson photosynthetic complex."
Context: Authors: Neill Lambert, Tarun Raheja, Simon Cross, Paul Menczel, Shahnawaz Ahmed, Alexander Pitchford, Daniel Burgarth, Franco Nori.
Confidence: high

### 7.2 QD3SET-1 Database -- QuTiP-Generated Benchmark Data

Claim: The QD3SET-1 database contains quantum dissipative dynamics datasets generated using QuTiP's HEOM implementation, including 1,000 spin-boson trajectories and 879 FMO complex trajectories, with individual datasets requiring up to 1,000 CPU hours.
Source: Frontiers in Physics 11, 1223973 (2023)
URL: https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2023.1223973/full
Date: 2023-07-03
Excerpt: "For all combinations of these parameters, the system RDM was propagated using the HEOM approach implemented in the QuTiP software package... In total, 1,000 HEOM calculations, 500 for symmetric and 500 for asymmetric spin-boson Hamiltonians, were performed... Computations for the database were executed on high-performance computing (HPC) clusters, with individual data sets demanding up to 1000 CPU hours for propagation."
Context: The SB dataset used QuTiP HEOM. The FMO HEOM dataset used the PHI code. Other FMO datasets used LTLME (modified quantum_HEOM package).
Confidence: high

### 7.3 QuTiP Tutorials for Photosynthetic Systems

Claim: The QuTiP tutorials page includes dedicated HEOM tutorial notebooks, with specific examples for bosonic baths, fermionic baths, and the FMO complex.
Source: QuTiP Tutorials
URL: https://qutip.org/tutorials
Date: N/A
Excerpt: "Hierarchical Equations of Motion -- Example notebooks demonstrating how to use QuTiP's HEOM solver for both fermionic and bosonic baths. These examples have been explained in detail in our recent publication [QuTiP-BoFiN paper]."
Context: Tutorial categories include: HEOM 1a (Drude-Lorentz bosonic bath), HEOM 2 (FMO example), HEOM 3a-d (discrete boson and fermionic leads), HEOM 4 (fitting spectral densities).
Confidence: high

---

## 8. Published Papers Using QuTiP for ENAQT-Related Simulations

### 8.1 QuTiP-BoFiN Paper (2023)

Claim: The QuTiP-BoFiN paper (Phys. Rev. Research 5, 013181) is the primary reference for QuTiP HEOM applications to ENAQT, demonstrating FMO energy transfer with non-Markovian environment protection against pure dephasing.
Source: Physical Review Research
URL: https://arxiv.org/abs/2010.10806 (preprint)
Date: 2020-10-21 (v1), 2023-03-15 (published)
Excerpt: "For the bosonic case, our results include demonstrations of fitting arbitrary spectral densities with different approaches, and a study of the dynamics of energy transfer in the Fenna-Matthews-Olson photosynthetic complex... we both clarify how a suitable non-Markovian environment can protect against pure dephasing, and model recent experimental results demonstrating the suppression of electronic coherence."
Context: Section II.D specifically addresses FMO energy transfer and dephasing protection.
Confidence: high

### 8.2 QD3SET-1 Paper (2023)

Claim: The QD3SET-1 database paper (arXiv:2301.12096) used QuTiP HEOM to generate 1,000 spin-boson trajectories as benchmark data for machine learning quantum dynamics.
Source: arXiv:2301.12096
URL: https://arxiv.org/abs/2301.12096
Date: 2023-01-28
Excerpt: "The system RDM was propagated using HEOM approach implemented in QuTiP software package... The total propagation time was t_max*Delta = 20 and the HEOM integration time-step was set to dt*Delta = 0.05. In total, 1,000 of HEOM calculations... were performed."
Context: Database now hosted at Zenodo (DOI: 10.5281/zenodo.7557558) with Python extraction tools at https://github.com/Arif-PhyChem/QD3SET.
Confidence: high

### 8.3 Robustness of ENAQT (Shabani et al., 2014)

Claim: Shabani et al. (2014) provided numerical evidence for the robustness of ENAQT in the FMO complex, simulating excitonic energy transfer with respect to variations in environmental interactions, Hamiltonian parameters, disorders, and dipole orientations.
Source: Physical Review E 89, 042706 (2014)
URL: https://dspace.mit.edu/handle/1721.1/89154
Date: 2014-09-03
Excerpt: "We simulated excitonic energy transfer in Fenna-Matthews-Olson photosynthetic complex. We found that ENAQT is robust with respect to many relevant parameters of environmental interactions and Frenkel-exciton Hamiltonians, including reorganization energy, bath-frequency cutoff, temperature, initial excitations, dissipation rate, trapping rate, disorders, and dipole moments orientations."
Context: This was one of the first comprehensive robustness studies of ENAQT in the FMO complex.
Confidence: high

### 8.4 Non-equilibrium Spin-Boson Model (2024)

Claim: A 2024 paper on the non-equilibrium spin-boson model used QuTiP-BoFiN HEOM as the exact reference to validate cumulant equations, showing that steady-state coherences persist and have consequences for quantum machine efficiencies.
Source: arXiv:2403.04488v2
URL: https://arxiv.org/html/2403.04488v2
Date: 2024-09-14
Excerpt: "We do this by means of the hierarchical equations of motion (HEOM), using the recent Qutip BoFiN implementation... This translates into the fact that the steady state is not diagonal in the basis of the system Hamiltonian, a fact that has consequences for the efficiencies of quantum machines."
Context: The paper provides open-source code for cumulant equation simulations validated against QuTiP-BoFiN HEOM.
Confidence: high

### 8.5 ENAQT in Perovskite Nanocrystal Superlattices (2025)

Claim: A 2025 Nature Communications paper provided the first experimental visualization of ENAQT in perovskite nanocrystal superlattices, using Anderson Hamiltonian and Haken-Strobl-Reineker model simulations.
Source: Nature Communications (2025)
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025-02-02
Excerpt: "At intermediate temperatures, the interplay between coherence and dephasing collaboratively enhances transport. We observed a maximum for the steady-state diffusion constant at temperatures where static disorder and thermal fluctuations balance each other, providing experimental evidence of ENAQT."
Context: These simulations could be replicated in QuTiP using mesolve with dephasing operators or HEOM for non-Markovian environments.
Confidence: high

---

## 9. QuTiP Gallery and Examples Repository

### 9.1 Official Tutorials Repository

Claim: The QuTiP tutorials repository (qutip/qutip-tutorials) contains Jupyter notebooks for all QuTiP features, including dedicated HEOM tutorials for bosonic and fermionic baths.
Source: QuTiP Tutorials
URL: https://qutip.org/tutorials
Date: N/A
Excerpt: "These notebooks demonstrate and introduce specific functionality in QuTiP... Hierarchical Equations of Motion -- HEOM tutorials demonstrate how to use QuTiP's hierarchical equations of motion solver for both fermionic and bosonic baths."
Context: The tutorials are organized into: Python Introduction, Visualization, Quantum circuits, Time evolution, Optimal control, Tomography, Permutational invariant Lindblad dynamics, HEOM, and Miscellaneous.
Confidence: high

### 9.2 QuTiP Paper v5 Examples

Claim: The qutip-paper-v5-examples repository contains all code used to generate figures in the QuTiP 5 paper, including HEOM examples for spin-boson dynamics and the localization-delocalization transition.
Source: GitHub -- qutip/qutip-paper-v5-examples
URL: https://github.com/qutip/qutip-paper-v5-examples
Date: 2024-12-01
Excerpt: "heom_example.py: A basic example of how to use the new environment interface with the HEOM solver. Used to generate Fig. 10 and 11... heom_transition.py: A more involved example of how to use the fitting features to capture the spin-boson phase transition with the HEOM solver."
Context: Also includes JAX examples, Monte Carlo, Floquet, and QIP demonstrations.
Confidence: high

### 9.3 BoFiN Examples Repository

Claim: The BoFiN repository (tehruhn/bofin) contains an extensive set of example notebooks including the FMO complex, dynamical decoupling, non-equilibrium heat flow, arbitrary spectral density fitting, and single-molecule electronics.
Source: GitHub -- tehruhn/bofin
URL: https://github.com/tehruhn/bofin/tree/main/examples
Date: 2020
Excerpt: "BoFiN also came with an extensive set of example notebooks that are available at https://github.com/tehruhn/bofin/tree/main/examples."
Context: Notebooks include: example-1a (spin-boson Drude-Lorentz), example-1b (spin-boson UnderDamped), example-2 (FMO complex), example-3 (dynamical decoupling), example-4a/b (fermionic impurity and quantum dot), example-5 (fitting arbitrary spectral densities).
Confidence: high

### 9.4 Legacy qutip-notebooks Repository

Claim: The older qutip-notebooks repository contains many additional HEOM examples that have not yet been ported to v5 but remain a valuable resource.
Source: GitHub -- qutip/qutip-notebooks
URL: https://github.com/qutip/qutip-notebooks
Date: N/A
Excerpt: "This older repository contains many more examples and tutorials that have not yet been ported to v5. However, they remain an important resource, and overtime we hope to have all relevant notebooks updated and included in the official qutip-tutorials repository."
Context: The HEOM folder contains: heom-index.ipynb, heom-1a-spin-boson.ipynb, heom-2-fmo-example.ipynb, heom-3a-3d-fermionic-examples.ipynb.
Confidence: high

---

## 10. PIQS (Permutationally Invariant Quantum Solver)

### 10.1 Scalability to Hundreds of Qubits

Claim: PIQS enables exact Lindbladian dynamics of open quantum systems consisting of identical qubits by exploiting permutational symmetry, reducing the Hilbert space from exponential (2^N) to polynomial size, enabling simulation of hundreds of qubits.
Source: PIQS Documentation
URL: https://piqs.readthedocs.io/en/latest/intro.html
Date: 2020-07-17
Excerpt: "In the case where local processes are included in the model of a system's dynamics, numerical simulation requires dealing with density matrices of size 2^N. This becomes infeasible for a large number of qubits. We can simplify the calculations by exploiting the permutational invariance of indistinguishable quantum particles which allows the user to study hundreds of qubits."
Context: PIQS uses Cython for performance and sparse matrices. It is integrated in QuTiP as `qutip.piqs` since version 4.3.1 (2018).
Confidence: high

### 10.2 Limitations for ENAQT

Claim: PIQS is limited to Lindblad master equations and identical qubits with permutationally symmetric initial states. It does NOT support HEOM, non-Markovian baths, or site-specific Hamiltonians required for realistic photosynthetic ENAQT simulations.
Source: QuTiP 5 Paper (arXiv:2412.04705)
URL: https://arxiv.org/html/2412.04705v1
Date: 2024-12-06
Excerpt: "Note however that this currently does not override the behavior of every function which can create a Qobj in QuTiP (e.g., PIQs and HEOMSolver do not support this option [JAX data layer])."
Context: PIQS is applicable to collective spin dynamics (Dicke states, GHZ states, superradiance) but not to spatially extended systems with local disorder and site-specific couplings like FMO.
Confidence: high

### 10.3 Applications

Claim: PIQS applications include: time evolution for permutationally symmetric states, phase transitions of driven-dissipative systems, correlation functions, resonance fluorescence, steady-state superradiance, spin squeezing, and boundary time crystals.
Source: PIQS GitHub Repository
URL: https://github.com/nathanshammah/piqs
Date: 2017-09-22
Excerpt: "The time evolution of the total density matrix of quantum optics and cavity QED systems for permutationally symmetric initial states... Phase transitions of driven-dissipative out-of-equilibrium quantum systems... Spin squeezing for quantum metrology, long-range interaction in noisy spin models, decoherence in quantum information processing."
Context: While not directly applicable to FMO-type ENAQT, PIQS could model simplified collective dephasing models relevant to understanding the interplay between coherence and dissipation in symmetric systems.
Confidence: high

---

## 11. Practical Implementation Guide for ENAQT in QuTiP

### 11.1 Recommended Workflow for ENAQT Simulations

Based on the research findings, here is a practical workflow for simulating ENAQT using QuTiP:

**For Small Systems (2-8 sites) with Non-Markovian Baths:**
```python
from qutip.solver.heom import HEOMSolver, DrudeLorentzBath
# 1. Define system Hamiltonian
# 2. Create DrudeLorentzBath for each site
# 3. Include terminators for accuracy
# 4. Run HEOMSolver with max_depth=4-8
# 5. Extract populations and coherences from result
```

**For Larger Systems with Markovian Approximation:**
```python
from qutip import mesolve, lindblad_dissipator
# 1. Define H_Sys with site energies and couplings
# 2. Add local dephasing operators: sqrt(gamma_deph) * |i><i|
# 3. Add trapping/receiving Lindblad operators
# 4. Use mesolve or QuTiP-JAX for GPU acceleration
```

**For Permutationally Symmetric Systems:**
```python
from qutip.piqs import Dicke, dicke_basis, liouvillian
# 1. Use Dicke algebra for collective systems
# 2. Build Liouvillian with collective and local rates
# 3. Scale to hundreds of qubits
```

### 11.2 Key Parameters for ENAQT

| Parameter | Symbol | Typical Range | Description |
|-----------|--------|---------------|-------------|
| Reorganization energy | lambda | 10-100 cm^-1 | System-bath coupling strength |
| Bath cutoff frequency | gamma | 50-500 cm^-1 | Memory timescale of bath |
| Temperature | T | 77-300 K | Thermal energy |
| Hierarchy depth | N_c | 4-8 | Truncation level |
| Matsubara terms | N_k | 2-8 | Correlation expansion terms |
| Dephasing rate | gamma_phi | 0.1-10 ps^-1 | Pure dephasing (Lindblad) |

### 11.3 Performance Considerations

- **HEOM scaling**: Number of ADOs grows combinatorially: N_tot = (N + N_max)! / [N! * N_max!]
- **7-site FMO at N_max=16**: ~240,000 coupled density matrices
- **QuTiP HEOM**: Pure Python, similar speed to BoFiN-fast C++
- **Julia HierarchicalEOM.jl**: Significant speedup over QuTiP HEOM
- **GPU HEOM**: Not yet available in QuTiP; listed as future project
- **QD3SET generation**: Up to 1,000 CPU hours per dataset

---

## 12. Summary Table of Key Resources

| Resource | URL | Purpose |
|----------|-----|---------|
| QuTiP HEOM API | https://qutip.readthedocs.io/en/latest/apidoc/heom.html | Complete API reference |
| QuTiP HEOM Tutorials | https://qutip.org/tutorials | Jupyter notebook examples |
| BoFiN FMO Example | https://github.com/tehruhn/bofin/blob/main/examples/example-2-FMO-example.ipynb | FMO complex simulation |
| QuTiP 5 Paper Code | https://github.com/qutip/qutip-paper-v5-examples | All figures from paper |
| QD3SET Database | https://github.com/Arif-PhyChem/QD3SET | Benchmark datasets |
| QuTiP-BoFiN Paper | https://arxiv.org/abs/2010.10806 | HEOM library paper |
| QuTiP 5 Paper | https://arxiv.org/abs/2412.04705 | Framework overview |
| HierarchicalEOM.jl | https://github.com/qutip/HierarchicalEOM.jl | Faster Julia HEOM |
| DM-HEOM | https://github.com/noma/dm-heom | Distributed memory HEOM |
| PIQS | https://piqs.readthedocs.io | Permutational symmetry solver |
| QuTiP-JAX | https://github.com/qutip/qutip-jax | GPU acceleration |

---

*Research compiled from 15+ independent searches across documentation, academic papers, code repositories, and benchmarks.*
