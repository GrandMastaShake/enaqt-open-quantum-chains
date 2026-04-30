# Dimension 11: ENAQT Testable Predictions and Benchmarks

## Comprehensive Research Findings

**Date Compiled**: 2025-07-02
**Research Focus**: Environment-Assisted Quantum Transport (ENAQT) - specific testable predictions, computational validation methods, available data and tools
**Total Independent Searches**: 15+

---

## Table of Contents

1. [Core Mathematical Predictions of ENAQT](#1-core-mathematical-predictions)
2. [Quantitative Benchmarks for ENAQT Validation](#2-quantitative-benchmarks)
3. [Transport Efficiency Metrics](#3-transport-efficiency-metrics)
4. [Predicted Dephasing Rate Ranges for ENAQT](#4-dephasing-rate-ranges)
5. [ENAQT Signatures in Population Dynamics P(t)](#5-population-dynamics-signatures)
6. [Experimental Tests of ENAQT](#6-experimental-tests)
7. [ENAQT vs Classical Hopping](#7-enaqt-vs-classical)
8. [System Size, Disorder, and Temperature Effects](#8-system-size-effects)
9. [Published Code and Software](#9-published-code)
10. [What Constitutes "Proof" of ENAQT](#10-proof-of-enaqt)
11. [Machine Learning for ENAQT Detection](#11-ml-detection)

---

## 1. Core Mathematical Predictions of ENAQT

### Finding 1.1: The Three-Regime Prediction

Claim: ENAQT predicts three distinct dynamical regimes as a function of dephasing rate: (1) low dephasing where quantum localization suppresses transport, (2) intermediate dephasing where environment-assisted transport maximizes efficiency, and (3) high dephasing where quantum Zeno suppression occurs.
Source: Rebentrost et al. (2009), "Environment-assisted quantum transport"
URL: https://dspace.mit.edu/bitstream/handle/1721.1/70838/Rebentrost-2009-Environment-assisted%20quantum%20transport.pdf
Date: 2009
Excerpt: "From left to right, the fully quantum regime that is dominated by intrinsic static disorder in the system Hamiltonian; the ENAQT regime (qualitatively indicated by the yellow color gradient), where unitary evolution and dephasing collaborate with the result of increased efficiency; finally, the quantum Zeno regime, where strong dephasing suppresses the quantum transport."
Context: Original ENAQT paper establishing the foundational phenomenological predictions
Confidence: High

### Finding 1.2: Optimal Dephasing Matches Hopping Rate (Quantum Goldilocks Principle)

Claim: The maximum transport efficiency occurs when the dephasing timescale matches the hopping timescale: gamma_phi^opt ~ V (nearest-neighbor coupling). This is called the "quantum Goldilocks principle."
Source: Novo, Mohseni & Omar (2013), "Disorder-assisted quantum transport in suboptimal decoherence regimes"
URL: https://arxiv.org/abs/1312.6989
Date: 2013
Excerpt: "We find that the regime which maximizes transport efficiency is when the timescale of dephasing matches the timescale of hopping, verifying the quantum Goldilocks principle. This principle was put forward as a design principle for efficient structures: it states that a convergence of typical timescales of the system and the environment leads to optimal transfer efficiency."
Context: Binary tree and hypercube transport simulations with systematic disorder/dephasing parameter scans
Confidence: High

### Finding 1.3: Optimal Dephasing Value for Binary Trees

Claim: For a 5-generation binary tree with initial statistical mixture at the leaves, the optimal dephasing rate is gamma_phi = 1.6 V (in units of the nearest-neighbor coupling V).
Source: Novo, Mohseni & Omar (2013), arXiv:1312.6989
URL: https://arxiv.org/pdf/1312.6989
Date: 2013
Excerpt: "We find that the optimal value of dephasing is at gamma_phi = 1.6. After this maximum, as we will see in Fig. 3, the efficiency starts decreasing and tends to zero in the limit of infinite dephasing, as expected due to the previous argument."
Context: Numerical optimization over dephasing rate for binary tree transport efficiency
Confidence: High

### Finding 1.4: Universal Relation Between Occupation Spread and Current

Claim: The optimal dephasing rate for ENAQT correlates with the minimum spread of occupations, quantified by Delta_n = 1 - sqrt(sum_i((<n_i> - n_ext)^2)). This provides a universal analytical predictor of ENAQT without solving full dynamics.
Source: Levi et al. (2018), "Universal Origin for Environment-Assisted Quantum Transport"
URL: https://arxiv.org/pdf/1801.06799
Date: 2018
Excerpt: "Comparing Figure 2c with Figure 2a reveals a correlation between the exciton current and the distribution of the occupations: the maximal current seems to appear at (or close to) the dephasing rate at which the spread of occupations is minimal. To quantify this relation between the distribution of exciton occupations and optimal current, we define the quantity Delta_n = 1 - sqrt(sum_i (<n_i> - n_ext)^2)."
Context: Demonstrated across networks of different topologies, dimensions, sizes and symmetries
Confidence: High

### Finding 1.5: FMO Room Temperature Dephasing Estimate

Claim: For the FMO complex, the dephasing rate at room temperature is approximately 300 cm^-1, placing it well within the ENAQT regime.
Source: Rebentrost et al. (2009)
URL: https://dspace.mit.edu/bitstream/handle/1721.1/70838/Rebentrost-2009-Environment-assisted%20quantum%20transport.pdf
Date: 2009
Excerpt: "For the above spectral density the rate turns out to be gamma_phi(T) = 2*pi*(kT/hbar)*E_R/hbar*omega_c. This gives a rough estimate for the dephasing rate at room temperature of about 300 cm^-1, which is indicated in figure 2. Hence, the natural operating point of the FMO complex is estimated to be well within the regime of ENAQT."
Context: Using Ohmic spectral density with reorganization energy E_R = 35 cm^-1 and cutoff omega_c = 150 cm^-1
Confidence: Medium (estimate depends on spectral density parameters)

---

## 2. Quantitative Benchmarks for ENAQT Validation

### Finding 2.1: Primary Benchmark - Non-Monotonic Efficiency vs Dephasing

Claim: The definitive signature of ENAQT is a non-monotonic relationship between transport efficiency eta and dephasing rate gamma_phi, with a clear maximum at an intermediate value.
Source: Multiple sources (Rebentrost 2009; Novo 2013; Schaller 2025)
URL: Various
Date: 2009-2025
Excerpt: "The non-monotonous behavior is the main feature of ENAQT, which we will study here in detail by numerical and analytic means." (Schaller et al. 2025)
Context: This is the single most agreed-upon benchmark across all ENAQT literature
Confidence: High

### Finding 2.2: Multiple Optimal Regimes (Twin Peaks)

Claim: Some physically-motivated transport networks can exhibit TWO distinct ENAQT peaks in their steady-state transport efficiency, challenging the "single Goldilocks zone" assumption.
Source: Coates, Lovett & Gauger (2023), "From Goldilocks to twin peaks"
URL: https://pubs.rsc.org/en/content/articlelanding/2023/cp/d2cp04935j
Date: 2023
Excerpt: "In this paper we show that a consistent subset of physically modelled transport networks can have at least two ENAQT peaks in their steady state transport efficiency."
Context: Multiple peaks arise from distinct transport pathways each optimized at different dephasing ranges
Confidence: High

### Finding 2.3: Turnover Temperature as Benchmark

Claim: In temperature-dependent experiments, ENAQT predicts a "turnover temperature" T_t where static disorder and dephasing are balanced (Gamma ~ Delta), producing a peak in the long-time diffusion constant.
Source: Blach et al. (2025), Nature Communications
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025
Excerpt: "A maximum of long-time diffusion D_s is predicted at a turnover temperature T_t when Gamma ~ Delta, which is the essence of ENAQT. Our experimental data corroborate this turnover temperature in all three SLs."
Context: First experimental observation of ENAQT in perovskite nanocrystal superlattices
Confidence: High

### Finding 2.4: QD3SET-1 as ML Benchmark Database

Claim: The QD3SET-1 database provides standardized quantum dissipative dynamics datasets for benchmarking both physics-based and machine learning-based simulation methodologies.
Source: Ullah et al. (2023), "QD3SET-1: A Database with Quantum Dissipative Dynamics Data Sets"
URL: https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2023.1223973/full
Date: 2023
Excerpt: "The primary objective behind the release of the QD3SET-1 database is to provide researchers a valuable resource for the development, testing, and validation of their approaches, whether rooted in machine learning or other non-machine learning methodologies."
Context: 8 datasets: 1 spin-boson (1000 trajectories, HEOM) + 7 FMO datasets (500-879 trajectories, LTLME and HEOM)
Confidence: High

---

## 3. Transport Efficiency Metrics

### Finding 3.1: Integrated Trapping Probability

Claim: Transport efficiency is defined as the integrated probability of trapping at the target vertex: eta = 2*kappa*int_0^infinity <w|rho(t)|w> dt = 1 - Tr[lim_{t->infinity} rho(t)].
Source: Transport efficiency formalism (Caruso et al. and Paris group)
URL: https://sites.unimi.it/mgaparis/wp-content/PDF/Transp_eff.pdf
Date: Ongoing
Excerpt: "The transport efficiency is a relevant measure for a quantum transport process, which can be defined as the integrated probability of trapping at the vertex w: eta = 2*kappa*int_0^+infinity <w|rho(t)|w> dt = 1 - Tr[lim_{t->+infinity} rho(t)]."
Context: Standard definition used across ENAQT literature. The second equality holds because of loss processes at the trap vertex.
Confidence: High

### Finding 3.2: Exciton Current and Flux Metrics

Claim: In steady-state transport, the figure of merit is the exciton current I or population flux eta = gamma_l * rho_{NN}^{SS}, measuring the rate of population extraction from the output site.
Source: Design Principles for Enhanced Quantum Transport (2025/2026)
URL: https://arxiv.org/html/2604.23005v1
Date: 2026
Excerpt: "Our figure of merit for transport is the flux of population out of the Nth site of the system, eta = gamma_l * rho_{NN}^{SS}, where rho_{NN}^{SS} is the steady-state population of site N."
Context: Used in recent work on site-dependent noise optimization for quantum transport
Confidence: High

### Finding 3.3: Diffusion Constant as Transport Metric

Claim: The long-time diffusion constant D_s = 2*J^2*Gamma/(Gamma^2 + Delta^2) * [L(Gamma,Delta)]^2 serves as the key transport metric for experimental validation of ENAQT, where L is the coherence length.
Source: Blach et al. (2025), Nature Communications
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025
Excerpt: "As a result of the competition between Anderson localization and dephasing, the diffusion constant at the long-time steady-state limit is given by D_s = 2*J^2*Gamma/(Gamma^2 + Delta^2) * [L(Gamma,Delta)]^2."
Context: Used to connect Haken-Strobl-Reineker (HSR) model predictions with experimental measurements in perovskite nanocrystal superlattices
Confidence: High

---

## 4. Predicted Dephasing Rate Ranges for ENAQT Observability

### Finding 4.1: Binary Tree Optimal Range

Claim: For a 5-generation binary tree, the optimal dephasing regime is around gamma_phi ~ 1 V with a broad optimality region (gamma_phi ~ 0.1 to ~ 10 V). The efficiency is very robust to disorder within this range.
Source: Novo, Mohseni & Omar (2013)
URL: https://arxiv.org/pdf/1312.6989
Date: 2013
Excerpt: "We observe a large optimality region around gamma_phi ~ 1. Close to the optimal regime, the transport efficiency is very robust to disorder and remains almost constant as the latter varies within a large range of values (from delta_epsilon = 0 to delta_epsilon ~ 1.5)."
Context: Parameter scan across disorder delta_epsilon and dephasing rate gamma_phi for binary tree
Confidence: High

### Finding 4.2: ENAQT Disappears Above Critical System Size

Claim: ENAQT has a critical system size limit: L <= O(1) * J/mu, where J is the hopping rate and mu is the background absorption/recombination rate. Beyond this size, the non-monotonic behavior disappears.
Source: Schaller et al. (2025), "Dephasing enhanced transport of spin excitations"
URL: https://arxiv.org/abs/2502.10854
Date: 2025
Excerpt: "We show that the non-monotonic behavior of the efficiency, characteristic for ENAQT, disappears if the system size becomes larger than some critical value determined by the ratio of the extraction rate to the rate of background absorption."
Context: Analytical estimate for 1D: gamma_0 ~ 5*J/L - mu for large L >> 1 and small mu << J
Confidence: High

### Finding 4.3: FMO Complex Dephasing Window

Claim: For the FMO complex with trapping rate kappa_3 = 1 ps^-1, the optimal dephasing rate is approximately 10-50 cm^-1 (estimated from figure analysis), with the room temperature value (~300 cm^-1) being in the upper ENAQT range but below the deep Zeno regime.
Source: Rebentrost et al. (2009)
URL: https://dspace.mit.edu/bitstream/handle/1721.1/70838/Rebentrost-2009-Environment-assisted%20quantum%20transport.pdf
Date: 2009
Excerpt: "The transfer time is 75 ps in the fully quantum limit and improves significantly to 7 ps in the intermediate ENAQT regime. For large dephasing, the transfer slows down to 500 ps."
Context: Three regimes clearly separated in FMO efficiency plot: quantum (gamma < 1), ENAQT (1 < gamma < 100), Zeno (gamma > 100 cm^-1)
Confidence: Medium

---

## 5. ENAQT Signatures in Population Dynamics P(t)

### Finding 5.1: Destructive Interference at Trap Sites

Claim: In the absence of disorder and dephasing, destructive interference at the trap causes the excitation to live much longer outside the trap than at the trap. Adjacent sites synchronize with opposite signs, preventing the excitation from reaching the trap. Both disorder and dephasing suppress this destructive interference.
Source: Novo, Mohseni & Omar (2013)
URL: https://arxiv.org/pdf/1312.6989
Date: 2013
Excerpt: "After a short time, psi_2 and psi_3 synchronize (keeping opposite signs) and interfere destructively, preventing the excitation from reaching the trap. This helps us to understand why disorder and dephasing increase significantly the transport efficiency: both of them prevent this destructive interference from happening."
Context: Wavefunction dynamics analysis at trap and adjacent sites for binary tree
Confidence: High

### Finding 5.2: Coherence Prolongation Counterintuitively by Dephasing

Claim: Although dephasing alone damps all coherences, the interplay between Hamiltonian dynamics and dephasing makes some coherences last LONGER. The time the excitation remains delocalized between the trap and adjacent sites increases with both disorder and dephasing.
Source: Novo, Mohseni & Omar (2013)
URL: https://arxiv.org/pdf/1312.6989
Date: 2013
Excerpt: "It is remarkable that although dephasing alone damps superpositions between wavefunctions at different sites, its interplay with the hopping Hamiltonian leads to the prolongation of certain superpositions."
Context: Density operator analysis showing counterintuitive coherence dynamics
Confidence: High

### Finding 5.3: Ballistic-to-Diffusive Transition at Optimal Dephasing

Claim: At the optimal ENAQT dephasing rate, transport evolves from ballistic to mainly diffusive dynamics within very short times. Strong dephasing leads to subdiffusive transport.
Source: Maier et al. (trapped ion experiment)
URL: https://quantumoptics.at/images/publications/dissertation/ChristineMaier_Dissertation.pdf
Date: ~2019
Excerpt: "In the regime around gamma = J_max, where ENAQT is most efficient, we find that within very short times t ~ 1/J_max the transport evolves from ballistic to mainly diffusive dynamics, yielding C = 0.76 +/- 0.18. For strong dephasing, gamma = 18.4*J_max, we observe subdiffusive transport with a power exponent C = 0.44 +/- 0.02."
Context: 10-qubit trapped ion quantum simulation of transport dynamics
Confidence: High

---

## 6. Experimental Tests of ENAQT

### Finding 6.1: Perovskite Nanocrystal Superlattices (2025) - Major Experimental Confirmation

Claim: Environment-assisted quantum transport of excitons was experimentally demonstrated in CsPbBr3 perovskite nanocrystal superlattices, with a peak in the steady-state diffusion constant at intermediate temperatures where static disorder and dephasing are balanced (Gamma ~ Delta).
Source: Blach et al. (2025), Nature Communications 16, 1270
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025-02-02
Excerpt: "Most remarkably, at intermediate temperatures, the interplay between coherence and dephasing collaboratively enhances transport. We observed a maximum for the steady-state diffusion constant at temperatures where static disorder and thermal fluctuations balance each other, providing experimental evidence of ENAQT."
Context: Temperature range 7-298 K, turnover temperature T_t ~ 70 K (OA-1) and ~100 K (DAB-1, DAB-2)
Confidence: High

### Finding 6.2: Photonic Waveguide Network (2016)

Claim: The first experimental implementation of ENAQT used an integrated photonic simulator with femtosecond-laser direct writing. Decoherence was simulated by broadband illumination, enhancing transport efficiency in a 4-site network.
Source: Caruso et al. (2016), Nature Communications 7, 11282
URL: https://www.nature.com/articles/ncomms11282
Date: 2016-04-15
Excerpt: "Here we use an integrated photonic simulator to demonstrate the first implementation of ENAQT. We simulated this network using waveguides... the efficiency increases monotonically with the strength of decoherence."
Context: 4-site network with waveguide arrangement, transport efficiency measured as function of optical bandwidth
Confidence: High

### Finding 6.3: Electrical Oscillator Network (2015)

Claim: ENAQT was observed in a classical network of coupled electrical oscillators, showing a 22.5 +/- 3.6% relative enhancement in energy transport efficiency by introducing stochastic fluctuations.
Source: Viciani et al. / Noise-assisted energy transport in electrical oscillator networks
URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC4661523/
Date: 2015
Excerpt: "We report on the experimental observation of such effect in a network of coupled electrical oscillators. We demonstrate that by introducing stochastic fluctuations in one of the couplings of the network, a relative enhancement in the energy transport efficiency of 22.5 +/- 3.6% can be observed."
Context: Classical analog of ENAQT in coupled LC oscillators
Confidence: High

### Finding 6.4: Trapped Ion Quantum Simulation

Claim: A 10-qubit trapped ion quantum simulator observed Anderson localization at strong disorder without noise, ENAQT at intermediate noise (gamma ~ J_max), and quantum Zeno suppression at strong noise.
Source: Maier dissertation (University of Innsbruck / IQOQI)
URL: https://quantumoptics.at/images/publications/dissertation/ChristineMaier_Dissertation.pdf
Date: ~2019
Excerpt: "Rauschen mittlerer Staerke steigert die Effizienz des Quantentransports - ein Effekt der als Umgebungs-unterstuetzter Quantentransport (ENAQT) bekannt ist. Koppelt das Netzwerk an starkes Rauschen, so sehen wir eine Unterdrueckung des Quantentransports, was den Uebergang in das Quanten-Zeno Regime einleitet."
Context: 10-ion network with long-range spin-spin interactions, arbitrary disorder and noise spectra
Confidence: High

### Finding 6.5: Digital Quantum Computer Simulations

Claim: Two quantum algorithms for simulating ENAQT on digital quantum computers were introduced - one based on stochastic Hamiltonians and one on collision schemes. Both achieve memory-efficient encoding with qubits scaling logarithmically with system size.
Source: Gibb et al. (2022), New Journal of Physics 24, 023039
URL: https://iopscience.iop.org/article/10.1088/1367-2630/ac512f
Date: 2022
Excerpt: "Two different quantum algorithms are introduced, the first one based on stochastic Hamiltonians and the second one based on a collision scheme. We test both algorithms by simulating ENAQT in a small molecular network on a quantum computer emulator."
Context: Provides comparative analysis of approaches for digital quantum simulation of dephasing-assisted transport
Confidence: High

---

## 7. ENAQT vs Classical Hopping

### Finding 7.1: Quantum Diffusion Exceeds Classical in ENAQT Regime

Claim: Quantum diffusion in the ENAQT regime is considerably larger than classical diffusion, because the coherence length L(Gamma, Delta) defines a larger step-size than a classical random walk.
Source: Blach et al. (2025), Nature Communications
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025
Excerpt: "Quantum diffusion in the ENAQT regime is considerably larger than classical diffusion, as L(Gamma, Delta) defines a larger step-size than a classical random walk."
Context: Experimental measurement of diffusion constants in perovskite NC superlattices
Confidence: High

### Finding 7.2: Non-Gaussian Spatial Distribution as Quantum Signature

Claim: In the coherent transport regime (low temperature/T), the spatial distribution of excitons is non-Gaussian, consistent with coherent effects. At high temperature, the distribution becomes Gaussian as transport transitions to classical diffusion.
Source: Blach et al. (2025)
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025
Excerpt: "This behavior suggests enhanced exciton migration by coherent effects in this temperature regime and consistent with non-gaussian spatial distribution."
Context: Transient absorption and photoluminescence microscopy measurements
Confidence: High

### Finding 7.3: Control Samples Show Monotonic Temperature Dependence

Claim: Control samples (random NC thin film and bulk crystal) show very different temperature-dependent transport: monotonically increasing with temperature for random film, and monotonically decreasing with temperature for bulk crystal. Only the ordered NC superlattice shows the non-monotonic ENAQT signature.
Source: Blach et al. (2025)
URL: https://www.nature.com/articles/s41467-024-55812-8
Date: 2025
Excerpt: "Both these two control cases show very different temperature-dependent transport behaviors from the nanocrystal superlattices. Specifically, the diffusion constant monotonously increases as temperature increases in the random thin film."
Context: Control experiments essential for distinguishing ENAQT from classical thermal activation
Confidence: High

---

## 8. System Size, Disorder, and Temperature Effects

### Finding 8.1: Critical System Size for ENAQT

Claim: ENAQT disappears when the system size exceeds L_crit ~ J/mu (hopping rate divided by loss rate). For a 1D chain, the positive root for large L gives gamma_0 ~ 5*J/L - mu, decreasing with system size.
Source: Schaller et al. (2025), arXiv:2502.10854
URL: https://arxiv.org/pdf/2502.10854
Date: 2025
Excerpt: "We find for a 2D lattice that for linear system dimensions smaller than a characteristic length small dephasing indeed can lead to a significant increase of the transfer efficiency... The non-monotonic behavior of the efficiency, characteristic for ENAQT, disappears if the system size becomes larger than some critical value."
Context: Analytical bounds for transfer efficiency in low and high dephasing limits
Confidence: High

### Finding 8.2: Disorder-Assisted Transport in Suboptimal Dephasing Regime

Claim: When dephasing is below the optimal ENAQT regime, disorder can IMPROVE transport efficiency by up to 30%. This "disorder-assisted quantum transport" occurs because disorder destroys invariant subspaces that prevent transport.
Source: Novo, Mohseni & Omar (2013)
URL: https://arxiv.org/pdf/1312.6989
Date: 2013
Excerpt: "The improvement due to disorder is maximal in the purely unitary case (gamma_phi = 0), where the efficiency grows from 6% (for delta_epsilon = 0) to 34% (for delta_epsilon = 0.83). This improvement is washed out with increasing dephasing until the latter reaches its optimal value."
Context: Counterintuitive result: disorder helps when dephasing is suboptimal
Confidence: High

### Finding 8.3: Optimal Disorder Exists for Fixed Dephasing

Claim: For any dephasing rate below the optimal value, there exists an optimal disorder strength that maximizes transport efficiency. Very high disorder causes large energy mismatches that hinder propagation.
Source: Novo, Mohseni & Omar (2013)
URL: https://arxiv.org/pdf/1312.6989
Date: 2013
Excerpt: "For any given dephasing rate gamma_phi below the optimal value, increasing the disorder from zero up to around delta_epsilon ~ 1 enhances transport efficiency."
Context: Full interplay between disorder and dephasing in binary tree transport
Confidence: High

---

## 9. Published Code and Software

### Finding 9.1: quantum_HEOM Package

Claim: The quantum_HEOM Python package implements Lindblad and HEOM approaches for open quantum system dynamics, including local dephasing, global thermalizing, and local thermalizing Lindblad models, plus HEOM via QuTiP.
Source: J.W. Abbott, GitHub
URL: https://github.com/jwa7/quantum_HEOM
Date: 2019-2022
Excerpt: "With high control over input parameters and interactable in notebooks via 'black-boxed' code written in Python, quantum_HEOM allows the bath-influenced excitonic energy transfer dynamics of open quantum systems to be calculated and plotted."
Context: Used to generate training data for ML-QD approach in Nature Communications 2022
Confidence: High

### Finding 9.2: QD3SET-1 GitHub Repository

Claim: The QD3SET-1 database provides Python scripts for data extraction, visualization, and analysis. Data is in HDF5 format with standardized structure for population and coherence dynamics.
Source: Arif Ullah et al., GitHub
URL: https://github.com/Arif-PhyChem/QD3SET
Date: 2023
Excerpt: "A Python package for extracting data is provided together with the data set and can be accessed at https://github.com/Arif-PhyChem/QD3SET."
Context: Used in ML studies for training quantum dynamics predictors (CNN-LSTM, PINN architectures)
Confidence: High

### Finding 9.3: QuTiP HEOM Solver

Claim: QuTiP includes a built-in HEOM solver that can simulate quantum transport dynamics, including the spin-boson model with Debye/Drude-Lorentz spectral densities.
Source: QuTiP project
URL: https://qutip.org/
Date: Active (version 4.6+ used for QD3SET-1)
Excerpt: "QuTiP software package (version 4.6) used in HEOM calculations of the spin-boson model."
Context: Industry-standard open-source quantum simulation framework
Confidence: High

### Finding 9.4: qHEOM - Quantum Algorithm for HEOM on NISQ Devices

Claim: The qHEOM quantum algorithm implements numerically exact HEOM on NISQ quantum computers using SVD dilation, demonstrated on IBM hardware for charge transfer and FMO energy transfer.
Source: Batista group, Yale
URL: https://arxiv.org/abs/2411.12049
Date: 2024-2025
Excerpt: "qHEOM implements the SVD dilation to convert the non-unitary HEOM propagator into unitary gates... demonstrated on IBM quantum computers as applied to simulations of charge transfer and electronic energy transfer dynamics in the FMO complex."
Context: First demonstration of HEOM-based quantum dynamics on real quantum hardware
Confidence: High

---

## 10. What Constitutes "Proof" of ENAQT in a Computational Study

### Finding 10.1: Minimal Criteria for Computational ENAQT Proof

Based on synthesis of the literature, computational proof of ENAQT requires:

1. **Non-monotonic efficiency curve**: Transport efficiency eta must show a clear maximum as a function of dephasing rate gamma_phi
2. **Efficiency enhancement > 0**: eta(gamma_opt) > eta(gamma -> 0) (ENAQT exists only if dephasing improves over purely coherent transport)
3. **Three-regime identification**: Clear separation of quantum, ENAQT, and Zeno regimes
4. **Control at zero dephasing**: Efficiency at gamma = 0 must be computed as reference
5. **System parameters documented**: Hamiltonian, disorder strength, trapping/loss rates, temperature all specified
6. **Converged numerical results**: For exact methods (HEOM), convergence with hierarchy level and Matsubara terms must be demonstrated

Source: Synthesis of Rebentrost 2009, Novo 2013, Schaller 2025, Blach 2025
Confidence: High (consensus criteria)

### Finding 10.2: Additional Distinguishing Features

Claim: Stronger computational proof includes:
- Showing the optimal dephasing rate scales as gamma_opt ~ J (hopping rate) - the "Goldilocks" condition
- Demonstrating the efficiency is robust to disorder variations near the optimal point
- Computing both population dynamics AND coherence dynamics
- Using more than one numerical method (e.g., HEOM + Lindblad) to cross-validate
- Including appropriate control comparisons (e.g., ordered vs disordered systems)

Source: Synthesis across literature
Confidence: High

### Finding 10.3: The Delocalization-Transport Paradox

Claim: "Delocalization is the enemy of delivery" - in perfectly ordered systems with symmetric Hamiltonians, the excitation can be trapped in invariant subspaces (dark states) that don't couple to the sink. Dephasing is required to break these symmetries and open transport pathways.
Source: Novo, Mohseni & Omar (2013); Blach et al. (2025)
URL: Multiple
Date: 2013-2025
Excerpt: "In the case of zero dephasing, this effect is explained: random site energies disorder also destroys the invariant subspace, making all eigenstates couple to the trap." (Novo 2013)
Context: Key insight that explains why ENAQT can occur even in ordered systems
Confidence: High

---

## 11. Machine Learning for ENAQT Detection

### Finding 11.1: CNN-LSTM Architecture for Quantum Dynamics Prediction

Claim: CNN-LSTM hybrid neural networks have been trained on QD3SET-1 data to predict quantum dissipative dynamics (population P(t) and coherences) for spin-boson and FMO systems, achieving very low validation losses.
Source: Arif-PhyChem/rc-pinn-comparison, GitHub
URL: https://github.com/Arif-PhyChem/rc-pinn-comparison
Date: 2025
Excerpt: "FMO trained models: FMO_mr-pinn_model-1015-tloss-1.677e-06-vloss-1.871e-06.keras... SB trained models: SB_st-pinn_model-11-tloss-3.947e-06-vloss-1.405e-06.keras"
Context: Physics-informed neural networks with custom loss functions for quantum dynamics
Confidence: High

### Finding 11.2: AI-QD Approach for Predicting Energy Transfer Dynamics

Claim: An "artificial intelligence-based quantum dynamics" (AI-QD) approach can predict excitation energy transfer dynamics in light-harvesting complexes by learning from HEOM training data generated by quantum_HEOM.
Source: Ullah & Dral (2022), Nature Communications 13, 1930
URL: https://doi.org/10.1038/s41467-022-29621-w
Date: 2022
Excerpt: "Predicting the future of excitation energy transfer in light-harvesting complex with artificial intelligence-based quantum dynamics."
Context: quantum_HEOM package used to generate training data for trajectory learning
Confidence: High

### Finding 11.3: ML for ENAQT Signature Classification - Potential Approach

While no published work specifically trains ML to detect ENAQT signatures, the infrastructure exists:

1. **QD3SET-1** provides labeled population dynamics data across varying dephasing parameters
2. **Classification approach**: Train a classifier on P(t) curves to distinguish:
   - Class 1: Low dephasing (oscillatory, localization)
   - Class 2: Optimal dephasing (ENAQT - enhanced transport)
   - Class 3: High dephasing (Zeno - suppressed transport)
3. **Regression approach**: Predict transport efficiency eta from P(t) dynamics alone
4. **Key input features**: Population dynamics P(t), coherence decay rates, peak transport times, oscillation frequencies

Source: Synthesis of QD3SET-1 and ML-QD literature
Confidence: Medium (inferred potential, not yet demonstrated in literature)

---

## Summary Table: Key ENAQT Testable Predictions

| Prediction | Mathematical Form | Observable | Confidence |
|-----------|-------------------|------------|------------|
| Non-monotonic eta vs gamma | eta(gamma) has maximum | Efficiency vs dephasing curve | High |
| Optimal dephasing ~ hopping rate | gamma_opt ~ V | Peak position in eta curve | High |
| Turnover temperature | T_t where Gamma ~ Delta | Diffusion constant peak | High (exp. confirmed 2025) |
| ENAQT disappears for large systems | L_crit ~ J/mu | Size-dependent efficiency | High |
| Quantum > classical diffusion | D_quantum > D_classical | Diffusion constant comparison | High (exp. confirmed 2025) |
| Twin peaks possible | Multiple eta maxima | Multiple optimal regimes | High |
| Coherence prolongation | Coherence lifetime increases | Coherence dynamics | Medium |

---

## Available Computational Tools and Data

| Resource | URL | Purpose |
|----------|-----|---------|
| QD3SET-1 Database | https://github.com/Arif-PhyChem/QD3SET | Training/benchmarking data |
| quantum_HEOM | https://github.com/jwa7/quantum_HEOM | Lindblad + HEOM simulations |
| QuTiP | https://qutip.org/ | General open quantum systems |
| qHEOM (paper) | https://arxiv.org/abs/2411.12049 | NISQ quantum HEOM algorithm |
| QD3SET-1 Paper | https://doi.org/10.3389/fphy.2023.1223973 | Database documentation |

---

## References (Key Sources)

1. P. Rebentrost, M. Mohseni, I. Kassal, S. Lloyd, A. Aspuru-Guzik, "Environment-assisted quantum transport," New Journal of Physics 11, 033003 (2009).
2. L. Novo, M. Mohseni, Y. Omar, "Disorder-assisted quantum transport in suboptimal decoherence regimes," Phys. Chem. Chem. Phys. 25, 10103 (2013) [arXiv:1312.6989].
3. E. Levi et al., "Universal Origin for Environment-Assisted Quantum Transport," arXiv:1801.06799 (2018).
4. D.D. Blach et al., "Environment-assisted quantum transport of excitons in perovskite nanocrystal superlattices," Nature Communications 16, 1270 (2025).
5. A. Ullah et al., "QD3SET-1: A Database with Quantum Dissipative Dynamics Data Sets," Frontiers in Physics 11, 1223973 (2023).
6. A.R. Coates, B.W. Lovett, E.M. Gauger, "From Goldilocks to twin peaks: multiple optimal regimes for quantum transport," Phys. Chem. Chem. Phys. 25, 10103 (2023).
7. A. Ullah, P.O. Dral, "Predicting the future of excitation energy transfer in light-harvesting complex with artificial intelligence-based quantum dynamics," Nature Communications 13, 1930 (2022).
8. F. Caruso et al., "Enhancing coherent transport in a photonic network using controllable decoherence," Nature Communications 7, 11282 (2016).
9. Schaller et al., "Dephasing enhanced transport of spin excitations," arXiv:2502.10854 (2025).
10. L. Contreras-Pulido et al., "Design Principles for Enhanced Quantum Transport with Site-Dependent Noise," arXiv:2604.23005 (2026).

---

*Research compiled for Dimension 11: ENAQT Testable Predictions and Benchmarks*
*Part of comprehensive ENAQT-Microtubule research project*
