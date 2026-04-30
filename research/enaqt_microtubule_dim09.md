# Dimension 9: Open Quantum System Tools Beyond QuTiP

## Comprehensive Survey of Computational Tools for Simulating Open Quantum Systems in Biological Quantum Transport

**Date compiled:** 2025-07-28
**Searches performed:** 20+ independent queries across web search engines
**Scope:** Python, Julia, C++, and GPU-accelerated tools for open quantum dynamics with emphasis on biological quantum transport

---

## Table of Contents

1. [Python-Based Process Tensor & Non-Markovian Solvers](#1-process-tensor--non-markovian-solvers)
2. [GPU-Accelerated & Differentiable Quantum Simulators](#2-gpu-accelerated--differentiable-simulators)
3. [HEOM Solvers and Implementations](#3-heom-solvers)
4. [Julia Ecosystem for Open Quantum Systems](#4-julia-ecosystem)
5. [Machine Learning for Quantum Dynamics](#5-ml-for-quantum-dynamics)
6. [Neural Quantum States](#6-neural-quantum-states)
7. [Quantum Machine Learning Platforms](#7-quantum-machine-learning-platforms)
8. [Spectral Density Fitting Tools](#8-spectral-density-fitting-tools)
9. [2D Electronic Spectroscopy Simulation](#9-2d-electronic-spectroscopy)
10. [Photosynthetic Complex Simulation Packages](#10-photosynthetic-complex-simulators)
11. [Workflow Management Tools](#11-workflow-management-tools)
12. [Tools Assessment Summary](#12-assessment-summary)

---

## 1. Process Tensor & Non-Markovian Solvers

### 1.1 OQuPy (Open Quantum System in Python)

```
Claim: OQuPy is a major open-source Python package that uses the process tensor approach with matrix product operator (PT-MPO) representations to efficiently simulate non-Markovian open quantum systems. It provides numerically exact methods for dynamics, multi-time correlations, control optimization, interacting chains, and mean-field dynamics of ensembles.
Source: Journal of Chemical Physics / arXiv
URL: https://arxiv.org/abs/2406.16650
Date: 2024-06-24
Excerpt: "We present a major release of our open source software package, OQuPy (Open Quantum System in Python), which provides several recently developed numerical methods that address this challenging task. It utilizes the process tensor approach to open quantum systems in which a single map, the process tensor, captures all possible effects of an environment on the system."
Context: OQuPy was formerly known as TEMPO. The package implements PT-TEMPO for Gaussian bosonic environments and can import PT-MPOs from other methods like ACE. It is particularly well-suited for quantum chemistry and quantum sensing applications.
Confidence: high
```

```
Claim: OQuPy's PT-TEMPO method allows systematic removal of negligible correlations between system and environment, with computational cost controlled by bond dimension reflecting non-Markovianity rather than environment Hilbert space dimension or correlation time.
Source: OQuPy paper (JCP 2024)
URL: https://pubs.aip.org/aip/jcp/article/161/12/124108/3313792
Date: 2024-09-24
Excerpt: "The necessary bond dimension of a PT-MPO reflects the degree of non-Markovianity in the interaction, but does not necessarily scale with the correlation time or the dimension of the environment Hilbert space."
Context: This makes OQuPy particularly efficient for structured environments where memory effects are important but not arbitrarily complex, such as photosynthetic complexes.
Confidence: high
```

**Installation:** `pip install oqupy`
**Documentation:** https://oqupy.github.io/
**GitHub:** https://github.com/oqupy/oqupy

---

### 1.2 ACE (Automated Compression of Environments)

```
Claim: ACE is a numerically exact method for simulating open quantum systems coupled to arbitrary environments of independent degrees of freedom (bosonic, fermionic, or spin). It iteratively constructs and compresses process tensors using matrix product state techniques, automatically projecting the environment onto its most relevant degrees of freedom.
Source: Nature Physics / arXiv
URL: https://arxiv.org/abs/2101.01653
Date: 2021-01-05
Excerpt: "We present a numerically exact method for simulating open quantum systems with arbitrary environments which consist of a set of independent degrees of freedom. Our approach automatically reduces the large number of environmental degrees of freedom to those which are most relevant."
Context: ACE can handle non-Gaussian baths (e.g., anharmonic environments, spin baths) that are beyond the scope of methods like HEOM which assume Gaussian environments. The C++ toolkit was released in 2024.
Confidence: high
```

```
Claim: ACE has been demonstrated on electron transport in quantum dots, phonon effects, radiative decay, central spin dynamics, anharmonic environments (Morse potential), dispersive coupling to lossy cavities, and superradiance.
Source: Nature Physics
URL: https://research-repository.st-andrews.ac.uk/bitstream/10023/26069/1/ACE_accepted.pdf
Date: 2022
Excerpt: "The versatility and efficiency of our automated compression of environments (ACE) method provides a practical general-purpose tool for open quantum systems."
Context: ACE is written in C++ and fully controllable by configuration files. It can be combined with OQuPy through PT-MPO import/export.
Confidence: high
```

**GitHub:** https://github.com/ace-oqs/ace (C++ toolkit released 2024)
**Paper:** Cygorek et al., "Numerically exact open quantum systems simulations for arbitrary environments using automated compression of environments", Nature Physics 2021

---

## 2. GPU-Accelerated & Differentiable Simulators

### 2.1 Dynamiqs

```
Claim: Dynamiqs is an open-source Python library for GPU-accelerated and differentiable quantum simulations built on JAX. It provides solvers for the Schrodinger equation, Lindblad master equation, and stochastic master equation, with up to 60x speedup over QuTiP for large systems.
Source: Alice & Bob blog / arXiv
URL: https://alice-bob.com/blog/dynamiqs-gpu-opensource-quantum-simulation-library/
Date: 2024-11-20
Excerpt: "Dynamiqs is a Python library that enables the simulation of time-dependent quantum systems. It is mainly designed to solve differential equations of quantum dynamics such as the Schrodinger equation, Lindblad master equation, stochastic master equation and their variants."
Context: Development started in early 2023 (originally PyTorch), completely rewritten in JAX in early 2024. Sponsored by Alice & Bob (superconducting cat qubit startup). Benchmarks show up to 30x speedup on GPU vs QuTiP CPU for dissipative cat CNOT, and 60x with sparse representation in v0.3.0.
Confidence: high
```

```
Claim: Dynamiqs features end-to-end differentiability via JAX automatic differentiation, enabling quantum optimal control, parameter estimation, and state tomography with machine-precision gradient accuracy.
Source: Dynamiqs documentation
URL: https://www.dynamiqs.org
Date: 2025
Excerpt: "Differentiability is again enabled by the underlying JAX codebase and by Diffrax, a state-of-the-art library for the integration of differential equations. Diffrax features memory-efficient automatic differentiation methods for large matrices."
Context: Dynamiqs supports batching (running thousands of simulations concurrently), switching between CPU/GPU with `dq.set_device`, and is QuTiP-compatible (accepts QuTiP objects as arguments).
Confidence: high
```

**Installation:** `pip install dynamiqs`
**Documentation:** https://www.dynamiqs.org
**GitHub:** https://github.com/dynamiqs/dynamiqs (2,000+ stars)
**Citation:** P. Guilmin, A. Bocquet, E. Genois, D. Weiss, R. Gautier. *Dynamiqs: an open-source Python library for GPU-accelerated and differentiable simulation of quantum systems* (2025), in preparation.

---

### 2.2 jaxquantum

```
Claim: jaxquantum leverages JAX to enable auto-differentiable and (CPU, GPU, TPU) accelerated simulation of quantum dynamical systems, serving as a QuTiP drop-in replacement written entirely in JAX. It has absorbed the bosonic and qcsys packages.
Source: GitHub / jaxquantum.org
URL: https://github.com/EQuS/jaxquantum
Date: 2022-07-10
Excerpt: "jaxquantum leverages JAX to enable the auto differentiable and (CPU, GPU, TPU) accelerated simulation of quantum dynamical systems, including tooling such as operator construction, unitary evolution and master equation solving. As such, jaxquantum serves as a QuTiP drop-in replacement written entirely in JAX."
Context: Developed by the EQuS group. jaxquantum is now a unified toolkit for quantum circuit design, simulation and control.
Confidence: medium
```

**Documentation:** https://equs.github.io/jaxquantum

---

## 3. HEOM Solvers

### 3.1 QuTiP-BoFiN (Integrated into QuTiP 5)

```
Claim: QuTiP-BoFiN is a numerical HEOM library integrated with the QuTiP platform, implementing the hierarchical equations of motion for both bosonic and fermionic environments. It includes spectral density fitting capabilities and applications to light-harvesting complexes.
Source: Physical Review Research
URL: https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.013181
Date: 2023-03-15
Excerpt: "Here we present a numerical library in Python, integrated with the powerful QuTiP platform, which implements the HEOM for both bosonic and fermionic environments. We demonstrate its utility with a series of examples consisting of benchmarks against important known results and examples demonstrating insights gained with this library."
Context: BoFiN was originally a separate package (https://github.com/tehruhn/bofin) but has been fully integrated into QuTiP 5 with greater functionality. The paper includes FMO complex dynamics, spectral density fitting, and dynamical decoupling benchmarks.
Confidence: high
```

**Documentation:** https://qutip.readthedocs.io/en/latest/guide/guide-heom.html
**GitHub (original):** https://github.com/tehruhn/bofin (now unmaintained; use QuTiP 5)

---

### 3.2 GPU-HEOM

```
Claim: GPU-HEOM is a CUDA-based implementation of the hierarchical equations of motion that achieves tremendous speedup by exploiting graphics processing units. It was used for systematic studies of energy transfer efficiency in the FMO complex at physiological temperature with full memory effects.
Source: Journal of Chemical Theory and Computation
URL: https://pubs.acs.org/doi/10.1021/ct200126d
Date: 2011-07-01
Excerpt: "Here, we show that HEOM are also solvable for larger systems, since the underlying algorithm is ideally suited for the usage of graphics processing units (GPU). The tremendous reduction in computational time due to the GPU allows us to perform a systematic study of the energy-transfer efficiency in the Fenna-Matthews-Olson (FMO) light-harvesting complex at physiological temperature under full consideration of memory effects."
Context: GPU-HEOM was developed by Tobias Kramer et al. It is publicly available as a ready-to-run cloud computing tool hosted on the nanoHUB platform. The code can handle systems up to ~4096 energy states.
Confidence: high
```

**nanoHUB:** https://nanohub.org/resources/gpuheompop
**Reference:** Kramer et al., *J. Chem. Theory Comput.* **7**, 7 (2011)

---

### 3.3 PHI (Parallel Hierarchy Integrator)

```
Claim: PHI is a multi-threaded software implementation of HEOM for shared-memory parallel computers, capable of handling large photosynthetic systems like LH2 (B850 ring) and LH1 (B875 ring) combined 50-sites systems requiring up to 108 GB of memory.
Source: Journal of Chemical Theory and Computation / UIUC
URL: https://www.ks.uiuc.edu/Research/phi/
Date: 2012-07-23
Excerpt: "PHI is a multi-threaded software implementation of a hierarchy equation of motion (HEOM) integrator to run on shared memory parallel computers. The HEOM is a method to calculate the noise-averaged time evolution of density matrix that describes a quantum system in contact with a thermal environment."
Context: PHI uses adaptive timestep integration (RKF45) and a unique partitioning scheme that exploits the communication structure of the HEOM (d-dimensional Pascal's simplex). It is tied to structureless Drude-Lorentz spectral densities.
Confidence: high
```

**Paper:** Strumpfer & Schulten, *J. Chem. Theory Comput.* **8**, 2808 (2012)

---

### 3.4 QMaster

```
Claim: QMaster is a high-performance HEOM computation tool using OpenCL that runs on multiple accelerator architectures (multicore CPUs, GPUs, Intel Xeon Phi), maintaining a common code basis across all devices. It overcomes GPU memory limitations by leveraging the larger system memory of CPU-based accelerators.
Source: Journal of Chemical Theory and Computation
URL: https://pubs.acs.org/doi/10.1021/ct500629s
Date: 2015
Excerpt: "Here, we introduce the QMaster-tool for high-performance computations of open-system dynamics across a wide range of parallel accelerators, including large-memory devices."
Context: QMaster follows a different strategy than GPU-HEOM: instead of multi-device distribution, it uses OpenCL to run on any accelerator with larger memory (e.g., CPU workstations with 128+ GB RAM), enabling larger spectral densities with many vibrational peaks.
Confidence: high
```

---

### 3.5 pyheom / libheom

```
Claim: pyheom is an open-source Python 3 library for HEOM-based open quantum dynamics simulations, providing a binding of the C++17/CUDA libheom library. It supports generalization of HEOM for efficient calculations with arbitrary correlation functions.
Source: GitHub / J. Chem. Phys.
URL: https://github.com/tatsushi-ikeda/pyheom
Date: 2020-08-31
Excerpt: "pyheom is an open-source library that supports open quantum dynamics simulations based on the hierarchical equations of motion (HEOM) theory. This library provides a python 3 binding of libheom (pylibheom) and high-level APIs."
Context: Written by Tatsushi Ikeda and Gregory D. Scholes. libheom provides low-level C++17/CUDA APIs. The paper generalized HEOM for arbitrary correlation functions beyond standard Drude-Lorentz baths.
Confidence: high
```

**Citation:** Ikeda & Scholes, *J. Chem. Phys.* **152**, 204101 (2020)

---

### 3.6 fHEOM (Factorized HEOM)

```
Claim: fHEOM implements low-rank factorization methods for correlated baths in HEOM simulations, reducing spatially correlated multi-mode bath problems to a smaller number of effective modes through eigendecomposition.
Source: GitHub
URL: https://github.com/rihp/fHEOM
Date: 2025-10-20
Excerpt: "fHEOM v1.0 — Low-rank factorization method for correlated baths in quantum simulations. This framework is open-source to enable community validation and replication, extension to other quantum biology systems, performance optimization and improvements, transparent, reproducible science."
Context: Provides ready-to-use FMO complex Hamiltonians and bath operators. Built on top of QuTiP 5.2+ HEOM solver. Supports spatial correlation kernels (exponential, Gaussian, power-law).
Confidence: medium
```

---

### 3.7 DM-HEOM (Distributed Memory HEOM)

```
Claim: DM-HEOM is a portable and scalable solver framework for HEOM using distributed memory (MPI), enabling computations on large systems across multiple compute nodes.
Source: GitHub / IEEE IPDPSW
URL: https://github.com/noma/dm-heom
Date: 2018
Excerpt: "DM-HEOM: A portable and scalable solver-framework for the hierarchical equations of motion"
Context: Described in Noack et al., IEEE IPDPSW 2018. Designed for HPC clusters with MPI distributed memory parallelism.
Confidence: medium
```

---

## 4. Julia Ecosystem

### 4.1 QuantumOptics.jl

```
Claim: QuantumOptics.jl is a numerical framework written entirely in Julia for simulating open quantum systems, offering speed comparable to compiled languages. It solves Schrodinger equation, master equation, and Monte Carlo wave-function approaches.
Source: Computer Physics Communications / arXiv
URL: https://arxiv.org/abs/1707.01060
Date: 2017-07-04 / 2018-06-01
Excerpt: "We present an open source computational framework geared towards the efficient numerical investigation of open quantum systems written in the Julia programming language. Built exclusively in Julia and based on standard quantum optics notation, the toolbox offers speed comparable to low-level statically typed languages, without compromising on the accessibility and code readability found in dynamic languages."
Context: Inspired by QuTiP. Related packages include Correlation Expansion Package and CollectiveSpins library. MIT licensed.
Confidence: high
```

**Documentation:** https://docs.qojulia.org/
**GitHub:** https://github.com/qojulia/QuantumOptics.jl
**Citation:** Kramer et al., *Comput. Phys. Commun.* **227**, 109 (2018)

---

### 4.2 OpenQuantumSystems.jl

```
Claim: OpenQuantumSystems.jl is a Julia framework specifically focused on quantum biology in finite basis, inspired by QuantumOptics.jl. It supports exact dynamics, Schrodinger and Liouville-von Neumann equations, and various quantum master equation variants for molecular aggregates.
Source: GitHub
URL: https://github.com/detrin/OpenQuantumSystems.jl
Date: 2020-07-15
Excerpt: "OpenQuantumSystems.jl is a numerical framework written in Julia that makes it easy to simulate various kinds of open quantum systems with main focus on quantum biology in finite basis. It is inspired by the QuantumOptics.jl."
Context: Capabilities include: Hamiltonian construction for molecular aggregates, bath as LHOs with shifted potentials, dynamics via exact/SE/LvN/QME methods, iterative bath correction, memory kernel as superoperator, and convenience constructors for dimers/trimers/chains.
Confidence: medium
```

**Installation:** `] add OpenQuantumSystems`

---

### 4.3 QuantumToolbox.jl

```
Claim: QuantumToolbox.jl is an efficient Julia framework for simulating open quantum systems, including stochastic Schrodinger and master equations for continuous measurement simulations.
Source: arXiv
URL: https://arxiv.org/html/2504.21440v1
Date: 2025-04-30
Excerpt: "QuantumToolbox.jl: An efficient Julia framework for simulating open quantum systems"
Context: A newer Julia package that appears to be a modern alternative to QuantumOptics.jl with additional features for stochastic processes.
Confidence: medium
```

---

### 4.4 QuantumCumulants.jl

```
Claim: QuantumCumulants.jl is a Julia package for the symbolic derivation of mean-field equations for quantum mechanical operators, using generalized cumulant expansions to arbitrary order for open quantum systems.
Source: GitHub
URL: https://github.com/qojulia/QuantumCumulants.jl
Date: 2020-02-19
Excerpt: "QuantumCumulants.jl is a package for the symbolic derivation of mean-field equations for quantum mechanical operators in Julia. The equations are derived using fundamental commutation relations of operators."
Context: Can convert derived equations to ModelingToolkit.jl and solve with DifferentialEquations.jl. Uses Symbolics.jl for symbolic simplification.
Confidence: medium
```

---

## 5. Machine Learning for Quantum Dynamics

### 5.1 MLQD (Machine Learning Quantum Dynamics)

```
Claim: MLQD is an open-source Python package providing machine learning-based approaches for efficient propagation of quantum dissipative dynamics. It implements recursive (KRR-based) and non-recursive (CNN-based AIQD and OSTL) methods with hyperparameter optimization and visualization.
Source: Computer Physics Communications / ChemRxiv
URL: https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2023-0xkv1
Date: 2023
Excerpt: "We have developed MLQD as a comprehensive framework that streamlines and supports the implementation of our recently published machine learning-based approaches for efficient propagation of quantum dissipative dynamics. This framework encompasses: (1) the recursive dynamics with kernel ridge regression (KRR) method, as well as the non-recursive approaches utilizing convolutional neural networks (CNN), namely (2) artificial intelligence-based quantum dynamics (AIQD), and (3) one-shot trajectory learning (OSTL)."
Context: MLQD is available on the XACS cloud computing platform via MLatom interface, and as a pip package. Demonstrated on 2-state spin-boson model and 8-site FMO complex.
Confidence: high
```

**GitHub:** https://github.com/Arif-PhyChem/MLQD
**Documentation:** http://mlatom.com/tutorial/mlqd/
**License:** Apache 2.0

---

### 5.2 QD3SET-1 Database

```
Claim: QD3SET-1 is a database of 8 datasets containing time-evolved population and coherence dynamics for the spin-boson model and Fenna-Matthews-Olson (FMO) light-harvesting complex, generated with HEOM and LTLME methods, designed for benchmarking ML and physics-based approaches.
Source: Frontiers in Physics
URL: https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2023.1223973/full
Date: 2023-07-03
Excerpt: "Here, we present a QD3SET-1 database, a collection of eight datasets of time-evolved population dynamics of the two systems: the spin-boson (SB) model and the Fenna-Matthews-Olson (FMO) light-harvesting complex. The primary objective behind the release of the QD3SET-1 database is to provide researchers a valuable resource for the development, testing, and validation of their approaches."
Context: The SB dataset contains 1000 trajectories per parameter combination (epsilon, lambda, gamma, beta). FMO data uses 8-site Hamiltonian with HEOM in QuTiP. Data is freely available on GitHub.
Confidence: high
```

**GitHub:** https://github.com/Arif-PhyChem/QD3SET
**Citation:** Ullah et al., *Front. Phys.* **11**, 1223973 (2023)

---

### 5.3 Machine Learning Exciton Dynamics (FMO)

```
Claim: Multi-layer perceptrons (neural networks) have been used to predict TDDFT excited state energies of bacteriochlorophylls in the FMO complex, achieving 0.01 eV (0.5%) prediction errors and enabling rapid computation of spectral densities and exciton populations from MD trajectories.
Source: NIH/PMC - Proceedings of the Royal Society A
URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC6020119/
Date: 2018
Excerpt: "Once multi-layer perceptrons are trained, predicting excited state energies was found to be significantly faster than the corresponding QM/MM calculations. We showed that multi-layer perceptrons can successfully reproduce the energies of QM/MM calculations to a high degree of accuracy with prediction errors contained within 0.01 eV (0.5%)."
Context: Uses Coulomb matrices as rotation-invariant molecular representations. Training on correlation-sampled data correctly predicts spectral density shapes. Acceleration enables longer and larger exciton dynamics simulations.
Confidence: high
```

---

## 6. Neural Quantum States

### 6.1 NetKet

```
Claim: NetKet is an open-source Python library built on JAX for studying many-body quantum systems with neural quantum states and machine learning techniques. It supports ground state search, quantum state tomography, and both unitary and Markovian dissipative dynamics.
Source: netket.org / GitHub
URL: https://www.netket.org/
Date: Active development (founded 2018, NetKet 3 released 2021+)
Excerpt: "NetKet is an open-source project delivering cutting-edge methods for the study of many-body quantum systems with artificial neural networks and machine learning techniques. Neural Quantum States. NetKet provides state-of-the-art Neural-Network Quantum states, and advanced learning algorithms to study many-body quantum systems."
Context: Based on JAX (CPU, GPU, TPU). Multi-framework model support (Flax, Equinox, Haiku). Can interoperate with OpenFermion and QuTiP. Affiliated with NumFOCUS.
Confidence: high
```

**Installation:** `pip install netket` (Python 3.11+); GPU: `pip install 'netket[cuda]'`
**GitHub:** https://github.com/netket/netket (1,300+ stars)
**Citation:** Vicentini et al., *Phys. Rev. Research* (NetKet 3 paper)

---

### 6.2 jVMC

```
Claim: jVMC is a versatile Python codebase for variational Monte Carlo with neural quantum states, leveraging JAX for automatic differentiation, just-in-time compilation to accelerators (GPU/TPU), and distributed computing via MPI. It supports arbitrary NQS architectures and model Hamiltonians.
Source: SciPost Physics Codebases / arXiv
URL: https://arxiv.org/abs/2108.03409
Date: 2021-08-07
Excerpt: "Here, we present a Python codebase that supports arbitrary NQS architectures and model Hamiltonians. Additionally leveraging automatic differentiation, just-in-time compilation to accelerators, and distributed computing, it is designed to facilitate the composition of efficient NQS algorithms."
Context: Uses Flax Linen for neural network modules. Supports holomorphic and non-holomorphic NQS, autoregressive NQS for direct sampling. Benchmarked on up to 50-site TFIM with symmetrized RNN. Achieves near-ideal MPI scaling up to 64 GPUs.
Confidence: high
```

**Documentation:** https://jvmc.readthedocs.io/
**GitHub:** https://github.com/markusschmitt/jVMC
**Citation:** Schmitt & Reh, *SciPost Phys. Codebases* **2** (2022)

---

### 6.3 Mapalus (TensorFlow NQS)

```
Claim: Mapalus implements neural-network quantum states with TensorFlow 2 for GPU acceleration, with applications for finding extreme eigenvalues of Hermitian matrices and transfer learning protocols to improve NQS scalability.
Source: GitHub
URL: https://github.com/remmyzen/mapalus
Date: 2021-07-29
Excerpt: "This repository implements the neural-network quantum states with Python 3 and Tensorflow 2 library to speed-up the process with graphics processing units (GPU). It also has some applications for finding extreme eigenvalues of a Hermitian matrix."
Context: Inspired by early NetKet. Used in transfer learning papers to improve NQS scalability and efficiency.
Confidence: medium
```

---

## 7. Quantum Machine Learning Platforms

### 7.1 PennyLane

```
Claim: PennyLane is an open-source Python library for quantum machine learning, providing automatic differentiation, hybrid quantum-classical workflows, and integration with PyTorch, TensorFlow, and JAX. It includes simulators (default.qubit, default.mixed, lightning.qubit GPU) and hardware backends.
Source: arXiv
URL: https://arxiv.org/html/2511.14786v1
Date: 2025-11-13
Excerpt: "PennyLane is an open-source Python library for quantum machine learning and quantum computing, first released in 2018. It provides a unified interface for constructing quantum circuits, computing gradients through automatic differentiation, and optimizing variational quantum algorithms."
Context: While primarily designed for quantum circuits and QML, PennyLane's default.mixed simulator supports mixed quantum states and could be adapted for Lindbladian dynamics. However, it is NOT primarily designed for open quantum system dynamics simulation.
Confidence: high (for QML); low (for open systems dynamics)
```

**Note:** PennyLane is primarily a quantum computing/ML framework, not an open quantum dynamics simulator. Its relevance to biological quantum transport is indirect (through variational algorithms).

---

### 7.2 TensorFlow Quantum

```
Claim: TensorFlow Quantum (TFQ), announced by Google in March 2020, is a library for quantum machine learning that bridges quantum circuits with TensorFlow. As of 2020, it was primarily geared toward executing quantum circuits on classical simulators.
Source: Google / TensorFlow blog
URL: https://blog.tensorflow.org/2020/03/announcing-tensorflow-quantum-open.html
Date: 2020-03-09
Excerpt: "Today, TensorFlow Quantum is primarily geared towards executing quantum circuits on classical quantum circuit simulators. In the future, TFQ will be able to execute quantum circuits on actual quantum processors that are supported by Cirq."
Context: TFQ is NOT designed for open quantum system dynamics. It focuses on hybrid quantum-classical ML models. There is NO documented application of TFQ to simulating Lindblad dynamics, HEOM, or biological quantum transport.
Confidence: high (that it is NOT relevant for open systems dynamics)
```

---

### 7.3 Strawberry Fields

```
Claim: Strawberry Fields is Xanadu's full-stack Python library for continuous-variable (CV) quantum optical circuits, with simulators (Gaussian, Fock, Bosonic, TensorFlow backend) and applications in graph optimization, machine learning, and chemistry. It integrates with OpenFermion for bosonic systems.
Source: Quantum Journal / Xanadu
URL: https://quantum-journal.org/papers/q-2019-03-11-129/
Date: 2019-03-11
Excerpt: "We introduce Strawberry Fields, an open-source quantum programming architecture for light-based quantum computers, and detail its key features. Built in Python, Strawberry Fields is a full-stack library for design, simulation, optimization, and quantum machine learning of continuous-variable circuits."
Context: NOT directly relevant to open quantum system dynamics or biological quantum transport. Its primary domain is photonic quantum computing. The SFOpenBoson plugin allows OpenFermion bosonic Hamiltonians to be simulated, which could include some vibronic models.
Confidence: medium (indirect relevance through bosonic simulations)
```

---

## 8. Spectral Density Fitting Tools

### 8.1 QuTiP 5 Environment API

```
Claim: QuTiP 5 provides a comprehensive Environment API for open quantum systems that includes pre-defined DrudeLorentzEnvironment and UnderDampedEnvironment classes, with methods for spectral density fitting (approx_by_sd_fit), correlation function fitting (approx_by_cf_fit), Matsubara expansion, and Pade expansion.
Source: QuTiP 5.1 Documentation
URL: https://qutip.readthedocs.io/en/v5.1.1/guide/guide-environments.html
Date: 2025-01-15
Excerpt: "All bosonic environments can be approximated by directly fitting their correlation function with a multi-exponential ansatz (approx_by_cf_fit) or by fitting their spectral density with a sum of Lorentzians (approx_by_sd_fit), which correspond to underdamped environments with known multi-exponential decompositions."
Context: The approx_by_sd_fit method fits J(omega) with sum of underdamped terms: J(omega) = sum_k (2*a_k*b_k*omega)/(((omega+c_k)^2+b_k^2)*((omega-c_k)^2+b_k^2)). The number of terms is determined iteratively based on RMSE threshold.
Confidence: high
```

```
Claim: QuTiP-BoFiN demonstrated fitting arbitrary spectral densities with different approaches (Matsubara, Pade, direct fitting) and applied these to the FMO photosynthetic complex, showing how non-Markovian environments protect against pure dephasing.
Source: Physical Review Research
URL: https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.013181
Date: 2023-03-15
Excerpt: "For the bosonic case, our results include demonstrations of how to fit arbitrary spectral densities with different approaches, and a study of the dynamics of energy transfer in the Fenna-Matthews-Olson photosynthetic complex."
Context: The paper shows k_J = 3 and 4 underdamped spectral density fitting for FMO and clarifies how suitable non-Markovian environments can protect against pure dephasing.
Confidence: high
```

---

## 9. 2D Electronic Spectroscopy

### 9.1 NISE (Nonlinear Spectroscopy of Excitons)

```
Claim: NISE is a software package for nonlinear spectral simulations of excitonic systems, with tutorials including the LH2 light-harvesting system. It can construct Hamiltonian trajectories from PDB files, calculate linear absorption, luminescence, circular dichroism, and 2DUVvis spectra.
Source: GitHub
URL: https://github.com/GHlacour/NISE_Tutorials
Date: Active
Excerpt: "The folder LH2 contains the files needed to construct a Hamiltonian trajectory for the LH2 light harvesting system of purple bacteria... Submit the 2DUVvis calculation using the submit2D_MPI script. Find the 2D spectra from the calculated response functions with: 2DFFT input2D."
Context: NISE appears to be a mature C/Fortran code with Python preprocessing and analysis scripts. It uses MD-generated Hamiltonian trajectories and computes nonlinear optical response functions.
Confidence: medium
```

---

### 9.2 GPU-HEOM for 2D Spectra

```
Claim: GPU-HEOM has been used to compute numerically intense 2D echo spectra of the FMO complex, the only tool employed to date for 2D spectra with structured spectral densities.
Source: J. Chem. Theory Comput.
URL: https://pubs.acs.org/doi/10.1021/ct500629s
Date: 2015
Excerpt: "GPU-HEOM is more flexible in the calculation modes and utilizes the high compute-throughput of modern GPUs for the efficient computation of numerically intense 2D echo spectra of the FMO complex."
Context: The combination of HEOM accuracy with GPU acceleration makes 2D spectroscopy simulations feasible for realistic spectral densities with many vibrational peaks.
Confidence: high
```

---

## 10. Photosynthetic Complex Simulators

### 10.1 quantum_HEOM (Python Package)

```
Claim: quantum_HEOM is a Python package for calculating bath-influenced excitonic energy transfer dynamics in photosynthetic pigment-protein complexes, implementing three Lindblad master equation forms and interfacing with QuTiP's HEOM solver. It includes built-in 2-site spin-boson and 7-site FMO complex models.
Source: GitHub / Nature Communications (citing paper)
URL: https://github.com/jwa7/quantum_HEOM
Date: 2019-2022
Excerpt: "quantum_HEOM allows the bath-influenced excitonic energy transfer dynamics of open quantum systems to be calculated and plotted. Users can simulate EET for in-built model systems, namely a 2-site spin-boson system or the 7-site FMO complex."
Context: Developed as an MSci project at the University of Bristol. The citing paper (Ullah & Dral, Nature Communications 2022) used it to generate training data for AI-based quantum dynamics of FMO.
Confidence: medium
```

---

### 10.2 QMaster for Large Light-Harvesting Complexes

```
Claim: QMaster has been applied to simulate excitonic energy transfer in LHC II from spinach (14 chlorophylls per monomeric unit, with structured spectral densities containing 48+ vibrational peaks), using HEOM with OpenCL acceleration across multiple hardware architectures.
Source: J. Chem. Theory Comput.
URL: https://pubs.acs.org/doi/10.1021/ct500629s
Date: 2015
Excerpt: "We deploy the QMaster-tool to accurately model the time scale of excitonic energy transfer in LHC II found in spinach... LHC II consists of three monomeric units with 14 chlorophylls (Chl), further divided into two types Chla/Chlb."
Context: LHC II is significantly larger than FMO (14 vs 7 chlorophylls). The structured spectral density with >48 vibrational peaks requires the memory capacity of CPU-based OpenCL accelerators rather than GPU-only approaches.
Confidence: high
```

---

## 11. Workflow Management Tools

### 11.1 Snakemake for Physics Workflows

```
Claim: Snakemake is a Python-based workflow management system that has been applied to orchestrate computational physics analyses, including quantum Monte Carlo calculations and high energy physics data analysis. It enables reproducible, parallel execution of multi-stage computational pipelines.
Source: arXiv / GitHub
URL: https://github.com/PNNL-CompBio/CME-QM
Date: 2022
Excerpt: "CME pipeline perform quantum mechanical simulations using computational chemistry code called NWChem for geometry optimization, chemical property prediction and computing spectral properties. Pipeline (Snakemake Workflow)."
Context: While not specific to quantum dynamics, Snakemake has been used for: (1) QMC-SW for quantum Monte Carlo chemistry workflows, (2) CME for quantum mechanical simulations with NWChem, (3) HEP data analysis workflows on REANA. It provides a general framework for reproducible quantum simulation pipelines.
Confidence: medium
```

---

## 12. Assessment Summary

### 12.1 Tools Specifically Mentioned in Phase 1 — Verified Status

| Tool | Status | Relevance to Biological OQS |
|------|--------|---------------------------|
| **QuTiP** | Dominant Python package | Direct — HEOM, master equations, optimal control |
| **pyOMNIS** | **NOT FOUND** — No evidence of existence as quantum tool | N/A |
| **Quipucamayoc** | **NOT A QUANTUM TOOL** — Data visualization library for historical quipu data | N/A |
| **Strawberry Fields** | Active (Xanadu) — CV quantum computing | Indirect — bosonic systems via SFOpenBoson |
| **TensorFlow Quantum** | Active but focused on quantum circuits | **NOT applicable** to open systems dynamics |
| **PennyLane** | Active (Xanadu) — QML framework | Indirect — could simulate open systems via default.mixed |

### 12.2 Most Relevant Tools for Biological Quantum Transport (Ranked)

| Rank | Tool | Language | GPU? | Key Feature | Learning Curve |
|------|------|----------|------|-------------|----------------|
| 1 | **QuTiP 5 + BoFiN** | Python | Partial (CuPy backend) | Full HEOM, spectral density fitting, largest community | Moderate |
| 2 | **Dynamiqs** | Python (JAX) | Yes (JAX/GPU) | Fast, differentiable, QuTiP-like API | Moderate |
| 3 | **OQuPy** | Python | No | Process tensor, non-Markovian, exact correlations | Moderate |
| 4 | **MLQD + QD3SET** | Python | Via TensorFlow | Pre-trained ML models for FMO/spin-boson | Low |
| 5 | **GPU-HEOM** | CUDA C | Yes (CUDA) | Reference HEOM for 2D spectra, FMO studies | High |
| 6 | **pyheom** | Python/C++ | Yes (CUDA) | Generalized HEOM for arbitrary correlations | High |
| 7 | **PHI** | C (OpenMP) | No (multicore CPU) | Large systems (50 sites), shared memory | High |
| 8 | **QuantumOptics.jl** | Julia | No | Fast, clean syntax, similar to QuTiP | Moderate |
| 9 | **OpenQuantumSystems.jl** | Julia | No | Quantum biology focus, finite basis | Moderate |
| 10 | **ACE** | C++ | No | Arbitrary non-Gaussian environments | High |
| 11 | **NetKet** | Python (JAX) | Yes (JAX) | Neural quantum states for ground/dynamics | High |
| 12 | **fHEOM** | Python (QuTiP) | Via QuTiP | Correlated bath factorization | Low |

### 12.3 Key Gaps and Limitations

1. **No dedicated dephasing-assisted transport simulator**: No tool specifically targets the simulation of dephasing-assisted/environment-assisted quantum transport (ENAQT). General-purpose HEOM/master equation solvers must be configured manually for this purpose.

2. **pyOMNIS does not appear to exist**: Extensive searches found no evidence of a Python package called "pyOMNIS" for quantum dynamics. This may be a misremembered or conflated name.

3. **Quipucamayoc is unrelated to quantum**: The package is a data visualization library for Incan quipu and historical data digitization — completely unrelated to quantum physics.

4. **TensorFlow Quantum is not for open dynamics**: TFQ focuses on parameterized quantum circuits for quantum machine learning, not Lindbladian/HEOM dynamics.

5. **Workflow management is underutilized**: While Snakemake has been used for QMC and QM workflows, there is no dedicated workflow framework for quantum dynamics simulations. Most researchers write ad-hoc Python scripts.

6. **2D spectroscopy tools are limited**: NISE appears to be the main open-source option for nonlinear spectroscopy of excitons, but it is a C/Fortran code with limited Python integration. GPU-HEOM can compute 2D spectra but is specialized for HEOM-based approaches.

### 12.4 Recommended Stack for Biological Quantum Transport Research

**For beginners:** QuTiP 5 (HEOM + environments) → MLQD (pre-trained models)

**For performance:** Dynamiqs (GPU, differentiable) or GPU-HEOM (reference CUDA HEOM)

**For non-Markovian/non-Gaussian:** OQuPy (process tensor) → ACE (arbitrary environments)

**For machine learning:** NetKet (neural states) + MLQD (dynamics prediction)

**For large complexes:** PHI (multicore CPU, 50 sites) or QMaster (OpenCL, multi-accelerator)

**For spectral density fitting:** QuTiP 5 Environment API (approx_by_sd_fit, approx_by_cf_fit)

**For 2D spectroscopy:** NISE (nonlinear exciton simulations) + GPU-HEOM (HEOM-based 2D echo)

---

## References

1. Fux et al., "OQuPy: A Python package to efficiently simulate non-Markovian open quantum systems with process tensors", *J. Chem. Phys.* **161**, 124108 (2024)
2. Cygorek et al., "Numerically exact open quantum systems simulations for arbitrary environments using automated compression of environments", *Nature Physics* **17**, 172 (2021)
3. Guilmin et al., "Dynamiqs: an open-source Python library for GPU-accelerated and differentiable simulation of quantum systems" (2025), in preparation
4. Lambert et al., "QuTiP-BoFiN: A bosonic and fermionic numerical hierarchical-equations-of-motion library", *Phys. Rev. Research* **5**, 013181 (2023)
5. Kramer et al., "QuantumOptics.jl: A Julia framework for simulating open quantum systems", *Comput. Phys. Commun.* **227**, 109 (2018)
6. Ullah et al., "MLQD: A package for machine learning-based quantum dissipative dynamics", *Comput. Phys. Commun.* **294**, 108940 (2024)
7. Ullah et al., "QD3SET-1: A database with quantum dissipative dynamics datasets", *Front. Phys.* **11**, 1223973 (2023)
8. Vicentini et al., "NetKet 3: Machine learning toolbox for quantum physics" (2021+)
9. Schmitt & Reh, "jVMC: Versatile and performant variational Monte Carlo", *SciPost Phys. Codebases* **2** (2022)
10. Kramer et al., "GPU-HEOM: High-Performance Solution of Hierarchical Equations of Motion", *J. Chem. Theory Comput.* **7**, 7 (2011)
11. Strumpfer & Schulten, "PHI: Parallel Hierarchy Equations of Motion Integrator", *J. Chem. Theory Comput.* **8**, 2808 (2012)
12. Ikeda & Scholes, "Generalization of the hierarchical equations of motion theory", *J. Chem. Phys.* **152**, 204101 (2020)
13. Johansson et al., "QuTiP 5: The Quantum Toolbox in Python", arXiv:2412.04705 (2024)
14. Abbott, "quantum_HEOM" (2022), https://github.com/jwa7/quantum_HEOM
15. Killoran et al., "Strawberry Fields: A Software Platform for Photonic Quantum Computing", *Quantum* **3**, 129 (2019)
16. Google AI Quantum, "Announcing TensorFlow Quantum" (2020)

---

*Compiled from 20+ independent web searches across Google Scholar, arXiv, GitHub, journal websites, and documentation portals. All claims include inline citations with URLs.*
