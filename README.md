# Structural Deformation and Modal Analysis of Truss Bridges

This numerical modeling project was developed as part of the GM (Mathematical & Mechanical Engineering) curriculum at **INSA Rouen Normandie** (Academic Year 2024-2025). Supervised by **Mr. A. Tonnoir**, the work was carried out by **Mouad Sheradj Drissi** and **Amr Zaki Salih**.

---

## Overview

The goal of this project is to simulate and analyze the static deformation and dynamic vibration behavior of elastic truss bridge structures under mechanical loads (uniform pressure, shear).

The physical system is modeled as a 3D network of nodes (mass points) and arcs (elastic bars acting as springs). Newton's second law leads to:
1. **Static Equilibrium Formulation:** Solving the high-dimensional linear system $\mathbb{K}\underline{d} = \underline{F}$, where $\mathbb{K}$ is the symmetric positive semi-definite global stiffness matrix, $\underline{d}$ is the displacement vector, and $\underline{F}$ is the external force vector.
2. **Dynamic Modal Formulation:** Solving the generalized eigenvalue problem $\mathbb{K}\underline{d} = \omega^2 \mathbb{M}\underline{d}$ to extract structural natural frequencies $\omega$ and vibrational mode shapes $\underline{d}$, where $\mathbb{M}$ is the diagonal nodal mass matrix.

---

## Key Numerical & Engineering Findings

- **Iterative Solvers Comparison (Jacobi, Gauss-Seidel, Gradient, Conjugate Gradient):**
  - **Conjugate Gradient (CG):** Demonstrates superior convergence speed and accuracy on well-conditioned problems without requiring parameter tuning[cite: 17]. However, on ill-conditioned topologies (e.g., the detailed Sydney Harbour Bridge model with $\text{cond}(\mathbb{K}) \approx 6.48 \times 10^{17}$), vanilla CG stalls, underscoring the necessity of preconditioning.
  - **Gauss-Seidel:** Displays robust, steady convergence across various geometries and relaxation factors.
  - **Jacobi & Gradient Descent:** Highly sensitive to the relaxation parameter $\theta$; Jacobi diverges when $\theta$ exceeds structural spectral limits ($\theta > 0.735$ on the inverted truss).
- **Modal Analysis & Resonance:**
  - Natural modes are extracted using an **Inverse Power Method with QR Factorization and Deflation**.
  - Lower eigenvalues $\lambda = \omega^2$ correspond to low-frequency, large-amplitude global swaying/stretching modes that pose the highest risk of structural failure and resonance under environmental loads.

---

## Repository Architecture

- **`main.py`**: Interactive CLI driver for bridge visualization, static load simulation, parameter relaxation benchmarking, and modal animations.
- **`algo.py`**: Numerical routines implementing splitting methods (Jacobi, Gauss-Seidel, Gradient), Conjugate Gradient (`GC`), QR factorization, back-substitution, and the deflation-based inverse power iteration (`puInv`).
- **`data.py`**: Assembly routines for the $3n \times 3n$ stiffness matrix $\mathbb{K}$, mass matrix $\mathbb{M}$, boundary condition projection (supports/fixations), and load vector generation.
- **`affichage.py`**: 3D Matplotlib visualizers for truss wireframes, static deformation vector fields, residual convergence plots, and dynamic harmonic mode oscillations.
- **`Pont*.data`**: Geometric topologies, member connectivity, nodal masses, and fixed support constraints for:
  - Benchmark bridge (`choixData = 1`)
  - Sydney Harbour Bridge model (`choixData = 2`)
  - Inverted truss bridge (`choixData = 3`)
- **`TP2_compte_rendu_Mouad_SD_Amr_ZS.pdf`**: Complete theoretical and experimental project report.

---

## Requirements & Execution

### Prerequisites
Make sure Python 3 is installed with standard scientific packages:
  
    pip install numpy matplotlib scipy

---

## Running the Program

Launch the interactive terminal application:

    python3 main.py

---
### The menu prompts you to choose an action:

  **Option 1**: Visualize the 3D bridge geometry and compute matrix properties (symmetry, condition number, positive semi-definiteness).   

  **Option 2**: Apply static load cases (Uniform Pressure / Shear) and solve deformations using Jacobi, Gauss-Seidel, Gradient, or Conjugate Gradient with real-time 3D deformation plotting.

  **Option 3**: Compute N natural vibration modes using inverse iteration with Rayleigh quotient refinement, followed by dynamic 3D modal oscillation animations.   

  **Option 4**: Benchmark residual histories across relaxation values θ∈[0,2] to evaluate solver stability domains.

---
## Authors & Acknowledgments

  **Mouad Sheradj Drissi**

  **Amr Zaki Salih**

  Academic Advisor: **Mr. Antoine Tonnoir** (INSA Rouen Normandie)

---


