# Shercliff-Hartmann MHD Velocity Profiles and Pressure Drop Characteristics in Liquid Metal Divertor Channels

**Paper IV Open-Science Computational Repository**  
**Author:** Jacob James Rasmussen (ORCID: 0009-0009-7084-6523)  
**Target Journal:** *Fusion Engineering and Design* (Elsevier — Traditional No-Fee Subscription Publishing Model)

---

## 1. Executive Summary & Physics Overview

This repository provides open-source, high-fidelity numerical solvers, analytical derivations, and benchmark datasets for magnetohydrodynamic (MHD) liquid metal flows in high-field fusion reactor divertor channels (specifically Pb-17Li eutectic fluid).

Under extreme transverse magnetic fields ($B_0 = 10\text{ T}$), Lorentz forces ($\mathbf{J} \times \mathbf{B}$) reshape standard parabolic velocity profiles into flattened core structures bounded by thin Hartmann boundary layers.

### Key MHD Governing Equations

1. **Hartmann Number ($Ha$):**
   $$Ha = B_0 a \sqrt{\frac{\sigma}{\mu}}$$
   *Measures the ratio of magnetic Lorentz forces to viscous forces.*

2. **Interaction Parameter / Stuart Number ($N$):**
   $$N = \frac{Ha^2}{Re} = \frac{\sigma B_0^2 a}{\rho u_0}$$
   *Quantifies magnetic force dominance relative to inertial forces.*

3. **Hartmann Boundary Layer Thickness ($\delta_H$):**
   $$\delta_H = \frac{a}{Ha}$$
   *Characterizes boundary layer thinning along channel walls perpendicular to the magnetic field.*

4. **MHD Pressure Drop Gradient ($dP/dz$):**
   $$\frac{dP}{dz} = - \sigma u_0 B_0^2 \left( \frac{c_w}{1 + c_w} \right)$$
   *where $c_w = \frac{\sigma_w t_w}{\sigma a}$ is the wall conductance ratio.*

---

## 2. Benchmark Operational Parameters

| Parameter / Metric | Symbol / Formula | Benchmark Value | Physical Unit |
| :--- | :--- | :--- | :--- |
| **Hartmann Number** | $Ha$ | **500.0** | Dimensionless |
| **Interaction Parameter** | $N$ | **1250.0** | Dimensionless |
| **Half-Channel Width** | $a$ | **0.1** | $\text{m}$ |
| **Magnetic Field Strength** | $B_0$ | **10.0** | $\text{T}$ |
| **Bulk Fluid Velocity** | $u_0$ | **1.0** | $\text{m/s}$ |
| **Hartmann Boundary Layer** | $\delta_H = a / Ha$ | **0.20** | $\text{mm}$ |
| **MHD Pressure Drop Gradient**| $dP/dz$ | **0.45** | $\text{MPa/m}$ |
| **Fluid Density (Pb-17Li)** | $\rho$ | **9300.0** | $\text{kg/m}^3$ |
| **Electrical Conductivity** | $\sigma$ | **$1.0 \times 10^6$** | $\text{S/m}$ |
| **Dynamic Viscosity** | $\mu$ | **$1.0 \times 10^{-3}$** | $\text{Pa}\cdot\text{s}$ |

---

## 3. Directory Layout

```
paper-IV-github-repository/
├── LICENSE                          # MIT Open Source License
├── README.md                        # Primary documentation & benchmark summary
├── CITATION.cff                     # Citation Metadata Format (ORCID grounded)
├── requirements.txt                 # Python dependencies
├── .github/workflows/
│   └── python-test.yml              # CI/CD automated validation workflow
├── docs/
│   ├── Main_Manuscript.pdf          # Full paper prepared for Fusion Engineering & Design
│   ├── Cover_Letter.pdf             # Official cover letter to Prof. Mohamed A. Abdou
│   ├── Highlights.pdf               # Core submission highlights (<= 85 chars/bullet)
│   └── Supplementary_Information.pdf # ESI: Derivations & Buckingham Pi scaling
├── src/
│   ├── __init__.py                  # Package initialization
│   ├── utils.py                     # MHD formula scaling functions
│   └── solver.py                    # Velocity profile numerical solver & plot generator
├── data/
│   ├── raw_simulation_data.csv      # 500-point sweep data across Ha = 10 to 1000
│   └── benchmark_summary.json       # Structured benchmark JSON metrics
└── figures/
    ├── figure1_300dpi.png           # High-resolution (300 DPI) velocity profile plot
    └── figure1_vector.svg           # Vector SVG plot asset
```

---

## 4. Quickstart & Usage

To replicate the simulation results, install dependencies and execute the solver:

```bash
# Install requirements
pip install -r requirements.txt

# Run MHD solver and generate data/figures
python src/solver.py .
```

---

## 5. Citation

If you utilize this repository or numerical dataset in your research, please cite:

```bibtex
@article{rasmussen2026shercliff,
  author    = {Rasmussen, Jacob James},
  title     = {Shercliff-Hartmann Magnetohydrodynamic Velocity Profiles and Pressure Drop Characteristics in Liquid Metal Divertor Channels},
  journal   = {Fusion Engineering and Design},
  year      = {2026},
  note      = {ORCID: 0009-0009-7084-6523}
}
```
