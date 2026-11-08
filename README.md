# A Mathematical Study of COVID-19 Using a SEIQR Epidemic Model

B.Sc. (Hons.) Mathematics final-semester project (January–March 2022)
Asansol Girls' College · Kazi Nazrul University
Author: Shrestha Sarkar

This repository contains the project report, the Excel workbooks, the code and the charts for a deterministic **SEIQR** model of COVID-19 (Susceptible, Exposed, Infected, Isolated, Recovered).

## What the project does

- Builds a five-compartment model with births, natural deaths and an isolation (quarantine) class.
- Derives the basic reproduction number **R₀** with the next-generation matrix and shows that it is a threshold:
  - R₀ < 1: the disease-free equilibrium is locally asymptotically stable and the infection dies out.
  - R₀ > 1: it is unstable, and a unique, locally stable endemic equilibrium appears.
- Simulates two cases: **R₀ = 2.35** (epidemic) and **R₀ = 0.59** (the infection disappears).
- Studies how isolation rates change R₀ and finds which parameters matter most (normalised sensitivity indices).
- Compares the parameter values with published estimates and runs the model with literature-based values.
- Models a lockdown as a step change in the transmission rate.
- Extends the model with a disease death rate and vaccination, and gives the vaccination threshold R_v = R₀(1 − εp).
- Fits the model to real daily COVID-19 data for India (4 March 2020 to 10 January 2021).
- Compares Euler's method with Runge–Kutta 4 (RK4).

## Model

```
dS/dt = A − μS − βS(E + I)
dE/dt = βS(E + I) − (α + μ + γ)E
dI/dt = αE − (σ + μ)I
dQ/dt = γE + σI − (θ + μ)Q
dR/dt = θQ − μR
```

| Symbol | Meaning |
|---|---|
| A, μ | recruitment (birth) rate, natural death rate |
| β | transmission rate |
| α | exposed → infected rate |
| γ, σ | isolation rates of exposed and infected people |
| θ | isolated → recovered rate |

R₀ = (βA/μ) · [1/k₁ + α/(k₁k₂)], with k₁ = α + μ + γ and k₂ = σ + μ.

## Main results

- For the illustrative parameters, R₀ = 2.35 gives an epidemic wave with a peak of about **698 infected people on day 63.5**, which then settles at an endemic level of 338. With R₀ = 0.59 the infection dies out.
- Isolating infected people alone may not be enough to bring R₀ below 1, because exposed people also transmit.
- A 60% lockdown from day 30 cuts the infected peak from 698 to 146. If it is lifted on day 120, a second wave of about 636 follows.
- With vaccination the infection dies out when the vaccinated share exceeds (1 − 1/R₀)/ε (61.6% for ε = 0.9 and a disease death rate of 0.01).
- Fitted to the real India data, the effective reproduction number falls from **2.66** before the lockdown to **1.43** during it, **1.15** during reopening and **0.89** after the September 2020 peak. The median error is about 14%.
- With the step 0.1 day used in Excel, Euler's largest error in I(t) is about 5 people (0.7% of the peak); the RK4 error is about 10⁻⁷.

## Selected charts

All 26 charts are in the `charts/` folder.

### The SEIQR model


<img width="1600" height="640" alt="Fig01_Flow_Diagram" src="https://github.com/user-attachments/assets/af26bdc2-4281-47c0-bffa-afe608a3cd79" />


*Flow diagram: people move from susceptible (S) to exposed (E), then to infected (I) or isolated (Q), and finally to recovered (R).*

### Case 1: R₀ = 2.35 (epidemic)


<img width="1440" height="700" alt="Fig08_Case1_E_I_Q_together" src="https://github.com/user-attachments/assets/2092f4ba-14d9-435b-a50d-b615be2b5ac0" />


*Exposed, infected and isolated people rise in one wave (E first, then I, then Q) and settle at the endemic equilibrium.*

### Case 2: R₀ = 0.59 (the infection dies out)

<img width="1440" height="700" alt="Fig14_Case2_E_I_Q_together" src="https://github.com/user-attachments/assets/1abcb028-d4b8-4782-940e-cafcfe89c51d" />

*All three infection classes fall to zero within about 70 days.*

### Phase plane, R₀ > 1

<img width="1440" height="760" alt="Fig15_PhasePlane_Case1" src="https://github.com/user-attachments/assets/4e30d7b8-e65e-4986-843e-bcfdaff8dbc2" />

*Solutions from different starting points spiral into the endemic point.*

### Bifurcation at R₀ = 1

<img width="1200" height="720" alt="Fig02_Bifurcation" src="https://github.com/user-attachments/assets/80e02f2a-0350-456b-a9cb-ed0ae1d34793" />

*The disease-free state loses stability at R₀ = 1 and the endemic equilibrium appears (forward bifurcation).*

### Which parameters matter most for R₀?


<img width="1440" height="720" alt="Fig18_Sensitivity_Indices" src="https://github.com/user-attachments/assets/5483b9ef-dffd-451a-a3d2-df2c920808d5" />


*Normalised sensitivity indices of R₀: β and A act in proportion, σ is more effective than γ, and θ has no effect.*

### Lockdown scenarios

<img width="1440" height="760" alt="Fig19_Lockdown_Scenarios" src="https://github.com/user-attachments/assets/18908988-d0ba-4311-94bb-7175be609d63" />

*An early, strong lockdown suppresses the epidemic; lifting it too soon causes a second wave; a late lockdown barely lowers the peak.*

### Vaccination threshold

<img width="1440" height="720" alt="Fig23_Rv_vs_Coverage" src="https://github.com/user-attachments/assets/f2a4cd99-8df4-4bfa-a76f-e36e8cdb8125" />

*R_v = R₀(1 − εp) against vaccination coverage for two vaccine efficacies.*

### Fit to real data for India

<img width="1600" height="840" alt="Fig26_India_Data_vs_Model" src="https://github.com/user-attachments/assets/fb224477-a6de-4de0-b6e1-44763a198d44" />

*The SEIQR model with four transmission phases (black dots = real 7-day average, red line = model). Median error about 14%.*

### Euler's method against RK4


<img width="1440" height="760" alt="Fig24_Euler_vs_RK4_Error" src="https://github.com/user-attachments/assets/51f16b13-deeb-409c-9f07-5d9324e3f652" />


*Euler's error shrinks in proportion to the step size (first order); RK4's error shrinks about 16 times when the step is halved (fourth order).*

## Repository structure

```
.
├── README.md
├── report/
│   └── COVID19_SEIQR_Project_Corrected.docx
├── excel/
│   ├── SEIQR_COVID19_Simulation.xlsx
│   └── SEIQR_COVID19_India_Fit.xlsx
├── code/
│   ├── python/
│   │   └── seiqr_charts.py
│   └── cpp/
│       ├── seiqr.cpp          (Euler's method)
│       └── seiqr_rk4.cpp      (Runge–Kutta 4)
└── charts/
    └── Fig01 … Fig26 (PNG)
```

### Excel workbooks

**SEIQR_COVID19_Simulation.xlsx** (the model, Case 1 and Case 2)

| Sheet | Contents |
|---|---|
| Read_Me | purpose, method, colour legend, assumptions |
| Parameters | inputs, published estimates with sources, R₀, equilibria, lockdown inputs |
| Sim_R0_gt1 / Sim_R0_lt1 | calculation only (Euler method, 3,650 steps) |
| Sim_R0_gt1_Graphs / Sim_R0_lt1_Graphs | dashboards with six charts each |
| R0_Sensitivity | R₀ against isolation rates and β; sensitivity indices |
| Sim_Lockdown / Sim_Lockdown_Graphs | lockdown scenario and its dashboard |

Change only the blue cells on the Parameters sheet; everything else recalculates.

**SEIQR_COVID19_India_Fit.xlsx** (fit to real data): Read_Me, Parameters, Data, Model, Fit_Calc and Fit_Graphs. The Read_Me sheet explains how to re-fit with Excel Solver.

## How to run the code

**Python** (needs `numpy`, `scipy`, `matplotlib`; also runs on Google Colab):

```bash
pip install numpy scipy matplotlib
python code/python/seiqr_charts.py
```

For Case 2 change `beta = 4e-5` to `beta = 1e-5` in the script.

**C++** (writes a CSV file that Excel can plot):

```bash
g++ -O2 -o seiqr code/cpp/seiqr.cpp && ./seiqr
g++ -O2 -o seiqr_rk4 code/cpp/seiqr_rk4.cpp && ./seiqr_rk4
```

Both programs print R₀ and the peak of I. For Case 1 the RK4 program gives a peak of 697.875 on day 63.5.

## Data source

The India fit uses the series *Daily new confirmed COVID-19 cases per million people (7-day rolling average)*, 4 March 2020 to 10 January 2021. It is the Our World in Data series, as reproduced in the Immigration and Refugee Board of Canada, *COVID-19 bulletin No. 802 – India* (15 January 2021). The values are typed into the Data sheet of the workbook as published.

## Limitations

- The parameter values in the main simulations are **illustrative**, not estimated from data. The workbook and report compare them with published estimates.
- The stability results are local, not global.
- The India fit has five fitted numbers and hand-chosen change dates. The third change date was set at the observed peak, so the decline after the peak is built into the fit. Reported cases are only a fraction of all infections. A good fit shows that the model can reproduce the curve; it does not prove why cases rose and fell.
- The lockdown is a sudden change in β, whereas real lockdowns change gradually.
- The extended model has no waning of vaccine protection and no age groups.

## References

- Kermack, W. O. and McKendrick, A. G. (1927). *Proceedings of the Royal Society A*, 115, 700–721.
- van den Driessche, P. and Watmough, J. (2002). *Mathematical Biosciences*, 180, 29–48.
- Diekmann, O., Heesterbeek, J. A. P. and Roberts, M. G. (2010). *Journal of the Royal Society Interface*, 7, 873–885.
- Lauer, S. A. et al. (2020). *Annals of Internal Medicine*, 172(9), 577–582.
- Li, Q. et al. (2020). *New England Journal of Medicine*, 382(13), 1199–1207.
- He, X. et al. (2020). *Nature Medicine*, 26, 672–675.

The full reference list is in the report.
