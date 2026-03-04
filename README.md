# High-Fidelity Earth Orbit Simulator

A Python-based orbital mechanics simulator that propagates satellite trajectories in three dimensions under Earth's non-spherical gravity field, modelled via zonal harmonic perturbations up to fourth order (J2, J3, J4). Supports Earth, Mars, and Moon as central bodies, orbital element input, impulsive maneuvers, ECI→ECEF frame transformation, ground track generation, and a Tkinter GUI.

> Built by a 1st-year B.Tech Mechanical Engineering student as a self-driven deep-dive into high-fidelity astrodynamics — bridging the gap between two-body Keplerian models and professional tools like GMAT and STK.

---

## Results

### 3D Orbit — ISS-Like LEO (a = 7000 km, i = 51.6°)

Static view showing one complete orbit with Earth rendered in ECI frame:

![3D Orbit Static](docs/fig2_3d_orbit.png)

---

### Live Animation Frame — Inclined LEO Propagation

Frame from the real-time orbit animation showing the satellite position mid-propagation:

![Animation Frame](docs/fig3_anim_frame.png)

---

### Orbital Element Drift — J2 + J3 + J4 Perturbations (6 hours, a = 7000 km, i = 63.4°)

The most physically meaningful output: secular RAAN regression and short-period oscillations in all five elements driven by zonal harmonic perturbations.

![Orbital Elements](docs/fig4_orbital_elements.png)

Key observations from this run:
- **RAAN (Ω):** Secular regression of −0.133 deg/hr, consistent with J2 analytical prediction (see validation below)
- **Semi-major axis & eccentricity:** Short-period oscillations with no secular drift — confirming energy conservation under conservative perturbations
- **Inclination:** Bounded oscillation with ~0.03° amplitude — correctly driven by J3 asymmetry
- **Argument of perigee:** Oscillates rapidly due to near-circular orbit (ω is poorly defined at low eccentricity — expected numerical behaviour)

---

### Ground Track — 6-Hour Propagation (Earth, i = 63.4°)

ECI trajectory converted to ECEF frame accounting for Earth's sidereal rotation (ω⊕ = 7.2921 × 10⁻⁵ rad/s):

![Ground Track](docs/fig5_ground_track.png)

---

## Validation

The J2-induced RAAN precession rate from the simulation is compared against the closed-form first-order secular formula:

$$\dot{\Omega}_{J2} = -\frac{3}{2} \frac{n J_2 R_E^2}{a^2 (1-e^2)^2} \cos i$$

| Parameter | Value |
|---|---|
| Semi-major axis *a* | 7000 km |
| Eccentricity *e* | 0.001 |
| Inclination *i* | 63.4° |
| Orbital period | 97.14 min |
| **Analytical RAAN drift** | **−0.1342 deg/hr** |
| **Simulated RAAN drift** | **−0.1333 deg/hr** |
| **Error** | **< 0.7%** |

Sub-1% agreement with the analytical formula validates the J2 force model implementation. The small residual is expected — it arises from J3/J4 cross-terms that are absent in the first-order secular formula.

---

## Physics Model

Spacecraft motion is governed by Newton's second law in an Earth-Centered Inertial (ECI) frame:

$$\ddot{\mathbf{r}} = -\frac{\mu}{r^3}\mathbf{r} + \mathbf{a}_{J2} + \mathbf{a}_{J3} + \mathbf{a}_{J4}$$

The gravitational potential expanded in zonal harmonics:

$$U(r, \phi) = \frac{\mu}{r} \left[ 1 - \sum_{n=2}^{4} J_n \left( \frac{R_E}{r} \right)^n P_n(\sin \phi) \right]$$

Where $P_n$ are Legendre polynomials, $J_n$ are zonal harmonic coefficients, and $\phi$ is the geocentric latitude.

### Force Model Coefficients (Earth)

| Constant | Value | Source |
|---|---|---|
| μ | 3.986004418 × 10⁵ km³/s² | EGM96 |
| Rₑ | 6378.137 km | WGS-84 |
| J2 | 1.08262668 × 10⁻³ | EGM96 |
| J3 | −2.5324105 × 10⁻⁶ | EGM96 |
| J4 | −1.6198976 × 10⁻⁶ | EGM96 |
| ω⊕ | 7.2921150 × 10⁻⁵ rad/s | IAU |

### Physical effects modelled

| Effect | Modelled By |
|---|---|
| Nodal regression (orbit plane rotation) | J2 secular term |
| Argument of perigee precession | J2 secular term |
| Pear-shape asymmetry (N/S asymmetry) | J3 |
| Higher-order flattening correction | J4 |
| Earth surface collision detection | Terminal event in integrator |
| Frame rotation (ECI → ECEF) | Sidereal rotation matrix |

---

## Numerical Integration

The equations of motion are integrated using **SciPy's DOP853** solver — an 8th-order Dormand-Prince explicit Runge-Kutta method with adaptive step-size control.

| Setting | Value |
|---|---|
| Method | DOP853 (explicit RK8) |
| Relative tolerance | 1 × 10⁻¹⁰ |
| Absolute tolerance | 1 × 10⁻¹² |
| Event detection | Surface collision (terminal) |

DOP853 is the same integration family used in professional astrodynamics tools for its favourable balance of accuracy and computational cost at moderate force-model complexity.

---

## Features

- **Multi-body support:** Earth, Mars, Moon (each with correct μ, Rₑ, Jₙ, ω)
- **Perturbation control:** J2, J3, J4 individually toggleable
- **Dual input modes:** Classical orbital elements or Cartesian state vector
- **Impulsive maneuvers:** Δv events injected at arbitrary mission times
- **Orbit propagation:** Full 3D ECI trajectory via DOP853
- **Orbital element history:** Real-time extraction and 5-panel drift plot
- **Ground track:** ECI → ECEF conversion with sidereal rotation, lat/lon plot
- **Visualization:** Static 3D orbit plot, real-time animation, element drift panels
- **GUI:** Tkinter interface for parameter entry without editing code
- **Modular architecture:** Force models as separate functions — straightforward to extend

---

## Code Structure

```
orbital_core.py      # Physics engine
│
├── planets{}             # Central body constants (Earth / Mars / Moon)
├── gravtwobody()         # Central gravity acceleration
├── gravj2/j3/j4()        # Zonal harmonic perturbation accelerations
├── eom()                 # Full equations of motion (all forces)
├── colsn_evt()           # Surface collision terminal event
├── propOrbManv()         # Orbit propagator with maneuver injection
├── elemToState()         # Classical elements → Cartesian state
├── stateToElem()         # Cartesian state → Classical elements
├── compElemHist()        # Element history extraction + RAAN unwrapping
├── plElemDft()           # 5-panel orbital element drift plot
├── eciToEcef()           # ECI → ECEF frame rotation
├── ecefToLL()            # ECEF → geodetic latitude / longitude
├── PtGrdTrk()            # Ground track visualization
├── PtOrbStatic()         # Static 3D orbit plot
├── animOrb()             # Real-time 3D orbit animation
└── run_cli()             # Interactive command-line interface

gui_app.py           # Tkinter GUI wrapping orbital_core
```

---

## Quickstart

```bash
git clone https://github.com/rindanisamarthya-del/high-fidelity-orbit-simulator
cd high-fidelity-orbit-simulator
pip install numpy scipy matplotlib
```

**CLI mode:**
```bash
python orbital_core.py
```

**GUI mode:**
```bash
python gui_app.py
```

**Minimal API usage:**
```python
import orbital_core as oc
import numpy as np

planet = oc.planets["Earth"]

# ISS-like orbit: a=7000km, e=0.001, i=51.6°
r0, v0 = oc.elemToState(7000, 0.001, 51.6, 0, 0, 0, planet["mu"])

# Propagate 24 hours with full J2/J3/J4
sol = oc.propOrbManv(r0, v0, planet, durhrs=24, use_J2=True)

# Visualise
oc.PtOrbStatic(sol.y[:3].T, planet["radius"], "ISS Orbit")

# Orbital element drift
time_h, a, e, inc, raan, argp = oc.compElemHist(sol, planet["mu"])
oc.plElemDft(time_h, a, e, inc, raan, argp)
```

---

## Planned Extensions

- [ ] Atmospheric drag (Jacchia / NRLMSISE-00 density model)
- [ ] Third-body perturbations (Sun, Moon)
- [ ] Hohmann transfer with bi-impulsive Δv computation
- [ ] Continuous thrust / low-thrust trajectory
- [ ] Osculating → mean element conversion
- [ ] Ground station contact window computation
- [ ] Export to CCSDS OEM format

---

## References

- Vallado, D. A. — *Fundamentals of Astrodynamics and Applications*, 4th ed.
- Montenbruck, O. & Gill, E. — *Satellite Orbits: Models, Methods and Applications*
- Curtis, H. D. — *Orbital Mechanics for Engineering Students*, 3rd ed.
- EGM96 Geopotential Model — NASA GSFC

---

## License

MIT License — free to use, modify, and extend with attribution.

## Disclaimer

This software is intended for educational and research purposes. It is not a replacement for certified mission-analysis tools such as GMAT, STK, or Orekit.
