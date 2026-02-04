# High-Fidelity Earth Orbit Simulator

A Python-based geo-centered orbital mechanics simulator that numerically propagates satellite motion in three dimensions while accounting for Earth's non-spherical gravity field using zonal harmonic perturbations (J2, J3, J4).

This project is designed to bridge the gap between introductory two-body orbit models and high-fidelity astrodynamics software, with a strong emphasis on physical correctness, visualization, and extensibility.

Note: J3 and J4 are implemented using zonal harmonic acceleration models assuming an oblate central body aligned with the inertial Z-axis.

---

## Features

- Complete 3D inertial orbit propagation
- Earth gravity modeled using zonal harmonics:
  - J2 (oblateness)
  - J3 (pear-shape asymmetry)
  - J4 (higher-order flattening)
- Numerical integration of equations of motion
- Visualization of geo-centered orbits
- Observable orbital precession effects
- Modular force-model architecture (easily extensible)

---

## Physics Model

The spacecraft motion is governed by:

- Newton’s second law in an Earth-centered inertial (ECI) frame
- Central gravitational acceleration
- Perturbative accelerations derived from Earth’s gravitational potential expanded in zonal harmonics

The gravity potential is modeled as:

\[
U(r, \phi) = \frac{\mu}{r} \left[ 1 - \sum_{n=2}^{N} J_n \left( \frac{R_E}{r} \right)^n P_n(\sin \phi) \right]
\]

Where:
- \( \mu \) is Earth’s gravitational parameter
- \( R_E \) is Earth’s mean radius
- \( J_n \) are zonal harmonic coefficients
- \( P_n \) are Legendre polynomials

This allows realistic modeling of:
- Nodal regression
- Argument of perigee precession
- Long-term orbital evolution

---

## Visualization

- 3D Earth rendering
- Orbit trajectory plotting
- Start point marking
- Visual confirmation of perturbation-induced precession

---

## Example Test Case

Typical low Earth orbit test demonstrating J2-induced precession:

- Initial position: ~7000 km from Earth center
- Initial velocity: circular LEO velocity
- Inclined orbit (non-equatorial)
- Simulation time: several hours to days

Expected behavior:
- Orbit plane slowly rotates
- Perigee shifts over time
- Orbit remains bounded and physically stable

---

## Code Structure

- Force model functions (central gravity, J2, J3, J4)
- Numerical propagator
- Event handling (Earth impact detection)
- Visualization module

---

## Future Work

Planned and potential extensions include:

- Mission planning interface
- Orbital element input/output
- Impulsive and continuous maneuvers
- Atmospheric drag (NRLMSISE-00)
- Third-body perturbations (Sun/Moon)
- Ground track visualization
- Maneuver optimization and targeting

---

## References

- Vallado, D. A. — *Fundamentals of Astrodynamics and Applications*
- Montenbruck & Gill — *Satellite Orbits*
- Curtis, H. — *Orbital Mechanics for Engineering Students*

These references are used strictly for theoretical validation of classical equations.

---

## License

MIT License. Free to use, modify, and extend.

---

## Disclaimer

This software is intended for educational and exploratory purposes.
It is **not** a replacement for certified mission-analysis tools such as GMAT, STK, or Orekit.

