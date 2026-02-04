# Advanced Orbital Mechanics Simulator
# Earth-centered, J2 perturbation, element drift, ground track, maneuvers, animation
# Author: Samarthya Rindani

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Consts & planets
SecPerHour = 3600.0
EarthSideRealDaySec = 86164.0905
EarthOmg = 2 * np.pi / EarthSideRealDaySec  # rad/s

planets = {
    "Earth": {
        "mu": 3.986004418e5,      # km³/s²
        "radius": 6378.137,       # km
        "J2": 1.08262668e-3,
        "J3": -2.5324105e-6,
        "J4": -1.6198976e-6,
        "omega": EarthOmg
    },
    "Mars": {
        "mu": 4.282837e4,
        "radius": 3396.19,
        "J2": 1.96045e-3,
        "J3": 3.145e-5,
        "J4": -1.537e-5,
        "omega": 2 * np.pi / (24.6229 * 3600)
    },
    "Moon": {
        "mu": 4.9048695e3,
        "radius": 1737.4,
        "J2": 2.032e-4,
        "J3": 0.0,
        "J4": 0.0,
        "omega": 2 * np.pi / (27.321661 * 86400)
    }
}


# Force models
def gravtwobody(r, mu):
    r_norm = np.linalg.norm(r)
    return -mu * r / r_norm**3

def gravj2(r, mu, R, J2):
    x, y, z = r
    r_norm = np.linalg.norm(r)
    r2 = r_norm * r_norm
    z2 = z * z
    ftr = 1.5 * J2 * mu * R**2 / r_norm**5
    ax = ftr * x * (5 * z2 / r2 - 1)
    ay = ftr * y * (5 * z2 / r2 - 1)
    az = ftr * z * (5 * z2 / r2 - 3)
    return np.array([ax, ay, az])

def gravj3(r, mu, R, J3):
    x, y, z = r
    r_norm = np.linalg.norm(r)
    r2 = r_norm**2
    z2 = z**2

    ftr = 0.5 * J3 * mu * R**3 / r_norm**7
    common = 5 * z * (7 * z2 / r2 - 3)

    ax = ftr * x * common
    ay = ftr * y * common
    az = ftr * (6 * z2 - 7 * z2**2 / r2 - 3 * r2)

    return np.array([ax, ay, az])


def gravj4(r, mu, R, J4):
    x, y, z = r
    r_norm = np.linalg.norm(r)
    r2 = r_norm**2
    z2 = z**2

    ftr = (5/8) * J4 * mu * R**4 / r_norm**7
    term = 35 * z2**2 / r2**2 - 30 * z2 / r2 + 3

    ax = ftr * x * term
    ay = ftr * y * term
    az = ftr * z * (35 * z2**2 / r2**2 - 42 * z2 / r2 + 9)

    return np.array([ax, ay, az])


# Eqns of motion
def eom(t, state, planet, use_J2=True, use_J3=False, use_J4=False):
    r = state[:3]
    v = state[3:]

    a = gravtwobody(r, planet["mu"])

    if use_J2 and planet.get("J2", 0.0) != 0.0:
        a += gravj2(r, planet["mu"], planet["radius"], planet["J2"])

    if use_J3 and planet.get("J3", 0.0) != 0.0:
        a += gravj3(r, planet["mu"], planet["radius"], planet["J3"])

    if use_J4 and planet.get("J4", 0.0) != 0.0:
        a += gravj4(r, planet["mu"], planet["radius"], planet["J4"])

    return np.hstack((v, a))


# Collision / surface event
def colsn_evt(t, state, radius):
    return np.linalg.norm(state[:3]) - radius
colsn_evt.terminal = True
colsn_evt.direction = -1
 
# Orbit propgn with maneuvers
def propOrbManv(r0, v0, planet, durhrs, use_J2=True, maneuvers=None):
    if maneuvers is None:
        maneuvers = []

    y0 = np.hstack((r0, v0))
    t_total = durhrs * SecPerHour

    manvtimes = sorted([t_sec for t_sec, _ in maneuvers])
    segEnds = [0.0] + manvtimes + [t_total]
    segEnds = np.unique(segEnds)

    all_t = []
    all_y = []
    current_y = y0.copy()
    t_start = 0.0

    for t_end in segEnds[1:]:
        sol = solve_ivp(
            fun=lambda t, y: eom(t, y, planet, use_J2, use_J3=True, use_J4=True),
            t_span=(t_start, t_end),
            y0=current_y,
            method="DOP853",
            rtol=1e-10,
            atol=1e-12,
            events=lambda t, y: colsn_evt(t, y, planet["radius"])
        )


        all_t.append(sol.t)
        all_y.append(sol.y)

        if sol.status == 1:
            print(f"Collision with the surface at t ≈ {sol.t_events[0][0]/3600:.2f} hours")
            break

        for tm, dv in maneuvers:
            if tm >= sol.t[0] and tm <= sol.t[-1]:
                current_y = sol.sol(tm)
                current_y[3:] += dv  # small tolerance
                print(f"Δv applied at t = {tm/3600:.2f} h: {dv} km/s")
                break
        else:
            current_y = sol.y[:,-1]

        t_start = t_end

    t_all = np.concatenate(all_t)
    y_all = np.hstack(all_y)

    class Result:
        def __init__(self, t, y):
            self.t = t
            self.y = y

    return Result(t_all, y_all)

# Orbital element conversions
def elemToState(a, e, i, raan, argp, nu, mu):
    i, raan, argp, nu = np.radians([i, raan, argp, nu])
    p = a * (1 - e**2)
    r_pf = p / (1 + e * np.cos(nu)) * np.array([np.cos(nu), np.sin(nu), 0])
    v_pf = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0])

    def R3(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def R1(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    Q = R3(raan) @ R1(i) @ R3(argp)
    return Q @ r_pf, Q @ v_pf

def stateToElem(r, v, mu):
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    h_vec = np.cross(r, v)
    h = np.linalg.norm(h_vec)
    if h < 1e-8:
        return 0, 0, 0, 0, 0

    n_vec = np.cross([0,0,1], h_vec)
    n = np.linalg.norm(n_vec)
    e_vec = (np.cross(v, h_vec) / mu) - (r / r_norm)
    e = np.linalg.norm(e_vec)

    energy = 0.5 * v_norm**2 - mu / r_norm
    a = -mu / (2 * energy) if energy < 0 else np.inf

    i = np.degrees(np.arccos(h_vec[2] / h))

    if n > 1e-8:
        raan = np.degrees(np.arccos(n_vec[0] / n))
        if n_vec[1] < 0:
            raan = 360 - raan
    else:
        raan = 0.0

    if n > 1e-8 and e > 1e-8:
        argp = np.degrees(np.arccos(np.dot(n_vec, e_vec) / (n * e)))
        if e_vec[2] < 0:
            argp = 360 - argp
    else:
        argp = 0.0

    return a, e, i, raan, argp

# Element history
def compElemHist(sol, mu):
    r = sol.y[:3].T
    v = sol.y[3:].T
    a_list, e_list, i_list, raan_list, argp_list = [], [], [], [], []

    for ri, vi in zip(r, v):
        ai, ei, ii, oi, wi = stateToElem(ri, vi, mu)
        a_list.append(ai)
        e_list.append(ei)
        i_list.append(ii)
        raan_list.append(oi)
        argp_list.append(wi)

    raan_unw = np.unwrap(np.radians(raan_list)) * 180 / np.pi
    argp_unw = np.unwrap(np.radians(argp_list)) * 180 / np.pi

    return (
        sol.t / 3600.0,
        np.array(a_list),
        np.array(e_list),
        np.array(i_list),
        raan_unw,
        argp_unw
    )

def plElemDft(time, a, e, i, raan, argp):
    fig, axs = plt.subplots(3, 2, figsize=(12, 10))
    axs = axs.flatten()
    axs[0].plot(time, a); axs[0].set_title("Semi-major axis a (km)")
    axs[1].plot(time, e); axs[1].set_title("Eccentricity e")
    axs[2].plot(time, i); axs[2].set_title("Inclination i (deg)")
    axs[3].plot(time, raan); axs[3].set_title("RAAN Ω (deg)")
    axs[4].plot(time, argp); axs[4].set_title("Argument of perigee ω (deg)")
    for ax in axs[:5]:
        ax.grid(True)
    axs[5].axis("off")
    plt.tight_layout()
    plt.show()

# Ground track
def eciToEcef(r_eci, t_sec, omega):
    theta = omega * t_sec
    if r_eci.ndim == 1:
        c, s = np.cos(theta), np.sin(theta)
        x =  c * r_eci[0] + s * r_eci[1]
        y = -s * r_eci[0] + c * r_eci[1]
        z = r_eci[2]
        return np.array([x, y, z])
    else:
        c = np.cos(theta)
        s = np.sin(theta)
        x =  c * r_eci[:,0] + s * r_eci[:,1]
        y = -s * r_eci[:,0] + c * r_eci[:,1]
        z = r_eci[:,2]
        return np.column_stack((x, y, z))

def ecefToLL(r_ecef, R):
    lon = np.degrees(np.arctan2(r_ecef[:,1], r_ecef[:,0]))
    lat = np.degrees(np.arcsin(r_ecef[:,2] / np.linalg.norm(r_ecef, axis=1)))
    return lat, lon

def PtGrdTrk(lat, lon, title="Ground Track"):
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111)
    ax.plot(lon, lat, 'b-', lw=1.2, label="Ground track")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.axhline(0, color='k', lw=0.8, alpha=0.4)
    ax.axvline(0, color='k', lw=0.8, alpha=0.4)
    plt.tight_layout()
    plt.show()

# 3D Orbit plot & animation
def PtOrbStatic(r, radius, title):
    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(r[:,0], r[:,1], r[:,2], 'b-', lw=1.2)
    ax.scatter(r[0,0], r[0,1], r[0,2], color='green', s=80, label='Start')
    
    u = np.linspace(0, 2*np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.15)
    
    ax.set_box_aspect([1,1,1])
    ax.set_title(title)
    ax.legend()
    plt.show()

def animOrb(r, radius, interval=40):
    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1,1,1])

    line, = ax.plot([], [], [], 'b-', lw=1.5)
    point, = ax.plot([], [], [], 'ro', ms=8, label='Satellite')

    u = np.linspace(0, 2*np.pi, 40)
    v = np.linspace(0, np.pi, 40)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.15)

    lim = 1.2 * np.max(np.abs(r))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    def init():
        line.set_data([], [])
        line.set_3d_properties([])
        point.set_data([], [])
        point.set_3d_properties([])
        return line, point

    def update(frame):
        line.set_data(r[:frame,0], r[:frame,1])
        line.set_3d_properties(r[:frame,2])
        point.set_data([r[frame,0]], [r[frame,1]])
        point.set_3d_properties([r[frame,2]])
        return line, point

    ani = FuncAnimation(fig, update, frames=len(r), init_func=init,
                        interval=interval, blit=False)
    ax.legend()
    plt.show()

# Main method
def run_cli():
    print("\n=== Advanced Orbital Mechanics Simulator ===\n")
    
    names = list(planets.keys())
    for i, name in enumerate(names, 1):
        print(f"{i}) {name}")
    print("4) Custom")
    
    c = int(input("Select central body (1-4): "))
    if c <= 3:
        planet = planets[names[c-1]]
        body_name = names[c-1]
    else:
        planet = {
            "mu": float(input("μ (km³/s²): ")),
            "radius": float(input("Radius (km): ")),
            "J2": float(input("J2 (0 if none): ")),
            "omega": float(input("Rotation rate ω (rad/s): "))
        }
        body_name = "Custom"

    use_J2 = input("Enable J2 perturbation? (y/n): ").lower().startswith('y')

    print("\nOrbit input mode:")
    print("1) State vector (r, v)")
    print("2) Classical orbital elements")
    mode = int(input("Choice (1 or 2): "))

    if mode == 1:
        r0 = np.array([float(x) for x in input("Position x y z (km): ").split()])
        v0 = np.array([float(x) for x in input("Velocity vx vy vz (km/s): ").split()])
    else:
        a    = float(input("Semi-major axis a (km): "))
        e    = float(input("Eccentricity e: "))
        i    = float(input("Inclination i (deg): "))
        raan = float(input("RAAN Ω (deg): "))
        argp = float(input("Argument of perigee ω (deg): "))
        nu   = float(input("True anomaly ν (deg): "))
        r0, v0 = elemToState(a, e, i, raan, argp, nu, planet["mu"])

    duration = float(input("Simulation duration (hours): "))

    # Maneuvers
    maneuvers = []
    if input("Add impulsive maneuvers? (y/n): ").lower().startswith('y'):
        print("Format: time(hours) dvx dvy dvz (km/s)   [one per line]")
        print("Press Enter twice to finish.")
        while True:
            line = input("> ").strip()
            if not line:
                break
            try:
                parts = list(map(float, line.split()))
                if len(parts) != 4:
                    print("Need 4 values: t_h dvx dvy dvz")
                    continue
                t_h, dvx, dvy, dvz = parts
                maneuvers.append((t_h * 3600.0, np.array([dvx, dvy, dvz])))
            except:
                print("Invalid input — skipped")

    # Propagation
    sol = propOrbManv(r0, v0, planet, duration, use_J2, maneuvers)

    r_history = sol.y[:3].T

    # Plots
    PtOrbStatic(r_history, planet["radius"], f"{body_name} - Orbit")
    animOrb(r_history, planet["radius"], interval=30)

    time_h, a, e, i_deg, raan_deg, argp_deg = compElemHist(sol, planet["mu"])
    plElemDft(time_h, a, e, i_deg, raan_deg, argp_deg)

    # Ground track (only meaningful for rotating bodies like Earth)
    if body_name == "Earth" or input("Attempt ground track? (y/n): ").lower().startswith('y'):
        print("Generating ground track...")
        t_sec = sol.t
        r_eci = r_history
        r_ecef = eciToEcef(r_eci, t_sec, planet["omega"])
        lat, lon = ecefToLL(r_ecef, planet["radius"])

        # Simple longitude unwrapping to reduce jumps
        dlon = np.diff(lon)
        jumps = np.where(np.abs(dlon) > 180)[0]
        for j in jumps:
            lon[j+1:] -= np.sign(dlon[j]) * 360

        PtGrdTrk(lat, lon, f"{body_name} Ground Track – {duration:.1f} h")
if __name__ == "__main__":
    run_cli()