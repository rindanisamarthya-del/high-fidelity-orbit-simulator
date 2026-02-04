import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import orbital_core as oc

class OrbitalSimulatorGUI:

    def __init__(self, root):
        self.root = root
        self.root.title("Advanced Orbital Mechanics Simulator")
        self.root.geometry("520x720")

        self.create_widgets()

    def create_widgets(self):
        frm = ttk.Frame(self.root, padding=10)
        frm.pack(fill="both", expand=True)

        # ---- Planet Selection ----
        ttk.Label(frm, text="Central Body", font=("Arial", 11, "bold")).pack(anchor="w")

        self.planet_var = tk.StringVar(value="Earth")
        ttk.Combobox(frm, textvariable=self.planet_var,
                     values=list(oc.planets.keys()),
                     state="readonly").pack(fill="x")

        # ---- Orbit Input Mode ----
        ttk.Label(frm, text="\nOrbit Input Mode", font=("Arial", 11, "bold")).pack(anchor="w")
        self.mode_var = tk.IntVar(value=2)

        ttk.Radiobutton(frm, text="State Vector (r, v)", variable=self.mode_var, value=1).pack(anchor="w")
        ttk.Radiobutton(frm, text="Orbital Elements", variable=self.mode_var, value=2).pack(anchor="w")

        # ---- Orbital Elements ----
        self.entries = {}
        grid = ttk.Frame(frm)
        grid.pack(fill="x", pady=5)

        fields = [
            ("a (km)", "7000"),
            ("e", "0.001"),
            ("i (deg)", "51.6"),
            ("RAAN (deg)", "0"),
            ("ω (deg)", "0"),
            ("ν (deg)", "0"),
        ]

        for i, (label, default) in enumerate(fields):
            ttk.Label(grid, text=label).grid(row=i, column=0, sticky="w")
            ent = ttk.Entry(grid)
            ent.insert(0, default)
            ent.grid(row=i, column=1, sticky="ew")
            self.entries[label] = ent

        grid.columnconfigure(1, weight=1)

        # ---- Simulation Settings ----
        ttk.Label(frm, text="\nPerturbations", font=("Arial", 11, "bold")).pack(anchor="w")
        self.j2 = tk.BooleanVar(value=True)
        self.j3 = tk.BooleanVar(value=True)
        self.j4 = tk.BooleanVar(value=True)

        ttk.Checkbutton(frm, text="J2", variable=self.j2).pack(anchor="w")
        ttk.Checkbutton(frm, text="J3", variable=self.j3).pack(anchor="w")
        ttk.Checkbutton(frm, text="J4", variable=self.j4).pack(anchor="w")

        ttk.Label(frm, text="\nSimulation Duration (hours)").pack(anchor="w")
        self.duration = ttk.Entry(frm)
        self.duration.insert(0, "12")
        self.duration.pack(fill="x")

        # ---- Run Button ----
        ttk.Button(frm, text="Run Simulation", command=self.run_sim).pack(pady=15)

    def run_sim(self):
        try:
            planet = oc.planets[self.planet_var.get()]

            # ---- Initial Conditions ----
            if self.mode_var.get() == 2:
                a = float(self.entries["a (km)"].get())
                e = float(self.entries["e"].get())
                i = float(self.entries["i (deg)"].get())
                raan = float(self.entries["RAAN (deg)"].get())
                argp = float(self.entries["ω (deg)"].get())
                nu = float(self.entries["ν (deg)"].get())

                r0, v0 = oc.elemToState(a, e, i, raan, argp, nu, planet["mu"])
            else:
                messagebox.showerror("Not Implemented", "State vector input coming next.")
                return

            duration = float(self.duration.get())

            # ---- Propagation ----
            sol = oc.propOrbManv(
                r0, v0,
                planet,
                duration,
                use_J2=self.j2.get(),
                maneuvers=[]
            )

            r_hist = sol.y[:3].T

            # ---- Plots ----
            oc.PtOrbStatic(r_hist, planet["radius"],
                           f"{self.planet_var.get()} Orbit")
            oc.animOrb(r_hist, planet["radius"])

            time_h, a, e, inc, raan, argp = oc.compElemHist(sol, planet["mu"])
            oc.plElemDft(time_h, a, e, inc, raan, argp)

            if "omega" in planet:
                r_ecef = oc.eciToEcef(r_hist, sol.t, planet["omega"])
                lat, lon = oc.ecefToLL(r_ecef, planet["radius"])
                oc.PtGrdTrk(lat, lon,
                            f"{self.planet_var.get()} Ground Track")

        except Exception as ex:
            messagebox.showerror("Simulation Error", str(ex))


if __name__ == "__main__":
    root = tk.Tk()
    app = OrbitalSimulatorGUI(root)
    root.mainloop()
