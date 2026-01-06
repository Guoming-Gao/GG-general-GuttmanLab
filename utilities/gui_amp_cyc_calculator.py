import tkinter as tk
from tkinter import messagebox
import math

AVOGADRO = 6.022e23  # mol^-1

TARGET_CONC_NM = 10.0  # target concentration in nM


def compute_cycles():
    try:
        # Read inputs
        n_cells = float(entry_cells.get())
        cdna_per_cell = float(entry_cdna.get())
        volume_ul = float(entry_volume.get())
        efficiency = float(entry_eff.get())

        # Basic checks
        if n_cells <= 0 or cdna_per_cell <= 0 or volume_ul <= 0 or efficiency <= 1.0:
            raise ValueError("Inputs must be positive and efficiency > 1.0")

        # Total starting molecules
        total_molecules = n_cells * cdna_per_cell

        # Convert molecules to moles
        n0_mol = total_molecules / AVOGADRO

        # Volume in liters
        volume_L = volume_ul * 1e-6

        # Initial concentration (M)
        C0 = n0_mol / volume_L

        # Target concentration in M
        Cf = TARGET_CONC_NM * 1e-9

        # Fold increase needed
        fold = Cf / C0

        if fold <= 1:
            result = "Already ≥ 10 nM; 0 cycles needed (or check inputs)."
        else:
            # cycles = ln(fold) / ln(efficiency)
            cycles = math.log(fold) / math.log(efficiency)
            result = f"Cycles needed to reach 10 nM: {cycles:.2f}"

        lbl_result.config(text=result)

    except ValueError as e:
        messagebox.showerror("Input error", str(e))


# GUI setup
root = tk.Tk()
root.title("PCR Cycle Calculator")

# Labels and entries
tk.Label(root, text="N cells per tube:").grid(row=0, column=0, sticky="e")
entry_cells = tk.Entry(root, width=15)
entry_cells.grid(row=0, column=1)

tk.Label(root, text="cDNA molecules per cell:").grid(row=1, column=0, sticky="e")
entry_cdna = tk.Entry(root, width=15)
entry_cdna.grid(row=1, column=1)

tk.Label(root, text="PCR volume (µL):").grid(row=2, column=0, sticky="e")
entry_volume = tk.Entry(root, width=15)
entry_volume.grid(row=2, column=1)

tk.Label(root, text="PCR efficiency (per cycle):").grid(row=3, column=0, sticky="e")
entry_eff = tk.Entry(root, width=15)
entry_eff.grid(row=3, column=1)

# Default example values (your earlier scenario)
entry_cells.insert(0, "30000")
entry_cdna.insert(0, "10000")
entry_volume.insert(0, "50")
entry_eff.insert(0, "1.7")

btn_compute = tk.Button(root, text="Compute cycles", command=compute_cycles)
btn_compute.grid(row=4, column=0, columnspan=2, pady=5)

lbl_result = tk.Label(root, text="Cycles needed to reach 10 nM:")
lbl_result.grid(row=5, column=0, columnspan=2, pady=5)

root.mainloop()
