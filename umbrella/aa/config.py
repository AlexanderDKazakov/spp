__author__ = "Alexander D. Kazakov (aleksandr[dot]kazakov_at_uv[dot]es)"
__version__ = "0.1.0"

from pathlib import Path
import openmm.unit as unit
import numpy as np
from enum import Enum

OUTPUT_PATH = Path("../out")
OUTPUT_PATH.mkdir(exist_ok=True)
GROMACS_TOP = Path("~/code/gromacs/installed/share/gromacs/top").expanduser()
ENERGY_DECOMPOSITION_FNAME = "energy_decompose.txt"

class MOTIONS(Enum):
    TRANSLATION = 0
    ROTATION = 1

CONFIG = {

    # GENERAL
    "box_l":           [80. * unit.angstroms for _ in range(3)],
    "pressure":        1.0*unit.bar,
    "temperature":     300.0*unit.kelvin,
    'frictionCoeff':   1/unit.picosecond,
    "step_size":       0.002*unit.picosecond,
    "nonBondedCutoff": 1.0*unit.nanometer,
    "switchDistance":  0.9*unit.nanometer,

    # WINDOWS SETUP
    "increment_steps": 1, # delta t
    "motion_type":     MOTIONS.ROTATION,
    "motion_type_params": {
        MOTIONS.ROTATION : {
            "variable":               "theta0",
            "initial_variable_value": 0*unit.radian,
            "atom_indices":           [(107, 96, 24, 13), ],                             # idx tortion
            "cv_formula":             "0.5 * fc_pull * (cv-theta0)^2",                   # formula
            "fc_pull":                1000.0 * unit.kilojoule_per_mole/unit.radians**2,  # kJ/rad**2 # force constant
            "velocity_pulling":       0.001 * unit.radians / unit.picosecond,            # rad/ps
            "windows":                np.linspace(-0.3, 0.6*np.pi, 32),                     # desired values
        },
        MOTIONS.TRANSLATION : {
            "variable":               "r0",
            "initial_variable_value": 0.3*unit.nanometer,
            "atom_indices":           [(12, 95), (5, 88), (6, 89), ],                      # idx
            "cv_formula":             "0.5 * fc_pull * (cv-r0)^2",                         # formula
            "fc_pull":                1000.0 * unit.kilojoule_per_mole/unit.nanometers**2, # kJ/nm**2 # force constant
            "velocity_pulling":       0.001 * unit.nanometers / unit.picosecond,           # nm/ps
            "windows":                np.linspace(0.2, 0.6, 26),                           # desired values
        }
    },

    # PER WINDOW SETUP
    "total_steps":                 480_000_000, # delta t
    "cv_record_steps":             10_000,      # delta t| CV
    "report_steps":                100_000,     # delta t| State|PDB
    "windows_short_equilib_steps": 10_000_000,  # delta t | 500_000 too short
}
