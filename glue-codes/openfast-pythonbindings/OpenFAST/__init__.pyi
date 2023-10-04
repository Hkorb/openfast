from __future__ import annotations

from enum import Enum, auto
import numpy.typing as npt
import numpy as np

array = npt.NDArray[np.float64]

class GlobTurbineDataType:
    def __init__(self, TurbID: int, FASTInputFileName: str, FASTRestartFileName: str,
                 TurbineBasePos: list[float], TurbineHubPos: list[float],
                 forcePtsBladeDistributionType: str,
                 numForcePtsBlade: int, numForcePtsTwr: int, 
                 nacelle_cd: float=0.0, nacelle_area: float=0.0,
                 air_density: float=0.0): ...

class SimStartType(Enum):
    INIT = auto()
    TRUE_RESTART = auto()
    RESTART_DRIVER_INIT_FAST = auto()

class FastInputs:
    def __init__(self, nTurbinesGlob: int, dryRun: bool, debug: bool, tStart: float,
                 simStart: SimStartType, 
                 nEveryCheckPoint: int, tMax: float, dtFAST: float, scStatus: bool,
                 scLibFile: str, globTurbineData: list[GlobTurbineDataType]): ...
    
class OpenFAST:
    def __init__(): ...
    def set_inputs(self, fast_inputs: FastInputs): ...
    def init(self): ...
    def solution0(self): ...
    def step(self): ...
    def step_no_write(self): ...
    def end(self): ...
    def calc_nacelle_force(self, u, v, w, cd, area, rho, fx, fy, fz): ...
    def set_velocities(self,
                       nacelle_velocities: array,
                       blade_velocities_x: array,
                       blade_velocities_y: array,
                       blade_velocities_z: array,
                       tower_velocities_x: array,
                       tower_velocities_y: array,
                       tower_velocities_z: array): ...

    def get_forces(self,
                   nacelle_forces: array,
                   blade_forces_x: array,
                   blade_forces_y: array,
                   blade_forces_z: array,
                   tower_forces_x: array,
                   tower_forces_y: array,
                   tower_forces_z: array): ...
    
    def get_coordinates(self,
                nacelle_coordinates_x: array,
                nacelle_coordinates_y: array,
                nacelle_coordinates_z: array,
                blade_coordinates_x: array,
                blade_coordinates_y: array,
                blade_coordinates_z: array,
                tower_coordinates_x: array,
                tower_coordinates_y: array,
                tower_coordinates_z: array): ...
    
    def allocate_turbines_to_procs_simple(self): ...
    
    def compute_torque_thrust(self, torque: array, thrust: array): ...
    def get_azimuths(self, azimuths: array): ...
    def get_rotor_speeds(self, rotor_speeds: array): ...