# Transfer solver — Lambert, porkchop sweep, Hohmann comparison

import math
import numpy as np
from .constants import GM_SUN, AU, SMA_EARTH, SMA_MARS

# Stumpff functions
def _stumpff_c(z):
    if abs(z) < 1e-6: return 0.5 - z/24.0 + z*z/720.0
    if z > 0: return (1.0 - math.cos(math.sqrt(z))) / z
    return (math.cosh(math.sqrt(-z)) - 1.0) / (-z)

def _stumpff_s(z):
    if abs(z) < 1e-6: return 1.0/6.0 - z/120.0 + z*z/5040.0
    if z > 0:
        sz = math.sqrt(z)
        return (sz - math.sin(sz)) / (sz**3)
    sz = math.sqrt(-z)
    return (math.sinh(sz) - sz) / (sz**3)


def lambert_solve(r1, r2, tof, mu, prograde=True):
    """Solve Lambert's problem (universal variable formulation).
    
    Args:
        r1: departure position vector [x, y, z] (m)
        r2: arrival position vector [x, y, z] (m)
        tof: time of flight (seconds)
        mu:  gravitational parameter (m³/s²)
        prograde: True for prograde (short way)
    
    Returns (v1, v2) velocity vectors, or (None, None) on failure.
    """
    r1, r2 = np.array(r1), np.array(r2)
    R1 = np.linalg.norm(r1)
    R2 = np.linalg.norm(r2)
    
    cross = np.cross(r1, r2)
    cos_dnu = np.dot(r1, r2) / (R1 * R2)
    cos_dnu = max(-1.0, min(1.0, cos_dnu))
    dnu = math.acos(cos_dnu)
    
    if prograde:
        if cross[2] < 0: dnu = 2*math.pi - dnu
    else:
        if cross[2] >= 0: dnu = 2*math.pi - dnu
    
    sin_dnu = math.sin(dnu)
    if abs(sin_dnu) < 1e-15:
        return None, None  # 0° or 180° degenerate
    
    A = sin_dnu * math.sqrt(R1 * R2 / (1.0 - cos_dnu))
    
    # Newton-Raphson
    z = 0.0
    for _ in range(100):
        C = _stumpff_c(z)
        S = _stumpff_s(z)
        y = R1 + R2 + A * (z * S - 1.0) / math.sqrt(C)
        if y < 0: z += 0.1; continue
        
        x = math.sqrt(y / C)
        tof_guess = (x**3 * S + A * math.sqrt(y)) / math.sqrt(mu)
        
        # Numerical derivative
        dz = 0.01
        C2 = _stumpff_c(z + dz)
        S2 = _stumpff_s(z + dz)
        y2 = R1 + R2 + A * ((z+dz) * S2 - 1.0) / math.sqrt(C2)
        if y2 < 0: y2 = abs(y2)
        x2 = math.sqrt(y2 / C2)
        tof2 = (x2**3 * S2 + A * math.sqrt(y2)) / math.sqrt(mu)
        dtdz = (tof2 - tof_guess) / dz
        
        if abs(dtdz) < 1e-20: break
        z_new = z - (tof_guess - tof) / dtdz
        if abs(z_new - z) < 1e-10: z = z_new; break
        z = z_new
        if z > 4*math.pi*math.pi: z = 4*math.pi*math.pi - 0.1
        if z < -1000: z = -1000
    
    C = _stumpff_c(z)
    S = _stumpff_s(z)
    y = R1 + R2 + A * (z * S - 1.0) / math.sqrt(C)
    if y < 0: return None, None
    
    f = 1.0 - y / R1
    g = A * math.sqrt(y / mu)
    gdot = 1.0 - y / R2
    if abs(g) < 1e-20: return None, None
    
    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    return v1, v2


def porkchop_sweep(r1_vals, rv1_vals, r2_vals, rv2_vals, tof_days,
                    departure_years, mu_sun=GM_SUN):
    """Sweep departure date × TOF, solving Lambert for each.
    
    Uses circular-coplanar approximation for positions.
    Returns numpy arrays: C3_grid, vinf_grid.
    """
    r_earth = SMA_EARTH
    r_mars  = SMA_MARS
    v_earth = math.sqrt(mu_sun / r_earth)
    v_mars  = math.sqrt(mu_sun / r_mars)
    omega_e = v_earth / r_earth * 86400  # rad/day
    omega_m = v_mars  / r_mars  * 86400
    
    n_dep = len(departure_years)
    n_tof = len(tof_days)
    c3_grid = np.full((n_tof, n_dep), np.nan)
    vinf_grid = np.full((n_tof, n_dep), np.nan)
    
    jd_start = 2451545.0 + (departure_years[0] - 2000) * 365.25
    
    for i, dep_year in enumerate(departure_years):
        jd_d = 2451545.0 + (dep_year - 2000) * 365.25
        theta_e = omega_e * (jd_d - jd_start)
        r1 = np.array([r_earth * math.cos(theta_e), r_earth * math.sin(theta_e), 0.0])
        rv1 = np.array([-v_earth * math.sin(theta_e), v_earth * math.cos(theta_e), 0.0])
        
        for j, tof in enumerate(tof_days):
            jd_a = jd_d + tof
            theta_m = omega_m * (jd_a - jd_start) + math.pi  # Mars offset
            r2 = np.array([r_mars * math.cos(theta_m), r_mars * math.sin(theta_m), 0.0])
            rv2 = np.array([-v_mars * math.sin(theta_m), v_mars * math.cos(theta_m), 0.0])
            
            v1, v2 = lambert_solve(r1, r2, tof * 86400, mu_sun, True)
            if v1 is None: continue
            
            vinf_d = np.linalg.norm(v1 - rv1)
            vinf_a = np.linalg.norm(v2 - rv2)
            if vinf_d > 50e3 or vinf_a > 50e3: continue
            
            c3_grid[j, i] = (vinf_d / 1000)**2  # km²/s²
            vinf_grid[j, i] = vinf_d / 1000      # km/s
    
    return c3_grid, vinf_grid
