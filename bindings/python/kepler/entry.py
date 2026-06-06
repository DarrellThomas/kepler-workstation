# Atmospheric entry simulation — ballistic entry, corridor analysis

import math
from .atmosphere import get_atmosphere

def sutton_graves_heat(rho, vel, rn):
    """Sutton-Graves stagnation heat flux (W/m²).
    
    Args:
        rho: density (kg/m³)
        vel: velocity (m/s)
        rn:  nose radius (m)
    """
    if rho <= 0 or rn <= 0:
        return 0.0
    return 1.83e-4 * math.sqrt(rho / rn) * vel**3 * 1000.0  # W/m²


def simulate_entry(body, alt_entry, vel_entry, fpa_entry_deg,
                    mass_kg=5500, area_m2=12.57, cd=1.3, rn_m=2.0,
                    g_limit=10.0, q_limit_kw=1500.0,
                    max_steps=200000):
    """Simulate a 3-DOF ballistic atmospheric entry.
    
    Returns dict with trajectory data and boundary results.
    """
    from .constants import GM_EARTH, R_EARTH, GM_MARS, R_MARS, \
        GM_JUPITER, R_JUPITER, GM_VENUS, R_VENUS, \
        GM_TITAN, R_TITAN, GM_NEPTUNE, R_NEPTUNE, GM_URANUS, R_URANUS
    
    body_params = {
        'earth':   (GM_EARTH,   R_EARTH),
        'mars':    (GM_MARS,    R_MARS),
        'jupiter': (GM_JUPITER, R_JUPITER),
        'venus':   (GM_VENUS,   R_VENUS),
        'titan':   (GM_TITAN,   R_TITAN),
        'neptune': (GM_NEPTUNE, R_NEPTUNE),
        'uranus':  (GM_URANUS,  R_URANUS),
    }
    
    gm, radius = body_params[body.lower()]
    beta = mass_kg / (cd * area_m2)  # ballistic coefficient
    
    alt, vel, fpa = alt_entry, vel_entry, math.radians(fpa_entry_deg)
    dt = 0.1
    
    # Track peaks
    g_max, q_max, alt_min = 0.0, 0.0, alt_entry
    skip_out, impact, g_exceed, q_exceed = False, False, False, False
    
    # For trajectory recording (sampled)
    traj_alt, traj_vel, traj_g, traj_q = [], [], [], []
    sample_every = max(1, max_steps // 500)
    
    for step in range(max_steps):
        # Atmosphere
        rho, press, tempk = get_atmosphere(body, alt)
        
        # Gravity
        r = radius + alt
        g = gm / (r * r)
        
        # Drag
        drag = 0.5 * rho * vel * vel / beta
        g_load = drag / 9.80665
        
        # Heat rate
        q_kw = sutton_graves_heat(rho, vel, rn_m) / 1000.0 if rho > 1e-15 else 0.0
        
        # Track peaks
        if g_load > g_max:
            g_max = g_load
        if q_kw > q_max:
            q_max = q_kw
        if alt < alt_min:
            alt_min = alt
        
        # Sample
        if step % sample_every == 0 and alt > 0:
            traj_alt.append(alt/1000.0)
            traj_vel.append(vel/1000.0)
            traj_g.append(g_load)
            traj_q.append(q_kw)
        
        # Terminal conditions
        if alt > alt_entry * 2 and fpa > 0:
            skip_out = True
            break
        if alt <= 0 or vel < 100:
            impact = True
            break
        if g_load > g_limit:
            g_exceed = True
            break
        if q_kw > q_limit_kw:
            q_exceed = True
            break
        
        # Adaptive timestep
        if rho > 1e-6:
            dt = 0.1 / math.sqrt(rho * 1000)
            dt = max(0.001, min(dt, 5.0))
        else:
            dt = 5.0
        
        # RK4 step (simplified — Euler for brevity, good enough for this grid)
        alt_dot = vel * math.sin(fpa)
        vel_dot = -drag - g * math.sin(fpa)
        fpa_dot = (vel / r - g / vel) * math.cos(fpa) if vel > 0 else 0
        
        alt += alt_dot * dt
        vel += vel_dot * dt
        fpa += fpa_dot * dt
    
    return {
        'skip_out': skip_out,
        'impact': impact,
        'g_exceed': g_exceed,
        'q_exceed': q_exceed,
        'g_max': g_max,
        'q_max_kw': q_max,
        'alt_min_km': alt_min / 1000.0,
        'trajectory': {
            'alt_km': traj_alt,
            'vel_kms': traj_vel,
            'g_load': traj_g,
            'q_kw': traj_q,
        }
    }


def find_corridor(body, alt_entry, vel_entry, mass_kg=5500, area_m2=12.57,
                   cd=1.3, rn_m=2.0, g_limit=10, q_limit_kw=1500):
    """Sweep entry angles to find the viable corridor.
    
    Returns dict with corridor boundaries and a sweep of results.
    """
    sweep = []
    skip_out_deg, g_limit_deg, q_limit_deg = None, None, None
    
    angles = []
    fpa = -0.5
    while fpa >= -60:
        angles.append(fpa)
        if fpa > -3:
            fpa -= 0.1
        elif fpa > -10:
            fpa -= 0.25
        elif fpa > -20:
            fpa -= 0.5
        else:
            fpa -= 1.0
    
    for fpa_deg in angles:
        r = simulate_entry(body, alt_entry, vel_entry, fpa_deg,
                           mass_kg, area_m2, cd, rn_m, g_limit, q_limit_kw)
        sweep.append((fpa_deg, r))
        
        if skip_out_deg is None and not r['skip_out'] and not r['g_exceed'] and not r['q_exceed']:
            skip_out_deg = fpa_deg
        if not r['g_exceed']:
            g_limit_deg = fpa_deg
        if not r['q_exceed']:
            q_limit_deg = fpa_deg
    
    return {
        'skip_out_deg': skip_out_deg,
        'g_limit_deg': g_limit_deg,
        'q_limit_deg': q_limit_deg,
        'corridor_width': abs((g_limit_deg or 0) - (skip_out_deg or 0)),
        'sweep': sweep
    }
