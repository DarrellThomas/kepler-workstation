# Orbital mechanics — vis-viva, circular orbits, Kepler's laws

import math
from .constants import G0

def circular_velocity(gm, r):
    """Orbital velocity for a circular orbit (m/s).
    
    Args:
        gm: gravitational parameter (m³/s²)
        r:  orbital radius from center of body (m)
    """
    return math.sqrt(gm / r)

def circular_period(gm, r):
    """Orbital period for a circular orbit (seconds)."""
    return 2.0 * math.pi * math.sqrt(r**3 / gm)

def vis_viva(gm, r, a):
    """Velocity at distance r on an orbit with semi-major axis a.
    
    Args:
        gm: gravitational parameter (m³/s²)
        r:  current distance from center (m)
        a:  semi-major axis (m), negative for hyperbolic
    """
    return math.sqrt(gm * (2.0/r - 1.0/a))

def hohmann_transfer(r1, r2, gm_central, gm_depart, gm_arrive,
                      r_depart_body, r_arrive_body, alt_park=300e3):
    """Hohmann transfer between two circular coplanar orbits.
    
    Returns dict with C3, dv_depart, dv_arrive, dv_total, tof_days.
    """
    a_trans = (r1 + r2) / 2.0
    v_p1 = math.sqrt(gm_central / r1)
    v_p2 = math.sqrt(gm_central / r2)
    v_t1 = math.sqrt(gm_central * (2.0/r1 - 1.0/a_trans))
    v_t2 = math.sqrt(gm_central * (2.0/r2 - 1.0/a_trans))
    
    v_inf_d = v_t1 - v_p1
    v_inf_a = v_p2 - v_t2
    
    r_park1 = r_depart_body + alt_park
    r_park2 = r_arrive_body + alt_park
    v_circ1 = math.sqrt(gm_depart / r_park1)
    v_circ2 = math.sqrt(gm_arrive / r_park2)
    
    dv_dep = math.sqrt(v_inf_d**2 + 2*gm_depart/r_park1) - v_circ1
    dv_arr = math.sqrt(v_inf_a**2 + 2*gm_arrive/r_park2) - v_circ2
    
    period = 2.0 * math.pi * math.sqrt(a_trans**3 / gm_central)
    tof_days = period / (2.0 * 86400.0)
    
    return {
        'c3_depart': v_inf_d**2,
        'v_inf_depart': v_inf_d,
        'v_inf_arrive': v_inf_a,
        'dv_depart': dv_dep,
        'dv_arrive': dv_arr,
        'dv_total': dv_dep + dv_arr,
        'tof_days': tof_days,
        'a_transfer_au': a_trans / 149597870700.0
    }

def rocket_equation(dv, isp, payload_kg, inert_frac=0.12):
    """Mass budget from Tsiolkovsky rocket equation.
    
    Args:
        dv:         delta-V (m/s)
        isp:        specific impulse (seconds)
        payload_kg: final mass delivered (kg)
        inert_frac: tank + structure fraction of propellant mass
    
    Returns dict with propellant_kg, initial_kg, mass_ratio.
    """
    ve = isp * G0
    mr = math.exp(dv / ve)
    propellant = payload_kg * (mr - 1.0) * (1.0 + inert_frac)
    initial = payload_kg * mr * (1.0 + inert_frac * (mr - 1.0))
    return {
        'mass_ratio': mr,
        'propellant_kg': propellant,
        'initial_kg': initial,
        'inert_kg': propellant * inert_frac
    }
