# Atmosphere models — density, pressure, temperature for 8 planetary bodies
# Mirrors planetary_environment.hpp

import math

def atmosphere_earth(alt_m):
    """US Standard Atmosphere 1976, 7-layer. Returns (rho, press, tempk)."""
    h = alt_m / 1000.0
    if h < 0: h = 0
    if h > 84.852: h = 84.852
    
    H  = [0, 11, 20, 32, 47, 51, 71, 84.852]
    LR = [-6.5, 0, 1.0, 2.8, 0, -2.8, -2.0]
    Tb = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
    R, g0, P0 = 287.058, 9.80665, 101325.0
    
    i = 0
    for k in range(1, 7):
        if h >= H[k]: i = k
    
    lapse, T_base = LR[i], Tb[i]
    tempk = T_base if abs(lapse) < 1e-9 else T_base + lapse * (h - H[i])
    
    Pb = P0
    for k in range(i):
        dh = H[k+1] - H[k]
        if abs(LR[k]) < 1e-9:
            Pb *= math.exp(-g0 * dh * 1000 / (R * Tb[k]))
        else:
            Tt = Tb[k] + LR[k] * dh
            Pb *= (Tt / Tb[k]) ** (-g0 / (R * LR[k] * 0.001))
    
    dh = h - H[i]
    if abs(lapse) < 1e-9:
        press = Pb * math.exp(-g0 * dh * 1000 / (R * T_base))
    else:
        Tt = T_base + lapse * dh
        press = Pb * (Tt / T_base) ** (-g0 / (R * lapse * 0.001))
    
    rho = press / (R * tempk)
    return max(rho, 0), max(press, 0), tempk


def atmosphere_mars(alt_m):
    """NASA GRC Mars atmosphere. Returns (rho, press, tempk)."""
    T_c = -31.0 - 0.000998 * alt_m if alt_m < 7000 else -23.4 - 0.00222 * alt_m
    tempk = T_c + 273.1
    p_kpa = 0.699 * math.exp(-0.00009 * alt_m)
    press = p_kpa * 1000.0
    rho = p_kpa / (0.1921 * tempk)
    return max(rho, 0), max(press, 0), tempk


def atmosphere_jupiter(alt_m):
    """Galileo Probe ASI data. Altitude above 1-bar level. Returns (rho, press, tempk)."""
    h = alt_m / 1000.0
    if h > 50:        tempk = 110.0 + 0.05 * (h - 50)
    elif h > 0:       tempk = 166.0 - 1.12 * h
    elif h > -60:     tempk = 166.0 - 2.0 * h
    else:             tempk = 286.0 - 2.2 * (h + 60)
    if tempk < 100: tempk = 100
    
    scale_h = 3745.0 * tempk / 24.79
    press = 1.0e5 * math.exp(-alt_m / scale_h)
    rho = press / (3745.0 * tempk)
    return max(rho, 0), max(press, 0), tempk


def atmosphere_venus(alt_m):
    """Venera/Pioneer Venus data. Returns (rho, press, tempk)."""
    h = alt_m / 1000.0
    if h < 0: h = 0
    if h > 100: h = 100
    at = [0, 10, 20, 30, 40, 50, 60, 70, 80, 100]
    tc = [462, 385, 306, 222, 143, 75, -10, -43, -76, -112]
    pa = [92.1, 47.39, 22.52, 9.851, 3.501, 1.066, 0.2357, 0.0369, 0.00476, 0.0000266]
    
    i = 0
    for k in range(9):
        if h >= at[k]: i = k
    if i >= 9: i = 8
    
    frac = (h - at[i]) / (at[i+1] - at[i])
    tempk = (tc[i] + frac * (tc[i+1] - tc[i])) + 273.15
    lp = math.log(pa[i]) + frac * (math.log(pa[i+1]) - math.log(pa[i]))
    press = math.exp(lp) * 101325.0
    rho = press / (191.0 * tempk)
    return max(rho, 0), max(press, 0), tempk


def atmosphere_titan(alt_m):
    """Cassini-Huygens data. Returns (rho, press, tempk)."""
    h = alt_m / 1000.0
    if h > 250:       tempk = 170.0 + 0.1 * (h - 250)
    elif h > 40:      tempk = 70.0 + (170 - 70) / (250 - 40) * (h - 40)
    elif h > 0:       tempk = 94.0 - 0.6 * h
    else:             tempk = 94.0
    
    scale_h = 296.8 * tempk / 1.352
    press = 1.47e5 * math.exp(-alt_m / scale_h)
    rho = press / (296.8 * tempk)
    return max(rho, 0), max(press, 0), tempk


def atmosphere_neptune(alt_m):
    """Voyager 2 RSS. Altitude above 1-bar. Returns (rho, press, tempk)."""
    h = alt_m / 1000.0
    if h > 50:        tempk = 55.0 + 0.4 * (h - 50)
    elif h > 0:       tempk = 72.0 - 0.34 * h
    else:             tempk = 72.0 - 0.9 * h
    
    scale_h = 3615.0 * tempk / 11.15
    press = 1.0e5 * math.exp(-alt_m / scale_h)
    rho = press / (3615.0 * tempk)
    return max(rho, 0), max(press, 0), tempk


def atmosphere_uranus(alt_m):
    """Voyager 2 RSS. Altitude above 1-bar. Returns (rho, press, tempk)."""
    h = alt_m / 1000.0
    if h > 50:        tempk = 53.0 + 0.3 * (h - 50)
    elif h > 0:       tempk = 76.0 - 0.46 * h
    else:             tempk = 76.0 - 0.8 * h
    
    scale_h = 3641.0 * tempk / 8.87
    press = 1.0e5 * math.exp(-alt_m / scale_h)
    rho = press / (3641.0 * tempk)
    return max(rho, 0), max(press, 0), tempk


# Map body name to atmosphere function
ATMOSPHERES = {
    'earth':   atmosphere_earth,
    'mars':    atmosphere_mars,
    'jupiter': atmosphere_jupiter,
    'venus':   atmosphere_venus,
    'titan':   atmosphere_titan,
    'neptune': atmosphere_neptune,
    'uranus':  atmosphere_uranus,
}

def get_atmosphere(body, alt_m):
    """Get (rho, press, tempk) for a named body at altitude (m)."""
    fn = ATMOSPHERES.get(body.lower())
    if fn is None:
        raise ValueError(f"Unknown body: {body}. Options: {list(ATMOSPHERES.keys())}")
    return fn(alt_m)
