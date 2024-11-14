# Author: Md Fahim Hasan
# PhD Candidate
# Colorado State University
# Fahim.Hasan@colostate.edu

# This script is a simple replication of part of surface aerodynamic temperature – energy balance (SAT-EB) model.
# Net radiation and soil heat flux estimates are provided, so their estimation is not implemented.
# The script can be modified by integrating other components of SEBAL model to make it a fully functional ET estimation model.


import numpy as np


def calc_heights(hc):
    """
    Calculates zero plane displacement height (d), aerodynamic roughness height for momentum transfer (z_om),
    aerodynamic roughness height for heat and vapor transfer (z_oh).

    :param hc: canopy height (m).

    :return: d, z_om, and z_oh (unit in m).
    """
    d = 0.67 * hc       # m
    z_om = 0.123 * hc    # m
    z_oh = 0.1 * z_om     # m

    return d, z_om, z_oh


def calc_friction_vel(u, d, z_m, z_om, corr_momentum_1, corr_momentum_2, k=0.41):
    """
    Calculates friction velocity (u*).

    :param u: Wind speed at measurement height (m/s).
    :param d: Zero plane displacement height (m).
    :param z_m: Height of wind measurements (m).
    :param z_om: Roughness length governing momentum transfer (m)).
    :param corr_momentum_1: Atmospheric stability correction factor for momentum transfer (m).
    :param corr_momentum_2: Atmospheric stability correction factor for momentum transfer (m).
    :param k: Von Karman's constant. Default set to 0.41.

    :return: friction velocity (m/s).
    """
    friction_vel = (u * k) / (np.log((z_m - d) / z_om) - corr_momentum_1 + corr_momentum_2)

    return friction_vel


def calc_rah(friction_vel, z_m, d, z_oh, corr_heat_1, corr_heat_2, k=0.41):
    """
    Calculates Surface aerodynamic resistance (rah).

    :param friction_vel: friction velocity (m/s).
    :param z_m: Height of wind measurements (m).
    :param d: Zero plane displacement height (m).
    :param z_oh: Aerodynamic roughness length governing transfer of heat and vapor (m).
    :param corr_heat_1: Atmospheric stability correction factor for heat transfer (m).
    :param corr_heat_2: Atmospheric stability correction factor for heat transfer (m).
    :param k: Von Karman's constant. Default set to 0.41.

    :return: Surface aerodynamic resistance (s/m).
    """
    rah = (np.log((z_m - d) / z_oh) - corr_heat_1 + corr_heat_2) / (friction_vel * k)

    return rah


def calc_To(Ts, Ta, LAI, u):
    """
    Calculates aerodynamic temperature using site-specific calibrated empirical equation.

    :param Ts: Radiometric temperature (deg C).
    :param Ta: Air temperature (deg C).
    :param LAI: LAI (m2/m2).
    :param u: Wins speed at measurement height of crop (m/s).

    :return: Aerodynamic temperature (deg C).
    """
    To = 0.534*Ts + 0.39*Ta + 0.224*LAI - 0.192*u + 1.68

    return To


def calc_H_with_To(To, Ta, rah, air_density=1.12, Cp=1004):
    """
    Calculates sensible heat flux (H) for a crop pixel.

    :param To: Aerodynamic temperature (deg C).
    :param Ta: Air temperature (deg C).
    :param rah: Aerodynamic resistance for the crop pixel (s/m).
    :param air_density: Air density (kg/m^3), default is 1.12 kg/m^3.
    :param Cp: Specific heat capacity of air (J/kg/K), default is 1004 J/kg/K.

    :return: Sensible heat flux (H) for the crop pixel (W/m^2).
    """
    To = To + 273.15
    Ta = Ta + 273.15
    H = (air_density * Cp * (To - Ta)) / rah

    return H


def calc_instant_LE(Rn, G, H):
    """
    Calculates instantaneous latent heat flux (LE) for a crop pixel.

    :param Rn: Net radiation (W/m^2).
    :param G: Soil heat flux (W/m^2).
    :param H: Sensible heat flux (W/m^2).

    :return: Instantaneous latent heat flux (LE) at the time of data acquisition (W/m^2).
    """
    LE = Rn - G - H

    return LE


def calc_LE_vap(Ts):
    """
    Calculates latent heat of vaporization based on the temperature.

    :param Ts: Temperature (K).

    :return: Latent heat of vaporization (LE_vap) (J/kg).
    """
    LE_vap = (2.501 - 0.00236 * (Ts - 273.15)) * (10**6)  # unit J/Kg

    return LE_vap


def calc_ET_inst(LE, Ts, t=3600):
    """
    Calculate instantaneous evapotranspiration (ET) in mm/day.

    :param LE: Instantaneous latent heat flux (W/m^2).
    :param Ts: Radiometric temperature (K).
    :param t: Time period (seconds), default is 3600s (1 hour).

    :return: Instantaneous ET (ET_inst) in mm/day.
    """
    LE_vap = calc_LE_vap(Ts)

    ET_inst = (t * LE / (LE_vap * 1000)) * 1000  # unit in mm/hr

    return ET_inst


def calc_ET_fraction(ET_inst, ETref_hour):
    """
    Calculate the evapotranspiration fraction as the ratio of instantaneous ET to ETref.

    :param ET_inst: Instantaneous ET (mm/hr).
    :param ETref_hour: Hourly alfalfa or grass reference ET (mm/hr).

    :return: Evapotranspiration fraction.
    """
    ET_frac = ET_inst / ETref_hour

    return ET_frac


def calc_ET_daily(ET_frac, ETref_daily):
    """
    Calculates daily ET using evaporative f

    :param ET_frac: Evapotranspiration fraction.
    :param ETref_daily: Daily alfalfa or grass reference ET (mm/d).

    :return: daily ET (mm/day).
    """
    ET_daily = ET_frac * ETref_daily

    return ET_daily


if __name__ == '__main__':
    Rn_sat = 700  # Net radiation at the time of the satellite overpass (W/m2)
    G_sat = 115.5  # Soil heat flux at the time of the satellite overpass (W/m2)
    Ts_sat = 32  # Radiometric surface temperature (Ts) at the time of the satellite overpass (deg C)
    Ta = 29.5  # Air temperature (Ta) at 3 m height (measured on site)
    LAI = 3.9  # Leaf area index (m2/m2)
    Cp = 1004  # Specific heat capacity of dry air (J/kg/K)
    air_density = 1.12  # Density of air (kg/m3)
    u = 1.8  # Wind speed measured at a height of 3.5 m above the surface (m/s); correction for height not needed as it is at crop height (not at weather station)
    hc = 1.49  # canopy height (m)

    corr_heat_1 = 0.3       # Ψh((Zm – d)/L) (unit m)
    corr_heat_2 = 0.06      # Ψh(Zoh/L) (unit m)
    corr_momentum_1 = 1.25  # Ψm((Zm – d)/L) (unit m)
    corr_momentum_2 = 0.11  # Ψm(Zom/L) (unit m)
    ETref_hourly_alfalfa = 1.05  # mm/hr
    ETref_daily_alfalfa = 8.9  # mm/d

    z_m = 3.5  # wind speed measurement height (m)

    # Step 1: calculate heights
    d, z_om, z_oh = calc_heights(hc)

    # Step 2: friction velocity (u*)
    fric_vel = calc_friction_vel(u, d, z_m, z_om, corr_momentum_1, corr_momentum_2, k=0.41)

    # Step 3: aerodynamic resistance (rah)
    rah = calc_rah(fric_vel, z_m, d, z_oh, corr_heat_1, corr_heat_2, k=0.41)

    # Step 4: aerodynamic temperature (To)
    To = calc_To(Ts_sat, Ta, LAI, u)

    # Step 5: Sensible heat (H)
    H_sat = calc_H_with_To(To, Ta, rah, air_density, Cp)

    # Step 6: Latent heat flux (LE)
    LE_inst = calc_instant_LE(Rn_sat, G_sat, H_sat)

    # Step 7: Hourly ET
    ET_hourly = calc_ET_inst(LE_inst, (Ts_sat + 273.15), t=3600)

    # Step 8: Evapotranspiration fraction
    ET_frac = calc_ET_fraction(ET_hourly, ETref_hourly_alfalfa)

    # Step 9: Daily ET
    ET_daily = calc_ET_daily(ET_frac, ETref_daily_alfalfa)

    print(f'1) Friction velocity (atmospheric stability corrected): {fric_vel: .2f} m/s')
    print(f'2) Surface aerodynamic resistance: {rah:.2f} s/m')
    print(f'3) Surface aerodynamic temperature estimated: {To:.2f} deg C')
    print(f'4) Sensible heat flux: {H_sat:.2f} W/m2')
    print(f'5) Latent heat flux: {LE_inst:.2f} W/m2')
    print(f'6) Hourly corn actual evapotranspiration (ETi): {ET_hourly:.2f} mm/hr')
    print(f'7) Alfalfa reference evapotranspiration fraction (ETrF): {ET_frac:.2f}')
    print(f'8) Daily corn actual evapotranspiration (ETd): {ET_daily:.2f} mm/d')