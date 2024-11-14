# Author: Md Fahim Hasan
# PhD Candidate
# Colorado State University
# Fahim.Hasan@colostate.edu

# This script is a simple replication of part of SEBAL model on how to use cold-hot pixel approach to estimate ETa.
# For simplicity, this script don't require estimation of rah for crop pixel and many crucial parts. Rather they are directly provided.

# The script can be modified by integrating other components of SEBAL model to make it a fully functional ET estimation model.
# For example: Surface Albedo, NDVI, Surface Emissivity, Net Radiation (Rn), Soil Heat Flux (G), Aerodynamic Resistance (rah),
# Roughness Length for Momentum (Z_om), Friction Velocity (u*).

import matplotlib.pyplot as plt

def cold_hot_dT(Rn_hot, G_hot, rah_hot, air_density=1.1, Cp=1004):
    """
    Calculate dT values for cold and hot pixels.

    # SEBAL fixes dT at two “anchor” pixels (cold and hot) and uses a linear relationship between Ts and dT.
    # Using this linear relationship, dT can be estimated using Ts for any pixel in a scene.

    :param Rn_hot: Net radiation at the hot pixel (W/m2).
    :param G_hot: Soil heat flux at the hot pixel (W/m2).
    :param rah_hot: Aerodynamic resistance to heat transport at the hot pixel (s/m) (corrected for atmospheric stability).
    :param air_density: Air density (kg/m3), default is 1.1 kg/m3.
    :param Cp: Specific heat capacity of air (J/kg/K), default is 1004 J/kg/K.

    :return: dT and H values for both cold and hot pixels (dT_cold, dT_hot, H_cold, H_hot).
    """
    # cold pixel
    H_cold = 0                      # considering neutral condition (recently wetted pixel)
    dT_cold = 0

    # hot pixel
    H_hot = Rn_hot - G_hot          # considering that there is crop but not transpiring (very water stressed crop)
    dT_hot = H_hot * rah_hot / (air_density * Cp)

    return dT_cold, dT_hot, H_cold, H_hot


def Ts_dT_linear(Ts_cold, dT_cold, Ts_hot, dT_hot):
    """
    Calculates the linear relationship (slope and intercept) between Ts and dT.

    :param Ts_cold: Radiometric temperature at the cold pixel (K).
    :param dT_cold: Temperature difference at the cold pixel (assumed to be 0).
    :param Ts_hot: Radiometric temperature at the hot pixel (K).
    :param dT_hot: Temperature difference at the hot pixel (K).

    :return: Intercept and slope for the Ts-dT relationship.
    """
    slope = (dT_hot - dT_cold) / (Ts_hot - Ts_cold)
    intercept = - slope * Ts_cold

    return intercept, slope


def plot_Ts_dT_eqn(Ts_cold, dT_cold, Ts_hot, dT_hot, slope, intercept):
    """
    Plots Ts-dT linear equation.
    """
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(Ts_cold, dT_cold)
    ax.scatter(Ts_hot, dT_hot)
    ax.plot([Ts_cold, Ts_hot], [dT_cold, dT_hot], 'k-')
    ax.text(0.2, 0.7, transform=ax.transAxes, s=f'dT = {round(slope, 2)} * Ts {round(intercept, 2)}')
    plt.xlabel('Ts (K)')
    plt.ylabel('dT (K)')
    plt.savefig('figs/Ts_dT_plot_SEBAL.png')


def calc_dT_crop(Ts, intercept, slope):
    """
    Calculate dT for a crop pixel using the Ts-dT linear relationship.

    :param Ts: Radiometric temperature of the crop pixel (K).
    :param intercept: Intercept from Ts-dT linear relationship.
    :param slope: Slope from Ts-dT linear relationship.

    :return: dT for the crop pixel.
    """
    dT_crop = intercept + slope * Ts

    return dT_crop


def calc_H(dT_crop, rah_crop, air_density=1.1, Cp=1004):
    """
    Calculate sensible heat flux (H) for a crop pixel during satellite overpass.

    :param dT_crop: Temperature difference for the crop pixel (K).
    :param rah_crop: Aerodynamic resistance for the crop pixel (s/m).
    :param air_density: Air density (kg/m3), default is 1.1 kg/m3.
    :param Cp: Specific heat capacity of air (J/kg/K), default is 1004 J/kg/K.

    :return: Sensible heat flux (H) for the crop pixel (W/m2).
    """
    H = (air_density * Cp * dT_crop) / rah_crop

    return H


def calc_instant_LE(Rn, G, H):
    """
    Calculate instantaneous latent heat flux (LE) for a crop pixel.

    :param Rn: Net radiation at the crop pixel (W/m2).
    :param G: Soil heat flux at the crop pixel (W/m2).
    :param H: Sensible heat flux at the crop pixel (W/m2).

    :return: Instantaneous latent heat flux (LE) for the crop pixel (W/m2).
    """
    LE = Rn - G - H

    return LE


def calc_LE_vap(Ts_cold):
    """
    Calculate latent heat of vaporization based on the temperature at the cold pixel.

    :param Ts_cold: Radiometric temperature at the cold pixel (K).

    :return: Latent heat of vaporization (LE_vap) (J/kg).
    """
    LE_vap = (2.501 - 0.00236 * (Ts_cold - 273.15)) * (10**6)  # unit J/Kg

    return LE_vap


def calc_ET_inst(LE, Ts_cold, t=3600):
    """
    Calculate instantaneous evapotranspiration (ET) in mm/day.

    :param LE: Instantaneous latent heat flux (W/m2).
    :param Ts_cold: Radiometric temperature at the cold pixel (K).
    :param t: Time period (seconds), default is 3600s (1 hour).

    :return: Instantaneous ET (ET_inst) in mm/day.
    """
    LE_vap = calc_LE_vap(Ts_cold)

    ET_inst = (t * LE / (LE_vap * 1000)) * 1000  # unit in mm/hr

    return ET_inst


def calc_EF(LE, Rn_crop, G_crop):
    """
    Calculate the evaporative fraction (EF), which is the fraction of available energy used for evapotranspiration.
    This approach is based on the self-preservation theory of daytime preservation theory of daytime
    fluxes which states that the ratio between the latent heat flux (LE)
    and the available energy (Rn –G) remains constant during the day) remains constant during the day.

    :param LE: Instantaneous latent heat flux (W/m2).
    :param Rn_crop: Net radiation at the crop pixel (W/m2).
    :param G_crop: Soil heat flux at the crop pixel (W/m2).

    :return: Evaporative fraction.
    """
    EF = LE / (Rn_crop - G_crop)

    return EF


def calc_ET_daily(Ts_cold, EF, Rn_daily_avg):
    """
    Calculate daily evapotranspiration (ET) for a crop pixel.

    :param Ts_cold: Radiometric temperature at the cold pixel (K).
    :param EF: Evaporative fraction, the fraction of available energy used for evapotranspiration.
    :param Rn_daily_avg: Average daily net radiation at the crop pixel (W/m2).

    :return: Daily evapotranspiration (ET_daily) in mm/day.
    """
    LE_vap = calc_LE_vap(Ts_cold)

    # the total soil heat flux for an entire day is assumed to tend to zero for vegetation and most bare soils
    ET_daily = (86400 * EF * Rn_daily_avg / (LE_vap * 1000)) * 1000  # unit in mm/day

    return ET_daily


if __name__ == '__main__':
    # Provided inputs
    Ts_cold = 300.95  # Cold pixel radiometric surface temperature (K)
    Ts_hot = 317.75  # Hot pixel radiometric surface temperature (K)
    Rn_cold = 800  # Cold pixel net radiation (W/m2)
    Rn_hot = 580  # Hot pixel net radiation (W/m2)
    G_cold = 40  # Cold pixel soil heat flux (W/m2)
    G_hot = 261  # Hot pixel soil heat flux (W/m2)
    rah_hot = 17  # Hot pixel surface aerodynamic resistance (s/m)
    rah_crop = 20  # Crop pixel surface aerodynamic resistance (s/m)
    Cp = 1004  # Specific heat capacity of dry air (J/kg/K)
    air_density = 1.1  # Density of air (kg/m3)
    Ts_crop = 305.5  # Crop radiometric surface temperature (K)
    Rn_sat = 725  # Crop pixel net radiation at satellite overpass (W/m2)
    G_sat = 54.4  # Crop pixel soil heat flux at satellite overpass (W/m2)
    Rn_daily_avg = 215  # Daily average daily net radiation (W/m2)

    # Step 1: cold and hot pixel dT values
    dT_cold, dT_hot, H_cold, H_hot = cold_hot_dT(Rn_hot, G_hot, rah_hot, air_density, Cp)

    # Step 2: linear relationship (slope and intercept) between Ts and dT
    intercept, slope = Ts_dT_linear(Ts_cold, dT_cold, Ts_hot, dT_hot)

    # Step 3: plot Ts-dT relationship plot
    plot_Ts_dT_eqn(Ts_cold, dT_cold, Ts_hot, dT_hot, slope, intercept)

    # Step 4: dT for the crop pixel
    dT_crop = calc_dT_crop(Ts_crop, intercept, slope)

    # Step 5: sensible heat flux (H) for the crop pixel
    H_sat = calc_H(dT_crop, rah_crop, air_density, Cp)

    # Step 6: instantaneous latent heat flux (LE) for the crop pixel
    LE_sat = calc_instant_LE(Rn_sat, G_sat, H_sat)

    # Step 7: hourly evapotranspiration rate for the crop pixel
    ET_hourly = calc_ET_inst(LE_sat, Ts_cold, t=3600)

    # Step 8: evaporative fraction (EF) for the crop pixel
    EF = calc_EF(LE_sat, Rn_sat, G_sat)

    # Step 9: Calculate daily evapotranspiration rate for the crop pixel
    ET_daily = calc_ET_daily(Ts_cold, EF, Rn_daily_avg)

    # Print results using f-strings
    print(f'1. Cold pixel dT value: {dT_cold:.2f} K')
    print(f'2. Hot pixel dT value: {dT_hot:.2f} K')
    print(f'3. Cold pixel sensible heat flux: {H_cold:.2f} W/m2')
    print(f'4. Hot pixel sensible heat flux: {H_hot:.2f} W/m2')
    print(f'5. dT equation (slope and intercept): slope = {slope:.2f}, intercept = {intercept:.2f}')
    print(f'6. Crop pixel sensible heat flux: {H_sat:.2f} W/m2')
    print(f'7. Crop pixel instantaneous latent heat flux: {LE_sat:.2f} W/m2')
    print(f'8. Crop pixel hourly evapotranspiration rate: {ET_hourly:.2f} mm/hr')
    print(f'9. Crop pixel evaporative fraction: {EF:.2f}')
    print(f'10. Crop pixel daily evapotranspiration rate: {ET_daily:.2f} mm/day')

