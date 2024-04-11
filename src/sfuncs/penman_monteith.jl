# 
function penman_monteith(ETo, G, T, Td, u2, es, ea, Ra)
        """
        Calculates the potential evapotranspiration (PET) using the Penman-Monteith equation.

        Parameters
        ----------
        ETo : Float64
            Reference evapotranspiration (mm/day).
        G : Float64
            Soil heat flux density (mm/day).
        T : Float64
            Air temperature (°C).
        Td : Float64
            Dew point temperature (°C).
        u2 : Float64
            Wind speed at 2 m height (m/s).
        es : Float64
            Saturation vapor pressure (kPa).
        ea : Float64
            Actual vapor pressure (kPa).
        Ra : Float64
            Aerodynamic resistance (s/m).

        Returns
        -------
        PET : Float64
            Potential evapotranspiration (mm/day).
        """
        # Constants
        R = 8.314 # J/mol/K
        cp = 1.013e-3 # kJ/g/K

        # Latent heat of vaporization (MJ/kg)
        Lambda = 2.501 - 0.002361 * T

        # Psychrometric constant (kPa/°C)
        gamma = cp * P / (0.622 * Lambda)

        # Slope of the saturation vapor pressure curve (kPa/°C)
        delta = 4098 * es / (T + 237.3) ^ 2

        # Net radiation (MJ/m2/day)
        Rn = (1 - 0.23) * ETo

        # Air density (kg/m3)
        rho = P * 1000 / (R * (T + 273.15))

        # Specific heat of air (kJ/kg/K)
        cpa = 1.013 * rho ^ -0.0065 * 1000

        # Delta term (MJ/m2/day/°C)
        delta_term = (delta / (delta + gamma)) * (Rn - G)

        # Psi term (MJ/m2/day)
        psi_term = (gamma / (delta + gamma)) * rho * cp * (es - ea) / Ra * u2

        # Potential evapotranspiration (mm/day)
        PET = (delta_term + psi_term) / Lambda

        return PET
    end

    