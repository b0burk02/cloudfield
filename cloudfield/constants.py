"""Physical and unit constants shared across modules."""
PC_CM = 3.086e18             # cm per pc
ALPHA_PREF = 59.0            # alphaG = 59 * I_UV / n  (I_UV Draine-normalized, n total H nuclei)
ISRF_BUNIT = "Draine78 ISRF" # BUNIT of the LightCube combined 912-2000 A cube
SOLAR_I_UV = 0.633           # Bish+2026 Sec III.1: 63.3% of the Draine Standard ISRF at the Sun
SOLAR_RADIUS_PC = 25.0       # averaging radius used by Bish+2026 for the solar value
R_NO_DUST = 69.0             # pc; Edenhofer+2024 has no dust inside this radius (cells are 0)
