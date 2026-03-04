
from galaxies import Galaxy

import astropy.units as u
from astropy.coordinates import SkyCoord

# Define galaxy parameters for our set.

# ngc6822
gal_ngc6822 = Galaxy("ncg6822")

# Namumba+2017 (https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.3761N/abstract)
gal_ngc6822.vsys = -55 * u.km/u.s
gal_ngc6822.inclination = 66 * u.deg
gal_ngc6822.position_angle = 118 * u.deg
gal_ngc6822.center_position = SkyCoord("19h44m58.0s −14d48m11.9s")
gal_ngc6822.distance = 0.526 * u.Mpc

# Vrot = 60

# Namumba+2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3365N/abstract)
gal_ic10 = Galaxy("ic10")
gal_ic10.vsys = -351 * u.km/u.s
gal_ic10.inclination = 47 * u.deg
gal_ic10.position_angle = 65 * u.deg
gal_ic10.center_position = SkyCoord("00h20m24.6s +59d17m30.0s")
gal_ic10.distance = 0.77 * u.Mpc

# Vrot = 30

# From OH+2015: Vrot = 36.5

# Oh+2015 (https://ui.adsabs.harvard.edu/abs/2015AJ....149..180O/abstract)
gal_ic1613 = Galaxy("ic1613")
gal_ic1613.vsys = -232 * u.km/u.s
gal_ic1613.inclination = 48 * u.deg
gal_ic1613.position_angle = 73.7 * u.deg
gal_ic1613.center_position = SkyCoord("01h04m49.6s +02d08m14.1s")
gal_ic1613.distance = 0.760 * u.Mpc

# Vrot = 17

# Oh+2015 (https://ui.adsabs.harvard.edu/abs/2015AJ....149..180O/abstract)
gal_wlm = Galaxy("wlm")
gal_wlm.vsys = -122 * u.km/u.s
gal_wlm.inclination = 74 * u.deg
gal_wlm.position_angle = 174.5 * u.deg
gal_wlm.center_position = SkyCoord("00h01m59.9s -15d27m57.2s")
gal_wlm.distance = 0.984 * u.Mpc

# Vrot = 37

# m33
gal_m33 = Galaxy("M33")

# Koch+2018
gal_m33.vsys = -180 * u.km/u.s
gal_m33.inclination = 55 * u.deg
gal_m33.position_angle = 201 * u.deg
# gal_m33.center_position = SkyCoord("19h44m58.0s −14d48m11.9s")  # Trust default.
gal_m33.distance = 0.869 * u.Mpc  #

# Vrot = 110


# m31
gal_m31 = Galaxy("M31")

# Chemin+2009
gal_m31.vsys = -300 * u.km/u.s
gal_m31.inclination = 77 * u.deg
gal_m31.position_angle = 35 * u.deg
# gal_m31.center_position = SkyCoord("19h44m58.0s −14d48m11.9s")  # Trust default.
gal_m31.distance = 0.776 * u.Mpc  #

# Vrot = 240

