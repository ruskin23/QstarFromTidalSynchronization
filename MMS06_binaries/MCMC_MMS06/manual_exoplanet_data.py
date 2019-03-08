"""Manually collected from the literature exoplanet system data."""

from astropy.constants import R_sun, M_sun
from math import pi

rho_sun = (M_sun / (4.0 * pi /3.0 * R_sun**3)).to('g/cm3').value

data = (
    dict(
        pl_hostname   = 'CoRoT-1',
        st_dens       = 0.660 * rho_sun,
        st_denserr1   = 0.019 * rho_sun,
        st_denserr2   = -0.019 * rho_sun
    ),
    dict(
        pl_hostname   = 'CoRoT-2',
        st_metfe      = 0.04,
        st_metfeerr1  = 0.05,
        st_metfeerr2  = -0.05
    ),
    dict(
        pl_hostname   = 'HAT-P-13',
        st_dens       = 0.244 * rho_sun,
        st_denserr1   = 0.013 * rho_sun,
        st_denserr2   = -0.013 * rho_sun
    ),
    dict(
        pl_hostname   = 'HAT-P-36',
        st_logg       = 4.370,
        st_loggerr1   = 0.04,
        st_loggerr2   = -0.04,
        st_metfe      = 0.22,
        st_metfeerr1  = 0.04,
        st_metfeerr2  = -0.04,
        pl_rvamp      = 343.1,
        pl_rvamperr1  = 21.3,
        pl_rvamperr2  = -21.3
    ),
    dict(
        pl_hostname   = 'HAT-P-37',
        st_dens       = 1.33 * rho_sun,
        st_denserr1   = 0.13 * rho_sun,
        st_denserr2   = -0.11 * rho_sun
    ),
    dict(
        pl_hostname   = 'KELT-1',
        st_logg       = 4.228,
        st_loggerr1   = 0.014,
        st_loggerr2   = -0.021,
        st_metfe      = 0.052,
        st_metfeerr1  = 0.079,
        st_metfeerr2  = -0.079,
        st_dens       = 0.594,
        st_denserr1   = 0.027,
        st_denserr2   = -0.042,
        pl_ratdor     = 3.619,
        pl_ratdorerr1 = 0.055,
        pl_ratdorerr2 = -0.087
    ),
    dict(
        pl_hostname   = 'OGLE-TR-056',
        st_logg       = 4.258,
        st_loggerr1   = 0.043,
        st_loggerr2   = -0.043,
        st_metfe      = 0.22,
        st_metfeerr1  = 0.1,
        st_metfeerr2  = -0.1,
        st_dens       = 0.62 * rho_sun,
        st_denserr1   = 0.21 * rho_sun,
        st_denserr2   = -0.21 * rho_sun,
        st_vsini      = 3.2,
        st_vsinierr1  = 1.0,
        st_vsinierr2  = -1.0
    ),
    dict(
        pl_hostname   = 'OGLE-TR-113',
        st_logg       = 4.552 ,
        st_loggerr1   = 0.009,
        st_loggerr2   = -0.017,
        st_metfe      = 0.090,
        st_metfeerr1  = 0.08,
        st_metfeerr2  = -0.08,
        st_dens       = 1.679 * rho_sun,
        st_denserr1   = 0.064 * rho_sun,
        st_denserr2   = -0.062 * rho_sun,
        pl_ratdor     = 6.37,
        pl_ratdorerr1 = 0.166,
        pl_ratdorerr2 = -0.166,
        pl_rvamp      = 267.0,
        pl_rvamperr1  = 34.0,
        pl_rvamperr2  = -34.0
    ),
    dict(
        pl_hostname   = 'Qatar-1',
        st_teff       = 4860.0,
        st_tefferr1   = 125.0,
        st_tefferr2   = 125.0
    ),
    dict(
        pl_hostname   = 'Qatar-2',
        pl_ratdor     = 6.50,
        pl_ratdorerr1 = 0.197,
        pl_ratdorerr2 = -0.197
    ),
    dict(
        pl_hostname   = 'TrES-2',
        st_dens       = 1.105 * rho_sun,
        st_denserr1   = 0.011 * rho_sun,
        st_denserr2   = -0.011 * rho_sun
    ),
    dict(
        pl_hostname   = 'TrES-3',
        st_logg       = 4.581,
        st_loggerr1   = 0.017,
        st_loggerr2   = -0.012,
        st_metfe      = -0.190,
        st_metfeerr1  = 0.08,
        st_metfeerr2  = -0.08,
        st_dens       = 1.648 * rho_sun,
        st_denserr1   = 0.041 * rho_sun,
        st_denserr2   = -0.041 * rho_sun,
        pl_ratdor     = 6.03,
        pl_ratdorerr1 = 0.176,
        pl_ratdorerr2 = -0.176,
        pl_rvamp      = 369.0,
        pl_rvamperr1  = 11.0,
        pl_rvamperr2  = -11.0
    ),
    dict(
        pl_hostname   = 'WASP-1',
        st_teff       = 6110.0,
        st_tefferr1   = 50.0,
        st_tefferr2   = 50.0
    ),
    dict(
        pl_hostname   = 'WASP-12',
        st_logg       = 4.38,
        st_loggerr1   = 0.1,
        st_loggerr2   = -0.1,
        st_metfe      = 0.30,
        st_metfeerr1  = 0.05,
        st_metfeerr2  = -0.15,
        st_dens       = 0.325 * rho_sun,
        st_denserr1   = 0.016 * rho_sun,
        st_denserr2   = -0.016 * rho_sun,
        pl_ratdor     = 2.98,
        pl_ratdorerr1 = 0.154,
        pl_ratdorerr2 = -0.154
    ),
    dict(
        pl_hostname   = 'WASP-18',
        st_logg       = 4.367,
        st_loggerr1   = 0.028,
        st_loggerr2   = -0.042,
        st_metfe      = 0.0,
        st_metfeerr1  = 0.09,
        st_metfeerr2  = -0.09,
        st_dens       = 0.687 * rho_sun,
        st_denserr1   = 0.062 * rho_sun,
        st_denserr2   = -0.062 * rho_sun,
        pl_ratdor     = 3.57,
        pl_ratdorerr1 = 0.187,
        pl_ratdorerr2 = -0.187
    ),
    dict(
        pl_hostname   = 'WASP-19',
        pl_ratdor     = 3.552,
        pl_ratdorerr1 = 0.093,
        pl_ratdorerr2 = -0.093,
        st_dens       = 0.893 * rho_sun,
        st_denserr1   = 0.015 * rho_sun,
        st_denserr2   = -0.015 * rho_sun,
        pl_rvamp      = 257.0,
        pl_rvamperr1  = 3.0,
        pl_rvamperr2  = -3.0
    ),
    dict(
        pl_hostname   = 'WASP-23',
        pl_rvamp      = 145.8,
        pl_rvamperr1  = 1.5,
        pl_rvamperr2  = -2.1,
        st_dens       = 1.843 * rho_sun,
        st_denserr1   = 0.025 * rho_sun,
        st_denserr2   = -0.027 * rho_sun,
        st_metfe      = 0.05,
        st_metfeerr1  = 0.1, #bullshit value
        st_metfeerr2  = -0.1, #bullshit value
    ),
    dict(
        pl_hostname   = 'WASP-24',
        st_metfe      = 0.07,
        st_metfeerr1  = 0.1,
        st_metfeerr2  = -0.1,
        st_dens       = 0.502 * rho_sun,
        st_denserr1   = 0.03 * rho_sun,
        st_denserr2   = -0.03 * rho_sun
    ),
    dict(
        pl_hostname   = 'WASP-26',
        st_metfe      = -0.02,
        st_metfeerr1  = 0.09,
        st_metfeerr2  = -0.09,
        st_dens       = 0.46 * rho_sun,
        st_denserr1   = 0.06 * rho_sun,
        st_denserr2   = -0.06 * rho_sun,
        pl_rvamp      = 135.5,
        pl_rvamperr1  = 3.5,
        pl_rvamperr2  = -3.5
    ),
    dict(
        pl_hostname   = 'WASP-3',
        st_logg       = 4.250,
        st_loggerr1   = 0.05,
        st_loggerr2   = -0.05,
        st_metfe      = 0.0,
        st_metfeerr1  = 0.2,
        st_metfeerr2  = -0.2,
        st_dens       = 0.495 * rho_sun,
        st_denserr1   = 0.024 * rho_sun,
        st_denserr2   = -0.024 * rho_sun,
        pl_ratdor     = 5.18,
        pl_ratdorerr1 = 0.35,
        pl_ratdorerr2 = -0.35
    ),
    dict(
        pl_hostname   = 'WASP-32',
        st_metfe      = -0.13,
        st_metfeerr1  = 0.1,
        st_metfeerr2  = -0.1,
        st_dens       = 0.8 * rho_sun,
        st_denserr1   = 0.1 * rho_sun,
        st_denserr2   = -0.1 * rho_sun,
        pl_rvamp      = 487.0,
        pl_rvamperr1  = 5.0,
        pl_rvamperr2  = -5.0
    ),
    dict(
        pl_hostname   = 'WASP-33',
        st_logg       = 4.30,
        st_loggerr1   = 0.2,
        st_loggerr2   = -0.2,
        st_metfe      = 0.1,
        st_metfeerr1  = 0.2,
        st_metfeerr2  = -0.2,
        st_dens       = 0.700,
        st_denserr1   = 0.034,
        st_denserr2   = -0.034
    ),
    dict(
        pl_hostname   = 'WASP-36',
        st_logg       = 4.498,
        st_loggerr1   = 0.012,
        st_loggerr2   = -0.012,
        st_metfe      = -0.31,
        st_metfeerr1  = 0.12,
        st_metfeerr2  = -0.12,
        st_dens       = 1.715,
        st_denserr1   = 0.075,
        st_denserr2   = -0.068,
        pl_ratdor     = 6.00,
        pl_ratdorerr1 = 0.157,
        pl_ratdorerr2 = -0.157
    ),
    dict(
        pl_hostname   = 'WASP-4',
        st_dens       = 1.23 * rho_sun,
        st_denserr1   = 0.022 * rho_sun,
        st_denserr2   = -0.022 * rho_sun,
        pl_ratdor     = 5.43,
        pl_ratdorerr1 = 0.162,
        pl_ratdorerr2 = -0.162,
        pl_rvamp      = 234.6,
        pl_rvamperr1  = 2.2,
        pl_rvamperr2  = -2.3
    ),
    dict(
        pl_hostname   = 'WASP-43',
        st_logg       = 4.646,
        st_loggerr1   = 0.059,
        st_loggerr2   = -0.044,
        st_metfe      = -0.05,
        st_metfeerr1  = 0.17,
        st_metfeerr2  = -0.17,
        st_dens       = 2.70,
        st_denserr1   = 0.61,
        st_denserr2   = -0.36,
        pl_ratdor     = 5.13,
        pl_ratdorerr1 = 0.36,
        pl_ratdorerr2 = -0.36,
        pl_rvamp      = 550.3,
        pl_rvamperr1  = 6.7,
        pl_rvamperr2  = -6.7
    ),
    dict(
        pl_hostname   = 'WASP-44',
        st_dens       = 1.19 * rho_sun,
        st_denserr1   = 0.32 * rho_sun,
        st_denserr2   = -0.22 * rho_sun,
        st_metfe      = 0.06,
        st_metfeerr1  = 0.1,
        st_metfeerr2  = -0.1,
        pl_rvamp      = 138.8,
        pl_rvamperr1  = 9.0,
        pl_rvamperr2  = -9.0
    ),
    dict(
        pl_hostname   = 'WASP-46',
        st_metfe      = -0.37,
        st_metfeerr1  = 0.13,
        st_metfeerr2  = -0.13,
        pl_rvamp      = 387.0,
        pl_rvamperr1  = 10.0,
        pl_rvamperr2  = -10.0
    ),
    dict(
        pl_hostname   = 'WASP-49',
        st_dens       = 1.0098 * rho_sun,
        st_denserr1   = 0.06 * rho_sun,
        st_denserr2   = -0.06 * rho_sun,
        st_metfe      = -0.23,
        st_metfeerr1  = 0.07,
        st_metfeerr2  = -0.07,
        pl_rvamp      = 56.8,
        pl_rvamperr1  = 2.4,
        pl_rvamperr2  = -2.4
    ),
    dict(
        pl_hostname   = 'WASP-5',
        pl_ratdor     = 5.71,
        pl_ratdorerr1 = 0.34,
        pl_ratdorerr2 = -0.34
    ),
    dict(
        pl_hostname   = 'WASP-50',
        st_logg       = 4.537,
        st_loggerr1   = 0.022,
        st_loggerr2   = -0.022,
        st_metfe      = -0.120,
        st_metfeerr1  = 0.08,
        st_metfeerr2  = -0.08,
        st_dens       = 2.07,
        st_denserr1   = 0.14,
        st_denserr2   = -0.126,
        pl_ratdor     = 7.53,
        pl_ratdorerr1 = 0.35,
        pl_ratdorerr2 = -0.35,
        pl_rvamp      = 256.6,
        pl_rvamperr1  = 4.4,
        pl_rvamperr2  = -4.4
    ),
    dict(
        pl_hostname   = 'WASP-52',
        pl_rvamp      = 84.3,
        pl_rvamperr1  = 3.0,
        pl_rvamperr2  = -3.0
    ),
    dict(
        pl_hostname   = 'WASP-57',
        pl_rvamp      = 100.0,
        pl_rvamperr1  = 7.0,
        pl_rvamperr2  = -7.0
    ),
    dict(
        pl_hostname   = 'WASP-77 A',
        pl_ratdor     = 5.43,
        pl_ratdorerr1 = 0.124,
        pl_ratdorerr2 = -0.124,
        pl_rvamp      = 321.9,
        pl_rvamperr1  = 3.9,
        pl_rvamperr2  = -3.9
    ),
    dict(
        pl_hostname   = 'WASP-81',
        st_vsini      = 1.20,
        st_vsinierr1  = 0.73,
        st_vsinierr2  = -0.73
    ),
    dict(
        pl_hostname   = 'WASP-98',
        pl_rvamp      = 150.0,
        pl_rvamperr1  = 10.0,
        pl_rvamperr2  = -10.0
    ),
    dict(
        pl_hostname   = 'WTS-2',
        pl_ratdor     = 5.25,
        pl_ratdorerr1 = 0.29,
        pl_ratdorerr2 = -0.29
    ),
    dict(
        pl_hostname   = 'XO-2 N',
        pl_orbincl    = 88.8,
        pl_orbinclerr1= 1.2,
        pl_orbinclerr2= -1.2,
        st_dens       = 1.034 * rho_sun,
        st_denserr1   = 0.127 * rho_sun,
        st_denserr2   = -0.058 * rho_sun
    )
)
