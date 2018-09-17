class struct :

    def __init__(self,**kwargs):

        for name,value in kwargs.items():
            setattr(self,name,value)



class test:

    def disp(self):

        print (self.set)

    def ini(self):

        for i in range(10):

            self.set.append(i)


    def __init__(self,observables,proposed_step,set = None):


        self.observables = observables
        self.proposed_step = proposed_step
        if set is None:

            self.set = []




    observation_data = struct(
                        Teff = struct( value = 1, sigma = 0.1 ),
                        feh = struct(value=1, sigma=0.1),
                        rvk = struct(value=1, sigma=0.1),
                        inclination = struct(value=1, sigma=0.1)
    )

    fixed_parameters = struct(
                                disk_dissipatoin_age = 5e-3
    )

    proposed_step = struct(
                        Teff_proposed_step = 0.2,
                        feh_proposed_step = 0.2,
                        rvk_proposed_step = 0.3,
                        inclination_step = 0.2
    )




p = [2,3,4,5]

d = []

for i in p:
    print (i)
    t = i**4
    d.append(t)


print (d)




ICS_TEST:


wsun = 2.0 * numpy.pi / 25.34


def create_planet(mass=(constants.M_jup / constants.M_sun).to('')):
    """Return a configured planet to use in the evolution."""

    planet = LockedPlanet(
        mass=mass,
        radius=(constants.R_jup / constants.R_sun).to('')
    )
    return planet


def create_star(interpolator, convective_phase_lag):
    """Create the star to use in the evolution."""

    star = EvolvingStar(mass=1.0,
                        metallicity=0.0,
                        wind_strength=0.17,
                        wind_saturation_frequency=2.45,
                        diff_rot_coupling_timescale=5.0e-3,
                        interpolator=interpolator)
    star.select_interpolation_region(star.core_formation_age())
    star.set_dissipation(zone_index=0,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=convective_phase_lag)
    star.set_dissipation(zone_index=1,
                         tidal_frequency_breaks=None,
                         spin_frequency_breaks=None,
                         tidal_frequency_powers=numpy.array([0.0]),
                         spin_frequency_powers=numpy.array([0.0]),
                         reference_phase_lag=0.0)
    return star


def create_system(star, planet, disk_lock_frequency):
    """Create the system which to evolve from the given star and planet."""

    porb_initial = 3.5
    disk_dissipation_age = 4e-3
    binary = Binary(primary=star,
                    secondary=planet,
                    initial_orbital_period=porb_initial,
                    initial_eccentricity=0.0,
                    initial_inclination=0.0,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    secondary_formation_age=disk_dissipation_age)
    binary.configure(age=star.core_formation_age(),
                     semimajor=float('nan'),
                     eccentricity=float('nan'),
                     spin_angmom=numpy.array([0.0]),
                     inclination=None,
                     periapsis=None,
                     evolution_mode='LOCKED_SURFACE_SPIN')
    planet.configure(age=disk_dissipation_age,
                     companion_mass=star.mass,
                     semimajor=binary.semimajor(porb_initial),
                     eccentricity=0.0,
                     spin_angmom=numpy.array([0.0]),
                     inclination=None,
                     periapsis=None,
                     locked_surface=False,
                     zero_outer_inclination=True,
                     zero_outer_periapsis=True)
    star.detect_stellar_wind_saturation()
    return binary


