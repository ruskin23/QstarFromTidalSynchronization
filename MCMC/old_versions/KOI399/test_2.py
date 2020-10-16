#!/usr/bin/env python3

from basic_utils import Structure, semimajor
from SamplingUtils import \
    sample_float_with_errors, \
    calculate_logpdf, \
    bounded_normal_normalization
from orbital_evolution.transformations import phase_lag
from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from orbital_evolution.initial_condition_solver import InitialConditionSolver
from orbital_evolution.binary import Binary
from orbital_evolution.dissipating_body import DissipatingBody
from os.path import expanduser, exists, abspath, dirname
from os.path import join as join_paths
import scipy.stats
from astropy import constants, units
import os

import sys

sys.path.append(join_paths(abspath(dirname(__file__)), '..', 'shared'))

# Needs to be after os.path and sys to allow adding the seach path.
# pylint: disable=wrong-import-position
from calculate_planet_mass import calculate_planet_mass


# pylint: enable=wrong-import-position

class MetropolisHastings:
    """Run a Metropolis-Hastings MCMC to get Q* samples."""

    def _minimum_semimajor(self):
        """
        Return the minimum semimajor axis for the planet to survive.
        Must be called after self._update_unrejectable_param and
        self._create_star.
        Args: None
        Returns:
            The semimajor axis below which the planet is assumed destroyed in
            solar radii.
        """

        rroche = (
                2.44
                *
                self.param.planet_radius
                *
                (self.param.star_mass / self.param.planet_mass) ** (-1.0 / 3.0)
        )
        return max(self.star.radius(self.param.age), rroche)

    def _create_planet(self):
        """
        Return a planet to use in the evolution.
        Must be called after _update_unrejectable_param.
        Sets self.planet to an instance of LockedPlanet with the appropriate
        mass and radius.
        Args: None
        Returns: None.
        """

        if hasattr(self, 'planet'): self.planet.delete()
        self.planet = LockedPlanet(mass=self.param.planet_mass,
                                   radius=self.param.planet_radius)

    def _create_star(self):
        """
        Create the star to use in the evolution.
        Must be called after _update_unrejectable_param.
        Sets self.star to an instance of EvolvingStar without defining its
        dissipation
        Args: None
        Returns: None
        """

        if hasattr(self, 'star'): self.star.delete()
        self.star = EvolvingStar(
            mass=self.param.star_mass,
            metallicity=self.param.feh,
            wind_strength=self.param.wind_strength,
            wind_saturation_frequency=self.param.wind_saturation_frequency,
            diff_rot_coupling_timescale=(
                self.param.coreenv_coupling_timescale
            ),
            interpolator=self.interpolator
        )
        self.star.select_interpolation_region(self.star.core_formation_age())

    def _set_stellar_dissipation(self, lgQ):
        """Set the dissipation of self.star per the given lg(Q*')."""

        if (self._dissipation_powerlaw):
            lgQ_min = self.lgQ_range[0]
            tidal_frequency_breaks = scipy.array([
                self._nominal_wtide
                *
                10.0 ** (
                        (lgQ - lgQ_min)
                        /
                        self._dissipation_powerlaw
                )
            ])
            message = (
                    'Assuming Q* = '
                    +
                    repr(10.0 ** lgQ_min)
                    +
                    ' * (wtide / '
                    +
                    repr(tidal_frequency_breaks[0])
                    +
                    ')^'
                    +
                    repr(-self._dissipation_powerlaw)
                    +
                    ' for wtide '
            )
            if (self._dissipation_powerlaw > 0):
                tidal_frequency_powers = scipy.array([
                    self._dissipation_powerlaw,
                    0.0
                ])
                message += '< '
            else:
                tidal_frequency_powers = scipy.array([
                    0.0,
                    self._dissipation_powerlaw
                ])
                message += '> '

            reference_phase_lag = phase_lag(lgQ_min)

            message += repr(tidal_frequency_breaks[0])
            print(message)

        else:
            tidal_frequency_breaks = None
            tidal_frequency_powers = scipy.array([0.0])
            reference_phase_lag = phase_lag(lgQ)
        self.star.set_dissipation(
            zone_index=0,
            tidal_frequency_breaks=tidal_frequency_breaks,
            spin_frequency_breaks=None,
            tidal_frequency_powers=tidal_frequency_powers,
            spin_frequency_powers=scipy.array([0.0]),
            reference_phase_lag=reference_phase_lag
        )
        self.star.set_dissipation(
            zone_index=1,
            tidal_frequency_breaks=None,
            spin_frequency_breaks=None,
            tidal_frequency_powers=scipy.array([0.0]),
            spin_frequency_powers=scipy.array([0.0]),
            reference_phase_lag=0.0
        )

    def _evolve_binary(self,
                       disk_lock_frequency,
                       disk_dissipation_age,
                       lgQ):
        """
        Create and evolve the system with the given star, planet and orbit.
        Must be called after _create_star and _create_planet.
        Sets self.binary to an instance of Binary with an evolution
        reproducinng the specified configuration.
        Args:
            - disk_lock_frequency:
                The proposed frequency at which the surface spin of the star
                is initially locked.
            - disk_dissipation_age:
                The age at which the disk disspates, releasing the surface
                spin of the star.
            - lgQ:
                The proposed value for log10 of the tidal quality factor of
                the star.
        Returns:
            - initial_orbital_period:
                The orbital period at which the planet first appears in
                order to reproduce the specified final orbital period in
                days.
            - final_stellar_period:
                The spin period of the stellar surface in days for the
                evolution which reproduces the specified final orbital
                period.
        """

        self._set_stellar_dissipation(lgQ)
        initial_orbital_period, final_stellar_period = (
            self.find_initial_period(
                star=self.star,
                planet=self.planet,
                target=Structure(
                    age=self.param.age,
                    Porb=self.param.final_orbital_period,
                    Wdisk=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age
                )
            )
        )
        self.binary = self.find_initial_period.binary
        return initial_orbital_period, final_stellar_period

    def _draw_unrejectable_observables(self):
        """Select Teff*, rho*, [Fe/H] and RV K per their distributions."""

        mass_age_solutions = []
        while (
                len(mass_age_solutions) == 0
                or
                self.param.age > 13.7
                or
                self.param.age < 0.1
        ):
            proposed_values = dict()
            for quantity_name, quantity_bounds in \
                    self.unrejectable_bounds.items():
                observable = getattr(self.observables,
                                     quantity_name)
                step_size = getattr(self.step_size, quantity_name)
                current_value = getattr(self.param, quantity_name)
                proposed_values[quantity_name] = sample_float_with_errors(
                    observable,
                    current_value,
                    step_size,
                    quantity_bounds
                )

            print('Looking for star with [Fe/H] = %(feh)s, Teff = %(teff)s, '
                  'rho = %(rho_star)s'
                  %
                  proposed_values)

            mass_age_solutions = self.interpolator.change_variables(
                feh=proposed_values['feh'],
                teff=proposed_values['teff'],
                rho=proposed_values['rho_star']
            )
            print('Got: ' + repr(mass_age_solutions))
            if (len(mass_age_solutions) > 0):
                self.param.star_mass, self.param.age = mass_age_solutions[
                    scipy.stats.randint.rvs(0, len(mass_age_solutions))
                ]

        for quantity_name, quantity_value in proposed_values.items():
            setattr(self.param, quantity_name, quantity_value)

    def _update_unrejectable_param(self):
        """
        Select star mass, age, [Fe/H] & planet mass for next iteration.
        Creates self.star.
        Args: None
        Returns:
            - min_semimajor:
                The smallest value the semimajor axis at the present age is
                allowed to be.
        """

        self._draw_unrejectable_observables()
        self._create_star()
        min_semimajor = 1.0
        final_semimajor = 0.0
        while final_semimajor <= min_semimajor:
            self.param.final_orbital_period = sample_float_with_errors(
                self.observables.porb,
                scipy.nan,
                None
            )
            self.param.planet_mass = calculate_planet_mass(self.param)
            min_semimajor = self._minimum_semimajor()
            final_semimajor = semimajor(self.param.star_mass,
                                        self.param.planet_mass,
                                        self.param.final_orbital_period)
        return min_semimajor

    def _propose_rejectable_param(self):
        """
        Return the next proposed final porb, disk lock frequency and lg(Q*).
        Must be called after self._update_unrejectable_param and
        self._create_star.
        Args: None
        Returns:
            - disk_lock_frequency:
                The value of the disk locking angular velocity proposed at
                the next step.
            - lgQ:
                The value of lg(Q*') proposed at the next step.
        """

        min_disk_lock_frequency, max_disk_lock_frequency = (
            self.disk_lock_frequency_range
        )
        disk_lock_frequency = min_disk_lock_frequency - 1.0
        while (
                disk_lock_frequency < min_disk_lock_frequency
                or
                disk_lock_frequency > max_disk_lock_frequency
        ):
            disk_lock_frequency = scipy.stats.norm.rvs(
                self.param.disk_lock_frequency,
                self.step_size.disk_lock_frequency
            )
        disk_dissipation_age = (
                self.param.disk_lifetime_param[0]
                *
                disk_lock_frequency ** self.param.disk_lifetime_param[1]
        )

        min_lgQ, max_lgQ = self.lgQ_range
        lgQ = min_lgQ - 1.0
        while lgQ < min_lgQ or lgQ > max_lgQ:
            lgQ = scipy.stats.norm.rvs(
                self.param.lgQ,
                self.step_size.lgQ
            )

        return disk_lock_frequency, disk_dissipation_age, lgQ

    def _step_logprobability_ratio(self,
                                   disk_lock_frequency,
                                   lgQ):
        """
        Return the ratio of probabilities of proposing a step back to forth.
        Args:
            The parameters returned by _propose_initial_conditions_lgQ in the
            same exact order.
        Returns:
            The ratio of the probability to take the proposed step in reverse
            to the probability of proposing in the forward direction.
        """

        return scipy.log(
            bounded_normal_normalization(
                disk_lock_frequency,
                self.step_size.disk_lock_frequency,
                *self.disk_lock_frequency_range
            )
            /
            bounded_normal_normalization(
                self.param.disk_lock_frequency,
                self.step_size.disk_lock_frequency,
                *self.disk_lock_frequency_range
            )
            *
            bounded_normal_normalization(
                lgQ,
                self.step_size.lgQ,
                *self.lgQ_range
            )
            /
            bounded_normal_normalization(
                self.param.lgQ,
                self.step_size.lgQ,
                *self.lgQ_range
            )
        )

    def _accept_last_evolution(self,
                               disk_lock_frequency,
                               disk_dissipation_age,
                               lgQ,
                               iteration):
        """
        Accept the latest evolution with the appropriate probability.
        Args:
            The parameters returned by _propose_initial_conditions_lgQ in the
            same exact order.
            - iteration:
                How many steps were previously accepted.
        Returns:
            True iff the latest proposed step should be accepted.
        """

        self.final_state = self.binary.final_state()
        assert (self.final_state.age == self.param.age)
        self.final_state.porb = (
            self.binary.orbital_period(self.final_state.semimajor)
        )
        self.final_state.wstar = (
                self.final_state.envelope_angmom
                /
                self.star.envelope_inertia(self.final_state.age)
        )
        self.final_state.star_vsini = (
                units.Quantity(self.final_state.wstar, 1.0 / units.day)
                *
                units.Quantity(self.star.radius(self.param.age), units.R_sun)
        ).to(units.km / units.s).value
        self.final_state.pstar = 2.0 * scipy.pi / self.final_state.wstar

        if scipy.isnan(self.final_state.porb):
            observables_logprobability = -float('inf')
        else:
            if hasattr(self.observables, 'pstar'):
                observables_logprobability = calculate_logpdf(
                    self.final_state.pstar,
                    self.observables.pstar
                )
            else:
                observables_logprobability = calculate_logpdf(
                    self.final_state.star_vsini,
                    self.observables.star_vsini
                )
        if not hasattr(self.param, 'observables_logprobability'):
            if scipy.isnan(self.final_state.porb):
                proposal_logprobability_ratio = -scipy.inf
                acceptance_probability = 0.0
            else:
                proposal_logprobability_ratio = 0.0
                acceptance_probability = 1.0
        else:
            proposal_logprobability_ratio = self._step_logprobability_ratio(
                disk_lock_frequency,
                lgQ
            )
            if (
                    not scipy.isfinite(observables_logprobability)
                    or
                    self.final_state.porb < (
                    self.observables.porb
                    -
                    max(self.observables.porb.minus_error,
                        100.0 * self.precision)
            )
                    or
                    self.final_state.porb > (
                    self.observables.porb
                    +
                    max(self.observables.porb.plus_error,
                        100.0 * self.precision)
            )
            ):
                acceptance_probability = 0.0
            elif (
                    not scipy.isfinite(self.param.observables_logprobability)
                    or
                    (
                            observables_logprobability
                            +
                            proposal_logprobability_ratio
                    ) > self.param.observables_logprobability
            ):
                acceptance_probability = 1.0
            else:
                acceptance_probability = scipy.exp(
                    observables_logprobability
                    -
                    self.param.observables_logprobability
                    +
                    proposal_logprobability_ratio
                )
        print('Final state:\n' + self.final_state.format())
        print('Acceptance probability: ' + repr(acceptance_probability))
        if scipy.stats.uniform.rvs() <= acceptance_probability:
            self.param.observables_logprobability = (
                observables_logprobability
            )
            self.param.proposal_logprobability_ratio = (
                proposal_logprobability_ratio
            )
            return True
        else:
            return False

    def _add_rejected_step(self, **substitute_values):
        """
        Add another rejected step to the rejected output file.

        Keyword only arguments:
            Specify value to use instead of self.param members with identical
            names.
        Returns: None
        """

        def iterate_output_param(get_value):
            for trace_source, trace_type in [(self.param, 'in'),
                                             (self.final_state, 'out')]:
                for param_name in dir(trace_source):
                    if (
                            param_name != 'format'
                            and
                            param_name != 'disk_lifetime_param'
                            and
                            param_name[0] != '_'
                    ):
                        trace_param_name = '%s(%s)' % (trace_type,
                                                       param_name)
                        if get_value:
                            if (
                                    trace_type == 'in'
                                    and param_name in substitute_values
                            ):
                                value = substitute_values[param_name]
                            else:
                                value = getattr(trace_source, param_name)
                            yield (trace_param_name, value)
                        else:
                            yield trace_param_name

        if not self._wrote_rejected_header:
            column_number = 0
            self.rejected_file_columns = dict()
            for param_name in iterate_output_param(False):
                self.rejected_file.write('%35s ' % param_name)
                self.rejected_file_columns[param_name] = column_number
                column_number += 1
            self.rejected_file.write('\n')
            self._wrote_rejected_header = True
        rejected_file_data = scipy.empty(len(self.rejected_file_columns))

        for param_name, param_value in iterate_output_param(True):
            rejected_file_data[
                self.rejected_file_columns[param_name]
            ] = param_value

        self.rejected_file.write(' '.join(['%35s' % repr(value)
                                           for value in rejected_file_data])
                                 +
                                 '\n')
        self.rejected_file.flush()

    def _step(self, iteration):
        """Propose MCMC steps until one is accepted and update self.param."""

        accepted = False
        while not accepted:
            min_semimajor = self._update_unrejectable_param()
            disk_lock_frequency, disk_dissipation_age, lgQ = (
                self._propose_rejectable_param()
            )
            self._create_planet()
            print(
                    'Running evolution from\n'
                    +
                    self.param.format()
                    +
                    'wdisk = %s, lgQ = %s, a_min = %s'
                    %
                    (disk_lock_frequency, lgQ, min_semimajor)
            )
            try:
                initial_orbital_period, final_stellar_period = (
                    self._evolve_binary(disk_lock_frequency,
                                        disk_dissipation_age,
                                        lgQ)
                )
            except:
                continue
            accepted = self._accept_last_evolution(disk_lock_frequency,
                                                   disk_dissipation_age,
                                                   lgQ,
                                                   iteration)
            print('Accepted' if accepted else 'Rejected')
            initial_semimajor = semimajor(self.param.star_mass,
                                          self.param.planet_mass,
                                          initial_orbital_period)
            if not accepted:
                self._add_rejected_step(
                    initial_orbital_period=initial_orbital_period,
                    initial_semimajor=initial_semimajor,
                    disk_lock_frequency=disk_lock_frequency,
                    disk_dissipation_age=disk_dissipation_age,
                    lgQ=lgQ,
                    min_semimajor=min_semimajor
                )
                self.trace_file.write(self.last_trace_line)
                self.trace_file.flush()
        self.param.initial_semimajor = initial_semimajor
        self.param.disk_lock_frequency = disk_lock_frequency
        self.param.disk_dissipation_age = disk_dissipation_age
        self.param.lgQ = lgQ
        self.param.min_semimajor = min_semimajor

    def _add_successful_step(self, iteration):
        """Add another successful step to the trace and output file."""

        def iterate_output_param(get_value):
            """Generator that iterates over all parameters to output."""

            for trace_source, trace_type in [(self.param, 'in'),
                                             (self.final_state, 'out')]:
                for param_name in dir(trace_source):
                    if (
                            param_name != 'format'
                            and
                            param_name != 'disk_lifetime_param'
                            and
                            param_name[0] != '_'
                    ):
                        trace_param_name = '%s(%s)' % (trace_type,
                                                       param_name)
                        if get_value:
                            yield (trace_param_name,
                                   getattr(trace_source, param_name))
                        else:
                            yield trace_param_name

        if not hasattr(self, 'trace'):
            assert (iteration == 0)
            self.trace = Structure()
            self.trace_file_columns = dict()
            column_number = 0
            for param_name in iterate_output_param(False):
                setattr(self.trace, param_name, [])
                self.trace_file.write('%35s ' % param_name)
                self.trace_file_columns[param_name] = column_number
                column_number += 1
            self.trace_file.write('\n')
        trace_file_data = scipy.empty(len(self.trace_file_columns))
        print('Iteration ' + repr(iteration))
        for param_name, param_value in iterate_output_param(True):
            destination = getattr(self.trace, param_name)
            assert (len(destination) == iteration)
            destination.append(param_value)
            trace_file_data[
                self.trace_file_columns[param_name]
            ] = param_value

        self.last_trace_line = (' '.join(['%35s' % repr(value)
                                          for value in trace_file_data])
                                +
                                '\n')
        self.trace_file.write(self.last_trace_line)
        self.trace_file.flush()

    def _sample(self, starting_iteration, num_iterations):
        """Perform the sampling, starting from the current state."""

        iteration = starting_iteration
        while iteration < num_iterations:
            self._step(iteration)
            self._add_successful_step(iteration)
            iteration += 1
        return trace

    def __init__(self,
                 stellar_evolution_interpolator,
                 fixed_parameters,
                 observables,
                 step_size,
                 disk_lock_frequency_range,
                 lgQ_range,
                 max_timestep,
                 precision,
                 dissipation_powerlaw):
        """
        Set up an MCMC run for a given observed system.
        Args:
            - stellar_evolution_interpolator:
                An instance of stellar_evolution.VarChangingInterpolator
                defining the stellar evolution.
            - fixed_parameters:
                Objet with at least the following attributes:
                    - planet_formation_age:
                        The age at which the planet forms in Gyr.
                    - planet_radius:
                        The radius of the planet in Jovian radii.
                    - coreenv_coupling_timescale:
                        The timescale on which the stellar core and envelope
                        couple in Gyr.
                    - wind_strength:
                        The constant defining the efficiency of angular
                        momentum loss due to wind.
                    - wind_saturation_frequency:
                        The spin angular velocity at which the angular
                        momentum loss due to stellar wind saturates.
                    - disk_lifetime_param:
                        Scaling and powerlaw index for the depnedence of the
                        age at which the stellar surface spin is released
                        from its lock to a disk in Gyr on the angulare
                        velocity at which the star is locked during that
                        time.
            - observables:
                Object with at least the following attributes, each of which
                should have further attributes called value (the measured
                value for the corresponding observable) and uncertainty (the
                1-sigma uncertainty of the corresponding observable):
                    - rho_star:
                        The density of the star.
                    - feh:
                        The [Fe/H] of the star.
                    - teff:
                        The stellar effective temperature in K.
                    - rvk:
                        The radial velocity semi-amplitude.
                    - porb:
                        The orbital period of the planet.
                    - orbital_inclination:
                        The angel between the orbital angular momentum and
                        the line of sight in degrees.
                    Exactly one of the following two attributes must be
                    defined:
                    - pstar:
                        The spin period of the stellar surface in days.
                    - star_vsini:
                        The angular velocity of the stellar surface in
                        rad/day.
            - step_size:
                The standard deviation of the step along each variable. If
                None, the variable is sampled directly according to its
                observable distribution. The following attributes must exist:
                    - disk_lock_frequency:
                        The size of the step for the disk locking angular
                        velocity.
                    - lgQ:
                        The size of the step for lg(Q'*).
                    - feh:
                        The size of the step for the stellar [Fe/H].
                    - teff:
                        The size of the step for the stellar effective
                        temperature.
                    - rho_star:
                        The size of the step for the stellar density.
                    - rvk:
                        The size of the step for the radial velocity
                        semi-amplitude.
                    - orbital_inclination:
                        The size of the step for the orbital inclination.
            - disk_lock_frequency_range:
                The angular velocity at which the stellar surface rotation is
                at first locked in rad/day.
            - lgQ_range:
                The range of lgQ values to try.
            - max_timestep:
                The maximum time step the orbital evolution is allowed to
                make.
            - precision:
                The required presicion of the orbital evolution.
            - dissipation_powerlaw:
                The powerlaw index of the phase lag scaling with tidal
                frequency. If not zero, the MCMC chain samples only Q* at the
                nominal tidal frequency.
        Returns: None
        """

        self.interpolator = stellar_evolution_interpolator

        self.find_initial_period = InitialConditionSolver(
            planet_formation_age=fixed_parameters.planet_formation_age,
            evolution_max_time_step=max_timestep,
            evolution_precision=precision
        )

        self.param = Structure()
        for quantity in ['planet_formation_age',
                         'planet_radius',
                         'coreenv_coupling_timescale',
                         'wind_strength',
                         'wind_saturation_frequency',
                         'disk_lifetime_param']:
            setattr(self.param,
                    quantity,
                    getattr(fixed_parameters, quantity))

        self.observables = observables
        self.step_size = step_size
        self.disk_lock_frequency_range = disk_lock_frequency_range
        self.lgQ_range = lgQ_range
        self.max_timestep = max_timestep
        self.precision = precision

        self.unrejectable_bounds = dict(
            feh=(float(self.interpolator.track_feh[0]),
                 float(self.interpolator.track_feh[-1])),
            rho_star=(-scipy.inf, scipy.inf),
            teff=(-scipy.inf, scipy.inf),
            rvk=(-scipy.inf, scipy.inf),
            orbital_inclination=(-scipy.inf, scipy.inf)
        )

        self._nominal_wtide = 4.0 * scipy.pi * (1.0 / observables.porb
                                                -
                                                1.0 / observables.pstar)
        print('Nominal tidal frequency: ' + repr(self._nominal_wtide))
        self._dissipation_powerlaw = dissipation_powerlaw

    def start(self,
              num_iterations,
              trace_file,
              rejected_file=None):
        """
        Run the MCMC for the given number of iterations.

        Args:
            - num_iterations:
                Propose steps until this many are accepted.
            - trace_file:
                The filename under which to store accepted steps. Must not
                exist.
            - rejected_file:
                The filename under which to store rejected steps. If None,
                the rejected steps are not stored. Must not exist.
        Returns:
            None, but the object receives a trace attribute containing the
            MCMC trace.
        """

        assert (not exists(trace_file))
        self.trace_file = open(trace_file, 'w')
        if rejected_file is not None:
            assert (not exists(rejected_file))
            self.rejected_file = open(rejected_file, 'w')
            self._wrote_rejected_header = False

        self.param.disk_lock_frequency = scipy.stats.uniform.rvs(
            self.disk_lock_frequency_range[0],
            (
                    self.disk_lock_frequency_range[1]
                    -
                    self.disk_lock_frequency_range[0]
            )
        )
        self.param.lgQ = scipy.stats.uniform.rvs(
            self.lgQ_range[0],
            self.lgQ_range[1] - self.lgQ_range[0]
        )
        stellar_age = None
        while (
                stellar_age is None
                or
                stellar_age > 13.7
                or
                stellar_age < 0.1
        ):
            for param_name in ['feh',
                               'teff',
                               'rho_star',
                               'rvk',
                               'orbital_inclination']:
                setattr(
                    self.param,
                    param_name,
                    sample_float_with_errors(
                        getattr(self.observables, param_name),
                        scipy.nan,
                        None
                    )
                )
            mass_age_solutions = self.interpolator.change_variables(
                feh=self.param.feh,
                teff=self.param.teff,
                rho=self.param.rho_star
            )
            if len(mass_age_solutions) > 0:
                stellar_age = mass_age_solutions[
                    scipy.stats.randint.rvs(0, len(mass_age_solutions))
                ][1]

        self._update_unrejectable_param()

        self._sample(0, num_iterations)

    def resume(self,
               num_iterations,
               trace_file,
               rejected_file=None):
        """
        Continue a chain from the last accepted step.
        Args: same as start, but the trace file must exist and not be empty.
        Returns:
            None, but the object receives a trace attribute containing the
            MCMC trace, including the steps contained in the input trace
            file.
        """

        def get_file_columns(fname):
            """
            Return the list of columns in the file in the correct order.
            """

            file_columns = []
            with open(fname, 'r') as file_content:
                for line in file_content:
                    file_columns = line.split()
                    if file_columns != []: break
            return {
                c: i for i, c in enumerate(file_columns)
            }

        trace_file_columns = get_file_columns(trace_file)
        if trace_file_columns == {}:
            print('Removing ' + repr(trace_file))
            os.remove(trace_file)
            if rejected_file is not None and exists(rejected_file):
                print('Removing ' + repr(rejected_file))
                os.remove(rejected_file)
            self.start(num_iterations, trace_file, rejected_file)
        self.trace_file_columns = trace_file_columns
        assert (bool(self.trace_file_columns))
        old_trace = scipy.genfromtxt(trace_file,
                                     names=True,
                                     deletechars='')
        param_name_list = [field_name for field_name, field_type
                           in old_trace.dtype.fields.items()]
        assert (set(param_name_list) == set(self.trace_file_columns.keys()))
        if old_trace.size > 1:
            self.trace = Structure(**{param_name: list(old_trace[param_name])
                                      for param_name in param_name_list})
        else:
            self.trace = Structure(
                **{param_name: [float(old_trace[param_name])]
                   for param_name in param_name_list}
            )
        old_iterations = old_trace.size
        self.final_state = Structure()
        for param_name in param_name_list:
            assert (param_name[-1] == ')')
            attribute = (param_name[3:-1] if param_name.startswith('in(')
                         else param_name[4:-1])
            if attribute == 'disk_lifetime_param': continue
            if param_name.startswith('in('):
                setattr(self.param, attribute,
                        getattr(self.trace, param_name)[-1])
            else:
                assert (param_name.startswith('out('))
                setattr(self.final_state,
                        attribute,
                        getattr(self.trace, param_name)[-1])

        self.trace_file = open(trace_file, 'a')
        self.trace_file.write('\n')
        if rejected_file is not None:
            self._wrote_rejected_header = exists(rejected_file)
            if self._wrote_rejected_header:
                self.rejected_file_columns = get_file_columns(rejected_file)
                if not self.rejected_file_columns:
                    self._wrote_rejected_header = False
            self.rejected_file = open(rejected_file, 'a')
            self.rejected_file.write('\n')

        self._sample(old_iterations, num_iterations)


if __name__ == '__main__':
    from stellar_evolution.manager import StellarEvolutionManager
    from orbital_evolution.evolve_interface import library as \
        orbital_evolution_library

    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        b"eccentricity_expansion_coef.txt"
    )

    serialized_dir = expanduser(
        '~/projects/git/tidal-orbital-evolution/'
        'stellar_evolution_interpolators'
    )
    manager = StellarEvolutionManager(serialized_dir)
    interpolator = manager.get_interpolator_by_name('default')

    mcmc = MetropolisHastings(
        stellar_evolution_interpolator=interpolator,
        fixed_parameters=Structure(
            planet_formation_age=5e-3,
            planet_radius=constants.R_jup.to(units.R_sun).value,
            coreenv_coupling_timescale=5e-3,
            wind_strength=0.17,
            wind_saturation_frequency=2.54,
            disk_lifetime=5e-3
        ),
        observables=Structure(
            rho_star=Structure(value=1.41,
                               uncertainty=0.1),
            feh=Structure(value=0.0,
                          uncertainty=0.1),
            teff=Structure(value=5777.0,
                           uncertainty=50.0),
            rvk=Structure(
                value=415,
                uncertainty=10.0
            ),
            porb=Structure(
                value=1.0,
                uncertainty=1e-5
            ),
            pstar=Structure(
                value=5.0,
                uncertainty=0.5
            )
        ),
        step_size=Structure(
            disk_lock_frequency=0.5 * MetropolisHastings._wsun,
            lgQ=0.5
        ),
    )
    print(mcmc.run(scipy.inf).format())