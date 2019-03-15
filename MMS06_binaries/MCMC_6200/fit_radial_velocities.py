#!/usr/bin/env python3
"""Fit or the amplitude, offset and shift describing an RV orbit."""

from argparse import ArgumentParser

from matplotlib import pyplot
import scipy
import scipy.optimize

def parse_command_line():
    """Parse the command line to object with attributes."""

    parser = ArgumentParser(
        description='Fit for the radial velocity semi-amplitude from a set of '
        'RV points. Assumes circular orbit and that the x-axis is the orbital '
        'phase.'
    )
    parser.add_argument(
        'rv_fname',
        help='The name of the file containing the RV measurements to fit.'
    )
    return parser.parse_args()

class EccentricOrbit:
    """Provide methods implementing eccentric orbit relations."""

    def __init__(self,
                 *,
                 eccentricity,
                 systemic_rv,
                 argument_of_pericenter,
                 rv_semi_amplitude):
        """Set-up the orbit with a given set of parameters."""

        self.eccentricity = eccentricity
        self.systemic_rv = systemic_rv
        self.argument_of_pericenter = argument_of_pericenter
        self.rv_semi_amplitude = rv_semi_amplitude
        self.target_mean_anomaly = 0.0

    def calc_mean_anomaly(self, eccentric_anomaly):
        """Return the mean anomaly that corresponds to an eccentric anomaly."""

        return (
            eccentric_anomaly
            -
            self.eccentricity * scipy.sin(eccentric_anomaly)
            -
            self.target_mean_anomaly
        )

    def calc_mean_anomaly_deriv(self, eccentric_anomaly):
        """Return the first deriv. of the mean anomaly wrt. eccentric anomaly"""

        return 1.0 + self.eccentricity * scipy.cos(eccentric_anomaly)

    def calc_mean_anomaly_deriv2(self, eccentric_anomaly):
        """Return the second deriv. of the mean anomaly wrt. eccentric anomaly"""

        return -self.eccentricity * scipy.sin(eccentric_anomaly)

    def calc_eccentric_anomaly(self, mean_anomaly):
        """Solve Kepler's eq. for the eccentric anomaly from mean anomaly."""

        self.target_mean_anomaly = mean_anomaly
        root_results = scipy.optimize.root_scalar(f=self.calc_mean_anomaly,
                                                  fprime=self.calc_mean_anomaly_deriv,
                                                  fprime2=self.calc_mean_anomaly_deriv2,
                                                  bracket=(mean_anomaly - 1.0,
                                                           mean_anomaly + 1.0))
        self.target_mean_anomaly = 0.0
        assert root_results.converged
        return root_results.root

    def calc_true_anomaly(self, eccentric_anomaly):
        """Return the true anomaly corresponding to the given eccentric one."""

        return scipy.arctan2(
            (1.0 - self.eccentricity**2)**0.5 * scipy.sin(eccentric_anomaly),
            scipy.cos(eccentric_anomaly) - self.eccentricity
        )

    def calc_radial_velocity(self, mean_anomaly):
        """Return the radial velocity at the given mean anomaly."""

        true_anomaly = self.calc_true_anomaly(
            self.calc_eccentric_anomaly(mean_anomaly)
        )
        return (
            self.rv_semi_amplitude * (scipy.cos(self.argument_of_pericenter
                                                +
                                                true_anomaly)
                                      +
                                      self.eccentricity
                                      *
                                      scipy.cos(self.argument_of_pericenter))
            +
            self.systemic_rv
        )

def residuals(params, rv_data):
    """
    Return the residuals of the RV estimate with the given parameters.

    Args:
        params:    The orbital parameters being fit in the following order:

            * radial velocity semi-amplitude (K)
            * eccentricity
            * argument of pericenter
            * systemic (center of mass) radial velocity

        rv_data:    The available radial velocity measurements.

    Returns:
        scipy array:
            The residuals of the RV function per the given parameters and RV
            data.
    """

    if params[1] < 0 or params[1] > 1:
        return scipy.full(rv_data[:, 1].shape, scipy.inf)
    orbit = EccentricOrbit(rv_semi_amplitude=params[0],
                           eccentricity=params[1],
                           argument_of_pericenter=params[2],
                           systemic_rv=params[3])
    return (
        scipy.vectorize(orbit.calc_radial_velocity)(2.0 * scipy.pi * rv_data[:, 0])
        -
        rv_data[:, 1]
    )

def fit_rv_data(rv_data_fname):
    """Fit or the amplitude, offset and shift describing an RV orbit."""

    rv_data = scipy.genfromtxt(rv_data_fname)
    orbital_params, flag = scipy.optimize.leastsq(residuals,
                                                  [10.0, 0.1, 0.0, 0.0],
                                                  args=(rv_data,),
                                                  maxfev=100000)
    assert(flag in [1, 2, 3, 4])
    return orbital_params, rv_data

if __name__ == '__main__':
    cmdline_args = parse_command_line()
    fit_results, measured_rvs = fit_rv_data(cmdline_args.rv_fname)
    print('Amplitude: ' + repr(fit_results[0]))
    print('eccentricity: ' + repr(fit_results[1]))
    print('arg of pericenter: ' + repr(fit_results[2]))
    print('gamma: ' + repr(fit_results[3]))


    print(fit_results)
    pyplot.plot(measured_rvs[:, 0], measured_rvs[:, 1], 'ok')
    plot_x = scipy.linspace(0.0, 1.0, 100.0)
    best_fit_orbit = EccentricOrbit(rv_semi_amplitude=fit_results[0],
                                    eccentricity=fit_results[1],
                                    argument_of_pericenter=fit_results[2],
                                    systemic_rv=fit_results[3])

    calc_rv = scipy.vectorize(best_fit_orbit.calc_radial_velocity)


    pyplot.plot(plot_x, calc_rv(2.0 * scipy.pi * plot_x), '-', linewidth=3)
    pyplot.show()
