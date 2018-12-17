from matplotlib import pyplot
import numpy

from astropy.io import fits
from astropy.stats import LombScargle

if __name__ == '__main__':
    eb_lc = fits.open('kepler_EB_test1.fits')
    time = eb_lc[1].data['TIME']
    flux = eb_lc[1].data['LC_DETREND']
    period = 2.0 * eb_lc[1].header['TPERIOD']
    in_eclipse = numpy.logical_or(
        numpy.logical_and(
            time % period > 0.92,
            time % period < 1.24,
        ),
        numpy.logical_and(
            time % period > 3.13,
            time % period < 3.45,
        )
    )
    out_of_eclipse = numpy.logical_and(
        numpy.isfinite(flux),
        numpy.logical_not(in_eclipse)
    )
    non_nan = numpy.isfinite(flux)

    frequency_all, power_all = LombScargle(time[non_nan],
                                           flux[non_nan]).autopower()
#    pyplot.plot(frequency_all, power_all, '-g')

    frequency_ooe, power_ooe = LombScargle(time[out_of_eclipse],
                                           flux[out_of_eclipse]).autopower()
    pyplot.semilogx(1.0 / frequency_ooe, power_ooe, '-r')
    pyplot.axvline(eb_lc[1].header['TPERIOD'])
    pyplot.show()

    pyplot.plot(time[out_of_eclipse] % period, flux[out_of_eclipse], '.g')
    pyplot.show()
