from matplotlib import pyplot
import numpy

from astropy.io import fits
from astropy.stats import LombScargle

if __name__ == '__main__':
    eb_lc = fits.open('KIC7289157.fits')
    time = eb_lc[1].data['TIME'] - 1
    flux = eb_lc[1].data['PDCSAP_FLUX']
    period = 5.266273

    in_eclipse = numpy.logical_or(
        numpy.logical_and(
            time % period > 1.45,
            time % period < 1.70,
        ),
        numpy.logical_and(
            time % period > 4.15,
            time % period < 4.44,
        )
    )
    out_of_eclipse = numpy.logical_and(
        numpy.isfinite(flux),
        numpy.logical_not(in_eclipse)
    )
    non_nan = numpy.isfinite(flux)

    pyplot.plot(time[out_of_eclipse],
                flux[out_of_eclipse],
                '.g')

    pyplot.plot(time[in_eclipse],
                flux[in_eclipse],
                '.r')
    pyplot.show()

    frequency_all, power_all = LombScargle(time[non_nan],
                                           flux[non_nan]).autopower()
    pyplot.semilogx(1.0 / frequency_all, power_all, '-g')

    frequency_ooe, power_ooe = LombScargle(time[out_of_eclipse],
                                           flux[out_of_eclipse]).autopower()
    pyplot.semilogx(1.0 / frequency_ooe, power_ooe, '-r')
    pyplot.axvline(period)
    pyplot.show()

    period = 7.75

    num_in_chunk = out_of_eclipse.sum() // 5
    for chunk, color in zip(range(5), 'rgbcm'):
        selected = numpy.copy(out_of_eclipse)
        selected[0: num_in_chunk * chunk] = False
        selected[num_in_chunk * (chunk + 1):] = False
        pyplot.plot(time[selected] % period,
                    flux[selected],
                    '.' + color)
    pyplot.show()
