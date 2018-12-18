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



    frequency_ooe, power_ooe = LombScargle(time[out_of_eclipse],
                                           flux[out_of_eclipse]).autopower()

    freq_min = 0.06
    freq_max = 0.2

    frequency_range = numpy.logical_and(
        frequency_ooe > freq_min,
        frequency_ooe < freq_max
    )

    new_frequency, new_power = frequency_ooe[frequency_range],power_ooe[frequency_range]

    peak_frequency = new_frequency[numpy.argmax(new_power)]

    spin_period = 1.0/peak_frequency

    print (spin_period)