from matplotlib import pyplot
import numpy
import argparse
from astropy.io import fits
from astropy.stats import LombScargle

if __name__ == '__main__':

    parser=argparse.ArgumentParser()
    parser.add_argument('fname',help='fits file to display light curve')

    args=parser.parse_args()


    eb_lc = fits.open(args.fname)
    time = eb_lc[1].data['TIME']
    flux = eb_lc[1].data['PDCSAP_FLUX']
    period = 8.215

    #time=time%period
    t=time%period
    pyplot.plot(t,flux)
    pyplot.show()


    time = numpy.array(time)

    in_eclipse = numpy.logical_or(
        numpy.logical_and(
            time % period > 3.21,
            time % period < 3.72,
        ),
        numpy.logical_and(
            time % period > 6.27,
            time % period < 6.70,
        )
    )


    #for external events
    in_event1 = numpy.logical_and(time>1283.89, time<1285.09)
    in_event2 = numpy.logical_and(time>1297.63, time<1298.84)
    in_event3 = numpy.logical_and(time>1302.45, time<1303.54)
    in_event4 = numpy.logical_and(time>1360.56, time<1361.64)

    in_event = numpy.logical_or.reduce((in_event1,in_event2,in_event3,in_event4))

    not_required = numpy.logical_or(in_eclipse, in_event)
    #not_required = in_eclipse

    out_of_eclipse = numpy.logical_and(
        numpy.isfinite(flux),
        numpy.logical_not(not_required)
    )

    non_nan = numpy.isfinite(flux)

    pyplot.plot(time[out_of_eclipse],
                flux[out_of_eclipse],
                '.g')

    pyplot.plot(time[not_required],
                flux[not_required],
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



    #period = 7.75

    #num_in_chunk = out_of_eclipse.sum() // 5
    #for chunk, color in zip(range(5), 'rgbcm'):
    #    selected = numpy.copy(out_of_eclipse)
    #    selected[0: num_in_chunk * chunk] = False
    #    selected[num_in_chunk * (chunk + 1):] = False
    #    pyplot.plot(time[selected] % period,
    #                flux[selected],
    #                '.' + color)
    #pyplot.show()


    #checking or 1sigma error

    check= numpy.logical_and(1.0/frequency_ooe<8.7,1.0/frequency_ooe>5.7)
    new_freq =frequency_ooe[check]
    new_power = power_ooe[check]

    pyplot.semilogx(1.0/new_freq,new_power, '-r' )

    mean = 1.0/(new_freq[numpy.argmax(new_power)])

    print(mean)
    #for index,value in enumerate(new_power):
    #    print(index)
    #    print(value)
#
#
    print(new_freq[6])
    pyplot.show()
