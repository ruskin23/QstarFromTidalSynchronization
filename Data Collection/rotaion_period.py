
from astropy.io import fits
from matplotlib import pyplot


image_file = 'kepler_EB_test1.fits'


fits.info(image_file)
image_data = fits.getdata(image_file, ext=0)


plt.figure()
plt.imshow(image_data, cmap='gray')
