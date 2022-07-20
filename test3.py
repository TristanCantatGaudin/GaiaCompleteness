from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from gu import load_data, selectionFunctionRADEC


l = np.arange(-180, 180, 1)
b = np.arange(-90, 90, 1)
l, b = np.meshgrid(l, b)
c = SkyCoord(l=l * u.degree,
			 b=b * u.degree,
			 frame='galactic')

dfm10 = load_data()

#the map is in (ra,dec) so we need to convert:
ra = c.icrs.ra.degree.flatten()
dec = c.icrs.dec.degree.flatten()
completeness21 = selectionFunctionRADEC(ra, dec, 21., dfm10)


#and plot
import matplotlib.pyplot as plt
plt.figure(1)
plt.scatter(l.flatten(), b.flatten(), c=selectionFunctionRADEC(ra, dec, 21., dfm10), vmin=0, vmax=1, s=1)
plt.colorbar(label='completeness at G=21')
plt.xlim(max(l.flatten()), min(l.flatten()))
plt.ylim(min(b.flatten()), max(b.flatten()))
plt.xlabel('l (degrees)')
plt.ylabel('b (degrees)')
plt.show()