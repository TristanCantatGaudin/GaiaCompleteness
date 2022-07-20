import numpy as np
import pandas as pd
import healpy as hp


def load_data(mapFile: str = 'allsky_m10_hpx7.h5') -> pd.DataFrame:
    """
    Load data from a map file
    """
    df = pd.read_hdf(mapFile, "data")
    print(f'File {mapFile:s} loaded.')
    return df


def get_model_parameters():
    """ Returns the selection function model parameters that we obtained through hierarchical fitting."""
    return dict(zip(
            ('a', 'b', 'c', 'd', 'e', 'f', 'x', 'y', 'z', 'lim'),
            (1.0154179774831278, -0.008254847738351057, 0.6981959151433699, -0.07503539255843136,
            1.7491113533977052, 0.4541796235976577, -0.06817682843336803, 1.5712714454917935,
            -0.12236281756914291, 20.53035927443456)
            ))


def sigmoid(G: np.ndarray,
            G0: np.ndarray,
            invslope: float,
            shape: float) -> np.ndarray:
    """ Generalized sigmoid function

    Parameters
    ----------
    G: nd.array
        where to evaluate the function
    G0: float
        inflection point
    invslope: float
        steepness of the linear part. Shallower for larger values
    shape: float
        if shape=1, model is the classical logistic function,
        shape converges to zero, then the model is a Gompertz function.

    Returns
    -------
    f(G) evaluation of the model.

	FIXME: this function is not robust to numerical issues (but works within the range of values we feed it)
    """
    delta = G - G0
    return 1 - (0.5 * (np.tanh(delta / invslope) + 1)) ** shape



	
	
	
def selectionFunction(G, m10):
	""" Predicts the completeness at magnitude G, given a value of M_10 read from a precomputed map.

	Parameters
	----------
	G:   float or nd.array
			where to evaluate the function
	m10: float or nd.array
			the value of M_10 in a given region

	Returns
	-------
	sf(G) between 0 and 1.
	The shape of the output will match the input:
		if given an array (i.e. an array of positions) the output is an array
		if given an array of Gmag and either one position or a matching array of positions, the output is also an array
		if only given scalars, the output is one number.
	
	"""
	#These are the best-fit value of the free parameters we optimised in our model:
	a, b, c, d, e, f, x, y, z, lim = get_model_parameters().values()
	
	if isinstance(m10, float):
		m10 = np.array([m10])
		m10wasFloat = True
	else:
		m10 = np.array(m10)
		m10wasFloat = False


	predictedG0 = a * m10 + b
	predictedG0[m10>lim] = c*m10[m10>lim] + (a-c) * lim + b

	predictedInvslope = x * m10 + y
	predictedInvslope[m10>lim] = z * m10[m10>lim] + (x - z) * lim + y

	predictedShape = d * m10 + e
	predictedShape[m10>lim] = f * m10[m10>lim] + (d - f) * lim + e
		
	if m10wasFloat and isinstance(G,(int, float))==True:
		return sigmoid(G, predictedG0, predictedInvslope, predictedShape)[0]	
	else:	
		return sigmoid(G, predictedG0, predictedInvslope, predictedShape)
	
	
	
	
		


def selectionFunctionRADEC(ra: float,
                           dec: float,
                           G: np.ndarray,
						   dfm10: pd.DataFrame) -> np.ndarray:
	"""Evaluate the completeness at ra, dec, for given G array.

	Parameters
	----------
	ra : float
		right ascension in degrees
	dec: float
		declination in degrees
	G : np.ndarray
		which magnitude to evaluate the completeness
	dfm10: pd.DataFrame
		the M10 data frame
	healpix_level: int, optional
		defaults to 7, healpix_level of the M10 data.

	Returns
	-------
	Evaluated completeness
	"""
	#first identify the healpix order from the number of hpx in the map:
	nbHpx = len(dfm10)
	healpix_level = int(np.log(nbHpx/12)/np.log(4))
	#then find the m10 value in the hpx corresponding to (ra,dec):
	m10 = dfm10[ hp.ang2pix(2 ** healpix_level, ra, dec, lonlat=True, nest=True) ]
	return selectionFunction(G, m10)


if __name__ == "__main__":
	##################################
	# This is the part where we call the function at chosen coordinates.

	from astropy import units as u
	from astropy.coordinates import SkyCoord

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
	print(completeness21)
	#completeness21 = [ selectionFunctionRADEC(rr, dd, 21, dfm10) for rr, dd in zip(ra, dec) ]

	#and plot
	import matplotlib.pyplot as plt
	plt.figure(1)
	plt.scatter(ra , dec, c=completeness21, vmin=0, vmax=1, s=1)
	plt.colorbar(label='completeness at G=21')
	plt.xlabel('ra (degrees)')
	plt.ylabel('dec (degrees)')

	plt.figure(2)
	plt.scatter(l.flatten(), b.flatten(), c=completeness21, vmin=0, vmax=1, s=1)
	plt.colorbar(label='completeness at G=21')
	plt.xlim(180, -180)
	plt.ylim(-90, 90)
	plt.xlabel('l (degrees)')
	plt.ylabel('b (degrees)')
	plt.show()