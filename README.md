The functions are located in gu.py

An example can be found in test3.py. The two necessary imports are:

    from gu import load_data, selectionFunctionRADEC


## To do:

* add pre-computed maps in healpix level 8 and 9, and allow load_data() to take arguments to choose map

* if needed: a hybrid map with a resolution scaling with stellar density

* maps are currently saved in uncompressed HDF5 format (pandas seems to struggle with gzipped hdf5)

* function to compute and save a custom map for any region and resoltuion (they are just 2D histograms anyway)
