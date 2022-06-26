#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2019-2022

import os, sys
import numpy as np
import h5py
from mytools import print, stop
import argparse

# ============= get_shape =============
def get_shape(filename, datasetname):
    # open HDF5 file
    f = h5py.File(filename, "r")
    # get shape
    dset = f[datasetname]
    shape = dset.shape
    f.close()
    return shape
# ============= end: get_shape =============

# ============= read =============
def read(filename, datasetname, ind=()):
    # open HDF5 file
    f = h5py.File(filename, "r")
    # grab the dataset as a numpy array (with requested index ind, if provided, else all)
    dset = f[datasetname]
    data = np.array(dset[ind])
    f.close()
    return data
# ============= end: read =============

# ============= write =============
def write(data, filename, datasetname, overwrite_file=False, overwrite_dataset=False):
    # first test if file exists
    openflag = "w"
    if os.path.isfile(filename):
        if not overwrite_file: openflag = "a" # append/modify
    # open HDF5 file
    f = h5py.File(filename, openflag)
    # check if dataset name exists in hdf5 file
    if datasetname in f.keys():
        if overwrite_dataset:
            del f[datasetname] # if we overwrite, then delete first
        else:
            print("dataset '"+datasetname+"' already in file '"+filename+"'. Use option 'overwrite_dataset=True' to overwrite.", error=True)
    # create and write data
    dset = f.create_dataset(datasetname, data=data)
    f.close()
    print("'"+datasetname+"' written in file '"+filename+"'", highlight=3)
    return data
# ============= end: write =============

# ============= delete =============
def delete(filename, datasetname, quiet=False):
    # open HDF5 file
    f = h5py.File(filename, "a")
    # check if dataset name exists in hdf5 file
    if datasetname in f.keys():
        if not quiet: print("deleting dataset '"+datasetname+"' in file '"+filename+"'.", warn=True)
        del f[datasetname]
    f.close()
# ============= end: delete =============

# ============= get_dataset_names =============
def get_dataset_names(filename):
    # open HDF5 file
    f = h5py.File(filename, "r")
    dsets = list(f.keys())
    f.close()
    return dsets
# ============= end: get_dataset_names =============


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Manipulate HDF files.')
    parser.add_argument("-i", "--i", dest='filename', nargs='*', type=argparse.FileType('r'), help="HDF5 file(s)")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_usage()
        exit()

    # get list of files
    filename = sorted([x.name for x in list(args.filename)])

    for filen in filename:
        print("Working on file '"+filen+"'")
        stop()

# ======================= MAIN End ===========================
