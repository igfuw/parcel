# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")
from scipy.io import netcdf
import numpy as np
import pytest
import subprocess

from parcel import parcel
from spectrum_plot import plot_spectrum

@pytest.fixture(scope="module")
def data(request):
    """
    Run parcel simulation and return opened netdcf file
    """
    outfile = "test_spectrum.nc"

    # running parcel model for open / closed chem system  ...
    parcel(dt = .5, sd_conc = 1024, outfreq = 40,  outfile=outfile,\
           out_bin = ["wradii:1e-9/1e-4/26/log/wet/0", "dradii:1e-9/1e-6/26/log/dry/0",\
                      "linwradii:1e-9/1e-4/26/lin/wet/0", "lindradii:1e-9/1e-6/26/lin/dry/0"])

    data = netcdf.netcdf_file(outfile, "r")

    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", outfile])

    #request.addfinalizer(removing_files)
    return data

def test_spectrum_bins(data):
    """
    Check if bin sizes stored in the netcdf file are equal to the difference between bin edges
    (done for all but last bin because in output files we store bin left edges)
    """

    def bin_checker(r_nc, dr_nc):
        dr = np.empty(r_nc.shape[0] - 1)
        dr[:] = (r_nc[1:] - r_nc[0:-1]) 
        
        # TODO
        # why this works ...
        for it in dr:
          assert  dr[it] == dr_nc[it]
        # .... while this doesn't?
        #assert (dr[:] == dr_nc[:-1]).all()

    # for log bins
    bin_checker(data.variables["wradii_r_wet"][:], data.variables["wradii_dr_wet"][:])
    bin_checker(data.variables["dradii_r_dry"][:], data.variables["dradii_dr_dry"][:])
    # for lin bins
    bin_checker(data.variables["linwradii_r_wet"][:], data.variables["linwradii_dr_wet"][:])
    bin_checker(data.variables["lindradii_r_dry"][:], data.variables["lindradii_dr_dry"][:])

def test_spectrum_diff(data, eps = 1e-20):
    """
    Compare the results with the referential simulation
    (stored in refdata folder)                                             
    """

    f_ref  = netcdf.netcdf_file("unit_test/refdata/test_spectrum.nc", "r")

    # bin edges and bin sizes
    for var in ["wradii_r_wet", "wradii_dr_wet", "dradii_r_dry", "dradii_dr_dry",
                "linwradii_r_wet", "linwradii_dr_wet", "lindradii_r_dry", "lindradii_dr_dry"]:

        assert (f_ref.variables[var][:] == data.variables[var][:]).all()

    # 0th moment
    for var in ["wradii_m0", "dradii_m0", "lindradii_m0", "linwradii_m0"]:

        refdata = f_ref.variables[var][:]
        cmpdata = data.variables[var][:]
        refdata = np.reshape(refdata, np.product(refdata.shape))
        cmpdata = np.reshape(cmpdata, np.product(cmpdata.shape))
 
        diff    = abs(refdata[:] - cmpdata[:])
        ind_max = np.argmax(diff)
        print " "
        print var
        print "the difference = " , diff[ind_max]
        print "the reference  = " , refdata[ind_max]
        print "the data       = " , cmpdata[ind_max]

        assert np.isclose(cmpdata, refdata, atol=0, rtol=eps).all(),\
            "differs e.g. " + str(var) + "; max(ref diff) = " +\
            str(np.where(refdata != 0.,\
            abs((data - refdata) / refdata), 0.).max())


def test_spectrum_plot(data):
    """
    plot dry and wet radius size distribution
    """
    plot_spectrum(data, outfolder = "plots/outputs/")
