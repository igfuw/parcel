import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
from scipy.io import netcdf
from parcel import parcel
import pytest
import os, glob
import filecmp
import subprocess
import pdb


# chwilowo tylko, aby szybko sie liczylo   
Dt_list = [1., 1.5]#[1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 8e-3, 1e-2, 2e-2, 4e-2, 8e-2, 1e-1, 2e-1, 4e-1, 8e-1, 1.] 

# tu korzystam wlasnie z tych wspomnianych fixture, jesli scope z domyslnego(f-cja) zmienie na module, to wykona sie tylko raz
@pytest.fixture(scope="module")
def data_all_timesteps(request):
    # initial values                                        
    RH_init = .99999
    T_init  = 280.
    p_init  = 100000.
    r_init  = common.eps * RH_init * common.p_vs(T_init) / (p_init - RH_init * common.p_vs(T_init))

    # lists to store RH_max and N at the end of the simulation from each test run    
    RH_list = []
    N_list  = []

    for dt in Dt_list:
        print "\nt time step", dt
        outfile_nc = "test_dt=" + str(dt) + ".nc" 
        parcel(dt=dt, outfreq = 1,   outfile = outfile_nc,\
                w = 1., T_0 = T_init, p_0 = p_init, r_0 = r_init, z_max = 20, \
                mean_r = 5e-8, gstdev = 1.5, n_tot = 1e9, sd_conc = 1000.)

        f_out  = netcdf.netcdf_file(outfile_nc, "r")
        RH_max = f_out.variables["RH"][:].max()
        N_end  = sum(f_out.variables["conc"][-1, -9:]) # concentration of drops > 1e-6 m                                                                                  
        RH_list.append((RH_max - 1)*100)  # [%]                                      
        N_list.append(N_end / 1e6)        # [1/mg]           

    # fixture moze zwrocic tylko jedna rzecz, weic wlozylam w slownik
    data = {"RH" : RH_list, "N" : N_list}

    print "TWORZE data"

    # to jest metoda zabijania tego co stworzylam w fixture (nie jest wymagana przez py.tes, ale rozumiem, ze chcemy)
    def removing_files():
        print "\n ZABIJAM data"
        # cos mi * przy rm nie dzialala - TODO
        for file in glob.glob("test_dt*"):
            subprocess.call(["rm", file])
    request.addfinalizer(removing_files)
    return data


def test_timestep_eps(data_all_timesteps, eps=0.2):
    # testing if the values of variables do not differ from ref. more than eps times
    for var in data_all_timesteps.values():
        for val in var:
            assert abs(val - var[0]) <= eps * var[0], "see figures...TODO" #dopisalabym sugestie, aby sprawdzic rysunek z plotu, jesli jest cos nie tak
 
def test_timestep_conv(data_all_timesteps, eps=0.05):
    # szczerze mowiac teraz juz nie pamietam co ma jak zbiegac...
    #RH ma byc jak najmniejsze? zgubilam sie
    pass


# sprawdzam dane ze stworzonymi przeze mnie danymi w katalogu refdata 
# niespecjalnie wiem dlaczego h5diff daje mi "Some objects are not comparable", ale status jest 0, wiec jest ok...
# zrobilam jednak dla kazdego dt osobno
@pytest.mark.parametrize("dt", Dt_list)
def test_timestep_diff(data_all_timesteps, dt, eps=0.2):
        filename =  "test_dt=" + str(dt) + ".nc"
        filename_nc4 = filename.replace(".nc", "_nc4.nc")
        #subprocess.call(["nccopy", "-k", "4", filename, filename_nc4])
        #subprocess.check_call(["h5diff", "--delta=1e-18", os.path.join("test/refdata", filename_nc4), filename_nc4])


#sprawdzam, czy program rysujacy dziala
# czy powinnam dodawac sprawdzanie stworzenia svg? wydaje mi sie, ze nic nie wnosi, ale moge dodac
def test_timestep_plot(data_all_timesteps):
# zmienic plik gnuplota aby czytal tablice N i RH
    #subprocess.check_call(["gnuplot", "test/plots/plot_timestep.gp"])
    pass
