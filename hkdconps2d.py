''''Test script written by Hampus Karlsson, August 2021.
    e-mail: hamka@chalmers.se/hkarlsson914@gmail.com
    Script for deconvoluting a series of 1D spectra from
    a pseudo 2D experiment.

    Literature:
    Petrakis, L. Spectral Line Shapes. Journal of Chemical Education,
    1967, Vol.44, No.8 pp.432-436'''

import datetime
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt


def lorentzian(x0,x,I0,R2):
    # x0 = center freq
    # x = offset
    # I0 = amplitude at x0
    # R2 = half width
    return I0*(R2/((x-x0)**2+R2**2))


def gaussian(x0,x,I0,R2):
    # x0 = center freq
    # x = offset
    # I0 = amplitude at x0
    # R2 = half width
    return I0*((np.sqrt(np.log(2))/np.sqrt(R2))*np.exp(-(x-x0)**2*(np.log(2)/R2**2)))


def residual(popt,y_exp,fix_args):

    # unpack fix_args
    pks,fns,lim = fix_args

    # check how many peaks to be fitted
    n_peaks = len(pks)

    # create a copy of experimental spectrum
    y_exp_cp = np.copy(y_exp)

    # replace the regions to be fitted with zeros
    for k in range(0,len(lim),2):
        y_exp_cp[lim[k]:lim[k+1]] = 0.0

    # OBS!!! Reason for a copy spectrum containing zero regions
    # is to be able to return a residual as an array and not a scalar
    # as this is demanded by scipy least_squares lm. Regions that
    # are not zero are equal to experimental spectrum and will not
    # contribute to the residual.

    # loop through all peaks and fit the
    # calculated line shape with current popt set
    for i in range(0,2*n_peaks,2):
        # Lorentizian part
        lor_part = fns[int(i/2)]*lorentzian(pks[int(i/2)],np.arange(lim[i],lim[i+1],1),popt[i],popt[i+1])
        # Gaussian part
        gaus_part = (1.0-fns[int(i/2)])*gaussian(pks[int(i/2)],np.arange(lim[i],lim[i+1],1),popt[i],popt[i+1])
        # sum them up
        y_calc = lor_part+gaus_part
        # add the calculated lineshape for a certain region
        # to the corresponding region in y_exp_cp
        y_exp_cp[lim[i]:lim[i+1]] += y_calc

    sq_diff = (y_exp-y_exp_cp)**2

    # possibly print params during optimization if interested
    #print(str([x for x in popt]))

    return sq_diff

# -------------------------------------
# 1. Get input data
# -------------------------------------
data_dir = "C:/Users/hamka/Documents/data/"
mw_on_data_fldr = "211019/5"


# Read spectra from text files
on_data = np.genfromtxt(data_dir+mw_on_data_fldr+"/pdata/2/hk_spec_regs.txt",skip_header=6)

# get vdlist
vd_list = np.genfromtxt(data_dir+mw_on_data_fldr+"/vdlist")


# -------------------------------------
# 2. Set necessary parameters for fit
# -------------------------------------

# Plot one of the spectra once
# to get indices of peaks:
plt.plot(on_data[2,:])
plt.show()

# or look which ppm value that
# correspond to certain index
#for i in range(len(on_data[0,:])):
#   print(str(i)+"  "+str(on_data[0,i]))

# Fill list with indices for peak centers
peaks = [257,1192,2394]

# Choose proportion lorentzian (rest is gaussian)
# for each individual peak
func = [0.1,0.1,0.1]

# Indices, left and right limit for
# limits for regions to be deconvuluted
dcon_limits = [91,399, 1039,1329, 2282,2506]



# -------------------------------------
# 3. Start deconvultion part
# -------------------------------------

# plot colors
clrs = ["r","g","b","m","y","c"]

# start create a string with output
out_str = "Time: "+str(datetime.datetime.now())+"\n"
out_str+= "Dataset: "+data_dir+mw_on_data_fldr+"\n"
out_str+= "ppms = "+str([round(on_data[0,:][peaks[ind]],2) for ind in range(len(peaks))])+"\n"

ppm_resolution = on_data[0,0]-on_data[0,1]
sfo1 = 100.6473721

for i in range(1,on_data.shape[0]):

    # Experimental data set
    dset = on_data[i]

    # array for initial guesses for deconvulotion
    p0 = []

    # array for low and high boundaries if wanted
    # (i.e limits of fitted paramters)
    lb,hb = [],[]

    # fill arrays with guesses and boundaries
    for k in range(len(peaks)):
        # peak intensity as initial guess
        p0.append(dset[peaks[k]])
        # peak intensity +/- inf is allowed
        lb.append(-np.inf)
        hb.append(np.inf)
        # Init guess line widts in term of indices
        p0.append(0.5)
        # limits of line width in indices
        lb.append(-np.inf)
        hb.append(np.inf)

    # peak indexes, fit func type and dcon region need
    # to be sent to residual funcion
    fix_params = [peaks,func,dcon_limits]

    # fit lineshape to spectrum
    plsq = least_squares(residual, p0, args=(dset, fix_params), method="lm")#,bounds=(lb,hb))

    # create two subplots with the shared x and y axes
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True,figsize=(6,6))

    # plot experimental spectrum
    ppms = on_data[0,:]
    ax1.plot(ppms, dset, "k-")

    # add data to output string
    out_str += str(vd_list[i-1])+" "

    # -------------------------------------
    # 4. Plotting
    # -------------------------------------

    # plot individual deconvolutions
    # with fitted parameters
    for l in range(0,len(p0),2):

        lorprt = func[int(l/2)]*lorentzian(peaks[int(l/2)],np.arange(0,dset.shape[0],1),plsq.x[l],plsq.x[l+1])
        gauprt = (1-func[int(l/2)])*gaussian(peaks[int(l/2)],np.arange(0,dset.shape[0],1),plsq.x[l],plsq.x[l+1])
        sol = lorprt+gauprt

        # calculate integral or intensity of fitted solution
        ig = trapezoid(sol[dcon_limits[l]:dcon_limits[l+1]])
        #dcon_peak_int = sol[peaks[int(l/2)]]
        fwhm_hz = 2*plsq.x[l+1]*ppm_resolution*sfo1

        # add integral to output string
        out_str+=str(round(ig,1))+"\t"

        ig_str = str(func[int(l/2)]*100)+"%,   AUC = "+"{:.3g}".format(ig)+",   HW = "+"{:.4g}".format(fwhm_hz)+" Hz"
        ax1.plot(ppms[dcon_limits[l]:dcon_limits[l+1]],sol[dcon_limits[l]:dcon_limits[l+1]],color = clrs[int(l/2)],linestyle="-",label = ig_str)

    out_str+="\n"
    ax1.set_ylabel("Peak intensity")
    #ax1.set_xlabel("13C (ppm)")

    ax1.set_title("VD = "+str(vd_list[i-1])+" s, individual fitted line shapes")
    ax1.legend(bbox_to_anchor=(0., 1.2, 1., 0.2), loc=3,ncol=1, mode="expand", borderaxespad=0.0)

    tot_lineshape = np.zeros(dset.shape[0])


    # plot the total lineshape
    # -------------------------------------
    for m in range(0,len(p0),2):

        lorprt = func[int(m/2)]*lorentzian(peaks[int(m/2)], np.arange(0,dset.shape[0],1),plsq.x[m],plsq.x[m+1])
        gauprt = (1-func[int(m/2)])*gaussian(peaks[int(m/2)],np.arange(0,dset.shape[0],1),plsq.x[m],plsq.x[m+1])
        tot_lineshape += lorprt+gauprt

    sum_res = "{:.3g}".format(np.sum(residual(plsq.x,dset,fix_params)))

    ax2.plot(ppms, dset, "k-")
    ax2.plot(ppms,tot_lineshape,"r-",label="res = "+str(sum_res))
    ax2.set_title("Total fitted lineshape")

    plt.xlim(on_data[0,:].max(),on_data[0,:].min())
    ax2.set_xlabel("13C (ppm)")
    ax2.set_ylabel("Peak intensity")
    plt.tight_layout()
    plt.legend()
    plt.subplots_adjust(top=0.7)
    plt.show()
    fig.savefig(data_dir+mw_on_data_fldr+"/pdata/2/"+"dcon_vd_"+str(i)+".png",format="png",dpi=300)

print(out_str)

outpf = open(data_dir+mw_on_data_fldr+"/pdata/2/"+"hk_dcon_integrals.txt","w")
outpf.write(out_str)
outpf.close()




