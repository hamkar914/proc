''''Test script written by Hampus Karlsson, August 2021. e-mail:
    hamka@chalmers.se/hkarlsson914@gmail.com. Script for fit and plot
    of Lorentzian/Gaussion functions to 1D NMR spectra. Literature:
    1. Petrakis, L. Spectral Line Shapes. Journal of Chemical Education,
    1967, Vol.44, No.8 pp. 432-436'''


import datetime
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt


# -------------------------------------
# 1. Input files
# -------------------------------------

ON_inpf = "C:/Users/hamka/Documents/data/211019/4/pdata/1/hk_spec_regs.txt"
OFF_inpf = "C:/Users/hamka/Documents/data/211019/3/pdata/1/hk_spec_regs.txt"

# Read spectra from text files
on_data = np.genfromtxt(ON_inpf,skip_header=6)
off_data = np.genfromtxt(OFF_inpf,skip_header=6)

# ppm_values of x-axis
ppms = np.asarray([round(x,2) for x in on_data[0,:]])


# -------------------------------------
# 2. Set necessary parameters for fit
# -------------------------------------

# Plot one of the spectra once
# to get indices of peaks:
plt.plot(on_data[1,:])
plt.show()

# or look which ppm value that
# correspond to certain index
#for i in range(len(on_data[0,:])):
#   print(str(i)+"  "+str(on_data[0,i]))

# Fill list with indices for peak centers
peaks = [2060,9590,19100]

# Choose proportion lorentzian (rest is gaussian)
# for each individual peak
func = [0.1,0.1,0.1]

# Indices, left and right limit for
# limits for regions to be deconvuluted
dcon_limits = [1020,2860, 8500,10530, 18200,19950]


# -------------------------------------
# 3. Functions for fitting line shape
# -------------------------------------

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
# 4. Start deconvulution part
# -------------------------------------

# plot colors
clrs = ["r","g","b","m","y","c"]

# start create a string with output
out_str = "Time: "+str(datetime.datetime.now())+"\n"
out_str+= "Dataset: "+ON_inpf[:-16]+"\n"
out_str+= "ppms = "+str([round(on_data[0,:][peaks[ind]],2) for ind in range(len(peaks))])+"\n"

ppm_resolution = on_data[0,0]-on_data[0,1]
sfo1 = 100.6473721


# -------------------------------------
# 5. Deconvulute MW ON spectrum
# -------------------------------------

# Experimental data set
dset_on = on_data[1]

# array for initial guesses for deconvulotion
p0 = []

# array for low and high boundaries if wanted
# (i.e limits of fitted paramters)
lb,hb = [],[]

# fill arrays with guesses and boundaries
for k in range(len(peaks)):
    # peak intensity as initial guess
    p0.append(dset_on[peaks[k]])
    # peak intensity +/- inf is allowed
    lb.append(-np.inf)
    hb.append(np.inf)
    # Init guess line widts in term of indices
    p0.append(5)
    # limits of line width in indices
    lb.append(-np.inf)
    hb.append(np.inf)

# peak indexes, fit func type and dcon region need
# to be sent to residual funcion
fix_params = [peaks,func,dcon_limits]

# fit lineshape to spectrum
plsq = least_squares(residual, p0, args=(dset_on, fix_params), method="lm") #,bounds=(lb,hb))


# -------------------------------------
# 6. Start plotting
# -------------------------------------

# create two subplots with the shared x and y axes
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12,6))

# plot experimental spectrum
ppms = on_data[0,:]
ax1.plot(ppms, dset_on, "k-")

# plot individual deconvolutions
# with fitted parameters

for l in range(0,len(p0),2):

    lorprt = func[int(l/2)]*lorentzian(peaks[int(l/2)],np.arange(dcon_limits[l],dcon_limits[l+1],1),plsq.x[l],plsq.x[l+1])
    gauprt = (1-func[int(l/2)])*gaussian(peaks[int(l/2)],np.arange(dcon_limits[l],dcon_limits[l+1],1),plsq.x[l],plsq.x[l+1])
    sol = lorprt+gauprt

    # calculate integral or intensity of fitted solution
    ig = trapezoid(sol)
    #dcon_peak_int = sol[peaks[int(l/2)]]
    fwhm_hz = 2*plsq.x[l+1]*ppm_resolution*sfo1

    # add integral to output string
    out_str+=str(round(ig,1))+"\t"

    ig_str = str(func[int(l/2)]*100)+"%,   AUC = "+"{:.3g}".format(ig)+"   HW = "
    ig_str += "{:.4g}".format(fwhm_hz)+" Hz  Reg = "+"{:.4g}".format(ppms[dcon_limits[l]])+"-"+"{:.4g}".format(ppms[dcon_limits[l+1]])+" ppm"
    ax1.plot(ppms[dcon_limits[l]:dcon_limits[l+1]],sol,color = clrs[int(l/2)],linestyle="-",label = ig_str)

out_str+="\n"
ax1.set_ylabel("Peak intensity")
ax1.set_title("MW ON")
ax1.set_xlabel("13C (ppm)")
ax1.legend(bbox_to_anchor=(0., 1.1, 1., 0.2), loc=3, mode="expand", borderaxespad=0.0)



# -------------------------------------
# 7. Deconvolute MW OFF spectrum
# -------------------------------------

# Experimental data set
dset_off = off_data[1]

# array for initial guesses for deconvulotion
p0 = []

# array for low and high boundaries if wanted
# (i.e limits of fitted paramters)
lb,hb = [],[]

# fill arrays with guesses and boundaries
for k in range(len(peaks)):
    # peak intensity as initial guess
    p0.append(dset_off[peaks[k]])
    # peak intensity +/- inf is allowed
    lb.append(-np.inf)
    hb.append(np.inf)
    # Init guess line widts in term of indices
    p0.append(5)
    # limits of line width in indices
    lb.append(-np.inf)
    hb.append(np.inf)

# peak indexes, fit func type and dcon region need
# to be sent to residual funcion
fix_params = [peaks,func,dcon_limits]

# fit lineshape to spectrum
plsq = least_squares(residual, p0, args=(dset_off, fix_params), method="lm") #,bounds=(lb,hb))


# -------------------------------------
# 8. Plot MW OFF deconvulotions
# -------------------------------------

# plot experimental spectrum
ppms = off_data[0,:]
ax2.plot(ppms, dset_off, "k-")

# plot individual deconvolutions
# with fitted parameters

for l in range(0,len(p0),2):

    lorprt = func[int(l/2)]*lorentzian(peaks[int(l/2)],np.arange(dcon_limits[l],dcon_limits[l+1],1),plsq.x[l],plsq.x[l+1])
    gauprt = (1-func[int(l/2)])*gaussian(peaks[int(l/2)],np.arange(dcon_limits[l],dcon_limits[l+1],1),plsq.x[l],plsq.x[l+1])
    sol = lorprt+gauprt

    # calculate integral or intensity of fitted solution
    ig = trapezoid(sol)
    #dcon_peak_int = sol[peaks[int(l/2)]]
    fwhm_hz = 2*plsq.x[l+1]*ppm_resolution*sfo1

    # add integral to output string
    out_str+=str(round(ig,1))+"\t"

    ig_str = str(func[int(l/2)]*100)+"%,   AUC = "+"{:.3g}".format(ig)+"   HW = "
    ig_str += "{:.4g}".format(fwhm_hz)+" Hz  Reg = "+"{:.4g}".format(ppms[dcon_limits[l]])+"-"+"{:.4g}".format(ppms[dcon_limits[l+1]])+" ppm"
    ax2.plot(ppms[dcon_limits[l]:dcon_limits[l+1]],sol,color = clrs[int(l/2)],linestyle="-",label = ig_str)

out_str+="\n"
ax2.set_ylabel("Peak intensity")
ax2.set_title("MW OFF")
ax2.set_xlabel("13C (ppm)")
ax2.legend(bbox_to_anchor=(0., 1.1, 1., 0.2), loc=3, mode="expand", borderaxespad=0.0)

plt.tight_layout()
plt.subplots_adjust(top=0.6)
plt.xlim(on_data[0,:].max(),on_data[0,:].min())
plt.show()
fig.savefig("C:/Users/hamka/Desktop/dcon1dspectrum.png",format="png",dpi=1000)




