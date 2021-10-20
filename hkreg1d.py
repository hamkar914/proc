'''Written by by Hampus Karlsson, Summer 2021.
   e-mail: hamka@chalmers.se/hkarlsson914@gmail.com
   processing script to extract a 1D spectrum. To be executed 
   within Topspin 4, from: C:\Bruker\TopSpin4.1.1\exp\stan\nmr\py\user\ 
   directory.'''
   
# get folder name
curdat = CURDATA()

# get some proc pars
proc_par_str = "Exp = "+ curdat[0]+" "+curdat[1]+"\n"
proc_par_str += "Pulse program = "+GETPAR("PULPROG")+"\n"
proc_par_str += "NS = "+GETPAR("NS")+"\n" 
proc_par_str += "LB = "+GETPAR("LB")+"\n"
proc_par_str += "PHC0 = "+GETPAR("PHC0")+"\n"
proc_par_str += "PHC1 = "+GETPAR("PHC1")+"\n"

# retrieve region limits
lft_lim = float(GETPAR("F1P"))
rght_lim = float(GETPAR("F2P"))

regions = []

# read the 1D	
RE([curdat[0],curdat[1],curdat[2],"C:/Users/hamka/Documents/data/"],"n")

# get slice of this particular 1D
real_spectrum_slice = GETPROCDATA(lft_lim,rght_lim)
regions.append(real_spectrum_slice)

# calculate the ppm vlaues for this region
spectrum_rng = lft_lim-rght_lim
resolution = spectrum_rng/len(real_spectrum_slice)  

ppm_str = ""
start_ppm_calc = float(lft_lim)

for i in range(len(real_spectrum_slice)):
    start_ppm_calc -= resolution
    ppm_str += str(start_ppm_calc)+" "


# Now write all line by line to output file
f = open("C:/Users/hamka/Documents/data/"+curdat[0]+"/"+curdat[1]+"/pdata/1/hk_spec_regs.txt",'w')

f.write(proc_par_str)
f.write(ppm_str+"\n")

for l in range(len(regions)):
	
    tmp_str = ""
    cur_spectrum = regions[l]
    
    for m in range(len(cur_spectrum)):
        value = cur_spectrum[m]
        tmp_str += str(value)+" "
        
    f.write(tmp_str+"\n")

f.close()


