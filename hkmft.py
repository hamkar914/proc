'''Written by by Hampus Karlsson, Summer 2021.
   e-mail: hamka@chalmers.se/hkarlsson914@gmail.com
   processing script for fourier transforming a set of 
   1D spectra from a pseudo 2D nmr experiment. To be executed 
   within Topspin 4, from: C:\Bruker\TopSpin4.1.1\exp\stan\nmr\py\user\ 
   directory.'''

td = GETPAR("1 TD")

curdat = CURDATA()

out_string = ""
out_string2 = ""

for i in range(2,int(td)+2):
	  
	  out_string += "efp "+str(i-1)+" "+str(i)+" y n\n"
	  
for i in range(2,int(td)+2):
	  
	  out_string2 += "re "+curdat[1]+" "+str(i)+"\nabs\n"
	
f = open("C:/Bruker/TopSpin4.1.1/exp/stan/nmr/lists/mac/user/hk_ft_all",'w')
f.write(out_string)
f.write(out_string2)
f.close()