# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 14:48:29 2021

@author: aduell
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate, scipy.integrate, sys
assert sys.version_info >= (3,6), 'Requires Python 3.6+'



# One more package: A units-and-constants package I wrote: http://pypi.python.org/pypi/numericalunits
# 
# Example: <span style="color:blue">x = 5 * cm</span> means "x equals 5 centimeters".
# 
# Example: <span style="color:blue">y = x / mm</span> means "y is the numerical value of x in millimeters". 

# In[2]:


from numericalunits import W, K, nm, m, um, cm, s, eV, meV, V, hPlanck, pi, mA, c0, kB, e
#import numericalunits
#help(numericalunits)

####EVENTUALLY CHANGE TO RUN SEAMLESSLY WITH 4 BEM CODE##############
TRfRbNoUnit = pd.read_csv('/Users/aduell/Desktop/CodeThings/pv-window-bem-master/pv-window-bem-master/Output/TRfRb.txt')
#TRfRbNoUnit = 
  #Import["/Users/lwheeler/Research/Projects/SwitchGlaze/BTO_SwitchGlaze_FY20/BEM/Window_inputs/20200722_BEM_big_run/TRfRb_50pct_select1.txt", "CSV", 
   #HeaderLines -> 1];
EQENoUnit = pd.read_csv('/Users/aduell/Desktop/CodeThings/pv-window-bem-master/pv-window-bem-master/Output/EQE.txt')
#EQENoUnit = 
  #Import["/Users/lwheeler/Research/Projects/SwitchGlaze/BTO_SwitchGlaze_FY20/BEM/Window_inputs/20200722_BEM_big_run/EQE_50pct_select1.txt", "CSV", 
   #HeaderLines -> 1];

   
EQENoUnit = np.array(EQENoUnit)
TRfRbNoUnit = np.array(TRfRbNoUnit)


SpeedOfLight = 2.9979e08 #m/s
h = 6.626e-34 #J*s


#Apply wavelength units of um
EQEum = EQENoUnit #um
#EQEum[:,0] *= um

TRfRbum = TRfRbNoUnit
#TRfRbum[:,0] *= um

#Apply energy units of eV
EQEeV = np.array(EQENoUnit)
EQEeV [:,0] = (h*SpeedOfLight/EQENoUnit[:,0]) 

TRfRbeV =  np.array(TRfRbNoUnit)
TRfRbeV[:,0] = (h*SpeedOfLight/TRfRbNoUnit[:,0])   #J*s*m/(s*um)
#Something is wrong with this calculation. Off by a factor of 100


print(EQEum-EQEeV)
print(TRfRbum-TRfRbeV)




plt.plot(TRfRbum[:,0], TRfRbum[:,1], color='magenta',marker=None,label="$T$")
plt.plot(TRfRbum[:,0], TRfRbum[:,2],color='green',marker=None,label="$R_f$")
plt.plot(TRfRbum[:,0], TRfRbum[:,3],color='purple',marker=None,label="$R_b$")
plt.plot(EQEum[:,0], EQEum[:,1],color='black',marker=None,label="EQE")
plt.legend(loc = 'upper right')
plt.xlabel('Wavelength, micron')
plt.show()
#This plot accurately displays the raw data range from 0.3 to 2.5 microns.



plt.plot(TRfRbeV[:,0] / eV, TRfRbeV[:,1], color='magenta',marker=None,label="$T$")
plt.plot(TRfRbeV[:,0] / eV, TRfRbeV[:,2],color='green',marker=None,label="$R_f$")
plt.plot(TRfRbeV[:,0] / eV, TRfRbeV[:,3],color='purple',marker=None,label="$R_b$")
plt.plot(EQEeV[:,0] / eV, EQEeV[:,1],color='black',marker=None,label="EQE")
plt.legend(loc = 'upper right')
plt.xlabel('Energy, eV')
plt.show()
#Range should be 4.13 eV to 0.50


#__________________________MATHEMATICA_______________________________#











                 
                          
                    #ISOLATE INDIVIDUAL CURVES#
#DarkTmicronNoInterp = TRfRbNoUnit [:, [0, 1]]#isolate T data from set #no longer needed but make sure the numebrs still work   #Can simplify by removing this line and changing interp line to TRfRbNoUnit instead of the output from this line      
#DarkTmicronNoInterp = ({#[[1]] micron, #[[2]]} &) /@ TRfRbNoUnit;
DarkTmicron = scipy.interpolate.interp1d(TRfRbum[:,0], TRfRbum[:,1])
#DarkTmicron = Interpolation[DarkTmicronNoInterp, InterpolationOrder -> 1];
#DarkTeVNoInterp = DarkTmicronNoInterp                     #No longer needed but make sure still works
#DarkTeVNoInterp = DarkTmicronNoInterp;
#DarkTeVNoInterp = (hPlanck*SpeedOfLight/DarkTeVNoInterp) #No longer needed but make sure still works
#DarkTeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/DarkTeVNoInterp[[All, 1]];
DarkTeV = scipy.interpolate.interp1d(TRfRbeV[:,0], TRfRbeV[:,1])
#DarkTeV = Interpolation[DarkTeVNoInterp, InterpolationOrder -> 1];


#DarkRfmicronNoInterp = ({#[[1]] micron, #[[3]]} &) /@ TRfRbNoUnit;
DarkRfmicron = scipy.interpolate.interp1d(TRfRbum[:,0], TRfRbum[:,2])
#DarkRfmicron = Interpolation[DarkRfmicronNoInterp, InterpolationOrder -> 1];
#DarkRfeVNoInterp = DarkRfmicronNoInterp;
#DarkRfeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/DarkRfeVNoInterp[[All, 1]];
DarkRfeV = scipy.interpolate.interp1d(TRfRbeV[:,0], TRfRbeV[:,2])
#DarkRfeV = Interpolation[DarkRfeVNoInterp, InterpolationOrder -> 1];

#DarkRbmicronNoInterp = ({#[[1]] micron, #[[4]]} &) /@ TRfRbNoUnit;
DarkRbmicron = scipy.interpolate.interp1d(TRfRbum[:,0], TRfRbum[:,3])
#DarkRbmicron = Interpolation[DarkRbmicronNoInterp, InterpolationOrder -> 1];
#DarkRbeVNoInterp = DarkRbmicronNoInterp;
#DarkRbeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/DarkRbeVNoInterp[[All, 1]];
DarkRbeV = scipy.interpolate.interp1d(TRfRbeV[:,0], TRfRbeV[:,3])
#DarkRbeV = Interpolation[DarkRbeVNoInterp, InterpolationOrder -> 1];

#AbsDarkMicronNoInterp = EQENoUnit 
#AbsDarkMicronNoInterp = ({# [[1]] micron, #[[2]]} &) /@ EQENoUnit;
AbsDarkMicron = scipy.interpolate.interp1d(EQEum[:,0], EQEum[:,1])  #it worked?
#AbsDarkMicron = Interpolation[AbsDarkMicronNoInterp, InterpolationOrder -> 1];
#AbsDarkeVNoInterp = AbsDarkMicronNoInterp
#AbsDarkeVNoInterp = AbsDarkMicronNoInterp;
#AbsDarkeVNoInterp = hPlanck*SpeedOfLight/AbsDarkeVNoInterp
#AbsDarkeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/AbsDarkeVNoInterp[[All, 1]];
AbsDarkeV= scipy.interpolate.interp1d(EQEeV[:,0], EQEeV[:,1])
#AbsDarkeV = Interpolation[AbsDarkeVNoInterp, InterpolationOrder -> 1];



#DarkAeVnointerp = TRfRbum[:,[0,1]]
#DarkAeVnointerp[:,1] = 1 - TRfRbum[:,2] - TRfRbum[:,1]
#DarkAeV = DarkAeVnointerp
#DarkAeV = scipy.interpolate.interp1d(DarkAeV[:,0], DarkAeV[:,1])







# this thingy takes two functions and one argument, 
# sticks the argument into each function then returns the result
#def add2funcs(f1,f2,x):
#    return f1(x) + f2(x)
  

#Ephoton = np.linspace(E_min, E_max, 100)
#addem = add2funcs(DarkRfeV, DarkTeV, Ephoton)
#print(addem)

##Units debug
#5.668124528078765e-29
#2.5628962674293367e-30
#Moopsy
 



################THIS IS THE CURRENT PROBLEM####################

#DarkAeV(Ephoton) = 1 - DarkRfeV(Ephoton) - DarkTeV(Ephoton)
#DarkAeV[Ephoton_] := 1 - DarkRfeV[Ephoton] - DarkTeV[Ephoton]

########################################################





#Plot[{DarkTeV[Ephoton eV], DarkRfeV[Ephoton eV], DarkRbeV[Ephoton eV], 
   #DarkAeV[Ephoton eV], AbsDarkeV[Ephoton eV]}, {Ephoton, Emin/eV , Emax/eV },
   #PlotRange -> {{Emin/eV, Emax/eV}, {-0.1, 1.1}}, PlotStyle -> {Thick}, 
  #Filling -> Axis, ImageSize -> 600, Frame -> True, 
  #FrameLabel -> {Style["Energy (eV)", 16]}, 
  #LabelStyle -> Directive[FontSize -> 16], 
  #PlotLegends -> {"\!\(\*SubscriptBox[\(T\), \(Inc\)]\)(E)", 
    #"\!\(\*SubscriptBox[\(R\), \(f\)]\)(E)", 
    #"\!\(\*SubscriptBox[\(R\), \(b\)]\)(E)", 
    #"\!\(\*SubscriptBox[\(A\), \(Total\)]\)(E)", 
    #"\!\(\*SubscriptBox[\(A\), \(ST\)]\)"}] // Quiet
    
  




