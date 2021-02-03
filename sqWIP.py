#!/usr/bin/env python
# coding: utf-8

# # The Shockley-Queisser limit
# 
# By Steven J. Byrnes ([https://sjbyrnes.com/](https://sjbyrnes.com/)). This document lives at [https://github.com/sbyrnes321/SolarCellEfficiencyLimits](https://github.com/sbyrnes321/SolarCellEfficiencyLimits). Please email me any feedback: steven.byrnes@gmail.com
# 
# The Shockley-Queisser (SQ) limit is a famous limit on the maximal possible efficiency of solar cells, limited only by fundamental physics. It applies to most solar cell designs in the world, except for "tandem solar cells" and some additional obscure exceptions (discussed at the end of the document). The most important parameter in the SQ model is the bandgap of the semiconductor: If the gap is right, the efficiency can be up to 34%, if the gap is way off, the efficiency limit may be much smaller. [Here is the original SQ paper](http://dx.doi.org/10.1063/1.1736034), but it’s also covered in every solar-cell textbook.
# 
# I’m using NREL’s data for the solar spectrum (AM1.5G) and intensity (1000 W/m²). In the original SQ paper, they assumed that
# the sun had a 6000-kelvin blackbody spectrum. So my graphs and values are slightly different. However, other papers and books
# that use AM1.5G spectrum get the same results as I do, for example [link 1](http://www.opticsinfobase.org/abstract.cfm?URI=OSE-2010-SWA1), [link 2](http://www.opticsinfobase.org/abstract.cfm?URI=OSE-2010-SWC4), *Practical Handbook of Photovoltaics* p128-9,
# [link 3](http://dx.doi.org/10.1109/T-ED.1984.21594).
# 
# I copied many of these graphs into the Wikipedia article on this topic - [http://en.wikipedia.org/wiki/Shockley-Queisser_limit](http://en.wikipedia.org/wiki/Shockley-Queisser_limit)
# 
# In this document you will find:
# 
# * A plot of the SQ efficiency limit as a function of bandgap
# * A plot of the SQ limit on short-circuit current, on open-circuit voltage, and on fill-factor, as a function of bandgap
# * A breakdown of exactly which factors lower the SQ limit for which bandgaps
# * A list of some "loopholes" to exceed the SQ limit.
# 
# Enjoy!
# 
# <p style="font-size:80%">Pronunciation of "Queisser": Hans-Joachim Queisser was German, so a German-speaker helped me guess how the name is pronounced. He guesses that "Queisser" rhymes with "nicer". ("Qu" as in "quick", "ei" as in "Einstein", "ss" as in "kiss", "er" as in "teacher"). (Thanks Florian!)</p>
# 
# Note: If you run all the code in this file, including re-creating all the graphs, it may take a few hours. (I made no effort to write efficient code.)
# 
# ## General program setup
# 
# This document is a mix of text and Python code, written using [Jupyter Notebook](http://jupyter.org/) (You can install Jupyter notebook through [Anaconda](https://www.anaconda.com/distribution/).)

# In[1]:


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


from numericalunits import W, K, nm, m, um, cm, s, eV, meV, V, hPlanck, pi, mA, c0, e
#import numericalunits
#help(numericalunits)
#hPlanck = 6.626e-34
SpeedOfLight = 299792458 #in m/s bleh figure out how to change later


# ## Program inputs
# 
# Solar cell temperature is 300 kelvin:

# In[3]:


#Tcell = 300 * K


# The incident light intensity and spectrum is assumed to be the NREL AM1.5G spectrum, which approximates the light coming from the sun and sky at a typical latitude on a clear day. For more information go to https://www.nrel.gov/grid/solar-resource/spectra.html

# In[4]:


worksheet = pd.read_excel('https://www.nrel.gov/grid/solar-resource/assets/data/astmg173.xls')
downloaded_array = np.array(worksheet)

# Wavelength is in column 0, AM1.5G data is column 2
AM15 = downloaded_array[1:, [0,2]]

# The first line should be 280.0 , 4.7309E-23
# The last line should be 4000.0, 7.1043E-03
#print(AM15)


# Tack on the appropriate units:

# In[5]:


AM15[:,0] *= nm
AM15[:,1] *= W / m**2 / nm


# The NREL data spans the following spectral range (in terms of both photon-wavelength and photon-frequency):

# In[6]:


λ_min = 280 * nm
λ_max = 4000 * nm
E_min = hPlanck * c0 / λ_max
E_max = hPlanck * c0 / λ_min


# Interpolate to get a continuous function which I will be able to do integrals on:

# In[7]:


AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])


# Here’s the plot, it looks correct:

# In[8]:


λs = np.linspace(λ_min, λ_max, num=500)
y_values = np.array([AM15interp(x) for x in λs])

plt.plot(λs / nm , y_values / (W/m**2/nm))
plt.xlabel('Wavelength (nm)')
plt.ylabel("Spectral intensity (W/m$^2$/nm)")
plt.title("Light from the sun")
plt.show()

# ## Properties of incident sunlight
# ### Put solar spectrum data in more convenient form
# It’s a bit more convenient for me to change the units for the solar spectrum, so that I can easily do integrals over photon energy, rather than wavelength, and calculate the number of photons instead of their energy. Therefore, I’ll define the function "SPhotonsPerTEA" which stands for Solar Photons per unit Time, per unit photon Energy-range, per unit Area of the solar cell (assuming the cell is facing normal to the sun). To convert from the AM1.5 data to these new units, the formula is:
# 
# $\text{SPhotonsPerTEA} = \frac{d(\text{number of photons per unit time per unit area})}{dE} = \frac{d(\text{photon power per unit area})}{d\lambda} \; \frac{(\text{number of photons per unit time per unit area})}{(\text{photon power per unit area})} \left| \frac{d\lambda}{dE} \right| = $
# $ = (\text{AM1.5 spectrum}) \; \frac{1}{\text{photon energy}} \; \frac{hc}{E^2}$
# 
# (I used $\left| \frac{d\lambda}{dE} \right| = \left| \frac{d}{dE} (\frac{hc}{E}) \right| = \frac{hc}{E^2}$.)

# In[9]:


def SPhotonsPerTEA(Ephoton):
    λ = hPlanck * c0 / Ephoton
    #print('λ = ',λ)
    return AM15interp(λ) * (1 / Ephoton) * (hPlanck * c0 / Ephoton**2)


# Example: The following calculation means that there are $1.43 \times 10^{18}$ solar photons with energy between 2eV and 2.001eV that hit a 1-square-meter patch each second:

# In[10]:


print('# of photons per m^2/s is',SPhotonsPerTEA(2 * eV) * (1 * meV) * (1 * m**2) * (1 * s))


# Next: The "Solar constant" is the sun's total irradiance. If I did this right, it should be 1000 watts/meter$^2$, because that's how NREL normalized their data.

# In[11]:


PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
# the messages warning about poor accuracy in integrating.
solar_constant = scipy.integrate.quad(PowerPerTEA,E_min,E_max, full_output=1)[0]
print('Solar Constant should be ~ 1000. Calculated as',solar_constant / (W/m**2))

# Close enough!
# 




#_______________________________________________________________________________________
########STOP HERE###################
################ALL NEW############## ADAM DUELL##############
#____________________________________________________________________________________________________



#Original mathematica code is commented out below python equivalent


#########description from mathematica#### change to match what will actually happen when finished####

#VLT = 50% selective absorber (C60/chloroaluminium phthalocyanine)
#Here I import the optical data from the transfer matrix method (TMM) python software.
#Import Python data and add units (micron or eV). The first file is T, Rf, and Rb. The second file is called EQE, but it is really the photons absorbered in the absorber layer. It is converted into a "proper" EQElater in the code (EQEdark[] function)





####EVENTUALLY CHANGE TO RUN SEAMLESSLY WITH 4 BEM CODE##############
TRfRbum = pd.read_csv('/Users/aduell/Desktop/CodeThings/pv-window-bem-master/pv-window-bem-master/Output/TRfRb.txt')
EQEum = pd.read_csv('/Users/aduell/Desktop/CodeThings/pv-window-bem-master/pv-window-bem-master/Output/EQE.txt')

   
EQEum = np.array(EQEum)
TRfRbum = np.array(TRfRbum)

SpeedOfLight = 299792458 #m/s
h = 4.135667516e-15 #eV*s

#Convert from um to energy units of eV
EQEeV = np.array(EQEum)
EQEeV [:,0] = (h*SpeedOfLight*1e6/EQEeV[:,0]) 


#Variables
T = TRfRbum[:,1]
Rf = TRfRbum[:,2]
Rb = TRfRbum[:,3]
EQE = EQEum[:,1]
lams = EQEum[:,0]
Ephoton = EQEeV [:,0]
Emin = min(Ephoton)
Emax = max(Ephoton)
print('Emin = ',Emin, 'and Emax = ',Emax, 'eV')

eta = .9 #unitless
Tcell=300 #K
kB = 8.61733034e-5 #eV/K



#Calculate Absorbance
DarkAeVnointerp = np.array(EQEeV)#Array is arbitrary so long as 0 is Ephoton, and 1 is mutable
DarkAeVnointerp[:,1] = 1 - Rf - T
Abs = DarkAeVnointerp[:,1]
DarkAeV = scipy.interpolate.interp1d(Ephoton, Abs, fill_value="extrapolate")
#Equivalent mathematica code normally comes below the interpolations: 'DarkAeV[Ephoton_] := 1 - DarkRfeV[Ephoton] - DarkTeV[Ephoton]'

#plot of data in microns
plt.plot(lams, T, color='magenta',marker=None,label="$T$")
plt.plot(lams, Rf, color='green',marker=None,label="$R_f$")
plt.plot(lams, Rb, color='purple',marker=None,label="$R_b$")
plt.plot(lams, EQE, color='black',marker=None,label="EQE")
plt.legend(loc = 'upper right')
plt.xlabel('Wavelength, micron')
plt.show()
#This plot accurately displays the raw data range from 0.3 to 2.5 microns.

#Plot of data in eV
plt.plot(Ephoton, T, color='magenta',marker=None,label="$T$")
plt.plot(Ephoton, Rf,color='green',marker=None,label="$R_f$")
plt.plot(Ephoton, Rb,color='purple',marker=None,label="$R_b$")
plt.plot(Ephoton, EQE,color='black',marker=None,label="EQE")
plt.plot(Ephoton, DarkAeVnointerp[:,1], color='blue',marker=None,label="$Abs$")
plt.legend(loc = 'upper right')
plt.xlabel('Energy, eV')
plt.show()
#Range should be 4.13 eV to 0.50. Is now accurate

#__________________________MATHEMATICA_______________________________#











                 
                          
                    #ISOLATE AND INTERPOLATE INDIVIDUAL CURVES#
#For T
DarkTmicron = scipy.interpolate.interp1d(lams, T)
DarkTeV = scipy.interpolate.interp1d(Ephoton, T)
#DarkTmicronNoInterp = ({#[[1]] micron, #[[2]]} &) /@ TRfRbNoUnit;
#DarkTmicron = Interpolation[DarkTmicronNoInterp, InterpolationOrder -> 1];
#DarkTeVNoInterp = DarkTmicronNoInterp;
#DarkTeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/DarkTeVNoInterp[[All, 1]];
#DarkTeV = Interpolation[DarkTeVNoInterp, InterpolationOrder -> 1];

#For Rf
DarkRfmicron = scipy.interpolate.interp1d(lams, Rf)
DarkRfeV = scipy.interpolate.interp1d(Ephoton, Rf)
#DarkRfmicronNoInterp = ({#[[1]] micron, #[[3]]} &) /@ TRfRbNoUnit;
#DarkRfmicron = Interpolation[DarkRfmicronNoInterp, InterpolationOrder -> 1];
#DarkRfeVNoInterp = DarkRfmicronNoInterp;
#DarkRfeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/DarkRfeVNoInterp[[All, 1]];
#DarkRfeV = Interpolation[DarkRfeVNoInterp, InterpolationOrder -> 1];

#For Rb
DarkRbmicron = scipy.interpolate.interp1d(lams, Rb)
DarkRbeV = scipy.interpolate.interp1d(lams, Rb)
#DarkRbmicronNoInterp = ({#[[1]] micron, #[[4]]} &) /@ TRfRbNoUnit;
#DarkRbmicron = Interpolation[DarkRbmicronNoInterp, InterpolationOrder -> 1];
#DarkRbeVNoInterp = DarkRbmicronNoInterp;
#DarkRbeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/DarkRbeVNoInterp[[All, 1]];
#DarkRbeV = Interpolation[DarkRbeVNoInterp, InterpolationOrder -> 1];

#For EQE
AbsDarkMicron = scipy.interpolate.interp1d(lams, EQE)  #it worked?
AbsDarkeV= scipy.interpolate.interp1d(Ephoton, EQE, fill_value="extrapolate")
#AbsDarkMicronNoInterp = ({# [[1]] micron, #[[2]]} &) /@ EQENoUnit;
#AbsDarkMicron = Interpolation[AbsDarkMicronNoInterp, InterpolationOrder -> 1];
#AbsDarkeVNoInterp = AbsDarkMicronNoInterp;
#AbsDarkeVNoInterp[[All, 1]] = hPlanck SpeedOfLight/AbsDarkeVNoInterp[[All, 1]];
#AbsDarkeV = Interpolation[AbsDarkeVNoInterp, InterpolationOrder -> 1];




#TEST SPACE#
nums = 200
TestSpaceum = np.linspace(.3,2.48,nums)
TestSpaceeV = np.linspace(.4,2.4,nums)
TestArray = np.array(TestSpaceum)
TestPlot = DarkTmicron(TestSpaceum)
plt.plot(TestPlot)
plt.title("Interpolation Test Plot")
plt.show()




################        Actively dealing with    ############################





#Here I calculate the power conversion efficiency using the absorption of the absorber calculated from TMM 


EQEDarknointerp = eta*AbsDarkeV(Ephoton)
#EQEdark[Ephoton_, Eta_] := Eta AbsDarkeV[Ephoton]
#Eta is the electron-hole pair extraction efficiency. You could probably call this the internal quantum efficiency too...
EQEDarknointerpCALC = EQEDarknointerp * Ephoton**2/(np.exp((Ephoton)/(kB * Tcell)) - 1)
EQEDark = scipy.interpolate.interp1d(Ephoton, EQEDarknointerpCALC, fill_value="extrapolate")


#print(EQEDarknointerp)
plt.plot(AbsDarkeV(Ephoton))
plt.plot(EQEDark(Ephoton))
plt.show()

RR0darkHalf_Integ = scipy.integrate.quad(EQEDark, Emin, Emax) #*((2 * pi)/(SpeedOfLight**2 * h**3))
RR0darkHalf = np.multiply(RR0darkHalf_Integ, ((2 * pi)/(SpeedOfLight**2 * h**3)))


theMoops
#RR0darkHalf[Tcell_, Eta] := (2 Pi)/(SpeedOfLight^2 hPlanck^3)
   #NIntegrate[
   #EQEdark[Ephoton, Eta] Ephoton^2/(xp[(Ephoton)/(kB Tcell)] - 1), {Ephoton, Emin, Emax}, 
   #Method -> "AdaptiveQuasiMonteCarlo"]
theMoops2  
   


   
#RR0darkApprox[V_?NumericQ, Tcell_, Eta_] :=  Exp[(e V)/(kB Tcell)] RR0darkHalf[Tcell, Eta]






#GeneratedDark[Eta_] := 
 #NIntegrate[SPhotonsPerTEA[Ephoton] EQEdark[Ephoton, Eta], {Ephoton, Emin, Emax}, 
  #Method -> "AdaptiveQuasiMonteCarlo"]
 
 #CurrentDensityDarkApprox[V_?NumericQ, Tcell_, \[Eta]_] := 
 #e (GeneratedDark[\[Eta]] - RR0darkApprox[V, Tcell, \[Eta]])

#PowerDarkApprox[V_?NumericQ, Tcell_, \[Eta]_] := 
 #V* CurrentDensityDarkApprox[V, Tcell, \[Eta]]

#(* Plot[{PowerDarkApprox[V volt,298kelvin,1]/(watt meter^-2),
#PowerDarkApprox[V volt,298kelvin,0.8]/(watt meter^-2)
#},
#{V,0,1.6}, \
#Filling->Axis,ImageSize->600,PlotRange->{{0,1.6},{0,300}},\
#Frame->True,FrameLabel->{Style["Voltage (V)",16],Style["Power (W m^-2)",16]},LabelStyle->Directive[FontSize->16],PlotLegends->{\
#"10 nm", "100 nm","1000 nm (\[Eta] = 1)","1000 nm (\[Eta] = 0.8)"}]//Quiet *)
    
#MaxPowerDark[Tcell_, \[Eta]_] := 
 #FindMaximum[(V* CurrentDensityDarkApprox[V, Tcell, \[Eta]]), {V, 1 volt, 0, 
    #Emax/e}, Method -> "PrincipalAxis"][[1]]
#VAtMPPDark[Tcell_, \[Eta]_] := 
 #V /. FindMaximum[(V* CurrentDensityDarkApprox[V, Tcell, \[Eta]]), {V, 1 volt,
      #0, Emax/e}, Method -> "PrincipalAxis"][[2]]
#JAtMPPDark[Tcell_, \[Eta]_] := 
 #CurrentDensityDarkApprox[VAtMPPDark[Tcell, \[Eta]], Tcell, \[Eta]]
#MaxEfficiencyDark[Tcell_, \[Eta]_] := 
 #MaxPowerDark[Tcell, \[Eta]]/SolarConstant

#MaxPowerDark[298 kelvin, 1]/(watt meter^-2)
#VAtMPPDark[298 kelvin, 1]/volt
#JAtMPPDark[298 kelvin, 1]/(mA cm^-2)
#MaxEfficiencyDark[298 kelvin, 1]












































































