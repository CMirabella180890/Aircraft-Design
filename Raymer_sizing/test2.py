from raymer_sizing_class import raymer_sizing_class
import math 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
# ===================================================================
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# ===================================================================
stat1  = PdfPages('figura17.pdf')
stat2  = PdfPages('figura18.pdf')
stat3  = PdfPages('figura19.pdf')
stat4  = PdfPages('figura20.pdf')
stat5  = PdfPages('figura21.pdf')
# ===================================================================
# ==============================================================================================================
# =============================================== USEFUL FUNCTIONS =============================================
# ==============================================================================================================
def meters_to_feet(x):
 return x/0.3048
# ==============================================================================================================
# ==============================================================================================================
def feet_to_meters(x):
 return x*0.3048
# ==============================================================================================================
empty_A      = 0.96
empty_c      = -0.05
empty_m      = 0.92
Endurance    = (30.0)*(60.0)                                                             # [s]
Range        = 4.861E6                                                                   # [ft]
W_crew       = 881.849                                                                   # [lb]
W_payload    = 21605.1                                                                   # [lb]
W_0          = 50200.0                                                                   # [lb]
# V_cruise     = 557.849409449                                                           # [ft * s^-1]
V_cruise     = 550                                                                       # [ft * s^-1]
eta_prop     = 0.8
PSFC         = 0.41                                                                      # [lb * (hr * bhp)^-1]
TSFC_cruise  = PSFC/(3600*550)*(V_cruise/eta_prop)                                       # [s^-1] 
TSFC_loiter  = PSFC/(3600*550)*(0.75*V_cruise/eta_prop)                                  # [s^-1] 
D_fus        = meters_to_feet(2.74)                                                      # [ft]
l_fus        = meters_to_feet(27.2)                                                      # [ft]
slend_ratio  = l_fus/D_fus
S_wet_fus    = math.pi*D_fus*l_fus*((1 - 2/slend_ratio)**(2/3))*(1 + 1/(slend_ratio**2)) # [ft^2]
S_wet_wing   = 125.41/(0.3048**2)                                                        # [ft^2]
taper_htail  = 0.6
taper_vtail  = 0.6 
thick_root   = 0.12  
thick_tip    = 0.09
tau          = thick_tip/thick_root  
S_wet_htail  = (2.0*(125.41/(0.3048**2)))*(1.0 + 0.25*thick_root*((1.0 + tau*taper_htail)/(1.0 + taper_htail)))      
S_wet_vtail  = 0.6*(2.0*(125.41/(0.3048**2)))*(1.0 + 0.25*thick_root*((1.0 + tau*taper_vtail)/(1.0 + taper_vtail))) 
S_wet_total  = S_wet_fus + S_wet_wing + S_wet_htail + S_wet_vtail    
S_ratio      = S_wet_total/(S_wet_wing*0.5)                                              # S_wet/S_ref 
Aspect_ratio = 12                                                                        # span^2/S_ref
AR_wet       = Aspect_ratio/S_ratio                                                      # AR/S_ratio
Eff_LD       = 13.5
# ==============================================================================================================
W_guess      = 50200.0                                                                                    # [lb]
Range        = np.linspace(0.0, 4.861E6, 200)
l            = len(Range)
W_0_calc     = np.zeros(l)
# ==============================================================================================================
# ================================================= TRADE STUDIES ==============================================
# ==============================================================================================================
# First graph: Effect of the range requirements variations on the MTOW
# ==============================================================================================================
for j in range(l): 
 R            = Range[j]
 MyAircraft2  = raymer_sizing_class(W_crew,W_payload,W_0,empty_A,empty_c,empty_m,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,R,Endurance, Eff_LD)
 err          = abs(W_guess - MyAircraft2.W_0_calculated)
 x            = MyAircraft2.W_0_calculated
 i            = 1
# ==============================================================================================================
 while err > 1E-3:
  print("###### Iterazione numero: ", i, " ######")
  MyAircraft3 = raymer_sizing_class(W_crew,W_payload,x,empty_A,empty_c,empty_m,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,R,Endurance, Eff_LD)
  err           = abs(x - MyAircraft3.W_0_calculated) 
  x             = MyAircraft3.W_0_calculated 
  print("====================================")
  print("W_0 calcolato: ", x, "\n") 
  i = i +1

 W_0_calc[j] = x
# ==============================================================================================================
# ===================================================================
fig1  = plt.figure()
plt.plot(Range/1E6, W_0_calc, color="red")
plt.xlabel(r'\textsc{range} - $[ft]\times 10^6$')              # x-label to the axes.
plt.ylabel(r'\textsc{mtow} - $[lb]$')               # y-label to the axes.
plt.title(r'Trade study - Range variation effect') 
plt.grid(True, linestyle='-.')
plt.show()
stat1.savefig(fig1)
stat1.close()
# ===================================================================
# ==============================================================================================================
# ==============================================================================================================
# Second graph: Effect of the loiter requirements variations on the MTOW
# ==============================================================================================================
# ==============================================================================================================
W_guess      = 50200.0                                                                                    # [lb]
Endurance    = np.linspace(0.0, 60.0*60.0, 200)
k            = len(Endurance)
W_0_calc1    = np.zeros(k)
# ==============================================================================================================
for j in range(k): 
 E            = Endurance[j]
 MyAircraft4  = raymer_sizing_class(W_crew,W_payload,W_0,empty_A,empty_c,empty_m,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,Range[199],E, Eff_LD)
 err          = abs(W_guess - MyAircraft4.W_0_calculated)
 x            = MyAircraft4.W_0_calculated
 i            = 1
# ==============================================================================================================
 while err > 1E-3:
  print("###### Iterazione numero: ", i, " ######")
  MyAircraft5 = raymer_sizing_class(W_crew,W_payload,x,empty_A,empty_c,empty_m,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,Range[199],E, Eff_LD)
  err           = abs(x - MyAircraft5.W_0_calculated) 
  x             = MyAircraft5.W_0_calculated 
  print("====================================")
  print("W_0 calcolato: ", x, "\n") 
  i = i +1

 W_0_calc1[j] = x

print(W_0_calc1)
# ==============================================================================================================
# ===================================================================
fig2  = plt.figure()
plt.plot(Endurance/60, W_0_calc1, color="red")
plt.xlabel(r'\textsc{endurance} - $[min]$')   # x-label to the axes.
plt.ylabel(r'\textsc{mtow} - $[lb]$')       # y-label to the axes.
#plt.ylim(0, 60)
#plt.ylim(50000, 65300.5)
plt.title(r'Trade study - Endurance variation effect') 
plt.grid(True, linestyle='-.')
plt.show()
stat2.savefig(fig2)
stat2.close()
# ===================================================================
# ==============================================================================================================
# ==============================================================================================================
# Third graph: Effect of the payload requirements variations on the MTOW
# ==============================================================================================================
# ==============================================================================================================
W_guess      = 50200.0      
Endurance    = 30*60                                                                              # [lb]
Payload      = np.linspace(0.0, W_payload, 200)
k            = len(Payload)
W_0_calc2    = np.zeros(k)
# ==============================================================================================================
for s in range(k): 
 P            = Payload[s]
 MyAircraft6  = raymer_sizing_class(W_crew,P,W_0,empty_A,empty_c,empty_m,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,Range[199],Endurance, Eff_LD)
 err          = abs(W_guess - MyAircraft6.W_0_calculated)
 x            = MyAircraft6.W_0_calculated
 i            = 1
# ==============================================================================================================
 while err > 1E-3:
  print("###### Iterazione numero: ", i, " ######")
  MyAircraft7 = raymer_sizing_class(W_crew,P,x,empty_A,empty_c,empty_m,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,Range[199],Endurance, Eff_LD)
  err           = abs(x - MyAircraft7.W_0_calculated) 
  x             = MyAircraft7.W_0_calculated 
  print("====================================")
  print("W_0 calcolato: ", x, "\n") 
  i = i +1

 W_0_calc2[s] = x

print(W_0_calc2)
# ==============================================================================================================
# ===================================================================
fig3  = plt.figure()
plt.plot(Payload, W_0_calc2, color="red")
plt.xlabel(r'\textsc{payload weight} - $[lb]$')   # x-label to the axes.
plt.ylabel(r'\textsc{mtow} - $[lb]$')       # y-label to the axes.
#plt.ylim(0, 60)
#plt.ylim(50000, 65300.5)
plt.title(r'Trade study - Payload variation effect') 
plt.grid(True, linestyle='-.')
plt.show()
stat3.savefig(fig3)
stat3.close()
# ===================================================================
# ==============================================================================================================
# ==============================================================================================================
# ==============================================================================================================
# Fourth graph: Effect of the composites empty - weight fraction on the MTOW
# ==============================================================================================================
# ==============================================================================================================
W_guess      = 50200.0                                                                                   # [lb]
Composites   = np.linspace(0.75, empty_m, 200)
k            = len(Composites)
W_0_calc3    = np.zeros(k)
# ==============================================================================================================
for f in range(k): 
 C            = Composites[f]
 MyAircraft8  = raymer_sizing_class(W_crew,P,W_0,empty_A,empty_c,C,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,Range[199],Endurance, Eff_LD)
 err          = abs(W_guess - MyAircraft8.W_0_calculated)
 x            = MyAircraft8.W_0_calculated
 i            = 1
# ==============================================================================================================
 while err > 1E-3:
  print("###### Iterazione numero: ", i, " ######")
  MyAircraft9 = raymer_sizing_class(W_crew,P,x,empty_A,empty_c,C,AR_wet,V_cruise,TSFC_cruise,TSFC_loiter,Range[199],Endurance, Eff_LD)
  err           = abs(x - MyAircraft9.W_0_calculated) 
  x             = MyAircraft9.W_0_calculated 
  print("====================================")
  print("W_0 calcolato: ", x, "\n") 
  i = i +1

 W_0_calc3[f] = x

print(W_0_calc3)
# ==============================================================================================================
# ===================================================================
fig4  = plt.figure()
plt.plot(Composites, W_0_calc3, color="red")
plt.xlabel(r'\textsc{Composites percentage}')   # x-label to the axes.
plt.ylabel(r'\textsc{mtow} - $[lb]$')       # y-label to the axes.
#plt.ylim(0, 60)
#plt.ylim(50000, 65300.5)
plt.title(r'Trade study - Composites weight fraction variation effect') 
plt.grid(True, linestyle='-.')
plt.show()
stat4.savefig(fig4)
stat4.close()
# ===================================================================
# ==============================================================================================================
# ==============================================================================================================
# Fifth graph: Effect of cruise speed variation on the MTOW
# ==============================================================================================================
# ==============================================================================================================
W_guess      = 50200.0                                                                                   # [lb]
CruiseSpeed  = np.linspace(400.0, V_cruise, 200)
k            = len(CruiseSpeed)
W_0_calc4    = np.zeros(k)
# ==============================================================================================================
for r in range(k): 
 V            = CruiseSpeed[r]
 MyAircraft10  = raymer_sizing_class(W_crew,W_payload,W_0,empty_A,empty_c,empty_m,AR_wet,V,TSFC_cruise,TSFC_loiter,Range[199],Endurance, Eff_LD)
 err          = abs(W_guess - MyAircraft10.W_0_calculated)
 x            = MyAircraft10.W_0_calculated
 i            = 1
# ==============================================================================================================
 while err > 1E-3:
  print("###### Iterazione numero: ", i, " ######")
  MyAircraft11 = raymer_sizing_class(W_crew,W_payload,x,empty_A,empty_c,empty_m,AR_wet,V,TSFC_cruise,TSFC_loiter,Range[199],Endurance, Eff_LD)
  err           = abs(x - MyAircraft11.W_0_calculated) 
  x             = MyAircraft11.W_0_calculated 
  print("====================================")
  print("W_0 calcolato: ", x, "\n") 
  i = i +1

 W_0_calc4[r] = x

print(W_0_calc4)
# ==============================================================================================================
# ===================================================================
fig5  = plt.figure()
plt.plot(CruiseSpeed, W_0_calc4, color="red")
plt.xlabel(r'\textsc{Cruise speed}- $[ft\cdot s^{-1}]$')   # x-label to the axes.
plt.ylabel(r'\textsc{mtow} - $[lb]$')       # y-label to the axes.
#plt.ylim(0, 60)
#plt.ylim(50000, 65300.5)
plt.title(r'Trade study - Cruise speed variation effect') 
plt.grid(True, linestyle='-.')
plt.show()
stat5.savefig(fig5)
stat5.close()
# ===================================================================
# ==============================================================================================================
