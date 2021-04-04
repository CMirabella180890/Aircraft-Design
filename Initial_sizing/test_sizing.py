# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 12:14:39 2020

@author: claum
"""
# ============================================================================
# ============================================================================
def knots2feetperseconds(x):
    """
    Converts speed in [kt] in [ft/s]

    Parameters
    ----------
    x : FLOAT
        Velocity in [kt].

    Returns
    -------
    V : FLOAT
        Velocity in [ft/s].

    """
    return 1.68781*x
# ============================================================================
def fpm2fps(x):
    """
    Converts speed in [fpm] in [ft/s]

    Parameters
    ----------
    x : FLOAT
        Velocity in [fpm].

    Returns
    -------
    V : FLOAT
        Velocity in [ft/s].

    """
    return x/60.0
# ============================================================================
# ============================================================================
from initial_sizing import isa_atmosphere, initial_sizing, mattingly_turboprop,\
    ratio_atmosphere
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
# ============================================================================
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# ============================================================================
from matplotlib.backends.backend_pdf import PdfPages
pp1 = PdfPages('figura1-1.pdf')
pp2 = PdfPages('figura1-2.pdf')
pp3 = PdfPages('figura1-3.pdf')
# ============================================================================
# DATA ATMOSPHERE
# ============================================================================
gamma = 1.4
R     = 1718.0    # [(ft * lb )/(slug * F°)]
T0    = 59        # [F°]
p0    = 2116      # [lb/ft^2]
rho0  = 0.002378  # [slug/ft^3] = [lb * ft^-1 * s^2]
h     = 8000.0    # [ft]
SC    = 25000.0   # [ft]
# ============================================================================
# ATMOSPHERE CALCULATIONS
# ============================================================================
my_atmo   = isa_atmosphere(T0, p0, rho0, h, gamma, R)
my_atmoSC = isa_atmosphere(T0, p0, rho0, SC, gamma, R)
# ============================================================================
# AIRCRAFT DATA 
# ============================================================================
AR        = 9.0
CDmin     = 0.025
CDTO      = 0.04
CLTO      = 0.5
mu        = 0.04
Sg        = 900.0                       # Take off length [ft]
rho1      = my_atmo.rho[0]
rho2      = my_atmo.rho[-1]
rhoSC     = my_atmoSC.rho[-1]
g         = 32.17404856                 # [slug * ft * s^-2]
V_liftoff = knots2feetperseconds(65.0)
Vmin      = knots2feetperseconds(61.0)  # [ft/s]
Vmax      = knots2feetperseconds(250.0) # [ft/s]
Vturn     = knots2feetperseconds(80.0)  # [ft/s]
V_climb   = knots2feetperseconds(135.0)
V_design  = knots2feetperseconds(160.0)
maxWS     = 60.0                        # [lb/ft^2]
nMAX      = 2.0 
maxROC    = fpm2fps(1500.0)             # foot per minute [ft/min] to [ft/s]
# ============================================================================
my_aircraft1 = initial_sizing(AR, rho1, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                              maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb)
my_aircraft2 = initial_sizing(AR, rho2, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                              maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb)
my_aircraft3 = initial_sizing(AR, rhoSC, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                              maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb)    
# =================================================================== 
fig1  = plt.figure()
plt.plot(my_aircraft1.ws, my_aircraft1.TWturn,\
         label=r"$\frac{T}{W}$ - Thrust to weight ratio for constant g turn - SL")
plt.plot(my_aircraft2.ws, my_aircraft2.TWturn,\
         label=r"$\frac{T}{W}$ - Thrust to weight ratio for constant g turn - $h_{\tiny \textup{max}}$")
plt.plot(my_aircraft1.ws, my_aircraft1.TWclimb,\
         label=r"$\frac{T}{W}$ - Thrust to weight ratio for selected max climb - SL")
plt.plot(my_aircraft2.ws, my_aircraft2.TWclimb,\
         label=r"$\frac{T}{W}$ - Thrust to weight ratio for selected max climb - $h_{\tiny \textup{max}}$")
plt.plot(my_aircraft1.ws, my_aircraft1.TWTO,\
         label=r"$\frac{T}{W} \,\, [lb/lb]$ - Thrust to weight ratio for selected Takeoff")  
plt.plot(my_aircraft2.ws, my_aircraft2.TWCS,\
         label=r"$\frac{T}{W} \,\, [lb/lb]$ - Thrust to weight ratio for selected Cruise speed")    
plt.plot(my_aircraft3.ws, my_aircraft3.TWCS,\
         label=r"$\frac{T}{W} \,\, [lb/lb]$ - Thrust to weight ratio for selected Service Ceiling")    
plt.xlabel(r'Wing loading - $\frac{W}{S} \,\, [lb/ft^2]$')  # x-label to the axes.
plt.ylabel(r'Thrust-to-weight ratio - $\frac{T}{W}$')     # y-label to the axes.
plt.title(r'Constraint diagram')                          # Title to the axes.
plt.legend(loc="upper right", prop={"size" : 6})
plt.ylim(0, 0.4)
plt.xlim(0, 60)
plt.grid(True, linestyle='-.')
plt.show()  
pp1.savefig(fig1)
pp1.close()
# =================================================================== 
# ============================================================================
# PRESTAZIONI DEL MOTORE TURBOELICA
# In questa sezione del codice si va a calcolare le prestazioni del motore 
# turboelica in termini di spinta e potenza rappresentate come una funzione 
# del numero di Mach e parametrizzate rispetto alla quota. Il modello applicato
# è proposto da Mattingly [Element of Propulsion, Gas Turbine and Rocke, 
# Capitolo 8, pag. 455 - pag. 459] e consiste nell'applicazione del concetto 
# del THROTTLE_RATIO e del così detto THETA_BREAK_RATIO. I valori di spinta e
# potenza saranno utilizzati per normalizzare le curve T/W e P/W rispetto alla
# quota zero (SL = Sea Level).
# ============================================================================
h1    = [0.0, 5000.0, 10000.0, 15000.0, 20000.0, 25000.0]    # [ft]
# ============================================================================
M_max      = 0.6
M_min      = 0.0
M          = np.linspace(M_min, M_max, 1000)
gamma      = 1.4
eta_p      = 0.8
BHP        = 2750.0
# ============================================================================
# LIST VARIABLES TO STORE VARIOUS OBJECT VARIABLES
# ============================================================================
my_atmo1   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Atmosphere's data 
my_ratio1  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Temperature and pressure ratio 
my_engine1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Engine's data 
V1         = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Design speed(h) = M_max * a(h)
F_SL       = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Nominal Sea Level Thrust at V1
thetac     = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Theta to calculate Throttle_ratio
# ============================================================================
for k in range(len(h1)):
    my_atmo1[k]   = isa_atmosphere(T0, p0, rho0, h1[k], gamma, R)
    V1[k]         = my_atmo1[k].a[-1]*M_max
    F_SL[k]       = (eta_p*550.0*BHP)/V1[k]
    my_ratio1[k]  = ratio_atmosphere(my_atmo1[k].T[-1], my_atmo1[k].p[-1],\
                                     M, gamma)
    thetac[k]     = my_ratio1[k].theta[-1]
    # ========================================================================
    # CLASSE CHE IMPLEMENTA IL METODO DI MATTINGLY
    # ========================================================================
    my_engine1[k] = mattingly_turboprop(F_SL[k], my_ratio1[k].delta,\
                                        my_ratio1[k].theta, thetac[k], M)
# =================================================================== 
Thrust_ratio = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
for l in range(len(h1)):
    Thrust_ratio[l] = my_engine1[l].F/F_SL[l] 
# ===================================================================     
fig2  = plt.figure()
plt.plot(M, Thrust_ratio[0], label=r"$h = 0\,ft$")
plt.plot(M, Thrust_ratio[1], label=r"$h = 5000\,ft$")
plt.plot(M, Thrust_ratio[2], label=r"$h = 10000\,ft$")
plt.plot(M, Thrust_ratio[3], label=r"$h = 15000\,ft$")
plt.plot(M, Thrust_ratio[4], label=r"$h = 20000\,ft$")
plt.plot(M, Thrust_ratio[5], label=r"$h = 25000\,ft$")
plt.xlabel(r'$M$ - Mach number')                         # x-label to the axes.
plt.ylabel(r"$F/F_{\textup{SL}} \,\, - \,\, [lb_{f}/lb_{f}]$ - Engine Thrust Ratio") # y-label to the axes.
plt.title(r'Turboprop Engine thrust ratio')              # Title to the axes.
plt.legend(loc="upper right", prop={"size" : 6})
#plt.ylim(0, 0.4)
plt.xlim(0, 0.6)
plt.grid(True, linestyle='-.')
plt.show()  
pp2.savefig(fig2)
pp2.close()
# ===================================================================  
# ===================================================================     
# =================================================================== 
Pow_ratio = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
for m in range(len(h1)):
    Pow_ratio[m] = Thrust_ratio[m]*((V1[m])/(eta_p*550.0))
# ===================================================================  
fig3  = plt.figure()
plt.plot(M, Pow_ratio[0], label=r"$h = 0\,ft$")
plt.plot(M, Pow_ratio[1], label=r"$h = 5000\,ft$")
plt.plot(M, Pow_ratio[2], label=r"$h = 10000\,ft$")
plt.plot(M, Pow_ratio[3], label=r"$h = 15000\,ft$")
plt.plot(M, Pow_ratio[4], label=r"$h = 20000\,ft$")
plt.plot(M, Pow_ratio[5], label=r"$h = 25000\,ft$")
plt.xlabel(r'$M$ - Mach number')                         # x-label to the axes.
plt.ylabel(r"$P/P_{\textup{SL}}\,\, - \,\, [SHP/SHP]$ - Engine Power ratio")  # y-label to the axes.
plt.title(r'Turboprop Engine power ratio')               # Title to the axes.
plt.legend(loc="upper right", prop={"size" : 6})
#plt.ylim(0, 0.4)
plt.xlim(0, 0.6)
plt.grid(True, linestyle='-.')
plt.show()  
pp3.savefig(fig3)
pp3.close()
# ===================================================================   