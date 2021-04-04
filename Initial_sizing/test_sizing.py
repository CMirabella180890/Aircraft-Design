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
    ratio_atmosphere, maximum_lift_coefficient
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import json 
from types import SimpleNamespace
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
# =======================================
JSONFileName = "atmosphere_input.json"
with open(JSONFileName, "r") as f:
    # ===================================
    # ATMOSPHERE DATA INSIDE A DICTIONARY
    # ===================================
    dict_atmosphere_data = json.load(f)
    # ===================================
# ===================================================
#   ATMOSPHERE DATA INSIDE A SIMPLE OBJECT VARIABLE
# ===================================================    
atmosphere_data = SimpleNamespace(**dict_atmosphere_data) 
# ===================================================    
gamma = atmosphere_data.GAMMA
R     = atmosphere_data.R                 # [(ft * lb )/(slug * F°)]
T0    = atmosphere_data.T_0               # [F°]
p0    = atmosphere_data.P_0               # [lb/ft^2]
rho0  = atmosphere_data.RHO_0             # [slug/ft^3] = [lb * ft^-1 * s^2]
h     = atmosphere_data.SELECTED_ALTITUDE # [ft]
SC    = atmosphere_data.SERVICE_CEILING   # [ft]
# ============================================================================
# ATMOSPHERE CALCULATIONS
# ============================================================================
my_atmo   = isa_atmosphere(T0, p0, rho0, h, gamma, R)
my_atmoSC = isa_atmosphere(T0, p0, rho0, SC, gamma, R)
# ============================================================================
# AIRCRAFT DATA 
# ============================================================================
# ===================================================
JSONFileName1 = "aircraft_input.json"
with open(JSONFileName1, "r") as f:
    # =================================
    # AIRCRAFT DATA INSIDE A DICTIONARY
    # =================================
    dict_aircraft_data = json.load(f)
    # =================================
# ===================================================
# AIRCRAFT DATA INSIDE A SIMPLE OBJECT VARIABLE
# ===================================================    
aircraft_data = SimpleNamespace(**dict_aircraft_data) 
# ===================================================
AR        = aircraft_data.ASPECT_RATIO
CDmin     = aircraft_data.CD_MIN
CDTO      = aircraft_data.CD_TAKEOFF
CLTO      = aircraft_data.CL_TAKEOFF
mu        = aircraft_data.GROUND_FRICTION_COEFFICIENT
Sg        = aircraft_data.GROUND_TAKEOFF_DISTANCE     # Take off length [ft]
rho1      = my_atmo.rho[0]
rho2      = my_atmo.rho[-1]
rhoSC     = my_atmoSC.rho[-1]
g         = aircraft_data.g                               # [slug * ft * s^-2]
V_liftoff = knots2feetperseconds(aircraft_data.V_LIFTOFF) # [ft/s]
Vmin      = knots2feetperseconds(aircraft_data.V_MIN)     # [ft/s]
Vmax      = knots2feetperseconds(aircraft_data.V_MAX)     # [ft/s]
Vturn     = knots2feetperseconds(aircraft_data.V_TURN)    # [ft/s]
V_climb   = knots2feetperseconds(aircraft_data.V_CLIMB)   # [ft/s]
V_design  = knots2feetperseconds(aircraft_data.V_DESIGN)  # [ft/s]
maxWS     = aircraft_data.maxWS                           # [lb/ft^2]
nMAX      = aircraft_data.nMAX 
maxROC    = fpm2fps(aircraft_data.maxROC)             # foot per minute [ft/min] to [ft/s]
# ============================================================================
my_aircraft1 = initial_sizing(AR, rho1, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                              maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb)
my_aircraft2 = initial_sizing(AR, rho2, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                              maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb)
my_aircraft3 = initial_sizing(AR, rhoSC, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                              maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb)    
# ============================================================================ 
# Create some mock data
t = np.arange(0.01, 10.0, 0.01)
data1 = np.exp(t)
data2 = np.sin(2 * np.pi * t)
V_STALL = np.linspace(knots2feetperseconds(59),knots2feetperseconds(80),5)
Max_Lift_Coeff = maximum_lift_coefficient(rho1, V_STALL, maxWS)
# ============================================================================
fig, ax1 = plt.subplots()
#color = 'tab:red'
ax1.set_xlabel(r'Wing loading - $\frac{W}{S} \,\, [lb/ft^2]$')
ax1.set_ylabel(r'Thrust-to-weight ratio - $\frac{T}{W}$')
ax1.plot(my_aircraft1.ws, my_aircraft1.TWturn,\
         label=r"$\frac{T}{W}$ - Constant g turn - SL",
         linestyle="solid")
ax1.plot(my_aircraft2.ws, my_aircraft2.TWturn,\
         label=r"$\frac{T}{W}$ - Constant g turn - $h_{\tiny \textup{max}}$",
         linestyle="solid")
ax1.plot(my_aircraft1.ws, my_aircraft1.TWclimb,\
         label=r"$\frac{T}{W}$ - Selected max climb - SL",
         linestyle="dashed")
ax1.plot(my_aircraft2.ws, my_aircraft2.TWclimb,\
         label=r"$\frac{T}{W}$ - Selected max climb - $h_{\tiny \textup{max}}$",
         linestyle="dashed")
ax1.plot(my_aircraft1.ws, my_aircraft1.TWTO,\
         label=r"$\frac{T}{W} \,\, [lb/lb]$ - Selected Takeoff",
         linestyle="solid")  
ax1.plot(my_aircraft2.ws, my_aircraft2.TWCS,\
         label=r"$\frac{T}{W} \,\, [lb/lb]$ - Selected Cruise speed",
         linestyle="dashdot")    
ax1.plot(my_aircraft3.ws, my_aircraft3.TWCS,\
         label=r"$\frac{T}{W} \,\, [lb/lb]$ - Selected Service Ceiling",
         linestyle="dotted")
ax1.legend(loc="upper right", prop={"size" : 5})
ax1.grid(True, linestyle='-.')
ax1.set_ylim((0, 0.45))
ax1.tick_params(axis='y')
# ============================================================================
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#color = 'tab:black'
ax2.set_ylabel('$C_L Max$')  # we already handled the x-label with ax1
ax2.plot(Max_Lift_Coeff.ws, Max_Lift_Coeff.CL_MAX[0])
ax2.plot(Max_Lift_Coeff.ws, Max_Lift_Coeff.CL_MAX[1])
ax2.plot(Max_Lift_Coeff.ws, Max_Lift_Coeff.CL_MAX[2])
ax2.plot(Max_Lift_Coeff.ws, Max_Lift_Coeff.CL_MAX[3])
ax2.plot(Max_Lift_Coeff.ws, Max_Lift_Coeff.CL_MAX[4])
ax2.set_ylim((0, 2.0))
ax2.tick_params(axis='y')
# ============================================================================
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title(r'Constraint diagram')                          # Title to the axes.
plt.xlim(0, 60)
# plt.grid(True, linestyle='-.')
plt.show()
pp1.savefig(fig)
pp1.close()
# ============================================================================
# OLD VERSION CONSTRAINTS DIAGRAM 
# fig1  = plt.figure()
# plt.plot(my_aircraft1.ws, my_aircraft1.TWturn,\
#          label=r"$\frac{T}{W}$ - Thrust to weight ratio for constant g turn - SL",
#          linestyle="solid")
# plt.plot(my_aircraft2.ws, my_aircraft2.TWturn,\
#          label=r"$\frac{T}{W}$ - Thrust to weight ratio for constant g turn - $h_{\tiny \textup{max}}$",
#          linestyle="solid")
# plt.plot(my_aircraft1.ws, my_aircraft1.TWclimb,\
#          label=r"$\frac{T}{W}$ - Thrust to weight ratio for selected max climb - SL",
#          linestyle="dashed")
# plt.plot(my_aircraft2.ws, my_aircraft2.TWclimb,\
#          label=r"$\frac{T}{W}$ - Thrust to weight ratio for selected max climb - $h_{\tiny \textup{max}}$",
#          linestyle="dashed")
# plt.plot(my_aircraft1.ws, my_aircraft1.TWTO,\
#          label=r"$\frac{T}{W} \,\, [lb/lb]$ - Thrust to weight ratio for selected Takeoff",
#          linestyle="solid")  
# plt.plot(my_aircraft2.ws, my_aircraft2.TWCS,\
#          label=r"$\frac{T}{W} \,\, [lb/lb]$ - Thrust to weight ratio for selected Cruise speed",
#          linestyle="dashdot")    
# plt.plot(my_aircraft3.ws, my_aircraft3.TWCS,\
#          label=r"$\frac{T}{W} \,\, [lb/lb]$ - Thrust to weight ratio for selected Service Ceiling",
#          linestyle="dotted")    
# plt.xlabel(r'Wing loading - $\frac{W}{S} \,\, [lb/ft^2]$')  # x-label to the axes.
# plt.ylabel(r'Thrust-to-weight ratio - $\frac{T}{W}$')     # y-label to the axes.
# plt.title(r'Constraint diagram')                          # Title to the axes.
# plt.legend(loc="upper right", prop={"size" : 6})
# plt.ylim(0, 0.4)
# plt.xlim(0, 60)
# plt.grid(True, linestyle='-.')
# plt.show()  
# pp1.savefig(fig1)
# pp1.close()
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
# =======================================
JSONFileName3 = "engine_input.json"
with open(JSONFileName3, "r") as f:
    # ===================================
    #   ENGINE DATA INSIDE A DICTIONARY
    # ===================================
    dict_engine_data = json.load(f)
    # ===================================
# ===================================================
#     ENGINE DATA INSIDE A SIMPLE OBJECT VARIABLE
# ===================================================    
engine_data = SimpleNamespace(**dict_engine_data) 
# ===================================================
h1    = engine_data.h1    # [ft]
# ============================================================================
M_max = engine_data.M_MAX
M_min = engine_data.M_MIN
M     = np.linspace(M_min, M_max, 1000)
eta_p = engine_data.ETA_P
BHP   = engine_data.BRAKE_HORSE_POWER
# ============================================================================
# LIST VARIABLES TO STORE VARIOUS OBJECT VARIABLES
# ============================================================================
my_atmo1   = engine_data.ATMO1             # Atmosphere's data 
my_ratio1  = engine_data.RATIO1            # Temperature and pressure ratio 
my_engine1 = engine_data.ENGINE1           # Engine's data 
V1         = engine_data.V1                # Design speed(h) = M_max * a(h)
F_SL       = engine_data.SEA_LEVEL_THRUST  # Nominal Sea Level Thrust at V1
thetac     = engine_data.DESIGN_TEMP_RATIO # Theta to calculate Throttle_ratio
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
Thrust_ratio = engine_data.THRUST_RATIO
for l in range(len(h1)):
    Thrust_ratio[l] = my_engine1[l].F/F_SL[l] 
# ===================================================================     
fig2  = plt.figure()
plt.plot(M, Thrust_ratio[0], label=r"$h = 0\,ft$",     linestyle="solid")
plt.plot(M, Thrust_ratio[1], label=r"$h = 5000\,ft$",  linestyle="dashed")
plt.plot(M, Thrust_ratio[2], label=r"$h = 10000\,ft$", linestyle="dotted")
plt.plot(M, Thrust_ratio[3], label=r"$h = 15000\,ft$", linestyle="dashdot")
plt.plot(M, Thrust_ratio[4], label=r"$h = 20000\,ft$", linestyle="solid")
plt.plot(M, Thrust_ratio[5], label=r"$h = 25000\,ft$", linestyle="dashed")
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
Pow_ratio = engine_data.POWER_RATIO
for m in range(len(h1)):
    Pow_ratio[m] = Thrust_ratio[m]*((V1[m])/(eta_p*550.0))
# ===================================================================  
fig3  = plt.figure()
plt.plot(M, Pow_ratio[0], label=r"$h = 0\,ft$",     linestyle="solid")
plt.plot(M, Pow_ratio[1], label=r"$h = 5000\,ft$",  linestyle="dashed")
plt.plot(M, Pow_ratio[2], label=r"$h = 10000\,ft$", linestyle="dotted")
plt.plot(M, Pow_ratio[3], label=r"$h = 15000\,ft$", linestyle="dashdot")
plt.plot(M, Pow_ratio[4], label=r"$h = 20000\,ft$", linestyle="solid")
plt.plot(M, Pow_ratio[5], label=r"$h = 25000\,ft$", linestyle="dashed")
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