from raymer_sizing_class import raymer_sizing_class
import math 
import numpy as np
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
Endurance    = 30*60                                                                     # [s]
Range        = 4.861E6                                                                   # [ft]
W_crew       = 881.849                                                                   # [lb]
W_payload    = 21605.3                                                                   # [lb]
W_0          = 50200.0                                                                   # [lb]
# V_cruise     = 557.849409449                                                             # [ft * s^-1]
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
S_wet_vtail  = (2.0*(125.41/(0.3048**2)))*(1.0 + 0.25*thick_root*((1.0 + tau*taper_vtail)/(1.0 + taper_vtail))) 
S_wet_total  = S_wet_fus + S_wet_wing + S_wet_htail + S_wet_vtail    
S_ratio      = S_wet_total/(657)                                                         # S_wet/S_ref 
Aspect_ratio = 12                                                                        # span^2/S_ref
AR_wet       = Aspect_ratio/S_ratio                                                      # AR/S_ratio
W_guess      = 50200.0                                                                   # [lb]
# ==============================================================================================================
MyAircraft   = raymer_sizing_class(W_crew, W_payload, W_0, empty_A, empty_c, empty_m, AR_wet, V_cruise, TSFC_cruise, TSFC_loiter, Range, Endurance)
x            = W_guess
err          = abs(x - MyAircraft.W_0_calculated) 
i            = 1
# ==============================================================================================================
while err > 1E-3:
 print("###### Iterazione numero: ", i, " ######")
 MyAircraft1   = raymer_sizing_class(W_crew, W_payload, x, empty_A, empty_c, empty_m, AR_wet, V_cruise, TSFC_cruise, TSFC_loiter, Range, Endurance)
 err           = abs(x - MyAircraft1.W_0_calculated) 
 x             = MyAircraft1.W_0_calculated 
 print("====================================")
 print("W_0 calcolato: ", x, "\n") 
 i = i +1
# ==============================================================================================================
print("--------------------------------- ")
print("ARwet: ", AR_wet, "\n")
print("W_0: ", MyAircraft.W_0_calculated, "in [lb]\n")
print("S_wet: ", S_wet_total, " in [ft^2]\n")
print("(L/D)max: ", MyAircraft.lift_to_drag_r)
# ==============================================================================================================
