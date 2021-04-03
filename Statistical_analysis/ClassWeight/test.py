# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 14:13:17 2021

@author: claum
"""
from Weight_estimation_class import myaircraft_weight

S_wing = 657.0
W_fuel = 11000.0 
AR = 12.0
Lambda_c4 = 0.0
thick_to_chord_ratio = 0.14
q_cruise = 0.5*0.00106526*(464.15**2)
taper_ratio = 0.6 
n_ultimate = 3.5
MTOW = 50700.0
V_max = 464.15

A1 = myaircraft_weight(S_wing, W_fuel, AR, Lambda_c4, thick_to_chord_ratio,\
                 q_cruise, taper_ratio, n_ultimate, MTOW, V_max)

print("Wing weight estimate: \n", "Raymer ==> ", A1.W_wing_raymer, " in [lb]")
print("Nicolai ==> ", A1.W_wing_nicolai, " in [lb]")