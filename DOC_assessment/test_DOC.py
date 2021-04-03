# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 20:04:32 2021

@author: claum
"""
from DOC_calc import DOC
# SAMPLE CALCULATION
# ------------------------
MTOW = 78000
switch = 1
haul = 1
N_engine = 2
k_engine = 0.08
k_AFspares = 0.08
TOT = 13607.771
USL = 15
P_residual_ratio = 0.1
p = 0.08 # Fokker 1993
n_pay = 15 # Fokker 1993
k_n0_ratio = 0.1
k_insurance = 0.004
flights_per_year = 1200 
P_fuel = 0.22
fuel_mass = 22000
cruise_speed = 908 # km/hr
crew_number = 2
Range = 6000
passengers = 150
# ------------------------
myAircraftcosts = DOC(MTOW, switch, haul, N_engine, k_engine, k_AFspares, TOT,\
                      USL, P_residual_ratio, p, n_pay, k_n0_ratio, k_insurance,\
                          flights_per_year, P_fuel, fuel_mass, cruise_speed,\
                              crew_number, Range, passengers)