# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 19:21:33 2021

@author: claum
"""
# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
# ---------------------------------------------------------------------------
class DOC(object):
    def __init__(self, MTOW, switch, haul, N_engine, k_engine, k_AFspares, TOT,\
                 USL, P_residual_ratio, p, n_pay, k_n0_ratio, k_insurance,\
                     flights_per_year, P_fuel, fuel_mass, cruise_speed, crew_number,\
                         Range, passengers):    
        """
        A first Direct Operating Costs (DOC) evaluation tool for conceptual 
        aircraft design.

        Parameters
        ----------
        MTOW   : TYPE
            Maximum Take-off weights 
        switch : INT
            Switch variable to decide correction factor inside depraciation costs
            estimation function 
            switch == 0 ---> k_correct == 1 
            switch == 1 ---> k_correct == 1.10
        haul   : INT
            Switch variable to decide P_delivery ratio 
            haul == 0 ---> Short and medium haul 
            haul == 1 ---> Long haul       
        N_engine : INT
            Number of engine installed
        k_engine : FLOAT 
            Factor to take a proportion of engine cost
        k_AFspares : FLOAT
            Factor to take a proportion of airframe cost
        TOT : FLOAT
            Take-off Thrust 
        USL : FLOAT
            Useful service life of the aircraft
        P_residual_ratio : FLOAT
            P_residual/P_total where P_residual is the residual aircraft value
            after useful service life 
            
        Partial outputs
        ---------------
        C_dep : TYPE
            Depreciation costs.
        C_int : TYPE
            Interest costs.
        C_ins : TYPE
            Insurance costs.
        C_fuel : TYPE
            Fuel usage costs.
        C_m : TYPE
            Maintainance costs.
        C_c : TYPE
            Crew costs.
        C_fee : TYPE
            Fee and charge costs.

        Returns
        -------
        Direct Operating Costs (DOC).

        """
        self.MTOW, self.switch, self.haul, self.N_engine, self.k_engine, self.k_AFspares,\
            self.TOT, self.P_residual_ratio, self.p, self.n_pay, self.k_n0_ratio,\
                self.k_insurance, self.flights_per_year, self.P_fuel, self.fuel_mass,\
                    self.cruise_speed, self.crew_number, self.Range, self.passengers = MTOW,\
                        switch, haul, N_engine, k_engine, k_AFspares, TOT, P_residual_ratio,\
                            p, n_pay, k_n0_ratio, k_insurance, flights_per_year, P_fuel,\
                                fuel_mass, cruise_speed, crew_number, Range, passengers
        self.Price   = self.total_price(MTOW,switch, haul, N_engine, k_engine, k_AFspares,TOT)
        self.C_dep   = self.depreciation_cost(self.Price[0], USL, P_residual_ratio)   
        self.C_int   = self.interest_cost(p, n_pay, k_n0_ratio, USL, self.Price[0])
        self.C_ins   = self.insurance_cost(k_insurance, self.Price[1])
        self.C_fuel  = self.fuel_cost(flights_per_year, P_fuel, fuel_mass)
        self.C_m     = self.maintainance_cost(self.Price[0], self.Price[2], N_engine)
        self.C_c     = self.crew_cost(MTOW, cruise_speed, crew_number)
        self.C_fee   = self.fee_cost(MTOW, flights_per_year, Range, passengers)
        self.DOC     = self.C_dep+self.C_int+self.C_ins+self.C_fuel+self.C_m+self.C_c+self.C_fee
        self.plot    = self.DOC_pie_chart(self.C_dep,self.C_int, self.C_ins, self.C_fuel,\
                                       self.C_m, self.C_c, self.C_fee, self.DOC)
# ---------------------------------------------------------------------------            
    def total_price(self, MTOW, switch, haul, N_engine, k_engine, k_AFspares,\
                          TOT):
        k_correct = 0
        P_delivery_ratio = 0
        # --------------------
        if switch==0:
            k_correct = 1
        elif switch==1:
            k_correct = 1.10
        # --------------------
        if haul==0:
            P_delivery_ratio == 500
        elif haul==1:
            P_delivery_ratio == 350
        # --------------------
        P_engine   = 293*(TOT)**0.81
        P_delivery = k_correct*P_delivery_ratio*MTOW
        P_airframe = P_delivery - N_engine*P_engine
        P_spares   = k_AFspares*P_airframe + k_engine*N_engine*P_engine
        return [P_delivery + P_spares, P_delivery, P_engine]
# ---------------------------------------------------------------------------   
    def depreciation_cost(self, P_total, USL, P_residual_ratio):
        # DEPRACIATION COSTS ESTIMATION 
        return (P_total/USL)*(1 - P_residual_ratio)
# ---------------------------------------------------------------------------        
    def interest_cost(self, p, n_pay, k_n0_ratio, USL, P_total):
        q = 1+p
        x = q**(n_pay) - k_n0_ratio
        y = q**(n_pay) - 1
        z = 1 - k_n0_ratio
        return (((x*(q-1)*n_pay/y)-z)/USL)*P_total
# ---------------------------------------------------------------------------        
    def insurance_cost(self, k_insurance, P_delivery):
        return k_insurance*P_delivery
# ---------------------------------------------------------------------------        
    def fuel_cost(self, flights_per_year, P_fuel, fuel_mass):
        return flights_per_year*P_fuel*fuel_mass
# ---------------------------------------------------------------------------        
    def maintainance_cost(self, P_total, P_engine, N_engine):
        # MATERIAL COSTS
        X = 3.3*((P_total-P_engine)/10**6)+14.2+(58*(P_engine/10**6)-26.1)*N_engine
        # CYCLE COSTS
        Y = 4.0*((P_total-P_engine)/10**6)+9.3+(7.5*(P_engine/10**6)+5.6)*N_engine
        return X+Y
# ---------------------------------------------------------------------------        
    def fee_cost(self, MTOW, flights_per_year, Range, passengers):
        C_fee_landing = 0.0090*MTOW*flights_per_year*0.0078
        C_fee_navig   = 0.00716*Range*np.sqrt(MTOW)*flights_per_year*0.0078
        C_fee_ground  = 182 + 7*passengers*flights_per_year*0.0078
        return C_fee_landing + C_fee_navig + C_fee_ground
# ---------------------------------------------------------------------------        
    def crew_cost(self, MTOW, cruise_speed, crew_number):
        """
        

        Parameters
        ----------
        MTOW : FLOAT
            Maximum take off weight in kg.
        cruise_speed : FLOAT
            Cruise speed in km/hrs.
        crew_number : INT
            Crew needed to operate the aircraft [2 or 3 in the cockpit].

        Returns
        -------
        Crew Costs
            Crew's payment annual costs.

        """
        x = 0
        if crew_number == 2:
            x = 74.5*(cruise_speed*(MTOW/10**5))**0.3 + 168.8
        elif crew_number == 3:
            x = 101*(cruise_speed*(MTOW/10**5))**0.3 + 273.2
        return x
# ---------------------------------------------------------------------------                   
    def DOC_pie_chart(self, C_dep, C_int, C_ins, C_fuel, C_m, C_c, C_fee, DOC):
        labels = 'Depreciation','Interest','Insurance','Fuel','Maintainance',\
            'Crew','Fee'
        sizes  = [C_dep/DOC, C_int/DOC, C_ins/DOC, C_fuel/DOC, C_m/DOC, C_c/DOC,\
                  C_fee/DOC]
        #explode = (0, 0, 0, 0)
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%',\
                shadow=False, startangle=90)
        ax1.axis('equal')
        plt.show()
