# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 13:46:00 2021

@author: claum
"""
# ============================================================================
import numpy as np 
# ============================================================================
class myaircraft_weight(object):
    def __init__(self, S_wing, W_fuel, AR, Lambda_c4, thick_to_chord_ratio,\
                 q_cruise, taper_ratio, n_ultimate, MTOW, V_max, Lambda_htail,\
                     S_htail, l_htail, taper_htail, AR_htail, b_htail, thick_htail,\
                        Lambda_vtail, AR_vtail, F_tail, taper_vtail, S_vtail, thick_vtail,\
                           S_fuselage, l_fuselage, depth_fuselage, DeltaP, Press_vol,\
                              max_width, max_depth, MLW, landing_gear_length,\
                                 nose_landing_gear_lenght, W_without_engine, N_engine):
        """
        A class to quickly estimate aircraft's weight

        Parameters
        ----------
        S_wing : FLOAT
            Wing area in squared feet [ft**2].
        W_fuel : FLOAT
            Weight of fuel in pounds [lb_f].
        AR : FLOAT
            Wing aspect ratio.
        Lambda_c4 : FLOAT
            Wing sweep at 25% MAC.
        thick_to_chord_ratio : FLOAT
            Wing thickness to chord ratio.
        q_cruise : FLOAT
            Dynamic pressure at cruise.
        taper_ratio : FLOAT
            Wing taper ratio.
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f].
        V_max : FLOAT
            Maximum level airspeed at Sea Level in KEAS [kn].
        Lambda_htail : FLOAT
            Wing sweep of horizontal tail at 25% MGC.
        S_htail : FLOAT
            Horizontal area in squared feet [ft**2].
        l_htail : FLOAT
            Horizontal tail arm from wing 25% MAC to the 25% * chord_htail in feet [ft].
        taper_htail : FLOAT
            Horizontal tail taper ratio.
        AR_htail : FLOAT
            Horizontal tail Aspect ratio.  
        b_htail : FLOAT
            Horizontal tail span in feet [ft].
        thick_htail : FLOAT
            Max root chord thickness of Horizontal tail in inches [in].
        Lambda_vtail : FLOAT
            Vertical tail sweep at 25% MGC.
        AR_vtail : FLOAT
            Vertical tail aspect ratio.
        F_tail : FLOAT
            Vertical tail factor | F_tail = 0 for conventional tail| 
                                 | F_tail = 1 for T tail           | 
        taper_vtail : FLOAT
            Vertical tail taper ratio.
        S_vtail : FLOAT
            Vertical tail area in squared feet [ft**2].
        thick_vtail : FLOAT
            Max root chord thickness of vertical tail in inches [in].
        S_fuselage : FLOAT
            Fuselage wetted area in squared feet [ft**2].
        l_fuselage : FLOAT
            Length of fuselage structure in feet [ft].
        depth_fuselage : FLOAT
            Depth of fuselage structure in feet [ft].
        DeltaP : FLOAT
            Cabin pressure differential in psi (typical value: 8 psi).
        Press_vol : FLOAT
            Volume of pressurized cabin section in cubic feet [ft**3].
        max_width : FLOAT
            Fuselage max width in feet [ft].
        max_depth : FLOAT
            Fuselage max depth in feet [ft]. 
        MLW : FLOAT
            Maximum landing weight
        landing_gear_length : FLOAT
            Lenght of the main landing gear strut in inches [in].
        nose_landing_gear_length : FLOAT
            Lenght of the nose landing gear strut in inches [in].
        W_without_engine : FLOAT
            MTOW - W_engine in pounds [lb].
        N_engine : FLOAT
            Number of engine. 

        Returns
        -------
        Aircraft's weights collected inside a class.

        """
# ============================================================================
        self.Se_wing, self.W_fuel, self.AR, self.Lambda_c4, self.thick_to_chord_ratio,\
            self.q_cruise, self.taper_ratio, self.n_ultimate, self.MTOW,\
                self.V_max, self.Lambda_htail, self.S_htail, self.taper_htail,\
                    self.AR_htail, self.b_htail, self.thick_htail, self.Lambda_vtail,\
                        self.AR_vtail, self.F_tail, self.taper_vtail, self.S_vtail,\
                            self.thick_vtail, self.S_fuselage, self.l_fuselage, self.depth_fuselage,\
                                self.DeltaP, self.Press_vol, self.max_width, self.max_depth, self.MLW,\
                                    self.landing_gear_length, self.nose_landing_gear_lenght, self.W_without_engine,\
                                       self.N_engine = S_wing, W_fuel, AR, Lambda_c4, thick_to_chord_ratio, q_cruise,\
                                          taper_ratio, n_ultimate, MTOW, V_max, Lambda_htail, S_htail, taper_htail,\
                                             AR_htail, b_htail, thick_htail, S_fuselage, l_fuselage, depth_fuselage,\
                                                DeltaP, Press_vol, max_width, max_depth, MLW, landing_gear_length,\
                                                   nose_landing_gear_length, W_without_engine, N_engine
# ============================================================================
        self.W_wing_raymer = self.raymer_wing_weight(S_wing, W_fuel, AR, Lambda_c4,\
                                              thick_to_chord_ratio, q_cruise,\
                                                  taper_ratio, n_ultimate, MTOW)
        self.W_wing_nicolai = self.nicolai_wing_weight(S_wing, W_fuel, AR, Lambda_c4, thick_to_chord_ratio,\
                           q_cruise, taper_ratio, n_ultimate, MTOW, V_max)
# ============================================================================
# ============================================================================
        self.W_htail_raymer = self.raymer_htail_weight(n_ultimate, MTOW, q_cruise, S_htail,\
                            thick_to_chord_ratio, Lambda_htail, AR_htail, taper_htail)
        self.W_htail_nicolai = self.nicolai_htail_weight(n_ultimate, MTOW, S_htail,\
                             AR_htail,l_htail, thick_htail)
# ============================================================================
# ============================================================================
        self.W_vtail_raymer = self.raymer_vtail_weight(F_tail, n_ultimate, MTOW,\
                                                      thick_to_chord_ratio, Lambda_vtail, AR_vtail, \
                                                      taper_vtail, q_cruise)
        self.W_vtail_nicolai = self.nicolai_vtail_weight(n_ultimate, MTOW,\
                                                        S_vtail, AR_vtail, thick_vtail)
# ============================================================================
# ============================================================================
        self.W_fus_raymer = self.raymer_fus_weight(S_fuselage, n_ultimate, MTOW,\
                                                  l_fuselage, depth_fuselage, q_cruise,\
                                                  Press_vol, DeltaP, l_htail)
        self.W_fus_nicolai = self.nicolai_fus_weight(n_ultimate, MTOW, l_fuselage,\
                                                     max_depth, max_width, V_max)
# ============================================================================
# ============================================================================
        self.W_main_land_gear_raymer = self.raymer_main_landing_gear_weight(n_ultimate, MLW, landing_gear_length)
        self.W_main_land_gear_nicolai = self.nicolai_main_landing_gear_weight(n_ultimate, MLW, landing_gear_length)
# ============================================================================
# ============================================================================
        self.W_nose_land_gear_raymer = self.raymer_nose_landing_gear_weight(n_ultimate, MLW, nose_landing_gear_length)
# ============================================================================
# ============================================================================
        self.W_engine_raymer = self.raymer_engine_weight(W_without_engine, N_engine)
        self.W_engine_nicolai = self.nicolai_engine_weight(W_without_engine, N_engine)
# ============================================================================
    # ========================================================================
    # ++++++++++++++++++++++++ WING WEIGHT ESTIMATE ++++++++++++++++++++++++++
    # ========================================================================                
    def raymer_wing_weight(self, S_wing, W_fuel, AR, Lambda_c4, thick_to_chord_ratio,\
                           q_cruise, taper_ratio, n_ultimate, MTOW):
        """
        Function that estimates wing weight. 
        
        Parameters
        ----------
        S_wing : FLOAT
            Wing area in squared feet [ft**2].
        W_fuel : FLOAT
            Weight of fuel in pounds [lb_f].
        AR : FLOAT
            Wing aspect ratio. 
        Lambda_c4 : FLOAT
            Wing sweep at 25% MAC.
        thick_to_chord_ratio : FLOAT
            Wing thickness to chord ratio.
        q_cruise : FLOAT
            Dynamic pressure at cruise.
        taper_ratio : FLOAT
            Wing taper ratio.
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f]
        
        Returns
        ------- 
        Wing weight in pounds [lb].
        """
        x = 0.036*(S_wing**0.758)
        y = W_fuel**0.0035
        z = (AR/(np.cos(Lambda_c4))**2)**0.6
        k = q_cruise**0.006
        l = taper_ratio**0.04
        m = (100*thick_to_chord_ratio/np.cos(Lambda_c4))**(-0.3)
        return x*y*z*k*l*m*(n_ultimate*MTOW)**0.49
    # ========================================================================    
    def nicolai_wing_weight(self, S_wing, W_fuel, AR, Lambda_c4, thick_to_chord_ratio,\
                           q_cruise, taper_ratio, n_ultimate, MTOW, V_max):
        """
        Function that estimates wing weight. 
        
        Parameters
        ----------
        S_wing : FLOAT
            Wing area in squared feet [ft**2].
        W_fuel : FLOAT
            Weight of fuel in pounds [lb_f].
        AR : FLOAT
            Wing aspect ratio. 
        Lambda_c4 : FLOAT
            Wing sweep at 25% MAC.
        thick_to_chord_ratio : FLOAT
            Wing thickness to chord ratio.
        q_cruise : FLOAT
            Dynamic pressure at cruise.
        taper_ratio : FLOAT
            Wing taper ratio.
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f]
        V_max : FLOAT
            Maximum level airspeed at Sea Level in KEAS [kn]. 
        
        Returns
        ------- 
        Wing weight in pounds [lb].
        """
        x = (n_ultimate*MTOW/1E5)**0.65
        y = (AR/(np.cos(Lambda_c4))**2)**0.57
        z = (S_wing/1E2)**0.61
        k = ((1 + taper_ratio)/(2*thick_to_chord_ratio))**0.36
        l = np.sqrt(1 +(V_max/5E2))
        return 96.948*(x*y*z*k*l)**0.993
    # ========================================================================    
    # ========================================================================
    # ++++++++++++++++++ HORIZONTAL TAIL WEIGHT ESTIMATE +++++++++++++++++++++
    # ========================================================================                
    def raymer_htail_weight(self, n_ultimate, MTOW, q_cruise, S_htail,\
                            thick_to_chord_ratio, Lambda_htail, AR_htail, taper_htail):
        """
        Function that estimates horizontal tail weight. 
        
        Parameters
        ----------
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f].
        q_cruise : FLOAT
            Dynamic pressure at cruise.
        thick_to_chord_ratio : FLOAT
            Wing thickness to chord ratio. 
        S_htail : FLOAT
            Horizontal area in squared feet [ft**2].
        Lambda_htail : FLOAT
            Wing sweep of horizontal tail at 25% MGC.
        AR_htail : FLOAT
            Vertical tail aspect ratio.
        taper_htail : FLOAT
            Horizontal tail taper ratio.
        
        Returns
        ------- 
        Horizontal tail weight in pounds [lb].
        """
        x = 0.016*(n_ultimate*MTOW)**0.414
        y = q_cruise**0.168
        z = S_htail**0.896
        k = 1/(((100*thick_to_chord_ratio)/(np.cos(Lambda_htail)))**(0.12))
        l = ((AR_htail)/((np.cos(Lambda_htail))**2))**0.043
        m = 1/(taper_htail**0.02)
        return x*y*z*k*l*m
    # ========================================================================    
    def nicolai_htail_weight(self, n_ultimate, MTOW, S_htail,\
                             AR_htail, l_htail, thick_htail):
        """
        Function that estimates horizontal tail weight. 
        
        Parameters
        ----------
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f]. 
        S_htail : FLOAT
            Horizontal tail area in squared feet [ft**2].
        l_htail : FLOAT
            Horizontal tail arm from wing 25% MAC to the 25% * chord_htail in feet [ft].
        AR_htail : FLOAT
            Vertical tail aspect ratio.
        thick_htail : FLOAT
            Max root chord thickness of Horizontal tail in inches [in].
        
        Returns
        ------- 
        Horizontal tail weight in pounds [lb].
        """
        x = (n_ultimate*MTOW/1E5)**0.87
        y = (S_htail/1E2)**1.2
        z = (l_htail/1E1)**0.483
        b = np.sqrt(AR_htail*S_htail)
        s = np.sqrt(b/thick_htail)
        return 127*(x*y*z*s)**0.458
    # ========================================================================  
    # ========================================================================
    # +++++++++++++++++++ VERTICAL TAIL WEIGHT ESTIMATE ++++++++++++++++++++++
    # ========================================================================                
    def raymer_vtail_weight(self, F_tail, n_ultimate, MTOW,\
                           thick_to_chord_ratio, Lambda_vtail, AR_vtail,
                           taper_vtail, q_cruise):
        """
        Function that estimates vertical tail weight. 
        
        Parameters
        ----------
        F_tail : FLOAT
            Vertical tail factor | F_tail = 0 for conventional tail| 
                                 | F_tail = 1 for T tail           | 
        n_ultimate : FLOAT
            Ultimate load factor. 
        MTOW : FLOAT
            Design gross weight in pounds [lb_f]. 
        thick_to_chord_ratio : FLOAT
            Wing thickness to chord ratio.
        Lambda_vtail : FLOAT
            Vertical tail sweep at 25% MGC.
        AR_vtail : FLOAT
            Vertical tail aspect ratio.
        taper_vtail : FLOAT
            Vertical tail taper ratio.
        q_cruise : FLOAT
            Dynamic pressure at cruise.
        
        Returns
        ------- 
        Vertical tail weight in pounds [lb].
        """
        x = 0.073*(1 + 0.2*F_tail)
        y = (n_ultimate*MTOW)**0.376
        z = q_cruise**0.122
        s = S_vtail**0.873
        k = 1/(((100*thick_to_chord_ratio)/(np.cos(Lambda_htail)))**(0.49))
        l = ((AR_htail)/((np.cos(Lambda_htail))**2))**0.357
        m = taper_htail**0.039
        return x*y*z*s*k*l*m
    # ========================================================================    
    def nicolai_vtail_weight(self, n_ultimate, MTOW, S_vtail, AR_vtail,\
                            thick_vtail):
        """
        Function that estimates vertical tail weight. 
        
        Parameters
        ----------
        n_ultimate : FLOAT
            Ultimate load factor. 
        MTOW : FLOAT
            Design gross weight in pounds [lb_f]. 
        S_vtail : FLOAT
            Vertical tail area in squared feet [ft**2].
        AR_vtail : FLOAT
            Vertical tail aspect ratio.
        thick_vtail : FLOAT
            Max root chord thickness of vertical tail in inches [in].
        
        Returns
        ------- 
        Vertical tail weight in pounds [lb].
        """
        x = (n_ultimate*MTOW/1E5)**0.87
        y = (S_htail/1E2)**1.2
        b = np.sqrt(AR_vtail*S_vtail)
        s = np.sqrt(b/thick_vtail)
        return 98.5*(x*y*s)
    # ========================================================================        
    # ========================================================================
    # ++++++++++++++++++++++ FUSELAGE WEIGHT ESTIMATE ++++++++++++++++++++++++
    # ========================================================================                
    def raymer_fus_weight(self, S_fuselage, n_ultimate, MTOW, l_fuselage, depth_fuselage, q_cruise,\
                         Press_vol, DeltaP, l_htail):
        """
        Function that estimates fuselage weight. 
        
        Parameters
        ----------
        S_fuselage : FLOAT
            Fuselage wetted area in squared feet [ft**2]. 
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f].
        l_fuselage : FLOAT
            Length of fuselage structure in feet [ft].
        depth_fuselage : FLOAT
            Depth of fuselage structure in feet [ft].
        q_cruise : FLOAT
            Dynamic pressure at cruise.  
        Press_vol : FLOAT
            Volume of pressurized cabin section in cubic feet [ft**3]. 
        DeltaP : FLOAT
            Cabin pressure differential in psi (typical value: 8 psi).   
        l_htail : FLOAT
            Horizontal tail arm from wing 25% MAC to the 25% * chord_htail in feet [ft].
          
        Returns
        ------- 
        Fuselage weight in pounds [lb].
        """
        x = 0.052*(S_fuselage**1.086)
        y = (n_ultimate*MTOW)**0.177
        z = q_cruise**0.241
        s = 1/(l_htail**0.051)
        k = 1/((l_fuselage/depth_fuselage)**0.072)
        l = 11.9*(Press_vol*DeltaP)**0.271
        return x*y*z*s*k*l
    # ========================================================================    
    def nicolai_fus_weight(self, n_ultimate, MTOW, l_fuselage, max_depth, max_width, V_max):
        """
        Function that estimates fuselage weight. 
        
        Parameters
        ---------- 
        n_ultimate : FLOAT
            Ultimate load factor.
        MTOW : FLOAT
            Design gross weight in pounds [lb_f].
        l_fuselage : FLOAT
            Length of fuselage structure in feet [ft].
        max_depth : FLOAT
            Fuselage max depth in feet [ft].
        max_width : FLOAT
            Fuselage max width in feet [ft].
        V_max : FLOAT
            Maximum level airspeed at Sea Level in KEAS [kn].  
          
        Returns
        ------- 
        Fuselage weight in pounds [lb].
        """
        x = (n_ultimate*MTOW/1E5)**0.286
        y = (l_fuselage/1E1)**0.857
        b = (max_depth + max_width)/1E1
        s = (V_max/1E2)**0.338
        return 200*(x*y*b*s)**1.1
    # ========================================================================  
    # ========================================================================
    # ++++++++++++++++ MAIN LANDING GEAR WEIGHT ESTIMATE +++++++++++++++++++++
    # ========================================================================                
    def raymer_main_landing_gear_weight(self, n_ultimate, MLW, landing_gear_length):
        """
        Function that estimates main landing gear weight. 
        
        Parameters
        ---------- 
        n_ultimate : FLOAT
            Ultimate load factor.
        MLW : FLOAT
            Design gross weight in pounds [lb_f].
        landing_gear_length : FLOAT
            Length of fuselage structure in feet [ft].  
          
        Returns
        ------- 
        Main landing gear weight in pounds [lb].
        """
        x = (n_ultimate*MLW)**0.768
        y = (landing_gear_length/12)**0.409
        return 0.095*x*y
    # ========================================================================    
    def nicolai_main_landing_gear_weight(self, n_ultimate, MLW, landing_gear_length):
        """
        Function that estimates main landing gear weight. 
        
        Parameters
        ---------- 
        n_ultimate : FLOAT
            Ultimate load factor.
        MLW : FLOAT
            Maximum landing weight in pounds [lb_f].
        landing_gear_length : FLOAT
            Lenght of the main landing gear strut in inches [in].  
          
        Returns
        ------- 
        Main landing gear weight in pounds [lb].
        """
        x = (n_ultimate*MTOW)**0.684
        y = (landing_gear_length/12)**0.601
        return 0.054*x*y
    # ======================================================================== 
    # ========================================================================
    # ++++++++++++++++ NOSE LANDING GEAR WEIGHT ESTIMATE +++++++++++++++++++++
    # ========================================================================                
    def raymer_nose_landing_gear_weight(self, n_ultimate, MLW, nose_landing_gear_length):
        """
        Function that estimates nose landing gear weight. 
        
        Parameters
        ---------- 
        n_ultimate : FLOAT
            Ultimate load factor.
        MLW : FLOAT
            Maximum landing weight in pounds [lb_f].
        nose_landing_gear_length : FLOAT
            Lenght of the nose landing gear strut in inches [in].  
          
        Returns
        ------- 
        Nose landing gear weight in pounds [lb].
        """
        x = (n_ultimate*MLW)**0.566
        y = (nose_landing_gear_length/12)**0.845
        return 0.125*x*y
    # ========================================================================
    # ========================================================================
    # ++++++++++++++++++++++ ENGINE WEIGHT ESTIMATE ++++++++++++++++++++++++++
    # ========================================================================                
    def raymer_engine_weight(self, W_without_engine, N_engine):
        """
        Function that estimates nose landing gear weight. 
        
        Parameters
        ---------- 
        W_without_engine : FLOAT
            Maximum landing weight in pounds [lb_f].
        nose_landing_gear_length : FLOAT
            Lenght of the nose landing gear strut in inches [in].  
          
        Returns
        ------- 
        Nose landing gear weight in pounds [lb].
        """
        x = W_without_engine**0.922
        return 2.575*x*N_engine
    # ========================================================================    
    def nicolai_engine_weight(self, W_without_engine, N_engine):
        """
        Function that estimates nose landing gear weight. 
        
        Parameters
        ---------- 
        W_without_engine : FLOAT
            Maximum landing weight in pounds [lb_f].
        nose_landing_gear_length : FLOAT
            Lenght of the nose landing gear strut in inches [in].  
          
        Returns
        ------- 
        Nose landing gear weight in pounds [lb].
        """
        x = W_without_engine**0.922
        return 2.575*x*N_engine
    # ========================================================================         
