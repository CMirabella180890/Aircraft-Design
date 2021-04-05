# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 11:37:52 2020

A class for initial Aircraft sizing

@author: claum
"""
# ============================================================================
import numpy as np 
import math 
# ============================================================================
# STANDARD ATMOSPHERE
# ============================================================================
class isa_atmosphere(object):
    """
    A class that calculates standard atmosphere in Engineering units.
    """
    def __init__(self, T0, p0, rho0, h, gamma, R):
        """
        Inizialization of class isa_atmosphere. In this application
        engineering units will be used [ft - lb - s].

        Parameters
        ----------
        T0 : FLOAT
            ISA Temperature [F°] (at SL = Sea Level).
        p0 : FLOAT
            ISA Pressure [lb/sqft] (at SL = Sea Level).
        p0 : FLOAT
            ISA Density [slug/cubic ft] (at SL = Sea Level).
        h : FLOAT
            Array of altitudes.
        gamma : FLOAT
            Air specific heat ratio.
        R : FLOAT
            Perfect gas constant evaluated for air. 
            **** VALUE USED INSIDE EQUATION **** 
            R = 1718.0 [(ft * lbf)/(slug * R°)]
            When used inside Standard Atmosphere equations, temperature
            must be converted to Rankyne deg using the following equation:
            T[R°] = T[F°] + 459.7 

        Returns
        -------
        OBJECT VARIABLE. For example:
            
                 my_ratio = isa_atmosphere(T0, p0, rho0, h, gamma, R)
                     
        Inside this variable will be stored:
            
        T   ==> Temperature [F°]
        p   ==> Pressure [lb/sqft] 
        rho ==> Density [slug/cubic ft]
        a   ==> Speed-of-sound [ft/s^2]

        """
        self.T0, self.p0, self.rho0, self.h, self.gamma, self.R = T0, p0,\
            rho0, h, gamma, R
        self.x   = np.linspace(0, self.h, 1000)
        self.T   = self.temperature(self.T0, self.x)
        self.p   = self.pressure(self.p0, self.T, self.x)
        self.rho = self.density(self.p, self.T)
        self.a   = self.speed_of_sound(gamma, self.T, R)
    # ========================================================================        
    def temperature(self, T0, x):
        """
        ISA Temperature in [F°]

        Parameters
        ----------
        T0 : FLOAT
            ISA Temperature [F°] (at SL = Sea Level).
        x : FLOAT
            Array which contains altitude values.

        Returns
        -------
        T : FLOAT
            ISA Temperature in [F°].

        """
        return self.T0 - 0.00356*(self.x)
    # ========================================================================    
    def pressure(self, p0, T, x):
        """
        ISA Pressure in [lb/sqft]

        Parameters
        ----------
        p0 : FLOAT
            ISA Pressure [lb/sqft] (at SL = Sea Level).
        x : FLOAT
            Array which contains altitude values.

        Returns
        -------
        p : FLOAT
            ISA Pressure in [lb/sqft].

        """    
        return self.p0*((self.T + 459.7)/518.6)**(5.2561)
    # ========================================================================    
    def density(self, p, T):
        """
        ISA Density in [slug/cubic ft]

        Parameters
        ----------
        rho0 : FLOAT
            ISA Density [slug/cubic ft] (at SL = Sea Level).
        x : FLOAT
            Array which contains altitude values.

        Returns
        -------
        rho : FLOAT
            ISA Density in [slug/cubic ft].

        """        
        return self.p/(1718.0*(self.T + 459.7))
    # ========================================================================
    def speed_of_sound(self, gamma, T, R): 
        """
        ISA Speed-of-sound [ft/s^2].

        Parameters
        ----------
        gamma : FLOAT
            Air specific heat ratio.
        T : FLOAT
            ISA Temperature in [F°].
        R : FLOAT
            Perfect gas constant evaluated for air. 
            **** VALUE USED INSIDE EQUATION **** 
            R = 1718.0 [(ft * lbf)/(slug * R°)]
            When used inside Standard Atmosphere equations, temperature
            must be converted to Rankyne deg using the following equation:
            T[R°] = T[F°] + 459.7 

        Returns
        -------
        aa : FLOAT
            ISA Speed-of-sound [ft/s^2] calculated as:
                aa = sqrt(gamma * (T[F°] + 459.7) * R)

        """
        aa = np.zeros(len(T))
        for i in range(len(T)):
            aa[i] = np.sqrt(gamma*(T[i] + 459.7)*R)
        return aa
    # ========================================================================    
# ============================================================================
# STANDARD ATMOSPHERE
# ============================================================================
# ============================================================================
# ATMOSPHERE RATIO
# ============================================================================
class ratio_atmosphere(object):
    """
    A class that calculates standard atmosphere ratio as a function of 
    Mach number.
    """
    def __init__(self, T, p, M, gamma):
        """
        Initialization of ratio_atmosphere class.

        Parameters
        ----------
        T : FLOAT
            ISA Temperature [F°] (at h in [ft]).
        p : FLOAT
            ISA Pressure [lb/sqft] (at h in [ft]).
        M : FLOAT
            An array which contains achievable flight Mach numbers.
        gamma : FLOAT
            Air Specific heat ratio.

        Returns
        -------
        OBJECT VARIABLE. For example:
            
                     my_ratio = ratio_atmosphere(T, p, M, gamma)
                     
        Inside this variable will be stored:
            
        theta  ==> Temperature ratio: T/T_SL
        delta  ==> Pressure ratio: p/p_SL 
         
        """
        self.T, self.p, self.M, self.gamma = T, p, M, gamma
        self.x         = np.zeros(len(M))
        for i in range(len(M)):
            self.x[i]     = self.partial_p(gamma, M[i]) # 1 + 0.5*(gamma - 1)*M[i]**2 # 
        self.T0    = 59
        self.p0    = 2116
        self.theta = np.zeros(len(M))
        self.delta = np.zeros(len(M))
        for j in range(len(M)):
            self.theta[j] = abs(self.theta_ratio(T, self.T0, self.x[j]))
            self.delta[j] = abs(self.delta_ratio(p, self.p0, self.x[j], gamma))  
    # ====================================================================    
    def theta_ratio(self, T, T0, x):
        """
        A function that calculates Temperature ratio: T/T_SL

        Parameters
        ----------
        T : FLOAT
            ISA Temperature [F°] (at h in [ft]).
        T0 : FLOAT
            ISA Temperature [F°] (at SL = Sea Level).
        x : FLOAT
            Common factor related to gamma and Mach numbers:
                      1 + 0.5*(gamma - 1)*(M**2)

        Returns
        -------
        theta : FLOAT
            Temperature ratio: T/T_SL.

        """
        d = T/T0
        return x*d
    # ====================================================================
    def delta_ratio(self, p, p0, x, gamma):
        """
        A function that calculates Pressure ratio: p/p_SL

        Parameters
        ----------
        p : FLOAT
            ISA Pressure [lb/sqft] (at h in [ft]).
        p0 : FLOAT
            ISA Pressure [lb/sqft] (at SL = Sea Level).
        x : FLOAT
            Common factor related to gamma and Mach numbers:
                      1 + 0.5*(gamma - 1)*(M**2)
        gamma : FLOAT    
            Air Specific heat ratio. 

        Returns
        -------
        delta : FLOAT
            Pressure ratio: p/p_SL.

        """
        d = p/p0
        y = gamma/(gamma-1)
        z = x**y
        return z*d
    # ====================================================================
    def partial_p(self, gamma, M):
        """
            Common factor related to gamma and Mach numbers:
                      1 + 0.5*(gamma - 1)*(M**2)        

        Parameters
        ----------
        gamma : FLOAT
            Air Specific heat ratio.
        M : FLOAT
            Reference Mach numbers.

        Returns
        -------
        zz : FLOAT
            Common factor inside atmosphere ratios.

        """
        xx = M**2
        yy = gamma - 1
        zz = 1 + 0.5*yy*xx
        return zz
    # ====================================================================            
# ============================================================================
# ATMOSPHERE RATIO
# ============================================================================
# ============================================================================
# MAXIMUM LIFT COEFFICIENT
# ============================================================================
class maximum_lift_coefficient(object):
    """
    Class to evaluate aircraft maximum lift coefficient as a function of 
    wing loading W/S (expressed in engineering units [lb/sqft])
    """
    def __init__(self, rho, V_STALL, maxWS):
        """
        Inizialization of the maximum_lift_coefficient class.

        Parameters
        ----------
        rho : FLOAT
            Density in [slug/cubic ft].
        V_STALL : FLOAT
            Array or list which contains selected stall speed values.
        maxWS : FLOAT
            Maximum wing loading W/S to build an array.

        Returns
        -------
        Maximum lift coefficient CL_MAX.

        """
        # ======================================================
        self.rho, self.V_STALL, self.maxWS = rho, V_STALL, maxWS
        # ======================================================
        self.V_SELECTED   = [0.0]*len(V_STALL)
        for i in range(len(V_STALL)):
            self.V_SELECTED[i] = V_STALL[i]
        self.q_STALL      = [0.0]*len(V_STALL)    
        for j in range(len(V_STALL)):
            self.q_STALL[j] = self.dynamic_press(rho, V_STALL[j])
        self.ws      = self.wing_loading(maxWS)
        self.CL_MAX  = [0.0]*len(V_STALL)
        for k in range(len(V_STALL)):
            self.CL_MAX[k] = (1/self.q_STALL[k])*self.ws
        # ======================================================
    def dynamic_press(self, rho, V):
        """
        Function that calculates a vector containing dynamic pressure values

        Parameters
        ----------
        rho : FLOAT
            Density at prescribed flight conditions.
        V : FLOAT
            Speed at prescribed flight conditions.

        Returns
        -------
        q : FLOAT
            Dynamic pressure (array or scalar) [lb/ft^2].

        """
        return 0.5*rho*V**2
    # ========================================================================
    # ======================================================================== 
    def wing_loading(self, maxWS):
        """
        Function that calculates wing loading values for plot

        Parameters
        ----------
        maxWS : FLOAT
            Max wing loading W/S [lb/ft^2].

        Returns
        -------
        ws : FLOAT
            A vector with wing loading values from 0 to MAX(W/S).

        """
        return np.linspace(5, maxWS, 1000)
    # ========================================================================    
# ============================================================================
# MAXIMUM LIFT COEFFICIENT
# ============================================================================
# ============================================================================
# INITIAL SIZING
# ============================================================================
class initial_sizing(object):
    """
    A class to assess engine power and wing surface of an aircraft.
    """
    def __init__(self, AR, rho, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                 maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb, eta_p):
        """
        Initialization of INITIAL_SIZING class.

        Parameters
        ----------
        AR : FLOAT
            Aircraft aspect ratio.
        rho : FLOAT
            Density at prescribed flight conditions.
        Vmax : FLOAT
            Maximum possible speed.
        Vmin : FLOAT
            Minimum achievable speed; it can be the aircrat stall speed.
        maxWS : FLOAT
            Max wing loading W/S [lb/ft^2].
        nMAX : FLOAT
            Maximum applicable load factor.
        CDmin : FLOAT
            Minimum drag coefficient.
        CDTO : FLOAT
            Drag coefficient at Take-Off conditions.
        maxROC : FLOAT
            Maximum Rate-of-Climb.
        g : FLOAT
            Gravitational acceleration:
            g = 32.17404856 [slug * ft * s^-2].
        mu : FLOAT
            Ground surface friction coefficient.
        V_liftoff : FLOAT
            Take-Off speed.
        CLTO : FLOAT
            Lift coefficient at Take-Off conditions.
        Sg : FLOAT
            Ground Take-Off distance in feet.
        V_design : FLOAT
            Design speed to calculate turn performances.
        V_climb : FLOAT
            Reference Horizontal speed for climb performances calculations.
        eta_p  : FLOAT 
            Propeller efficiency.

        Returns
        -------
        OBJECT VARIABLE. For example:
            
                     my_aircraft = initial_sizing(AR, ..., V_climb)
                     
        Inside this variable will be stored:
            
        TWturn  ==> Thrust-to-Weight ratio for constant turn
        TWclimb ==> Thrust-to-Weight ratio to achieve a prescribed 
                    Rate-of-Climb
        TWTO    ==> Thrust-to-Weight ratio to achieve a prescribed 
                    Take-Off Ground distance  
        TWCS    ==> Thrust-to-Weight ratio to achieve a prescribed 
                    Cruise speed
        TWSC    ==> Thrust-to-Weight ratio to achieve a prescribed 
                    Service ceiling         
        """
        # ====================================================================
        self.AR, self.rho, self.Vmax, self.Vmin, self.maxWS, self.nMAX,\
        self.CDmin, self.CDTO, self.maxROC, self.g, self.mu, self.V_liftoff,\
        self.CLTO, self.Sg, self.V_design, self.V_climb, self.eta_p = AR, rho,\
        Vmax, Vmin, maxWS, nMAX, CDmin, CDTO, maxROC, g, mu, V_liftoff, CLTO,\
        Sg, V_design, V_climb, eta_p
        # ====================================================================
        # THRUST-TO-WEIGHT RATIO CALCULATIONS 
        # ====================================================================
        self.V       = np.linspace(Vmin, Vmax, 1000)
        self.e       = self.efficiency_factor(AR)
        self.MaxROC  = 1.66667 # ROC in [ft/s]
        self.q       = self.dynamic_press(rho, self.V)
        self.qLO     = self.dynamic_press(rho, V_liftoff/np.sqrt(2.0))
        self.qd      = self.dynamic_press(rho, V_design)
        self.qclimb  = self.dynamic_press(rho, V_climb)
        self.qcruise = self.dynamic_press(rho, self.V[-1])
        self.ws      = self.wing_loading(maxWS)
        self.k       = self.lift_induced_drag_constant(AR, self.e)
        self.TWturn  = self.constant_turn(self.qd, CDmin, self.ws,\
                                          self.k, nMAX)
        self.TWclimb = self.rate_of_climb_TW_ratio(self.qclimb, CDmin,\
                                                   self.ws, self.k, nMAX,\
                                                       maxROC, V_climb)
        self.TWTO    = self.take_off_TW_ratio(self.qLO, CDTO, CLTO, g, Sg, self.ws,\
                                              mu, V_liftoff)
        self.TWCS    = self.cruise_speed_TW_ratio(self.qcruise, CDmin, self.k,\
                                                  self.ws)
        self.TWSC    = self.service_ceiling_TW_ratio(self.MaxROC, CDmin,\
                                                     self.k, self.ws, rho)
        # ====================================================================
        # POWER-TO-WEIGHT RATIO CALCULATIONS 
        # ====================================================================
        self.PWturn  = (self.TWturn*V_design)/(eta_p*550.0)
        self.PWclimb = (self.TWclimb*V_climb)/(eta_p*550.0)
        self.PWTO    = (self.TWTO*(V_liftoff/(np.sqrt(2.0))))/(eta_p*550.0)
        self.PWCS    = (self.TWclimb*Vmax)/(eta_p*550.0)
        self.PWSC    = (self.TWclimb*V_climb)/(eta_p*550.0)
    # ========================================================================        
    def efficiency_factor(self, AR):
        """
        Function that calculates Oswald efficiency factor

        Parameters
        ----------
        AR : FLOAT
            Aircraft aspect ratio.

        Returns
        -------
        e : FLOAT
            Oswald's efficiency factor.

        """
        x = 1.78
        y = 1 - 0.045*((self.AR)**0.68)
        return x*y - 0.64
    # ======================================================================== 
    def lift_induced_drag_constant(self, AR, e):
        """
        Function that calculates the lift induced drag constant
               1
        k = -------
            pi*AR*e
            
        Parameters
        ----------
        AR : FLOAT
            Aspect ratio.
        e : TYPE
            Oswald efficiency factor.

        Returns
        -------
        k : FLOAT
            Lift - induced drag constant.

        """
        return 1/(math.pi*AR*e) 
    # ======================================================================== 
    def feet2meters(self, h):
        """
        Converts h in feet to h in meters

        Parameters
        ----------
        h : TYPE
            Heigth in feet.

        Returns
        -------
        TYPE
            Height in meters.

        """
        return (1200/3937)*h
    # ======================================================================== 
    def dynamic_press(self, rho, V):
        """
        Function that calculates a vector containing dynamic pressure values

        Parameters
        ----------
        rho : FLOAT
            Density at prescribed flight conditions.
        V : FLOAT
            Speed at prescribed flight conditions.

        Returns
        -------
        q : FLOAT
            Dynamic pressure (array or scalar) [lb/ft^2].

        """
        return 0.5*rho*V**2
    # ======================================================================== 
    def knots2feetperseconds(self, x):
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
    # ======================================================================== 
    def wing_loading(self, maxWS):
        """
        Function that calculates wing loading values for plot

        Parameters
        ----------
        maxWS : FLOAT
            Max wing loading W/S [lb/ft^2].

        Returns
        -------
        ws : FLOAT
            A vector with wing loading values from 0 to MAX(W/S).

        """
        return np.linspace(5, maxWS, 1000)
    # ======================================================================== 
    def constant_turn(self, q, CDmin, ws, k, nMAX):
        """
        Function that calculates Thrust-to-weight ratio T/W for constant 
        turn at a specified load factor nMAX.

        Parameters
        ----------
        q : FLOAT
            Dynamic pressure evaluated in flight at a prescribed 
            horizontal speed and altitude.
        CDmin : FLOAT
            Minimum drag coefficient.
        ws : FLOAT
            An array containing wing loading values [lb/sqft].
        k : FLOAT
            Induced drag efficiency factor [pi * AR * e]^-1.
        nMAX : FLOAT
            Maximum applicable load factor.

        Returns
        -------
        zz : FLOAT
            T/W required to perform a constant turn.

        """
        xx = np.zeros(len(ws))
        yy = np.zeros(len(ws))
        zz = np.zeros(len(ws))
        for i in range(len(ws)):
            xx[i] = (q)*(CDmin)*(1/ws[i])
            yy[i] = (k)*(q)*(ws[i])*((nMAX/q)**2)
            zz[i] = xx[i] + yy[i]
            
        return zz
    # ======================================================================== 
    def rate_of_climb_TW_ratio(self, q, CDmin, ws, k, nMAX, maxROC, V):
        """
        A function to calculate T/W ratio required for a certain ROC
        at a prescribed horizontal speed and altitude.        

        Parameters
        ----------
        q : FLOAT
            Dynamic pressure evaluated in flight at a prescribed 
            horizontal speed and altitude.
        CDmin : FLOAT
            Minimum drag coefficient.
        ws : FLOAT
            An array containing wing loading values [lb/sqft].
        k : FLOAT
            Induced drag efficiency factor [pi * AR * e]^-1.
        nMAX : FLOAT
            Maximum applicable load factor.
        maxROC : FLOAT
            Maximum Rate-of-Climb.
        V : FLOAT
            Reference Horizontal speed.

        Returns
        -------
        kk : FLOAT
            T/W to achieve a certain ROC at a prescribed horizontal
            speed and altitude.

        """
        xx = np.zeros(len(ws))
        yy = np.zeros(len(ws))
        zz = np.zeros(len(ws))
        kk = np.zeros(len(ws))
        for i in range(len(ws)):
            xx[i] = (q)*(CDmin)*(1/ws[i])
            yy[i] = (ws[i])*(k)*(1/q)
            zz[i] = maxROC/V
            kk[i] = xx[i] + yy[i] + zz[i]

        return kk
    # ======================================================================== 
    def take_off_TW_ratio(self, q, CDTO, CLTO, g, Sg, ws, mu, V_liftoff):
        """
        A function to calculate T/W ratio required for a certain Ground
        Take-Off distance.        

        Parameters
        ----------
        q : FLOAT
            Dynamic pressure evaluated at Take-Off conditions.
        CDTO : FLOAT
            Drag coefficient at Take-Off conditions.
        CLTO : FLOAT
            Lift coefficient at Take-Off conditions.
        g : FLOAT
            Gravitational acceleration:
            g = 32.17404856 [slug * ft * s^-2].
        Sg : FLOAT
            Ground Take-Off distance in feet.
        ws : FLOAT
            An array containing wing loading values [lb/sqft].
        mu : FLOAT
            Ground surface friction coefficient.
        V_liftoff : FLOAT
            Take-Off speed.

        Returns
        -------
        kk : FLOAT
            T/W ratio to achieve a prescribed ground Take-Off distance.

        """
        xx = ((V_liftoff)**2)/(2*g*Sg)
        yy = np.zeros(len(ws))
        zz = np.zeros(len(ws))
        kk = np.zeros(len(ws))
        for i in range(len(ws)):
            yy[i] = (1/ws[i])*(CDTO)*(q)
            zz[i] = 1 - (1/ws[i])*(CLTO)*(q)
            kk[i] = xx + yy[i] + mu*zz[i]

        return kk 
    # ========================================================================  
    def cruise_speed_TW_ratio(self, q, CDmin, k, ws):
        """
        A function to calculate T/W ratio required for a certain cruise
        speed and a prescribed cruise altitude.        

        Parameters
        ----------
        q : FLOAT
            Dynamic pressure evaluated at cruise speed and cruise
            altitude.
        CDmin : FLOAT
            Minimum drag coefficient.
        k : FLOAT
            Induced drag efficiency factor [pi * AR * e]^-1.
        ws : FLOAT
            An array containing wing loading values [lb/sqft].
            
        Returns
        -------
        kk : FLOAT
            T/W ratio to achieve a certain cruise speed at the
            prescribed cruise altitude.

        """
        xx = np.zeros(len(ws))
        yy = np.zeros(len(ws))
        zz = np.zeros(len(ws))
        kk = np.zeros(len(ws))
        for i in range(len(ws)):
            xx[i] = q*CDmin*(1/ws[i])
            yy[i] = k/q
            zz[i] = yy[i]*ws[i]
            kk[i] = xx[i] + zz[i]
        
        return kk 
    # ========================================================================
    # ========================================================================  
    def service_ceiling_TW_ratio(self, MaxROC, CDmin, k, ws, rho):
        """
        A function to calculate T/W ratio required for a certain service 
        ceiling.

        Parameters
        ----------
        MaxROC : FLOAT
            Maximum Rate-of-Climb at service ceiling as prescribed by 
            CS 23/FAR 23 which is equal to 100 fpm [1.66667 ft/s].
        CDmin : FLOAT
            Minimum drag coefficient.
        k : FLOAT
            Induced drag efficiency factor [pi * AR * e]^-1.
        ws : FLOAT
            An array containing wing loading values [lb/sqft].
        rho : FLOAT
            Density at service ceiling [slug/cubic ft].

        Returns
        -------
        kk : FLOAT
            T/W ratio required for a certain prescribed service ceiling.

        """
        xx = np.zeros(len(ws))
        yy = np.zeros(len(ws))
        zz = np.zeros(len(ws))
        kk = np.zeros(len(ws))
        for i in range(len(ws)):
            xx[i] = np.sqrt(k/(3*CDmin))
            yy[i] = np.sqrt((2/rho)*(ws[i])*(xx[i]))
            zz[i] = MaxROC/yy[i]
            kk[i] = zz[i] + 4*zz[i]
        
        return kk 
    # ========================================================================
# ============================================================================
class mattingly_turboprop(object):
    """
    A class to apply Mattingly method to calculate thrust of a turboprop
    """
    # ========================================================================
    def __init__(self, F_SL, delta0, theta0, thetac, M):
        """
        Initialization of the Mattingly's method

        Parameters
        ----------
        F_SL : FLOAT
            Nominal Sea Level Thrust at design speed V1.
        delta0 : FLOAT
            ISA Pressure ratio evaluated at h[ft].
        theta0 : FLOAT
            ISA Temperature ratio evaluated at h[ft].
        thetac : FLOAT
            Reference THETA to obtain TR = Throttle ratio.
        M : FLOAT
            An array of achievable flight Mach numbers.

        Returns
        -------
        F/F_SL : FLOAT
            Engine Thrust ratio F/F_SL as a function of Mach number. 

        """
        # ====================================================================
        self.F_SL, self.delta0, self.theta0, self.thetac, self.M, = F_SL,\
        delta0, theta0, thetac, M
        # ====================================================================    
        self.T0 = 59                            # ISA Temperature at Sea Level 
        self.TR = self.throttle_ratio(thetac, self.T0)
        TR      = self.TR
        # ====================================================================
        # THRUST CALCULATION
        # ====================================================================
        self.F  = np.zeros(len(M))
        for i in range(len(M)):
            if (M[i] <= 0.1):
                self.F[i] = self.func1(F_SL, delta0[i])
            elif (M[i] > 0.1 and theta0[i] <= TR):
                self.F[i] = self.func2(F_SL, delta0[i], M[i])
            elif (M[i] > 0.1 and theta0[i] > TR):
                self.F[i] = self.func3(F_SL, delta0[i], theta0[i], M[i], TR)
        # ====================================================================    
    def throttle_ratio(self, thetac, T0):
        """
        A function to calculate TR = Throttle ratio

        Parameters
        ----------
        thetac : FLOAT
            Reference THETA to obtain TR = Throttle ratio.
        T0 : FLOAT
            ISA Temperature at Sea Level in [F°].

        Returns
        -------
        TR : FLOAT
            Engine Throttle ratio.

        """
        T  = thetac*T0
        TR = T/T0
        return TR          
    # ========================================================================    
    def func1(self, F_SL, delta0):
        """
        Thrust selection 1. 

        Parameters
        ----------
        F_SL : FLOAT
            Nominal Sea Level Thrust at design speed V1.
        delta0 : FLOAT
            ISA Pressure ratio evaluated at h[ft].

        Returns
        -------
        F : FLOAT
            Engine thrust.

        """
        return F_SL*delta0
    # ========================================================================    
    def func2(self, F_SL, delta0, M):
        """
        Thrust selection 2.

        Parameters
        ----------
        F_SL : FLOAT
            Nominal Sea Level Thrust at design speed V1.
        delta0 : FLOAT
            ISA Pressure ratio evaluated at h[ft].
        M : FLOAT
            A scalar variable containing the Mach number.

        Returns
        -------
        F : FLOAT
            Engine thrust.

        """
        x = (M - 0.1)**0.25
        y = 1 - 0.96*x
        z = delta0*y
        return F_SL*z
    # ========================================================================    
    def func3(self, F_SL, delta0, theta0, M, TR):
        """
        Thrust selection 2.        

        Parameters
        ----------
        F_SL : FLOAT
            Nominal Sea Level Thrust at design speed V1.
        delta0 : FLOAT
            ISA Pressure ratio evaluated at h[ft].
        theta0 : FLOAT
             ISA Temperature ratio evaluated at h[ft].
        M : FLOAT
            A scalar variable containing the Mach number.
        TR : FLOAT
            Engine Throttle ratio.

        Returns
        -------
        F : FLOAT
            Engine thrust.

        """
        x = (M - 0.1)**0.25
        y = 3*(theta0 - TR)
        z = 8.13*(M - 0.1)
        t = 1 - 0.96*x - (y/z)
        return F_SL*delta0*t
    # ========================================================================        
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
