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
        self.T0, self.p0, self.rho0, self.h, self.gamma, self.R = T0, p0,\
            rho0, h, gamma, R
        self.x   = np.linspace(0, self.h, 1000)
        self.T   = self.temperature(self.T0, self.x)
        self.p   = self.pressure(self.p0, self.T, self.x)
        self.rho = self.density(self.p, self.T)
        self.a   = self.speed_of_sound(gamma, self.T, R)
    # ========================================================================        
    def temperature(self, T0, x):
        return self.T0 - 0.00356*(self.x)
    # ========================================================================    
    def pressure(self, p0, T, x):
        return self.p0*((self.T + 459.7)/518.6)**(5.2561)
    # ========================================================================    
    def density(self, p, T): 
        return self.p/(1718.0*(self.T + 459.7))
    # ========================================================================
    def speed_of_sound(self, gamma, T, R): 
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
    A class that calculates standard atmosphere in Engineering units.
    """
    def __init__(self, T, p, M, gamma):
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
        d = T/T0
        return x*d
    # ====================================================================
    def delta_ratio(self, p, p0, x, gamma):
        d = p/p0
        y = gamma/(gamma-1)
        z = x**y
        return z*d
    # ====================================================================
    def partial_p(self, gamma, M):
        xx = M**2
        yy = gamma - 1
        zz = 1 + 0.5*yy*xx
        return zz
    # ====================================================================            
# ============================================================================
# ATMOSPHERE RATIO
# ============================================================================
# ============================================================================
# INITIAL SIZING
# ============================================================================
class initial_sizing(object):
    def __init__(self, AR, rho, Vmax, Vmin, maxWS, nMAX, CDmin, CDTO,\
                 maxROC, g, mu, V_liftoff, CLTO, Sg, V_design, V_climb):
        # ====================================================================
        self.AR, self.rho, self.Vmax, self.Vmin, self.maxWS, self.nMAX,\
        self.CDmin, self.CDTO, self.maxROC, self.g, self.mu, self.V_liftoff,\
        self.CLTO, self.Sg, self.V_design, self.V_climb = AR, rho, Vmax, Vmin,\
        maxWS, nMAX, CDmin, CDTO, maxROC, g, mu, V_liftoff, CLTO, Sg, V_design,\
        V_climb
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
    # ========================================================================        
    def efficiency_factor(self, AR):
        """
        Function that calculates Oswald efficiency factor

        Parameters
        ----------
        AR : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

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
        rho : TYPE
            DESCRIPTION.
        V : TYPE
            DESCRIPTION.

        Returns
        -------
        q : FLOAT
            Dynamic pressure array [lb/ft^2].

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
        WS : FLOAT
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
        Function that calculates Thrust - to - weight ratio T/W for constant 
        turn at a specified load factor nMAX

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        CDmin : TYPE
            DESCRIPTION.
        ws : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        nMAX : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

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
        

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        CDmin : TYPE
            DESCRIPTION.
        ws : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        nMAX : TYPE
            DESCRIPTION.
        maxROC : TYPE
            DESCRIPTION.
        V : TYPE
            DESCRIPTION.

        Returns
        -------
        kk : TYPE
            DESCRIPTION.

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
        

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        CDTO : TYPE
            DESCRIPTION.
        CLTO : TYPE
            DESCRIPTION.
        g : TYPE
            DESCRIPTION.
        Sg : TYPE
            DESCRIPTION.
        ws : TYPE
            DESCRIPTION.
        mu : TYPE
            DESCRIPTION.
        V_liftoff : TYPE
            DESCRIPTION.

        Returns
        -------
        kk : TYPE
            DESCRIPTION.

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
        

        Parameters
        ----------
        q : TYPE
            DESCRIPTION.
        CDmin : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.

        Returns
        -------
        kk : TYPE
            DESCRIPTION.

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
        

        Parameters
        ----------
        MaxROC : TYPE
            DESCRIPTION.
        CDmin : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        ws : TYPE
            DESCRIPTION.
        rho : TYPE
            DESCRIPTION.

        Returns
        -------
        kk : TYPE
            DESCRIPTION.

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
    def __init__(self, F_SL, delta0, theta0, thetac, M):
        """

        Parameters
        ----------
        F_SL : FLOAT
            Thrust rated at sea level.
        delta0 : FLOAT
            Pressure ratio.
        theta0 : FLOAT
            Temperature ratio.
        M : FLOAT
            Mach number.
        TR : FLOAT 
            Throttle ratio: modeling the temperature effect on engine
            performance.

        Returns
        -------
        Thrust - M number.

        """
        self.F_SL, self.delta0, self.theta0, self.thetac, self.M, = F_SL,\
            delta0, theta0, thetac, M
        self.T0 = 59       
        self.TR = self.throttle_ratio(thetac, self.T0)
        TR      = self.TR
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
        T = thetac*T0
        TR = T/T0
        return TR          
    # ========================================================================    
    def func1(self, F_SL, delta0):
        return F_SL*delta0
    # ========================================================================    
    def func2(self, F_SL, delta0, M):
        x = (M - 0.1)**0.25
        y = 1 - 0.96*x
        z = delta0*y
        return F_SL*z
    # ========================================================================    
    def func3(self, F_SL, delta0, theta0, M, TR):
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
