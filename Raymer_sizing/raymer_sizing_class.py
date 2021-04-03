import math 
import numpy as np 
                    
# ==============================================================================================    
class raymer_sizing_class(object): 
    """
    INSERT DESCRIPTION HERE
    """
    def __init__(self, W_crew, W_payload, W_0, empty_A, empty_c, empty_m, AR_wet, V_cruise, TSFC_cruise, TSFC_loiter, Range, Endurance, Eff_LD):   
        self.W_crew, self.W_payload, self.W_0, self.empty_A, self.empty_c, self.empty_m, self.AR_wet, self.V_cruise,\
           self.TSFC_cruise,self.TSFC_loiter, self.Range, self.Endurance, self.Eff_LD = W_crew, W_payload, W_0, empty_A, empty_c, empty_m,\
             AR_wet, V_cruise, TSFC_cruise, TSFC_loiter, Range, Endurance, Eff_LD
    # ==========================================================================================
        self.W_empty_ratio  = self.empty_weight_fraction(W_0, empty_A, empty_c, empty_m)                     # W_empty/W_0
        self.lift_to_drag_r = self.lift_to_drag_ratio(Eff_LD, AR_wet)                                        # (L/D)max
        self.lift_to_drag_loi = self.lift_to_drag_ratio(0.866*Eff_LD, AR_wet)                                # (L/D) during loiting
        self.W_warmup       = 0.985                                                                          # W_1/W_0
        self.W_climb        = 0.97                                                                           # W_2/W_1
        self.W_cruise       = self.breguet_cruise(Range, self.lift_to_drag_r, TSFC_cruise, V_cruise)         # W_3/W_2
        self.W_loiter       = self.breguet_loiter(Endurance, self.lift_to_drag_loi, TSFC_loiter, V_cruise)   # W_4/W_3
        self.W_diversion    = self.breguet_cruise(Range/8, self.lift_to_drag_r, TSFC_cruise, V_cruise)       # W_5/W_4
        self.W_loiter1      = self.breguet_loiter(Endurance/2, self.lift_to_drag_loi, TSFC_loiter, V_cruise) # W_6/W_5
        self.W_landing      = 0.995                                                                          # W_7/W_6
        W_warmup, W_climb, W_cruise, W_loiter, W_diversion, W_loiter1, W_landing = self.W_warmup, self.W_climb, self.W_cruise, self.W_loiter, self.W_diversion, self.W_loiter1, self.W_landing                                                                       # Constant definition
        self.W_7            = W_warmup*W_climb*W_cruise*W_loiter*W_diversion*W_loiter1*W_landing             # W_x/W_0
        self.W_f            = self.fuel_to_gross_weight_ratio(self.W_7)                                      # W_f/W_0
        self.W_0_calculated = self.final_gross_weight(W_payload, W_crew, self.W_f, self.W_empty_ratio)       # W_0 estimated
        W_0_calculated      = self.W_0_calculated
    # ==========================================================================================
    def empty_weight_fraction(self, W_0, empty_A, empty_c, empty_m):
        """
        Function that calculates W_e/W_0
        """
        return empty_m*(W_0**empty_c)*empty_A
    # ==========================================================================================
    def lift_to_drag_ratio(self, E, AR_wet):
        """
        Function that calculates L/D_max
        """
        return E*np.sqrt(AR_wet)
    # ==========================================================================================
    def breguet_cruise(self, Range, lift_to_drag_ratio, TSFC, V):
        """
        Function that calculates fuel cruise fraction weight.
        """
        return 1/np.exp((Range*TSFC)/(V*lift_to_drag_ratio))
    # ==========================================================================================
    def breguet_loiter(self, Endurance, lift_to_drag_ratio, TSFC, V):
        """
        Function that calculates fuel loiter fraction weight.
        """
        return 1/np.exp((Endurance*TSFC)/((0.76*V)*lift_to_drag_ratio*0.866))
    # ==========================================================================================
    def fuel_to_gross_weight_ratio(self, Wx):
        """
        Function that calculates final fuel fraction weight.
        """
        return 1.06*(1 - Wx)
    # ==========================================================================================
    def final_gross_weight(self, W_payload, W_crew, W_fuel, W_empty):
        """
        Function that calculates final fuel fraction weight.
        """
        return (W_payload + W_crew)/(1 - W_fuel - W_empty) 
