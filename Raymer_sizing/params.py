# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 14:52:38 2021

@author: claum
"""
# maps label to attribute name and types
label_attr_map = {
# ==============================================================================================================
        "empty_A=": ["empty_A", float],
        "empty_c=": ["empty_c", float],
        "empty_m=": ["empty_m", float],
        "Endurance=": ["Endurance", float],    # [s]
        "Range=": ["Range", float],        # [ft]
        "W_crew=": ["W_crew", float],       # [lb]
        "W_payload=": ["W_payload", float],    # [lb]
        "W_0=": ["W_0", float],          # [lb]
      # V_cruise     = 557.849409449               # [ft * s^-1]
        "V_cruise=": ["V_cruise", float],     # [ft * s^-1]
        "eta_prop=": ["eta_prop", float],
        "PSFC=": ["PSFC", float],         # [lb * (hr * bhp)^-1]
        "TSFC_cruise=": ["TSFC_cruise", float],  # [s^-1] 
        "TSFC_loiter=": ["TSFC_loiter", float],  # [s^-1] 
        "D_fus=": ["D_fus", float],        # [ft]
        "l_fus=": ["l_fus", float],        # [ft]
        "slend_ratio=": ["slend_ratio", float],
        "S_wet_fus=": ["S_wet_fus", float],    # [ft^2]
        "S_wet_wing=": ["S_wet_wing", float],   # [ft^2]
        "taper_htail=": ["taper_htail", float],
        "taper_vtail=": ["taper_vtail", float], 
        "thick_root=": ["thick_root", float], 
        "thick_tip=": ["thick_tip", float],
        "tau=": ["tau", float],  
        "S_wet_htail=": ["S_wet_htail", float],  # [ft^2]   
        "S_wet_vtail=": ["S_wet_vtail", float],  # [ft^2]
        "S_wet_total=": ["S_wet_total", float],  # [ft^2]  
        "S_ratio=": ["S_ratio", float],      # S_wet/S_ref 
        "Aspect_ratio=": ["Aspect_ratio", float], # span^2/S_ref
        "AR_wet=": ["AR_wet ", float],      # AR/S_ratio
        "W_guess=": ["W_guess", float],      # [lb]
# ==============================================================================================================
}
# ============================================================================
class Params(object):
    def __init__(self, input_file_name):
        with open(input_file_name, 'r') as input_file:
            for line in input_file:
                row = line.split()
                label = row[0]
                data = row[1:]  # rest of row is data list

                attr = label_attr_map[label][0]
                datatypes = label_attr_map[label][1:]

                values = [(datatypes[i](data[i])) for i in range(len(data))]
                self.__dict__[attr] = values if len(values) > 1 else values[0]
# ============================================================================
