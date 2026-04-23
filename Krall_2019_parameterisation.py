# -*- coding: utf-8 -*-
"""
Created on Mon May 30 12:37:55 2022

This is an implementation of Krall et.al 2019 for seawater only.
Note that equation 11 is used and necessary coefficents are found in the paper supplement.

Krall, K. E., Smith, A. W., Takagaki, N., and Jähne, B.: 
Air–sea gas exchange at wind speeds up to 85 m s−1, 
Ocean Sci., 15, 1783–1799, https://doi.org/10.5194/os-15-1783-2019, 2019. 
https://os.copernicus.org/articles/15/1783/2019/

Important note: This parameterisation was developed in the Heidelberg wind wave
tank. It is not always clear which processes are being accurately simulated in
wind/wave tanks. For this reason it should not be assumed that results from 
tank experiments are completely representative of the open ocean. If users wish 
to use this parameterisation they should be aware of its potential limitations 
when applied in the open ocean.  

@author: rps207
"""

from fluxengine.core.rate_parameterisation import KCalculationBase
from fluxengine.core.datalayer import DataLayer

class k_Krall2019(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["friction_velocity", "scskin","solubility_skin"];
        
    def output_names(self):
        return ["k","ks","kc","kr"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Krall2019.__call__)";
        print("%s Using the Krall et al., 2019 (Kr19) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide" ;
            data["k"].longName="Krall et al., 2019 (Kr19) gas transfer velocity";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #determine the Krall et al., 2019 k relationship
        for i in arange(len(self.k)):   
            self.k[i] = DataLayer.missing_value
            if ( (self.friction_velocity[i] != DataLayer.missing_value) and (self.solubility_skin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
                
                #Krall Ks for Seawater
                # for friction velocity between 0.75 and 5.8 cms-1
                if (friction_velocity[i] >= 0.75) and (friction_velocity[i] < 5.8):
                    self.ks[i] = (3600/7.19)*(self.friction_velocity[i])*pow(600,-0.5); 
                # for friction velocity between 5.8 and 15 cms-1
                elif (friction_velocity[i] >= 5.8) and (friction_velocity[i] <= 15):
                    self.ks[i] = (0.605)*pow(self.friction_velocity[i],3); 
                else:
                    self.ks[i] = DataLayer.missing_value;
                    
                #Krall Kc for Seawater
                # for friction velocity between 0.75 and 5.8 cms-1
                if (friction_velocity[i] >= 0.75) and (friction_velocity[i] < 5.8):
                    self.kc[i] = 0; 
                # for friction velocity between 5.8 and 15 cms-1
                elif (friction_velocity[i] >= 5.8) and (friction_velocity[i] <= 15):
                    self.kc[i] = (51.5)*pow(self.friction_velocity[i]-5.8,1.82); 
                else: 
                    self.kc[i] = DataLayer.missing_value;
            
                #Krall Kr for Seawater
                # for friction velocity between 0.75 and 5.8 cms-1
                if (friction_velocity[i] >= 0.75) and (friction_velocity[i] < 5.8):
                    self.kr[i] = 0; 
                # for friction velocity between 5.8 and 15 cms-1
                elif (friction_velocity[i] >= 5.8) and (friction_velocity[i] <= 15):
                    self.kr[i] = (1.3)*pow(self.friction_velocity[i]-5.8,1.75); 
                else: 
                    self.kr[i] = DataLayer.missing_value;
                    
                # equation 11 from Krall 19
                self.k[i] = (self.ks[i] * sqrt(600.0/self.scskin[i]))+(self.kr[i]/self.solubility_skin[i])*(1-exp(((self.solubility_skin[i]*-1*self.kc[i])/self.kr[i])*(sqrt(600.0/self.scskin[i]))))
            else:
                self.k[i] = DataLayer.missing_value
        
        return True;
