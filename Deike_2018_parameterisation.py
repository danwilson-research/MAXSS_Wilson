# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 11:06:01 2022

This is an implementation of Deike et.al 2018.

Deike, L., & Melville, W. K. (2018).
Gas transfer by breaking waves,
Geophysical Research Letters, 45, 10,482–10,492.
https://doi.org/10.1029/2018GL078758

https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GL078758?casa_token=o10woxAFDCYAAAAA%3ADNcM2uCmE9CzDRnesz12Y7j_x7_R472sk9TzEAbcQLQewMnAmAMlA1kFwd78PyYMpUjYyrVtgEIwwww

Important note:this parameterisation is a parameterisation for higher solubility gases like CO2 and DMS.
Use with caution for less soluble gases.

Important note: This parameterisation raises friction velocity and significant
wave height to fractional powers e.g. 5/3 and 4/3.
This means that to capture the variability in these two parameters as we do by taking the
second moment of the wind e.g <u10>^2 for other parameterisations. It is necessary
to calculate the moment to this power for both these variables beforehand
and you must pass that as an input in the configuration file.
e.g. <u*>^1.666 <Hs>^1.333. Here these variables are crudely called
friction_velocity_1_666 and sig_wv_ht_1_333

Orignal script written by Richard Sims (2022), updated by Daniel Ford (2026) and further updated and vectorised by Daniel Wilson (2026)
"""

from fluxengine.core.rate_parameterisation import KCalculationBase
from fluxengine.core.datalayer import DataLayer
import numpy as np
from numpy import arange


#Deike 2018
class k_Deike2018(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;

    def input_names(self):
        return ["friction_velocity", "sig_wv_ht","solubility_skin", "scskin"];

    def output_names(self):
        return ["k","kb","kd"]; # Here k (total k), kd (direct k) and kb (bubble k) are output

    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Deike2018.__call__)";
        print("%s Using the Deike et al., 2018 (D18) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata)
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide"
            data["k"].longName="Deike et al., 2018 (D18) gas transfer velocity"
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function)
            print(type(e), e.args)
            return False


        # This is the implementation of Deike et.al 2018 (their equation 12)
        A_nb=1.55e-4 # Deike note this constant is the non breaking constant from Fairall et al., 2011
        A_b=1e-5
        g=9.81
        gamma=5/3
        zeta=4/3

        # Create valid mask
        valid_mask = (
            (self.friction_velocity != DataLayer.missing_value) &
            (self.sig_wv_ht != DataLayer.missing_value) &
            (self.solubility_skin != DataLayer.missing_value) &
            (self.scskin != DataLayer.missing_value) &
            (self.scskin > 0.0)
        )

        self.k[:] = DataLayer.missing_value
        self.kb[:] = DataLayer.missing_value # Setting up arrays
        self.kd[:] = DataLayer.missing_value

        # Vectorized calculations for valid points
        u_star = self.friction_velocity[valid_mask]
        hs = self.sig_wv_ht[valid_mask]
        sol_skin = self.solubility_skin[valid_mask]
        sc_skin = self.scskin[valid_mask]

        # Need the Ostwald solubility and not the solubility of CO2 as provided by fluxengine - Conversions below are from: Sander, R.: Compilation of Henry's law constants (version 5.0.0) for water as solvent, Atmos. Chem. Phys., 23, 10901–12440, https://doi.org/10.5194/acp-23-10901-2023, 2023.
        #Note 2478.96 is the molar volume of an ideal gas at 25deg expressed in Joules per mole
        ostwald = (sol_skin / 101.325) * 2478.96

        # 1. Non-breaking component (m/s)
        k_nb_ms = A_nb * u_star * np.sqrt(660.0 / sc_skin)

        # 2. Bubble component (m/s)
        k_b_ms = (A_b / ostwald) * (np.power(u_star, gamma) * np.power(np.sqrt(g * hs), zeta))

        # 3. Total transfer velocity converted to cm/h
        self.k[valid_mask] = (k_nb_ms + k_b_ms) * 360000.0
        self.kd[valid_mask] = (k_nb_ms) * 360000.0 # Adding the direct and bubble k to the outputs
        self.kb[valid_mask] = (k_b_ms) * 360000.0

        return True
