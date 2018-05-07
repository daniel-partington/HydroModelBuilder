class Radon_EC_simple(object):
    """TODO: Docs"""

    def __init__(self, df, Ini_cond, Rn_decay=0.181, units='m,d'):
        self.Ini_Cond = Ini_cond  # 3 item list containing Initial flow, radon and ec concentrations
        self.Rn_decay = Rn_decay  # constant for Radon decay
        self.units = units  # units being used ... not really used yet!

        # Dataframe variables
        self.HZ_Poro = df['HZ_poro'].tolist()  # Hyporheic zone porosity
        self.HZ_Prod_Rate = df['HZ_Prod_Rate'].tolist()  # Hyporheic zone radon production rate
        self.HZ_RTime = df['HZ_RTime'].tolist()  # Hyporheic zone residence time
        self.R_depth_HZ = df['R_depth_HZ'].tolist()
        self.GTV = df['GTV'].tolist()  # Gas transfer velocity
        self.GW_Rn_conc = df['GW_Rn_conc'].tolist()  # Groundwater radon concentration
        self.GW_EC = df['GW_EC'].tolist()  # Groundwater EC
        self.Trib_EC = df['Trib_EC'].tolist()
        self.Trib_Rn = df['Trib_Rn'].tolist()
        self.dx = df['dx'].tolist()
        # Variables inherited from the output of SFR package
        self.Evap_Rate = df['Qet'].tolist()  # Evaporation rate
        self.Trib_inflow = df['Qovr'].tolist()

        try:
            self.R_Pump = df['R_pump'].tolist()
        except:
            self.R_Pump = [0.] * df.shape[0]
        # end try
        self.R_depth = df['depth'].tolist()
        self.R_width = df['width'].tolist()
        self.GW_inflow = (df['Qaquifer'] * -1.0).tolist()
    # End __init__()

    def Fl_Rn_EC_simul(self):
        """TODO: Docs"""

        # Specify radon decay constant
        Rn_decay = self.Rn_decay

        # Initialise temp vars
        Flow_temp = self.Ini_Cond[0]
        Rn_temp = self.Ini_Cond[1]
        EC_temp = self.Ini_Cond[2]
        Fl_simul = [Flow_temp]
        Rn_simul = [Rn_temp]
        EC_simul = [EC_temp]

        # Calculate simulated Radon and EC for the
        for i in xrange(len(self.dx) - 1):
            # Steady state stream flow calculation ... check against modelled
            # as this won't consider rate of change in storage
            Flow_temp += self.GW_inflow[i] - self.Evap_Rate[i] - \
                self.R_Pump[i] + self.Trib_inflow[i]
            # Radon concentration calculation based on steady state stream flow
            if self.GW_inflow[i] < 0.0:
                GW_Rn_conc_adj = Rn_temp
            else:
                GW_Rn_conc_adj = self.GW_Rn_conc[i]
            # End if
            Rn_temp += 1.0 / (Flow_temp) * (
                self.GW_inflow[i] * (GW_Rn_conc_adj - Rn_temp) +
                self.Trib_inflow[i] * (self.Trib_Rn[i] - Rn_temp) +
                Rn_temp * (self.Evap_Rate[i] - self.dx[i] *
                           self.R_width[i] * (self.GTV[i] + self.R_depth[i] * Rn_decay)) +
                self.dx[i] * self.HZ_Poro[i] / (1 + Rn_decay * self.HZ_RTime[i]) *
                self.R_width[i] * self.R_depth_HZ[i] * (self.HZ_Prod_Rate[i] - Rn_temp))
            if Rn_temp < 0:
                Rn_temp = 0.0
            # EC calculation based on steady state stream flow
            EC_temp += 1 / (Flow_temp) * (
                self.GW_inflow[i] * (self.GW_EC[1] - EC_temp) +
                self.Trib_inflow[i] * (self.Trib_EC[i] - EC_temp) +
                EC_temp * self.Evap_Rate[i])

            Fl_simul.append(Flow_temp)
            Rn_simul.append(Rn_temp)
            EC_simul.append(EC_temp)

        return Fl_simul, Rn_simul, EC_simul
    # End Fl_Rn_EC_simul()

# End Radon_EC_simple()
