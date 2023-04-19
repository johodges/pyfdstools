7#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
# 
# DESCRIPTION:
# This software is part of a python library to assist in developing and
# analyzing simulation results from Fire Dynamics Simulator (FDS).
# FDS is an open source software package developed by NIST. The source
# code is available at: https://github.com/firemodels/fds
# 
# This script defines the data types available for each type of line in
# FDS.
#
#=======================================================================
# # IMPORTS
#=======================================================================

from collections import defaultdict

class fdsLineTypes(object):
    def __init__(self, version="6.7.4"):
        self.bndf = self.getBNDFtypes(version)
        self.clip = self.getCLIPtypes(version)
        self.comb = self.getCOMBtypes(version)
        self.ctrl = self.getCTRLtypes(version)
        self.devc = self.getDEVCtypes(version)
        self.dump = self.getDUMPtypes(version)
        self.geom = self.getGEOMtypes(version)
        self.head = self.getHEADtypes(version)
        self.hole = self.getHOLEtypes(version)
        self.hvac = self.getHVACtypes(version)
        self.init = self.getINITtypes(version)
        self.matl = self.getMATLtypes(version)
        self.mesh = self.getMESHtypes(version)
        self.misc = self.getMISCtypes(version)
        self.mult = self.getMULTtypes(version)
        self.obst = self.getOBSTtypes(version)
        self.part = self.getPARTtypes(version)
        self.pres = self.getPREStypes(version)
        self.prof = self.getPROFtypes(version)
        self.prop = self.getPROPtypes(version)
        self.radi = self.getRADItypes(version)
        self.ramp = self.getRAMPtypes(version)
        self.reac = self.getREACtypes(version)
        self.slcf = self.getSLCFtypes(version)
        self.spec = self.getSPECtypes(version)
        self.surf = self.getSURFtypes(version)
        self.tabl = self.getTABLtypes(version)
        self.time = self.getTIMEtypes(version)
        self.trn2 = self.getTRNXtypes(version)
        self.trny = self.getTRNYtypes(version)
        self.trnz = self.getTRNZtypes(version)
        self.vent = self.getVENTtypes(version)
        self.wind = self.getWINDtypes(version)
        self.zone = self.getZONEtypes(version)

    def getCLIPtypes(self, version="6.7.4"):
        clipTypes = defaultdict(bool)
        clipTypes['ID'] = 'string'
        clipTypes['MAXIMUM_DENSITY'] = 'float'
        clipTypes['MINIMUM_DENSITY'] = 'float'
        clipTypes['CLIP_DT_RESTRICTIONS_MAX'] = 'float'
        return clipTypes

    def getCOMBtypes(self, version="6.7.4"):
        combTypes = defaultdict(bool)
        combTypes['ID'] = 'string'
        return combTypes
        
    def getGEOMtypes(self, version="6.7.4"):
        geomTypes = defaultdict(bool)
        geomTypes['ID'] = 'string'
        geomTypes['N_LEVELS'] = 'int'
        geomTypes['SPHERE_ORIGIN'] = 'listfloat'
        geomTypes['SPHERE_RADIUS'] = 'float'
        geomTypes['SURF_ID'] = 'string'
        geomTypes['SURF_IDS'] = 'liststring'
        geomTypes['XB'] = 'listfloat'
        return geomTypes
        
    def getWINDtypes(self, version="6.7.4"):
        windTypes = defaultdict(bool)
        windTypes['GROUND_LEVEL'] = 'float'
        windTypes['L'] = 'float'
        windTypes['MEAN_FORCING'] = 'matrixbool'
        windTypes['SPEED'] = 'float'
        windTypes['STRATIFICATION'] = 'bool'
        windTypes['U0'] = 'float'
        windTypes['V0'] = 'float'
        windTypes['W0'] = 'float'
        windTypes['Z_0'] = 'float'
        windTypes['Z_REF'] = 'float'
        windTypes['newline'] = 'ignore'
        return windTypes
    
    def getPROFtypes(self, version="6.7.4"):
        profTypes = defaultdict(bool)
        profTypes['ID'] = 'string'
        profTypes['XYZ'] = 'listfloat'
        profTypes['QUANTITY'] = 'string'
        profTypes['IOR'] = 'int'
        profTypes['FORMAT_INDEX'] = 'int'
        profTypes['newline'] = 'ignore'
        return profTypes
        
    def getRAMPtypes(self, version="6.7.4"):
        rampTypes = defaultdict(bool)
        rampTypes['T'] = 'listrowfloat'
        rampTypes['F'] = 'listrowfloat'
        rampTypes['ID'] = 'string'
        rampTypes['CTRL_ID'] = 'string'
        rampTypes['DEVC_ID'] = 'string'
        rampTypes['newline'] = 'ignore'
        return rampTypes
        
    def getTABLtypes(self, version="6.7.4"):
        tablTypes = defaultdict(bool)
        tablTypes['ID'] = 'string'
        tablTypes['TABLE_DATA'] = 'listfloat'
        return tablTypes
        
    def getPROPtypes(self, version="6.7.4"):
        propTypes = defaultdict(bool)
        propTypes['ACTIVATION_OBSCURATION'] = 'float'
        propTypes['ACTIVATION_TEMPERATURE'] = 'float'
        propTypes['ALPHA_C'] = 'float'
        propTypes['ALPHA_E'] = 'float'
        propTypes['BETA_C'] = 'float'
        propTypes['BETA_E'] = 'float'
        propTypes['C_FACTOR'] = 'float'
        propTypes['FLOW_RATE'] = 'float'
        propTypes['FLOW_RAMP'] = 'string'
        propTypes['FLOW_TAU'] = 'float'
        propTypes['GAUGE_EMISSIVITY'] = 'float'
        propTypes['GAUGE_TEMPERATURE'] = 'float'
        propTypes['HISTOGRAM'] = 'bool'
        propTypes['HISTOGRAM_CUMULATIVE'] = 'bool'
        propTypes['HISTOGRAM_LIMITS'] = 'listfloat'
        propTypes['HISTOGRAM_NBINS'] = 'int'
        propTypes['ID'] = 'string'
        propTypes['LENGTH'] = 'float'
        propTypes['OFFSET'] = 'float'
        propTypes['PART_ID'] = 'string'
        propTypes['PARTICLE_VELOCITY'] = 'float'
        propTypes['PARTICLES_PER_SECOND'] = 'int'
        propTypes['PDPA_END'] = 'float'
        propTypes['PDPA_INTEGRATE'] = 'bool'
        propTypes['PDPA_M'] = 'float'
        propTypes['PDPA_N'] = 'float'
        propTypes['PDPA_NORMALIZE'] = 'bool'
        propTypes['PDPA_RADIUS'] = 'float'
        propTypes['PDPA_START'] = 'float'
        propTypes['QUANTITY'] = 'string'
        propTypes['RTI'] = 'float'
        propTypes['SMOKEVIEW_ID'] = 'string'
        propTypes['SMOKEVIEW_PARAMETERS'] = 'liststring'
        propTypes['SPRAY_ANGLE'] = 'listfloat'
        propTypes['SPRAY_PATTERN_TABLE'] = 'string'
        propTypes['newline'] = 'ignore'
        return propTypes
    
    def getSPECtypes(self, version="6.7.4"):
        specTypes = defaultdict(bool)
        specTypes['AEROSOL'] = 'bool'
        specTypes['BACKGROUND'] = 'bool'
        specTypes['CONDUCTIVITY'] = 'float'
        specTypes['CONDUCTIVITY_SOLID'] = 'float'
        specTypes['DENSITY_LIQUID'] = 'float'
        specTypes['DENSITY_SOLID'] = 'float'
        specTypes['DIFFUSIVITY'] = 'float'
        specTypes['ENTHALPY_OF_FORMATION'] = 'float'
        specTypes['FORMULA'] = 'string'
        specTypes['HEAT_OF_VAPORIZATION'] = 'float'
        specTypes['ID'] = 'string'
        specTypes['LUMPED_COMPONENT_ONLY'] = 'bool'
        specTypes['MASS_FRACTION_0'] = 'float'
        specTypes['MAX_DIAMETER'] = 'float'
        specTypes['MEAN_DIAMETER'] = 'float'
        specTypes['MELTING_TEMPERATURE'] = 'float'
        specTypes['MIN_DIAMETER'] = 'float'
        specTypes['MW'] = 'float'
        specTypes['N_BINS'] = 'int'
        specTypes['RAMP_CP'] = 'string'
        specTypes['RAMP_K'] = 'string'
        specTypes['SPEC_ID'] = 'liststring'
        specTypes['SPECIFIC_HEAT'] = 'float'
        specTypes['SPECIFIC_HEAT_LIQUID'] = 'float'
        specTypes['THERMOPHORETIC_DIAMETER'] = 'float'
        specTypes['VAPORIZATION_TEMPERATURE'] = 'float'
        specTypes['VISCOSITY'] = 'float'
        specTypes['VOLUME_FRACTION'] = 'listfloat'
        specTypes['newline'] = 'ignore'
        return specTypes
        
    def getMATLtypes(self, version="6.7.4"):
        matlTypes = defaultdict(bool)
        matlTypes['A'] = 'listfloat'
        matlTypes['ABSORPTION_COEFFICIENT'] = 'float'
        matlTypes['CONDUCTIVITY'] = 'float'
        matlTypes['CONDUCTIVITY_RAMP'] = 'string'
        matlTypes['DENSITY'] = 'float'
        matlTypes['E'] = 'listfloat'
        matlTypes['EMISSIVITY'] = 'float'
        matlTypes['FYI'] = 'string'
        matlTypes['HEAT_OF_COMBUSTION'] = 'float'
        matlTypes['HEAT_OF_REACTION'] = 'listfloat'
        matlTypes['HEATING_RATE'] = 'float'
        matlTypes['ID'] = 'string'
        matlTypes['MATL_ID'] = 'matrixstring'
        matlTypes['NU_MATL'] = 'matrixfloat'
        matlTypes['NU_SPEC'] = 'matrixfloat'
        matlTypes['N_O2'] = 'listfloat'
        matlTypes['N_REACTIONS'] = 'int'
        matlTypes['N_S'] = 'listfloat'
        matlTypes['N_T'] = 'listfloat'
        matlTypes['NU_SPEC'] = 'matrixfloat'
        matlTypes['PYROLYSIS_RANGE'] = 'listfloat'
        matlTypes['REFERENCE_TEMPERATURE'] = 'listfloat'
        matlTypes['SPEC_ID'] = 'matrixstring'
        matlTypes['SPECIFIC_HEAT'] = 'float'
        matlTypes['SPECIFIC_HEAT_RAMP'] = 'string'
        matlTypes['newline'] = 'ignore'
        return matlTypes
    
    def getMULTtypes(self, version="6.7.4"):
        multTypes = defaultdict(bool)
        multTypes['DX'] = 'float'
        multTypes['DY'] = 'float'
        multTypes['DZ'] = 'float'
        multTypes['I_UPPER'] = 'int'
        multTypes['ID'] = 'string'
        multTypes['J_UPPER'] = 'int'
        multTypes['K_UPPER'] = 'int'
        return multTypes
    
    def getREACtypes(self, version="6.7.4"):
        reacTypes = defaultdict(bool)
        reacTypes['AUTO_IGNITION_TEMPERATURE'] = 'float'
        reacTypes['C'] = 'float'
        reacTypes['CO_YIELD'] = 'float'
        reacTypes['CRITICAL_FLAME_TEMPERATURE'] = 'float'
        reacTypes['EPUMO2'] = 'float'
        reacTypes['FORMULA'] = 'string'
        reacTypes['FUEL'] = 'string'
        reacTypes['FYI'] = 'string'
        reacTypes['H'] = 'float'
        reacTypes['HEAT_OF_COMBUSTION'] = 'float'
        reacTypes['ID'] = 'string'
        reacTypes['IDEAL'] = 'bool'
        reacTypes['N'] = 'float'
        reacTypes['N_SIMPLE_CHEMISTRY_REACTIONS'] = 'int'
        reacTypes['O'] = 'float'
        reacTypes['RADIATIVE_FRACTION'] = 'float'
        reacTypes['RAMP_CHI_R'] = 'string'
        reacTypes['SOOT_YIELD'] = 'float'
        reacTypes['SOOT_H_FRACTION'] = 'float'
        reacTypes['SPEC_ID_NU'] = 'liststring'
        reacTypes['NU'] = 'listfloat'
        reacTypes['newline'] = 'ignore'
        return reacTypes

    def getRADItypes(self, version="6.7.4"):
        radiTypes = defaultdict(bool)
        radiTypes['ANGLE_INCREMENT'] = 'int'
        radiTypes['INITIAL_RADIATION_ITERATIONS'] = 'int'
        radiTypes['KAPPA0'] = 'float'
        radiTypes['MIE_MINIMUM_DIAMETER'] = 'float'
        radiTypes['MIE_MAXIMUM_DIAMETER'] = 'float'
        radiTypes['MIE_NDG'] = 'float'
        radiTypes['NMIEANG'] = 'int'
        radiTypes['NUMBER_RADIATION_ANGLES'] = 'int'
        radiTypes['PATH_LENGTH'] = 'float'
        radiTypes['RADIATION'] = 'bool'
        radiTypes['RADIATIVE_FRACTION'] = 'float'
        radiTypes['RADIATION_ITERATIONS'] = 'int'
        radiTypes['RADTMP'] = 'float'
        radiTypes['TIME_STEP_INCREMENT'] = 'int'
        radiTypes['WIDE_BAND_MDOEL'] = 'bool'
        radiTypes['newline'] = 'ignore'
        return radiTypes

    def getMESHtypes(self, version="6.7.4"):
        meshTypes = defaultdict(bool)
        meshTypes['CYLINDRICAL'] = 'bool'
        meshTypes['ID'] = 'string'
        meshTypes['IJK'] = 'listint'
        meshTypes['MPI_PROCESS'] = 'int'
        meshTypes['MULT_ID'] = 'string'
        meshTypes['TRNX_ID'] = 'string'
        meshTypes['TRNY_ID'] = 'string'
        meshTypes['TRNZ_ID'] = 'string'
        meshTypes['XB'] = 'listfloat'
        meshTypes['unknownCounter'] = 'ignore'
        meshTypes['newline'] = 'ignore'
        return meshTypes

    def getSURFtypes(self, version="6.7.4"):
        surfTypes = defaultdict(bool)
        surfTypes['ADIABATIC'] = 'bool'
        surfTypes['BACKING'] = 'string'
        surfTypes['BURN_AWAY'] = 'bool'
        surfTypes['C_FORCED_RE'] = 'float'
        surfTypes['C_FORCED_CONSTANT'] = 'float'
        surfTypes['C_FORCED_RE_EXP'] = 'float'
        surfTypes['C_FORCED_PR_EXP'] = 'float'
        surfTypes['C_HORIZONTAL'] = 'float'
        surfTypes['C_VERTICAL'] = 'float'
        surfTypes['CELL_SIZE_FACTOR'] = 'float'
        surfTypes['COLOR'] = 'string'
        surfTypes['CONE_HEAT_FLUX'] = 'float'
        surfTypes['CONVECTION_LENGTH_SCALE'] = 'float'
        surfTypes['DEFAULT'] = 'bool'
        surfTypes['DRAG_COEFFICIENT'] = 'float'
        surfTypes['DT_INSERT'] = 'float'
        surfTypes['E_COEFFICIENT'] = 'float'
        surfTypes['EMISSIVITY'] = 'float'
        surfTypes['EMISSIVITY_BACK'] = 'float'
        surfTypes['EXTERNAL_FLUX'] = 'float'
        surfTypes['EXTERNAL_FLUX_RAMP'] = 'string'
        surfTypes['EXTINCTION_TEMPERATURE'] = 'float'
        surfTypes['FREE_SLIP'] = 'bool'
        surfTypes['FYI'] = 'string'        
        surfTypes['GEOMETRY'] = 'string'
        surfTypes['HEAT_TRANSFER_COEFFICIENT'] = 'float'
        surfTypes['HEAT_TRANSFER_COEFFICIENT_BACK'] = 'float'
        surfTypes['HRRPUA'] = 'float'
        surfTypes['ID'] = 'string'
        surfTypes['IGNITION_TEMPERATURE'] = 'float'
        surfTypes['INNER_RADIUS'] = 'float'
        surfTypes['LEAK_PATH'] = 'listint'
        surfTypes['LENGTH'] = 'float'
        surfTypes['MASS_FLUX'] = 'float'
        surfTypes['MASS_FLUX_TOTAL'] = 'float'
        surfTypes['MASS_FRACTION'] = 'matrixfloat'
        surfTypes['MATL_ID'] = 'matrixstring'
        surfTypes['MATL_MASS_FRACTION'] = 'matrixfloat'
        surfTypes['MAXIMUM_SCALING_HEAT_FLUX'] = 'float'
        surfTypes['MINIMUM_SCALING_HEAT_FLUX'] = 'float'
        surfTypes['MLRPUA'] = 'float'
        surfTypes['NO_SLIP'] = 'bool'
        surfTypes['NPPC'] = 'float'
        surfTypes['PART_ID'] = 'string'
        surfTypes['PARTICLE_MASS_FLUX'] = 'float'
        surfTypes['PROFILE'] = 'string'
        surfTypes['RADIUS'] = 'float'
        surfTypes['RAMP_EF'] = 'string'
        surfTypes['RAMP_MF'] = 'string'
        surfTypes['RAMP_Q'] = 'string'
        surfTypes['RAMP_T'] = 'string'
        surfTypes['RAMP_TMP_GAS_FRONT'] = 'string'
        surfTypes['RAMP_V'] = 'string'
        surfTypes['REFERENCE_HEAT_FLUX'] = 'float'
        surfTypes['REFERENCE_HEAT_FLUX_TIME_INTERVAL'] = 'float'
        surfTypes['RGB'] = 'listfloat'
        surfTypes['SPEC_ID'] = 'string'
        surfTypes['STRETCH_FACTOR'] = 'listfloat'
        surfTypes['SURFACE_VOLUME_RATIO'] = 'float'
        surfTypes['TAU_MF'] = 'float'
        surfTypes['TAU_PART'] = 'float'
        surfTypes['TAU_Q'] = 'float'
        surfTypes['TAU_T'] = 'float'
        surfTypes['TAU_V'] = 'float'
        surfTypes['TEXTURE_MAP'] = 'string'
        surfTypes['TEXTURE_WIDTH'] = 'float'
        surfTypes['TEXTURE_HEIGHT'] = 'float'
        surfTypes['THICKNESS'] = 'listindfloat'
        surfTypes['TMP_BACK'] = 'float'
        surfTypes['TMP_FRONT'] = 'float'
        surfTypes['TMP_GAS_FRONT'] = 'float'
        surfTypes['TMP_INNER'] = 'float'
        surfTypes['TRANSPARENCY'] = 'float'
        surfTypes['VEL'] = 'float'
        surfTypes['VEL_GRAD'] = 'float'
        surfTypes['VEL_T'] = 'listfloat'
        surfTypes['VOLUME_FLOW'] = 'float'
        surfTypes['WIDTH'] = 'float'
        surfTypes['newline'] = 'ignore'
        return surfTypes

    def getTRNXtypes(self, version="6.7.4"):
        trnTypes = defaultdict(bool)
        trnTypes['CC'] = 'float'
        trnTypes['ID'] = 'string'
        trnTypes['PC'] = 'float'
        return trnTypes
    
    def getTRNYtypes(self, version="6.7.4"):
        trnTypes = defaultdict(bool)
        trnTypes['CC'] = 'float'
        trnTypes['ID'] = 'string'
        trnTypes['PC'] = 'float'
        return trnTypes

    def getTRNZtypes(self, version="6.7.4"):
        trnTypes = defaultdict(bool)
        trnTypes['CC'] = 'float'
        trnTypes['ID'] = 'string'
        trnTypes['PC'] = 'float'
        return trnTypes
    
    def getVENTtypes(self, version="6.7.4"):
        surfTypes = defaultdict(bool)
        surfTypes['COLOR'] = 'string'
        surfTypes['CTRL_ID'] = 'string'
        surfTypes['DEVC_ID'] = 'string'
        surfTypes['DYNAMIC_PRESSURE'] = 'float'
        surfTypes['ID'] = 'string'
        surfTypes['IOR'] = 'int'
        surfTypes['MB'] = 'string'
        surfTypes['PBX'] = 'float'
        surfTypes['PBY'] = 'float'
        surfTypes['PBZ'] = 'float'
        surfTypes['RADIUS'] = 'float'
        surfTypes['RGB'] = 'listint'
        surfTypes['SURF_ID'] = 'string'
        surfTypes['TRANSPARENCY'] = 'float'
        surfTypes['UVW'] = 'listfloat'
        surfTypes['WIND'] = 'bool'
        surfTypes['XB'] = 'listfloat'
        surfTypes['XYZ'] = 'listfloat'
        surfTypes['number'] = 'ignore'
        surfTypes['unknownCounter'] = 'ignore'
        surfTypes['newline'] = 'ignore'
        return surfTypes

    def getPREStypes(self, version="6.7.4"):
        presTypes = defaultdict(bool)
        presTypes['CHECK_POISSON'] = 'bool'
        presTypes['FISHPAK_BC'] = 'listfloat'
        presTypes['MAX_PRESSURE_ITERATIONS'] = 'int'
        presTypes['PRESSURE_RELAX_TIME'] = 'float'
        presTypes['PRESSURE_TOLERANCE'] = 'float'
        presTypes['SOLVER'] = 'string'
        presTypes['SUSPEND_PRESSURE_ITERATIONS'] = 'bool'
        presTypes['TUNNEL_PRECONDITIONER'] = 'bool'
        presTypes['VELOCITY_TOLERANCE'] = 'float'
        presTypes['newline'] = 'ignore'
        return presTypes

    def getOBSTtypes(self, version="6.7.4"):
        obstTypes = defaultdict(bool)
        obstTypes['ID'] = 'string'
        obstTypes['XB'] = 'listfloat'
        obstTypes['CTRL_ID'] = 'string'
        obstTypes['DEVC_ID'] = 'string'
        obstTypes['SURF_ID'] = 'string'
        obstTypes['SURF_ID6'] = 'liststring'
        obstTypes['SURF_IDS'] = 'liststring'
        obstTypes['BNDF_OBST'] = 'bool'
        obstTypes['PERMIT_HOLE'] = 'bool'
        obstTypes['TEXTURE_ORIGIN'] = 'listfloat'
        obstTypes['THICKEN'] = 'bool'
        obstTypes['TRANSPARENCY'] = 'float'
        obstTypes['COLOR'] = 'string'
        obstTypes['BULK_DENSITY'] = 'float'
        obstTypes['RGB'] = 'listint'
        obstTypes['TRANSPARENCY'] = 'float'
        obstTypes['number'] = 'ignore'
        obstTypes['unknownCounter'] = 'ignore'
        obstTypes['newline'] = 'ignore'
        return obstTypes
    
    def getPARTtypes(self, version="6.7.4"):
        partTypes = defaultdict(bool)
        partTypes['AGE'] = 'float'
        partTypes['CHECK_DISTRIBUTION'] = 'bool'
        partTypes['CNF_RAMP_ID'] = 'string'
        partTypes['COLOR'] = 'string'
        partTypes['COMPLEX_REFRACTIVE_INDEX'] = 'float'
        partTypes['DENSE_VOLUME_FRACTION'] = 'float'
        partTypes['DIAMETER'] = 'float'
        partTypes['DISTRIBUTION'] = 'string'
        partTypes['DRAG_COEFFICIENT'] = 'listfloat'
        partTypes['DRAG_LAW'] = 'string'
        partTypes['FREE_AREA_FRACTION'] = 'float'
        partTypes['HORIZONTAL_VELOCITY'] = 'float'
        partTypes['ID'] = 'string'
        partTypes['INITIAL_TEMPERATURE'] = 'float'
        partTypes['MASSLESS'] = 'bool'
        partTypes['MINIMUM_DIAMETER'] = 'float'
        partTypes['MONODISPERSE'] = 'bool'
        partTypes['ORIENTATION'] = 'matrixint'
        partTypes['PERMEABILITY'] = 'listfloat'
        partTypes['POROUS_VOLUME_FRACTION'] = 'float'
        partTypes['PROP_ID'] = 'string'
        partTypes['QUANTITIES'] = 'matrixstring'
        partTypes['RADIATIVE_PROPERTY_TABLE'] = 'string'
        partTypes['REAL_REFRACTIVE_INDEX'] = 'float'
        partTypes['SAMPLING_FACTOR'] = 'float'
        partTypes['SPEC_ID'] = 'string'
        partTypes['STATIC'] = 'bool'
        partTypes['SURF_ID'] = 'string'
        partTypes['newline'] = 'ignore'
        return partTypes
    
    def getHOLEtypes(self, version="6.7.4"):
        holeTypes = defaultdict(bool)
        holeTypes['BLOCK_WIND'] = 'bool'
        holeTypes['CTRL_ID'] = 'string'
        holeTypes['DEVC_ID'] = 'string'
        holeTypes['ID'] = 'string'
        holeTypes['XB'] = 'listfloat'
        holeTypes['number'] = 'ignore'
        holeTypes['unknownCounter'] = 'ignore'
        holeTypes['newline'] = 'ignore'
        return holeTypes
    
    def getHVACtypes(self, version="6.7.4"):
        hvacTypes = defaultdict(bool)
        hvacTypes['AREA'] = 'float'
        hvacTypes['DUCT_ID'] = 'string'
        hvacTypes['ID'] = 'string'
        hvacTypes['NODE_ID'] = 'liststring'
        hvacTypes['TAU_VF'] = 'float'
        hvacTypes['TYPE_ID'] = 'string'
        hvacTypes['VENT_ID'] = 'string'
        hvacTypes['VENT2_ID'] = 'string'
        hvacTypes['VOLUME_FLOW'] = 'float'
        hvacTypes['number'] = 'ignore'
        hvacTypes['unknownCounter'] = 'ignore'
        hvacTypes['newline'] = 'ignore'
        return hvacTypes

    def getINITtypes(self, version="6.7.4"):
        initTypes = defaultdict(bool)
        initTypes['CELL_CENTERED'] = 'bool'
        initTypes['DEVC_ID'] = 'string'
        initTypes['DRY'] = 'bool'
        initTypes['DT_INSERT'] = 'float'
        initTypes['DX'] = 'float'
        initTypes['HRRPUV'] = 'float'
        initTypes['ID'] = 'string'
        initTypes['MASS_FRACTION'] = 'matrixfloat'
        initTypes['MASS_PER_TIME'] = 'float'
        initTypes['MASS_PER_VOLUME'] = 'float'
        initTypes['N_PARTICLES'] = 'int'
        initTypes['N_PARTICLES_PER_CELL'] = 'int'
        initTypes['PART_ID'] = 'string'
        initTypes['PARTICLE_WEIGHT_FACTOR']='float'
        initTypes['PATH_RAMP'] = 'liststring'
        initTypes['SPEC_ID'] = 'matrixstring'
        initTypes['TEMPERATURE'] = 'float'
        initTypes['UVW'] = 'listfloat'
        initTypes['XB'] = 'listfloat'
        initTypes['XYZ'] = 'listfloat'
        initTypes['newline'] = 'ignore'
        return initTypes

    def getDEVCtypes(self, version="6.7.4"):
        devcTypes = defaultdict(bool)
        devcTypes['BYPASS_FLOWRATE'] = 'float'
        devcTypes['CONVERSION_FACTOR'] = 'float'
        devcTypes['CTRL_ID'] = 'string'
        devcTypes['DELAY'] = 'float'
        devcTypes['DEPTH'] = 'float'
        devcTypes['DEVC_ID'] = 'string'
        devcTypes['DUCT_ID'] = 'string'
        devcTypes['FLOWRATE'] = 'float'
        devcTypes['HIDE_COORDINATES'] = 'bool'
        devcTypes['ID'] = 'string'
        devcTypes['INITIAL_STATE'] = 'bool'
        devcTypes['INIT_ID'] = 'string'
        devcTypes['IOR'] = 'int'
        devcTypes['LATCH'] = 'bool'
        devcTypes['MATL_ID'] = 'string'
        devcTypes['NO_UPDATE_DEVC_ID'] = 'string'
        devcTypes['NO_UPDATE_CTRL_ID'] = 'string'
        devcTypes['ORIENTATION'] = 'listfloat'
        devcTypes['OUTPUT'] = 'bool'
        devcTypes['PART_ID'] = 'string'
        devcTypes['POINTS'] = 'int'
        devcTypes['PROP_ID'] = 'string'
        devcTypes['QUANTITY'] = 'string'
        devcTypes['RELATIVE'] = 'bool'
        devcTypes['SETPOINT'] = 'float'
        devcTypes['SPATIAL_STATISTIC'] = 'string'
        devcTypes['SPEC_ID'] = 'string'
        devcTypes['STATISTICS'] = 'string'
        devcTypes['STATISTICS_START'] = 'float'
        devcTypes['SURF_ID'] = 'string'
        devcTypes['TEMPORAL_STATISTIC'] = 'string'
        devcTypes['TIME_AVERAGED'] = 'bool'
        devcTypes['TIME_HISTORY'] = 'bool'
        devcTypes['X_ID'] = 'string'
        devcTypes['XB'] = 'listfloat'
        devcTypes['XYZ'] = 'listfloat'
        devcTypes['Z_ID'] = 'string'
        devcTypes['UNITS'] = 'string'
        devcTypes['unknownCounter'] = 'ignore'
        devcTypes['newline'] = 'ignore'
        return devcTypes

    def getBNDFtypes(self, version="6.7.4"):
        bndfTypes = defaultdict(bool)
        bndfTypes['CELL_CENTERED'] = 'bool'
        bndfTypes['ID'] = 'ignore'
        bndfTypes['PART_ID'] = 'string'
        bndfTypes['PROP_ID'] = 'string'
        bndfTypes['QUANTITY'] = 'string'
        bndfTypes['SPEC_ID'] = 'string'
        bndfTypes['unknownCounter'] = 'ignore'
        bndfTypes['newline'] = 'ignore'
        return bndfTypes
    
    def getSLCFtypes(self, version="6.7.4"):
        slcfTypes = defaultdict(bool)
        slcfTypes['CELL_CENTERED'] = 'bool'
        slcfTypes['FYI'] = 'string'
        slcfTypes['ID'] = 'ignore'
        slcfTypes['PBX'] = 'float'
        slcfTypes['PBY'] = 'float'
        slcfTypes['PBZ'] = 'float'
        slcfTypes['QUANTITY'] = 'string'
        slcfTypes['SPEC_ID'] = 'string'
        slcfTypes['VECTOR'] = 'bool'
        slcfTypes['XB'] = 'listfloat'
        slcfTypes['unknownCounter'] = 'ignore'
        slcfTypes['newline'] = 'ignore'
        return slcfTypes

    def getCTRLtypes(self, version="6.7.4"):
        ctrlTypes = defaultdict(bool)
        ctrlTypes['ID'] = 'string'
        ctrlTypes['FUNCTION_TYPE'] = 'string'
        ctrlTypes['INITIAL_STATE'] = 'bool'
        ctrlTypes['INPUT_ID'] = 'liststring'
        ctrlTypes['DELAY'] = 'float'
        ctrlTypes['CONSTANT'] = 'float'
        ctrlTypes['LATCH'] = 'bool'
        ctrlTypes['RAMP_ID'] = 'string'
        ctrlTypes['newline'] = 'ignore'
        return ctrlTypes
    
    def getZONEtypes(self, version="6.7.4"):
        zoneTypes = defaultdict(bool)
        zoneTypes['ID'] = 'string'
        zoneTypes['XB'] = 'listfloat'
        zoneTypes['XYZ'] = 'listfloat'
        zoneTypes['LEAK_AREA'] = 'listindfloat'
        zoneTypes['newline'] = 'ignore'
        return zoneTypes
    
    def getDUMPtypes(self, version="6.7.4"):
        dumpTypes = defaultdict(bool)
        dumpTypes['COLUMN_DUMP_LIMIT'] = 'bool'
        dumpTypes['DT_CTRL'] = 'float'
        dumpTypes['DT_HRR'] = 'float'
        dumpTypes['DT_DEVC'] = 'float'
        dumpTypes['DT_BNDF'] = 'float'
        dumpTypes['DT_MASS'] = 'float'
        dumpTypes['DT_PART'] = 'float'
        dumpTypes['DT_PROF'] = 'float'
        dumpTypes['DT_SLCF'] = 'float'
        dumpTypes['DT_SL3D'] = 'float'
        dumpTypes['DT_PL3D'] = 'float'
        dumpTypes['DT_RESTART'] = 'float'
        dumpTypes['ID'] = 'ignore'
        dumpTypes['MASS_FILE'] = 'bool'
        dumpTypes['MMS_TIMER'] = 'float'
        dumpTypes['NFRAMES'] = 'int'
        dumpTypes['PLOT3D_QUANTITY'] = 'liststring'
        dumpTypes['PLOT3D_SPEC_ID'] = 'liststring'
        dumpTypes['RENDER_FILE'] = 'string'
        dumpTypes['SIG_FIGS'] = 'int'
        dumpTypes['SMOKE3D'] = 'bool'
        dumpTypes['WRITE_XYZ'] = 'bool'
        dumpTypes['newline'] = 'ignore'
        return dumpTypes
    
    def getTIMEtypes(self, version="6.7.4"):
        timeTypes = defaultdict(bool)
        timeTypes['DT'] = 'float'
        timeTypes['ID'] = 'ignore'
        timeTypes['LIMITING_DT_RATIO'] = 'float'
        timeTypes['LOCK_TIME_STEP'] = 'bool'
        timeTypes['T_BEGIN'] = 'float'
        timeTypes['T_END'] = 'float'
        timeTypes['WALL_INCREMENT'] = 'float'
        timeTypes['newline'] = 'ignore'
        return timeTypes
    
    def getMISCtypes(self, version="6.7.4"):
        miscTypes = defaultdict(bool)
        miscTypes['AEROSOL_SCRUBBING'] = 'bool'
        miscTypes['BNDF_DEFAULT'] = 'bool'
        miscTypes['CFL_MAX'] = 'float'
        miscTypes['CFL_VELOCITY_NORM'] = 'float'
        miscTypes['CHECK_HT'] = 'bool'
        miscTypes['CONSTANT_SPECIFIC_HEAT_RATIO'] = 'bool'
        miscTypes['DEPOSITION'] = 'bool'
        miscTypes['FLUX_LIMITER'] = 'string'
        miscTypes['FREEZE_VELOCITY'] = 'bool'
        miscTypes['GRAVITATIONAL_DEPOSITION'] = 'bool'
        miscTypes['GRAVITATIONAL_SETTLING'] = 'bool'
        miscTypes['GVEC'] = 'matrixfloat'
        miscTypes['HUMIDITY'] = 'float'
        miscTypes['ID'] = 'ignore'
        miscTypes['NOISE'] = 'bool'
        miscTypes['P_INF'] = 'float'
        miscTypes['PARTICLE_CFL'] = 'bool'
        miscTypes['PERIODIC_TEST'] = 'float'
        miscTypes['POROUS_FLOOR'] = 'bool'
        miscTypes['RESTART'] = 'bool'
        miscTypes['SIMULATION_MODE'] = 'string'
        miscTypes['SOLID_PHASE_ONLY'] = 'bool'
        miscTypes['SOOT_OXIDATION'] = 'bool'
        miscTypes['STRATIFICATION'] = 'bool'
        miscTypes['SUPPRESSION'] = 'bool'
        miscTypes['THERMOPHORETIC_DEPOSITION'] = 'bool'
        miscTypes['TMPA'] = 'float'
        miscTypes['TURBULENT_DEPOSITION'] = 'bool'
        miscTypes['VISIBILITY_FACTOR'] = 'float'
        miscTypes['VN_MAX'] = 'float'
        miscTypes['Y_CO2_INFTY'] = 'float'
        miscTypes['Y_O2_INFTY'] = 'float'
        miscTypes['newline'] = 'ignore'
        return miscTypes
    
    def getHEADtypes(self, version="6.7.4"):
        headTypes = defaultdict(bool)
        headTypes['CHID'] = 'string'
        headTypes['TITLE'] = 'string'
        headTypes['ID'] = 'ignore'
        headTypes['newline'] = 'ignore'
        return headTypes
