#-----------------------------------------------------------------------
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
        self.matl = self.getMATLtypes(version)
        self.reac = self.getREACtypes(version)
        self.radi = self.getRADItypes(version)
        self.mesh = self.getMESHtypes(version)
        self.surf = self.getSURFtypes(version)
        self.vent = self.getVENTtypes(version)
        self.pres = self.getPREStypes(version)
        self.obst = self.getOBSTtypes(version)
        self.hole = self.getHOLEtypes(version)
        self.init = self.getINITtypes(version)
        self.devc = self.getDEVCtypes(version)
        self.bndf = self.getBNDFtypes(version)
        self.slcf = self.getSLCFtypes(version)
        self.ctrl = self.getCTRLtypes(version)
        self.zone = self.getZONEtypes(version)
        self.dump = self.getDUMPtypes(version)
        self.time = self.getTIMEtypes(version)
        self.misc = self.getMISCtypes(version)
        self.head = self.getHEADtypes(version)
        self.ramp = self.getRAMPtypes(version)
        self.part = self.getPARTtypes(version)
        
    def getRAMPtypes(self, version="6.7.4"):
        rampTypes = defaultdict(bool)
        rampTypes['T'] = 'listrowfloat'
        rampTypes['F'] = 'listrowfloat'
        rampTypes['ID'] = 'string'
        return rampTypes
        
    def getMATLtypes(self, version="6.7.4"):
        matlTypes = defaultdict(bool)
        matlTypes['ID'] = 'string'
        matlTypes['FYI'] = 'string'
        matlTypes['SPECIFIC_HEAT'] = 'float'
        matlTypes['CONDUCTIVITY'] = 'float'
        matlTypes['DENSITY'] = 'float'
        matlTypes['EMISSIVITY'] = 'float'
        matlTypes['SPECIFIC_HEAT_RAMP'] = 'string'
        matlTypes['CONDUCTIVITY_RAMP'] = 'string'
        matlTypes['N_REACTIONS'] = 'int'
        matlTypes['NU_SPEC'] = 'float'
        matlTypes['SPEC_ID'] = 'string'
        matlTypes['REFERENCE_TEMPERATURE'] = 'float'
        matlTypes['HEAT_OF_REACTION'] = 'float'
        matlTypes['HEAT_OF_COMBUSTION'] = 'float'
        return matlTypes
    
    def getREACtypes(self, version="6.7.4"):
        reacTypes = defaultdict(bool)
        reacTypes['FUEL'] = 'string'
        reacTypes['ID'] = 'string'
        reacTypes['FORMULA'] = 'string'
        reacTypes['AUTO_IGNITION_TEMPERATURE'] = 'float'
        reacTypes['SOOT_YIELD'] = 'float'
        reacTypes['CO_YIELD'] = 'float'
        reacTypes['HEAT_OF_COMBUSTION'] = 'float'
        reacTypes['C'] = 'float'
        reacTypes['H'] = 'float'
        reacTypes['O'] = 'float'
        reacTypes['N'] = 'float'
        reacTypes['FYI'] = 'string'
        reacTypes['RADIATIVE_FRACTION'] = 'float'
        reacTypes['IDEAL'] = 'bool'
        return reacTypes

    def getRADItypes(self, version="6.7.4"):
        radiTypes = defaultdict(bool)
        radiTypes['RADIATIVE_FRACTION'] = 'float'
        return radiTypes

    def getMESHtypes(self, version="6.7.4"):
        meshTypes = defaultdict(bool)
        meshTypes['ID'] = 'string'
        meshTypes['IJK'] = 'listint'
        meshTypes['XB'] = 'listfloat'
        meshTypes['MPI_PROCESS'] = 'int'
        meshTypes['unknownCounter'] = 'ignore'
        return meshTypes

    def getSURFtypes(self, version="6.7.4"):
        surfTypes = defaultdict(bool)
        surfTypes['ID'] = 'string'
        surfTypes['MATL_ID'] = 'matrixstring'
        surfTypes['MATL_MASS_FRACTION'] = 'matrixfloat'
        surfTypes['THICKNESS'] = 'listindfloat'
        surfTypes['RGB'] = 'listfloat'
        surfTypes['COLOR'] = 'string'
        surfTypes['BACKING'] = 'string'
        surfTypes['GEOMETRY'] = 'string'
        surfTypes['FYI'] = 'string'
        surfTypes['LENGTH'] = 'float'
        surfTypes['WIDTH'] = 'float'
        surfTypes['C_FORCED_RE'] = 'float'
        surfTypes['C_FORCED_CONSTANT'] = 'float'
        surfTypes['C_FORCED_RE_EXP'] = 'float'
        surfTypes['C_FORCED_PR_EXP'] = 'float'
        surfTypes['CONVECTION_LENGTH_SCALE'] = 'float'
        surfTypes['C_HORIZONTAL'] = 'float'
        surfTypes['C_VERTICAL'] = 'float'
        surfTypes['EMISSIVITY'] = 'float'
        surfTypes['LEAK_PATH'] = 'listint'
        surfTypes['HRRPUA'] = 'float'
        surfTypes['RAMP_Q'] = 'string'
        surfTypes['ADIABATIC'] = 'bool'
        surfTypes['DEFAULT'] = 'bool'
        surfTypes['VOLUME_FLOW'] = 'float'
        surfTypes['VEL_T'] = 'listfloat'
        surfTypes['BURN_AWAY'] = 'bool'
        surfTypes['TMP_FRONT'] = 'float'
        surfTypes['RADIUS'] = 'float'
        surfTypes['TRANSPARENCY'] = 'float'
        return surfTypes

    def getVENTtypes(self, version="6.7.4"):
        surfTypes = defaultdict(bool)
        surfTypes['ID'] = 'string'
        surfTypes['XB'] = 'listfloat'
        surfTypes['MB'] = 'string'
        surfTypes['CTRL_ID'] = 'string'
        surfTypes['SURF_ID'] = 'string'
        surfTypes['IOR'] = 'int'
        surfTypes['number'] = 'ignore'
        surfTypes['unknownCounter'] = 'ignore'
        return surfTypes

    def getPREStypes(self, version="6.7.4"):
        presTypes = defaultdict(bool)
        presTypes['SOLVER'] = 'string'
        presTypes['VELOCITY_TOLERANCE'] = 'float'
        presTypes['MAX_PRESSURE_ITERATIONS'] = 'int'
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
        obstTypes['THICKEN'] = 'bool'
        obstTypes['TRANSPARENCY'] = 'float'
        obstTypes['COLOR'] = 'string'
        obstTypes['BULK_DENSITY'] = 'float'
        obstTypes['number'] = 'ignore'
        obstTypes['unknownCounter'] = 'ignore'
        return obstTypes
    
    def getPARTtypes(self, version="6.7.4"):
        partTypes = defaultdict(bool)
        partTypes['ID'] = 'string'
        partTypes['SURF_ID'] = 'string'
        partTypes['STATIC'] = 'bool'
        return partTypes
    
    def getHOLEtypes(self, version="6.7.4"):
        holeTypes = defaultdict(bool)
        holeTypes['ID'] = 'string'
        holeTypes['XB'] = 'listfloat'
        holeTypes['CTRL_ID'] = 'string'
        holeTypes['number'] = 'ignore'
        holeTypes['unknownCounter'] = 'ignore'
        return holeTypes

    def getINITtypes(self, version="6.7.4"):
        initTypes = defaultdict(bool)
        initTypes['ID'] = 'string'
        initTypes['XB'] = 'listfloat'
        initTypes['XYZ'] = 'listfloat'
        initTypes['PART_ID'] = 'string'
        initTypes['N_PARTICLES'] = 'int'
        initTypes['DX'] = 'float'
        initTypes['PART_ID'] = 'string'
        return initTypes

    def getDEVCtypes(self, version="6.7.4"):
        devcTypes = defaultdict(bool)
        devcTypes['ID'] = 'string'
        devcTypes['XB'] = 'listfloat'
        devcTypes['XYZ'] = 'listfloat'
        devcTypes['QUANTITY'] = 'string'
        devcTypes['INITIAL_STATE'] = 'bool'
        devcTypes['INIT_ID'] = 'string'
        devcTypes['SETPOINT'] = 'float'
        devcTypes['IOR'] = 'int'
        devcTypes['TIME_AVERAGED'] = 'bool'
        devcTypes['SPATIAL_STATISTIC'] = 'string'
        devcTypes['STATISTICS'] = 'string'
        devcTypes['DUCT_ID'] = 'string'
        devcTypes['SPEC_ID'] = 'string'
        devcTypes['unknownCounter'] = 'ignore'
        return devcTypes

    def getBNDFtypes(self, version="6.7.4"):
        bndfTypes = defaultdict(bool)
        bndfTypes['QUANTITY'] = 'string'
        bndfTypes['CELL_CENTERED'] = 'bool'
        bndfTypes['ID'] = 'ignore'
        bndfTypes['unknownCounter'] = 'ignore'
        return bndfTypes
    
    def getSLCFtypes(self, version="6.7.4"):
        slcfTypes = defaultdict(bool)
        slcfTypes['QUANTITY'] = 'string'
        slcfTypes['PBX'] = 'float'
        slcfTypes['PBY'] = 'float'
        slcfTypes['PBZ'] = 'float'
        slcfTypes['XB'] = 'listfloat'
        slcfTypes['VECTOR'] = 'bool'
        slcfTypes['SPEC_ID'] = 'string'
        slcfTypes['CELL_CENTERED'] = 'bool'
        slcfTypes['FYI'] = 'string'
        slcfTypes['ID'] = 'ignore'
        slcfTypes['unknownCounter'] = 'ignore'
        return slcfTypes

    def getCTRLtypes(self, version="6.7.4"):
        ctrlTypes = defaultdict(bool)
        ctrlTypes['ID'] = 'string'
        ctrlTypes['FUNCTION_TYPE'] = 'string'
        ctrlTypes['INPUT_ID'] = 'liststring'
        ctrlTypes['DELAY'] = 'float'
        return ctrlTypes
    
    def getZONEtypes(self, version="6.7.4"):
        zoneTypes = defaultdict(bool)
        zoneTypes['ID'] = 'string'
        zoneTypes['XB'] = 'listfloat'
        zoneTypes['LEAK_AREA'] = 'listindfloat'
        return zoneTypes
    
    def getDUMPtypes(self, version="6.7.4"):
        dumpTypes = defaultdict(bool)
        dumpTypes['DT_CTRL'] = 'float'
        dumpTypes['DT_HRR'] = 'float'
        dumpTypes['DT_DEVC'] = 'float'
        dumpTypes['DT_BNDF'] = 'float'
        dumpTypes['DT_SLCF'] = 'float'
        dumpTypes['DT_SL3D'] = 'float'
        dumpTypes['DT_PL3D'] = 'float'
        dumpTypes['DT_RESTART'] = 'float'
        dumpTypes['RENDER_FILE'] = 'string'
        dumpTypes['COLUMN_DUMP_LIMIT'] = 'bool'
        dumpTypes['WRITE_XYZ'] = 'bool'
        dumpTypes['ID'] = 'ignore'
        return dumpTypes
    
    def getTIMEtypes(self, version="6.7.4"):
        timeTypes = defaultdict(bool)
        timeTypes['T_BEGIN'] = 'float'
        timeTypes['T_END'] = 'float'
        timeTypes['ID'] = 'ignore'
        return timeTypes
    
    def getMISCtypes(self, version="6.7.4"):
        miscTypes = defaultdict(bool)
        miscTypes['TMPA'] = 'float'
        miscTypes['HUMIDITY'] = 'float'
        miscTypes['SUPPRESSION'] = 'bool'
        miscTypes['BNDF_DEFAULT'] = 'bool'
        miscTypes['ID'] = 'ignore'
        return miscTypes
    
    def getHEADtypes(self, version="6.7.4"):
        headTypes = defaultdict(bool)
        headTypes['CHID'] = 'string'
        headTypes['TITLE'] = 'string'
        headTypes['ID'] = 'ignore'
        return headTypes
