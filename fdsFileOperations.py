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
# EXAMPLES:
# See the examples subroutine for example operation.
#
#=======================================================================
# # IMPORTS
#=======================================================================

import numpy as np
import pandas as pd
from collections import defaultdict
import datetime
import re
import scipy.spatial as scsp
import os
import zipfile
from .fdsTypes import fdsLineTypes

class fdsFileOperations(object):
    """
    A class used to represent an FDS input file

    ...

    Attributes
    ----------
    bndfs : defaultdict
        dictionary containing each key in the bndf namelist
    ctrls : defaultdict
        dictionary containing each key in the ctrl namelist
    customLines : defaultdict
        dictionary containing custom lines to be added to the input file
    devcs : defaultdict
        dictionary containing each key in the devc namelist
    dump : defaultdict
        dictionary containing each key in the dump namelist
    head : defaultdict
        dictionary containing each key in the head namelist
    holes : defaultdict
        dictionary containing each key in the hole namelist
    inits : defaultdict
        dictionary containing each key in the init namelist
    matls : defaultdict
        dictionary containing each key in the matl namelist
    meshes : defaultdict
        dictionary containing each key in the mesh namelist
    meshOrder : defaultdict
        dictionary containing the order meshes are to be defined in the
        input file. This is an intermediate variable used after
        assigning mpi processes to meshes.
    misc : defaultdict
        dictionary containing each key in the misc namelist
    mpiProcesses : int
        integer number of mpi processes to use when building the fds
        input file.
    obsts : defaultdict
        dictionary containing each key in the obst namelist
    pres : defaultdict
        dictionary containing each key in the pres namelist
    props : defaultdict
        dictionary containing each key in the prop namelist
    radis : defaultdict
        dictionary containing each key in the radi namelist
    ramps : defaultdict
        dictionary containing each key in the ramp namelist
    reacs : defaultdict
        dictionary containing each key in the reac namelist
    slcfs : defaultdict
        dictionary containing each key in the slcf namelist
    specs : defaultdict
        dictionary containing each key in the spec namelist
    surfs : defaultdict
        dictionary containing each key in the surf namelist
    time : defaultdict
        dictionary containing each key in the time namelist
    vents : defaultdict
        dictionary containing each key in the vent namelist
    version : str
        string containing the fds version for the input file.
        Syntax is '#.#.#'. Currently supports 6.7.1 and 6.7.4.
    zones : defaultdict
        dictionary containing each key in the zone namelist


    Methods
    -------
    addBNDF(Qty, CELL_CENTERED=None)
        Adds a bndf key to the bndfs namelist.
    addCTRL(ID, FUNCTION_TYPE, INPUT_ID, DELAY=None)
        Adds a ctrl key to the ctrls namelist.
    addDEVC(ID, QUANTITY, XYZ=None, XB=None, IOR=None, SPEC_ID=None,
            TIME_AVERAGED=None, SPATIAL_STATISTIC=None, STATISTICS=None,
            INITIAL_STATE=None, INIT_ID=None, SETPOINT=None,
            DUCT_ID=None, PROP_ID=None)
        Adds a devc key to the devcs namelist.
    addDUMP(RENDER_FILE=None, COLUMN_DUMP_LIMIT=False, WRITE_XYZ=False,
            DT_PL3D=None, DT_SL3D=None, DT_SLCF=None, DT_BNDF=None,
            DT_DEVC=None, DT_CTRL=None, DT_HRR=None, DT_RESTART=None)
        Adds a dump key to the dump namelist.
    addHEAD(chid, title=None)
        Adds a head key to the head namelist.
    addHOLE(ID, XB)
        Adds a hole key to the holes namelist.
    addMATL(ID, Emi=None, Den=None, Con=None, Spe=None, kramp=None,
            cpramp=None, fyi=None)
        Adds a matl key to the matls namelist.
    addMESH(ID, IJK, XB)
        Adds a mesh key to the meshes namelist.
    addMISC(BNDF_DEFAULT=None, TMPA=None)
        Adds a misc key to the misc namelist.
    addMPIprocesses(numberOfProcesses, allowMeshSplitting=True,
                    splitMultiplier=1.20,
                    meshSplitAxes=[True, True, False])
        Adds mpi processes to meshes. Can be used to automatically
        split meshes to balance load on mpi processes.
    addOBST(ID, XB, SURF_IDS=None, SURF_ID=None, SURF_ID6=None,
            BNDF_OBST=True, THICKEN=None, TRANSPARENCY=None, COLOR=None)
        Adds obst key to the obsts namelist.
    addPRES(VELOCITY_TOLERANCE=None, MAX_PRESSURE_ITERATIONS=None)
        Adds pres keys to the pres namelist.
    addRAMP(ID, T, F)
        Adds ramp keys to the ramps namelist.
    addREAC(ID, FUEL=None, FORMULA=None, AIT=None, SY=None, COY=None,
            HOC=None, C=None, H=None, O=None, N=None, FYI=None, RF=None)
        Adds reac keys to the reacs namelist.
    addSLCF(Qty, PBX=None, PBY=None, PBZ=None,
            Vec=False, XB=None, SPEC_ID=None)
        Adds slcf key to the slcfs namelist.
    addSURF(ID, Mid=None, Col=None, Thi=None, Bac=None, Geo=None,
            Fyi=None, Len=None, LeaPat=None, Hrrpua=None, qramp=None,
            Rgb=None, adiabatic=False, VOLUME_FLOW=None, VEL_T=None)
        Adds surf key to the surfs namelist.
    addTIME(T_END=0.0, T_BEGIN=0.0)
        Adds time key to the times namelist.
    addVENT(ID, SURF_ID, XB=None, CTRL_ID=None, MB=None, IOR=None)
        Adds vent key to the vents namelist.
    addZONE(ID, XB, LEAK_AREA=None)
        Adds zone key to the zones namelist.
    calculateMeshCells()
        Returns a list of mesh keys and number of cells in each mesh.
    checkOverlappingMESH()
        Returns True if any meshes are overlapping else False
    dictFromLine(line, lineType, types)
        Returns a dictionary with keys and values from a namelist line.
    dictMerge(template, master, path=None)
        Returns merged dictionary where keys in master overwrite keys
        in template.
    generateFDStext()
        Returns str of input file.
    getDefaultFields()
        Returns default field order.
    getLineType(line)
        Returns namelist key from str line.
    getMeshLimits()
        Returns a dictionary containing a key 'XB' with an array of the
        total extents defined in meshes.
    getNewlineFromTypes()
        Returns a dictionary containing default new line parameters.
    getPolygonNamesFromFdsFile()
        Returns a list of polygons defined in the fds input file.
    importFile(file=None, text=None, textList=None)
        Adds keys to each namelist from an input file, text, or text
        list.
    interpretKey(key, lineType, types)
        Intermediate function which processes a key from a namelist
        key pair to returns the keyID, keyType, and keyValue.
    keyAssist(text, types, dic, internalKeys=['counter'], newline=False)
        Returns a namelist text line based on an input dictionary and
        type dictionary.
    keyFromLineType(lineType)
        Returns internal attribute name from namelist type.
    makeFDSLines(textFDS)
        Returns a list of namelist lines.
    makeLinesFromDict(items, types, prefix, newline=False)
        Returns a str generated from a namelist dictionary.
    makeMESH(meshes, meshTypes, meshOrder=False)
        Returns a str generated from a meshes namelist dictionary.
    makeRAMP(ramps)
        Returns a str generated from a ramps namelist dictionary.
    mergeTypeFromLineType(lineType)
        Returns internal merge type based on namelist type.
    parseFDSLines(lines)
        Adds each line to internal attribute namelist dictionaries.
    parseLine(line, lineType, types, key)
        Adds one line to the corresponding internal attribute namelist
        dictionary.
    saveModel(mpiProcesses, location, allowMeshSplitting=True,
              splitMultiplier=1.2)
        Saves an fds input file based on internal attribute namelist
        dictionaries. Allows splitting of meshes to optimize mpi
        processes balance.
    splitLineIntoKeys(line2)
        Returns namelist key pairs from a line.
    splitMESHonce(mesh)
        Splits a mesh along its largest axis.
    zopen(file)
        Opens a file or zip archive for reading.
    """
    
    def __init__(self, version="6.7.4"):
        """
        Parameters
        ----------
        version : str
            string containing the fds version for the input file.
            Syntax is '#.#.#'. Currently supports 6.7.1 and 6.7.4.
        """
        
        self.head = defaultdict(bool)
        self.devcs = defaultdict(bool)
        self.inits = defaultdict(bool)
        self.obsts = defaultdict(bool)
        self.holes = defaultdict(bool)
        self.vents = defaultdict(bool)
        self.surfs = defaultdict(bool)
        self.ramps = defaultdict(bool)
        self.ctrls = defaultdict(bool)
        self.meshes = defaultdict(bool)
        self.slcfs = defaultdict(bool)
        self.bndfs = defaultdict(bool)
        self.time = defaultdict(bool)
        self.dump = defaultdict(bool)
        self.misc = defaultdict(bool)
        self.zones = defaultdict(bool)
        self.reacs = defaultdict(bool)
        self.matls = defaultdict(bool)
        self.radis = defaultdict(bool)
        self.pres = defaultdict(bool)
        self.parts = defaultdict(bool)
        self.profs = defaultdict(bool)
        self.props = defaultdict(bool)
        self.specs = defaultdict(bool)
        self.winds = defaultdict(bool)
        self.customLines = []
        
        self.devcs['unknownCounter'] = 0
        self.obsts['unknownCounter'] = 0
        self.holes['unknownCounter'] = 0
        self.vents['unknownCounter'] = 0
        self.meshes['unknownCounter'] = 0
        self.slcfs['unknownCounter'] = 0
        self.bndfs['unknownCounter'] = 0
        self.profs['unknownCounter'] = 0
        
        self.meshOrder = False
        self.version = version
    
    
    def addBNDF(self, QUANTITY, CELL_CENTERED=None):
        """Adds a bndf key to internal attribute bndfs
        
        Adds a bndf key to internal attribte bndfs. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        QUANTITY : str
            Quantity of the bndf
        CELL_CENTERED : bool, optional
            Flag specifying if the quantity should be exported at cell
            centers or at cell edges (default None).
        """
        
        bndf = defaultdict(bool)
        bndf['ID'] = "BNDF-%05.0f"%(self.bndfs['unknownCounter'])
        bndf['QUANTITY'] = QUANTITY
        if CELL_CENTERED != None: bndf['CELL_CENTERED'] = CELL_CENTERED
        self.bndfs['unknownCounter'] += 1
        self.bndfs[bndf['ID']] = bndf
    
    
    def addCTRL(self, ID, FUNCTION_TYPE, INPUT_ID, DELAY=None,
                CONSTANT=None, RAMP_ID=None):
        """Adds a ctrl key to internal attribute ctrls
        
        Adds a bndf key to internal attribte ctrls. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            Identifier for this control
        FUNCTION_TYPE : str
            Identify for type of control.
            Valid entries are: ANY, ALL
        INPUT_ID : str
            Identifier for device or control for logic.
        DELAY : float, optional
            Time delay for activation of control (default None)
        CONSTANT : float, optional
            Value for constant defined in input id
        RAMP_ID : str, optional
            Name of ramp to be used to map control output
        """
        
        ctrl = defaultdict(bool)
        ctrl['ID'] = ID
        ctrl['FUNCTION_TYPE'] = FUNCTION_TYPE
        ctrl['INPUT_ID'] = INPUT_ID
        if DELAY != None: ctrl['DELAY'] = DELAY
        if CONSTANT != None: ctrl['CONSTANT'] = CONSTANT
        if RAMP_ID != None: ctrl['RAMP_ID'] = RAMP_ID
        self.ctrls[ID] = ctrl
    
    
    def addDEVC(self, ID, QUANTITY, XYZ=None, XB=None, IOR=None,
                SPEC_ID=None, TIME_AVERAGED=None,
                SPATIAL_STATISTIC=None, STATISTICS=None,
                INITIAL_STATE=None, INIT_ID=None, SETPOINT=None,
                DUCT_ID=None, NO_UPDATE_DEVC_ID=None, CTRL_ID=None,
                PROP_ID=None):
        """Adds a devc key to internal attribute devcs
        
        Adds a devc key to internal attribte devcs. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            Identifier for this device
        QUANTITY : str
            Quantity of the device
        XYZ : float array(3), optional
            Three component array containing X, Y, Z coordinates
            (default None)
        XB : float array(6), optional
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates (default None)
        IOR : int, optional
            Integer specifying the orientation of the device
            (default None)
        SPEC_ID : str, optional
            String specifying the species of the device (default None)
        TIME_AVERAGED : bool, optional
            Flag specifying if the device is time averaged
            (default None)
        SPATIAL_STATISTIC : str, optional
            String specifying spatial statistic of the device
            (default None)
        STATISTICS : str, optional
            String specifying statistic type
        INITIAL_STATE : bool, optional
            Flag specifying if device is initially active (defualt None)
        INIT_ID : str, optional
            String specifying init namelist identifier
        SETPOINT : float, optional
            Flag value used to determine activation of device
            (default None)
        DUCT_ID : str, optional
            String identifier of duct containing device
        NO_UPDATE_DEVC_ID : str, optional
            String identifier of device activation to stop updating
        CTRL_ID : str, optional
            String identifier of control for device
        PROP_ID : str, optional
            String identifier of properties for device
        """
        
        devc = defaultdict(bool)
        devc['ID'] = ID
        devc['QUANTITY'] = QUANTITY
        if XYZ != None:
            if type(XYZ) is list: XYZ = np.array(XYZ)
            devc['XYZ'] = XYZ
        if XB != None:
            if type(XB) is list: XB = np.array(XB)
            devc['XB'] = XB
        if INITIAL_STATE != None: devc['INITIAL_STATE'] = INITIAL_STATE
        if INIT_ID != None: devc['INIT_ID'] = INIT_ID
        if SETPOINT != None: devc['SETPOINT'] = SETPOINT
        if IOR != None: devc['IOR'] = IOR
        if TIME_AVERAGED != None: devc['TIME_AVERAGED'] = TIME_AVERAGED
        if SPATIAL_STATISTIC != None:
            devc['SPATIAL_STATISTIC'] = SPATIAL_STATISTIC
        if STATISTICS != None: devc["STATISTICS"] = STATISTICS
        if DUCT_ID != None: devc['DUCT_ID'] = DUCT_ID
        if SPEC_ID != None: devc['SPEC_ID'] = SPEC_ID
        if NO_UPDATE_DEVC_ID != None: devc['NO_UPDATE_DEVC_ID'] = NO_UPDATE_DEVC_ID
        if CTRL_ID != None: devc['CTRL_ID'] = CTRL_ID
        if SETPOINT != None: devc['SETPOINT'] = SETPOINT
        if PROP_ID != None: devc['PROP_ID'] = PROP_ID
        self.devcs[ID] = devc
    
    
    def addDUMP(self, RENDER_FILE=None, COLUMN_DUMP_LIMIT=False,
                WRITE_XYZ=False, DT_PL3D=None, DT_SL3D=None,
                DT_SLCF=None, DT_BNDF=None, DT_DEVC=None, DT_CTRL=None,
                DT_HRR=None, DT_RESTART=None):
        """Adds a dump key to internal attribute dumps
        
        Adds a dump key to internal attribute dumps. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        RENDER_FILE : str, optional
            Filename for render file (default None)
        COLUMN_DUMP_LIMIT : bool, optional
            Flag specifying if number of columns in CSV file should be
            limited (default False)
        WRITE_XYZ : bool, optional
            Flag specifying if an XYZ file should be generated by FDS
            (default False)
        DT_PL3D : float, optional
            Time interval to output PL3D data (default None)
        DT_SL3D : float, optional
            Time interval to output SL3D data (default None)
        DT_SLCF : float, optional
            Time interval to output SLCF data (default None)
        DT_BNDF : float, optional
            Time interval to output BNDF data (default None)
        DT_DEVC : float, optional
            Time interval to output DEVC data (default None)
        DT_CTRL : float, optional
            Time interval to output CTRL data (default None)
        DT_HRR : float, optional
            Time interval to output HRR data (default None)
        DT_RESTART : float, optional
            Time interval to save restart files (default None)
        """
        
        dump = defaultdict(bool)
        if RENDER_FILE != None: dump['RENDER_FILE'] = RENDER_FILE
        if COLUMN_DUMP_LIMIT:
            dump['COLUMN_DUMP_LIMIT'] = COLUMN_DUMP_LIMIT
        if WRITE_XYZ: dump['WRITE_XYZ'] = WRITE_XYZ
        if DT_PL3D != None: dump['DT_PL3D'] = DT_PL3D
        if DT_SL3D != None: dump['DT_SL3D'] = DT_SL3D
        if DT_SLCF != None: dump['DT_SLCF'] = DT_SLCF
        if DT_BNDF != None: dump['DT_BNDF'] = DT_BNDF
        if DT_DEVC != None: dump['DT_DEVC'] = DT_DEVC
        if DT_CTRL != None: dump['DT_CTRL'] = DT_CTRL
        if DT_HRR != None: dump['DT_HRR'] = DT_HRR
        if DT_RESTART != None: dump['DT_RESTART'] = DT_RESTART
        self.dump['ID'] = dump
    
    
    def addHEAD(self, chid, title=None):
        """Adds a head key to internal attribute head
        
        Adds a head key to internal attribute head. Note if no title is
        specified, title will be set to the same as chid.
        
        Parameters
        ----------
        chid: str
            Chid for use in the input file
        title: str, optional
            Title for use in the input file (default None)
        """
        
        head = defaultdict(bool)
        head['CHID'] = chid
        if title != None:
            head['TITLE'] = title
        else:
            head['TITLE'] = chid
        self.head['ID'] = head
        
        
    def addHOLE(self, ID, XB):
        """Adds a hole key to internal attribute holes
        
        Adds a hole key to internal attribute holes. 
        
        Parameters
        ----------
        ID : str
            String identifier for the hole
        XB : float array(6)
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates
        """
        
        hole = defaultdict(bool)
        hole['XB'] = XB
        hole['ID'] = ID
        self.holes[ID] = hole
    
    
    def addMATL(self, ID, Emi=None, Den=None, Con=None, Spe=None,
                kramp=None, cpramp=None, fyi=None):
        """Adds a matl key to internal attribute matls
        
        Adds a matl key to internal attribute matls. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            String identifier for the material
        Emi : float, optional
            Emissivity of the material (default None)
        Den : float, optional
            Density of the material (default None)
        Con : float, optional
            Conductivity of the material (default None)
        Spe : float, optional
        kramp : str, optional
            String identifier of thermal conductivity ramp
            (default None)
        cpramp : str, optional
            String identifier of specific heat capacity ramp
            (default None)
        fyi : str, optional
            String containing comment field to be included in input file
            (default None)
        """
        
        matl = defaultdict(bool)
        matl['ID'] = ID
        if Emi != None: matl['EMISSIVITY'] = Emi
        if Den != None: matl['DENSITY'] = Den
        if Con != None: matl['CONDUCTIVITY'] = Con
        if Spe != None: matl['SPECIFIC_HEAT'] = Spe
        if kramp != None: matl['CONDUCTIVITY_RAMP'] = kramp
        if cpramp != None: matl['SPECIFIC_HEAT_RAMP'] = cpramp
        if fyi != None: matl['FYI'] = fyi
        self.matls[ID] = matl
        
        
    def addMESH(self, ID, IJK, XB):
        """Adds a mesh key to internal attribute meshes
        
        Adds a mesh key to internal attribute meshes.
        
        Parameters
        ----------
        ID : str
            String identifier for the mesh
        IJK : int array(3)
            Three component array containing number of grid cells in
            each axis
        XB : float array(6)
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates
        """
        
        mesh = defaultdict(bool)
        if type(IJK) is list: IJK = np.array(IJK)
        if type(XB) is list: XB = np.array(XB)
        mesh['ID'] = ID
        mesh['IJK'] = IJK
        mesh['XB'] = XB
        self.meshes[ID] = mesh
    
    
    def addMISC(self, BNDF_DEFAULT=None, TMPA=None):
        """Adds a misc key to internal attribute misc
        
        Adds a misc key to internal attribute misc. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        BNDF_DEFAULT : bool
            Flag specifying if boundary data is to be output for all
            boundary surfaces by default (default None)
        TMPA : float
            Ambient air temperature
        """
        
        misc = defaultdict(bool)
        if BNDF_DEFAULT != None: misc['BNDF_DEFAULT'] = BNDF_DEFAULT
        if TMPA != None: misc['TMPA'] = TMPA
        self.misc['ID'] = misc
    
    
    def calculateCellsPerProcess(self):
        """Calculates the number of cells per mpi process based on
        information stored in internal attributes
        """
        meshes, numCells = self.calculateMeshCells()
        numProcesses = self.mpiProcesses
        IdealCellsPerProcess = np.sum(numCells)/numProcesses
        
        cellsPerProcess = np.zeros((numProcesses,))
        for i, mesh in enumerate(meshes):
            process = int(self.meshes[mesh]['MPI_PROCESS'])
            cellsPerProcess[process] += numCells[i]
        return IdealCellsPerProcess, cellsPerProcess
    
    def addMPIprocesses(self, numberOfProcesses,
                        allowMeshSplitting=True, splitMultiplier=1.20,
                        meshSplitAxes=[True, True, False]):
        """Adds mpi processes to meshes stored in internal attributes
        
        Adds mpi processes to meshes stored in internal attributes.
        Can be used to automatically split meshes to balance load on
        mpi processes.
        
        Parameters
        ----------
        numberOfProcesses : int
            Number of mpi processes
        allowMeshSplitting : bool
            Flag specifying whether meshes can be split
        splitMultiplier : float
            Threshold used in splitting meshes
        meshSplitAxes : list of booleans
            Specifies along which axes the software is allowed to split
            meshes.
        """
        
        meshes, numCells = self.calculateMeshCells()
        cellsPerProcess = np.sum(numCells)/numberOfProcesses
        mpiConverged = False
        splitConverged = False
        assumedConverged = False
        while not mpiConverged and not assumedConverged:
            mpiConverged = True
            while not splitConverged and allowMeshSplitting:
                splitConverged = True
                meshes, numCells = self.calculateMeshCells()
                for mesh, numCell in zip(meshes, numCells):
                    if numCell > cellsPerProcess*splitMultiplier:
                        self.splitMESHonce(self.meshes[mesh], meshSplitAxes)
                        splitConverged = False
            
            meshes, numCells = self.calculateMeshCells()
            #print(len(meshes), numberOfProcesses)
            if len(meshes) / 10 > numberOfProcesses:
                print("Warning: Number of meshes 10x greater than number of requested processes (%0.0f, %0.0f)"%(len(meshes), numberOfProcesses))
                print("AssumingConvergence")
                assumedConverged = True
            mpiProcessInds = np.zeros((len(numCells),))-1
            mpiProcess = np.zeros((numberOfProcesses,))
            while np.argwhere(mpiProcessInds == -1).shape[0] > 0:
                ind = np.argmax(numCells)
                ind2 = np.argmin(mpiProcess)
                mpiProcessInds[ind] = ind2
                mpiProcess[ind2] += numCells[ind]
                numCells[ind] = 0
            if np.max(mpiProcess) > cellsPerProcess*splitMultiplier and allowMeshSplitting:
                mpiConverged = False
                splitConverged = False
                splitMultiplier = splitMultiplier*0.9
        for key, mp in zip(meshes, mpiProcessInds):
            self.meshes[key]['MPI_PROCESS'] = mp
        self.mpiProcesses = numberOfProcesses
        self.meshOrder = np.argsort(mpiProcessInds)
        
        
    def addOBST(self, ID, XB, SURF_IDS=None, SURF_ID=None,
                SURF_ID6=None, BNDF_OBST=True, THICKEN=None,
                TRANSPARENCY=None, COLOR=None):
        """Adds an obst key to internal attribute obsts
        
        Adds an obst key to internal attribute obsts. Optional
        parameters that are specified as None will not be explicitly
        specified in a generated input file. These values at runtime
        will default to current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            String identifier for the obstruction
        XB : float array(6)
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates
        SURF_IDS : str array(3), optional
            Three component array specifying surface definition
            (default None)
        SURF_ID : str, optional
            String specifing surface for all faces
        SURF_ID6 : str array(6), optional
            Six component array specifying surface definition for
            X-, X+, Y-, Y+, Z-, Z+ (default None)
        BNDF_OBST : bool
            Flag specifying if boundary data is to be output for all
            faces of this obstruction (default True)
        THICKEN : bool
            Flag specifying if obstruction is to be thickened to be at
            least one grid cell thick (default None)
        TRANSPARENCY : float
            Value specifying how transparent this obstruction should be
            in visualization (default None)
        COLOR : str
            String specifiying a color for the obstruction
        """
        
        obst = defaultdict(bool)
        obst['XB'] = XB
        obst['ID'] = ID
        if SURF_IDS != None: obst['SURF_IDS'] = SURF_IDS
        if SURF_ID != None: obst['SURF_ID'] = SURF_ID
        if SURF_ID6 != None: obst['SURF_ID6'] = SURF_ID6
        if BNDF_OBST: obst['BNDF_OBST'] = True
        if THICKEN != None: obst['THICKEN'] = THICKEN
        if TRANSPARENCY != None: obst['TRANSPARENCY'] = TRANSPARENCY
        if COLOR != None: obst['COLOR'] = COLOR
        if self.obsts[ID]:
            counter = self.obsts[ID]['counter']
            counter += 1
            self.obsts["%s-%0.0f"%(ID, counter)] = obst
            self.obsts[ID]['counter'] = counter
        else:
            obst['counter'] = 0
            self.obsts[ID] = obst
        
        
    def addPRES(self, VELOCITY_TOLERANCE=None,
                MAX_PRESSURE_ITERATIONS=None):
        """Adds a pres key to internal attribute pres
        
        Adds a pres key to internal attribute pres. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        VELOCITY_TOLERANCE : float
            Value for the velocity error tolerance
        MAX_PRESSURE_ITERATIONS : int
            Maxmium number of iterations allowed in the pressure solver
        """
        
        pres = defaultdict(bool)
        if VELOCITY_TOLERANCE != None:
            pres['VELOCITY_TOLERANCE'] = VELOCITY_TOLERANCE
        if MAX_PRESSURE_ITERATIONS != None:
            pres['MAX_PRESSURE_ITERATIONS'] = MAX_PRESSURE_ITERATIONS
        self.pres['ID'] = pres
        
        
    def addRAMP(self, ID, T, F, appendZero=False, appendTime=1.0):
        """Adds a ramp key to internal attribute ramps
        
        Adds a ramp key to internal attribute ramps.
        
        Parameters
        ----------
        ID : str
            String identifier for the obstruction
        T : float array(N)
            Array specifying the x-axis of the ramp
        F : float array(N)
            Array specifying the y-axis of the ramp
        """
        if type(T) == pd.core.frame.DataFrame: T = T.values
        if type(T) == pd.core.series.Series: T = T.values
        if type(T) == np.ndarray: T = list(T)
        if type(F) == pd.core.frame.DataFrame: F = F.values
        if type(F) == pd.core.series.Series: F = F.values
        if type(F) == np.ndarray: F = list(F)
        if appendZero: 
            T.append(T[-1] + appendTime)
            F.append(0)
        if self.ramps[ID]:
            Ts = self.ramps[ID]['T']
            Fs = self.ramps[ID]['F']
            for t, f in zip(T, F):
                Ts.append(t)
                Fs.append(f)
            self.ramps[ID]['T'] = Ts
            self.ramps[ID]['F'] = Fs
        else:
            self.ramps[ID] = defaultdict(bool)
            self.ramps[ID]['T'] = T
            self.ramps[ID]['F'] = F
            self.ramps[ID]['ID'] = ID
        
        
    def addREAC(self, ID, FUEL=None, FORMULA=None, AIT=None, SY=None, 
                COY=None, HOC=None,
                C=None, H=None, O=None, N=None, FYI=None, RF=None):
        """Adds a reac key to internal attribute reacs
        
        Adds a reac key to internal attribute reacs. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            String identifier for the reaction
        FUEL : str, optional
            String name of the fuel in the reaction (default None)
        FORMULA : str, optional
            String formula of the reaction (default None)
        AIT : float, optional
            Float auto ignition temperature of the reaction
            (default None)
        SY : float, optional
            Float soot yield of the reaction (default None)
        COY : float, optional
            Float carbon monoxide yield of the reaction (default None)
        HOC : float, optional
            Float heat of combustion of the reaction (default None)
        C : float, optional
            Float number of carbon atoms in the chemical formula of the
            reaction (default None)
        H : float, optional
            Float number of hydrogen atoms in the chemical formula of
            the reaction (default None)
        O : float, optional
            Float number of oxygen atoms in the chemical formula of the
            reaction (default None)
        N : float, optional
            Float number of nitrogen atoms in the chemical formula of
            the reaction (default None)
        FYI : string, optional
            String containing comment field to be included in input file
        RF : float, optional
            Float radiative fraction of the reaction (default None)
        """
        
        reac = defaultdict(bool)
        reac['ID'] = ID
        if FUEL != None: reac['FUEL'] = FUEL
        if FORMULA != None: reac['FORMULA'] = FORMULA
        if AIT != None: reac['AUTO_IGNITION_TEMPERATURE'] = AIT
        if SY != None: reac['SOOT_YIELD'] = SY
        if COY != None: reac['CO_YIELD'] = COY
        if HOC != None: reac['HEAT_OF_COMBUSTION'] = HOC
        if C != None: reac['C'] = C
        if H != None: reac['H'] = H
        if O != None: reac['O'] = O
        if N != None: reac['N'] = N
        if FYI != None: reac['FYI'] = FYI
        if RF != None: reac['RADIATIVE_FRACTION'] = RF
        self.reacs[ID] = reac
        
        
    def addSLCF(self, QUANTITY, PBX=None, PBY=None, PBZ=None, Vec=False, XB=None, SPEC_ID=None):
        """Adds a slcf key to internal attribute slcfs
        
        Adds a slcf key to internal attribute slcfs. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        QUANTITY : str
            Quantity of the slice
        PBX : float, optional
            Value along x-axis of the plane (default None)
        PBY : float, optional
            Value along y-axis of the plane (default None)
        PBZ : float, optional
            Value along z-axis of the plane (default None)
        Vec : bool, optional
            Flag specifying if the slice is a vector slice
            (default False)
        XB : float array(6), optional
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates
        SPEC_ID : str, optional
            String specifying the species of the slice
        """
        
        slcf = defaultdict(bool)
        slcf['ID'] = "SLCF-%05.0f"%(self.slcfs['unknownCounter'])
        slcf['QUANTITY'] = QUANTITY
        if PBX != None: slcf['PBX'] = PBX
        if PBY != None: slcf['PBY'] = PBY
        if PBZ != None: slcf['PBZ'] = PBZ
        if SPEC_ID != None: slcf['SPEC_ID'] = SPEC_ID
        if Vec: slcf['VECTOR'] = 'TRUE'
        if XB != None:
            if type(XB) is list: XB = np.array(XB)
            slcf['XB'] = XB
        self.slcfs['unknownCounter'] += 1
        self.slcfs[slcf['ID']] = slcf
        
        
    def addSURF(self, ID, Mid=None, Col=None, Thi=None, Bac=None,
                Geo=None, Fyi=None, Len=None, LeaPat=None, Hrrpua=None,
                qramp=None, Rgb=None, adiabatic=False, VOLUME_FLOW=None,
                VEL_T=None):
        """Adds a surf key to internal attribute surfs
        
        Adds a surf key to internal attribute surfs. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            String identifier for the surface
        Mid : str array(N), optional
            Array of material IDs in the surface (default None)
        Col : str, optional
            String specifying the color of the surface (default None)
        Thi : float array(N), optional
            Array of floats specifying the thickness of each material
            in the surface (default None)
        Bac : str, optional
            String specifying the type of back boundary condition
            (default None)
        Geo : str, optional
            String specifying the type of geometry to use for the
            surface (default None)
        Fyi : str, optional
            String containing comment field to be included in input file
        Len : float, optional
            Value of length to be used in heat transfer calculation
            (default None)
        LeaPat : array(2), optional
            Array specifying leak path for the surface
        HRRPUA : float, optional
            Value of heat release rate per unit area of the surface
            (default None)
        qramp : str, optional
            String identifier of ramp for the heat release rate per unit
            area (default None)
        Rgb : float array(3), optional
            Array specifying the color of the surface (default None)
        adiabatic : bool, optional
            Flag specifying if the surface is adiabatic (default False)
        VOLUME_FLOW : float, optional
            Value of specified volume flow from the surface
            (default None)
        VEL_T : float, optional
            Value of specified tangential velocity from the surface
            (default None)
        """
        
        surf = defaultdict(bool)
        surf['ID'] = ID
        if Mid != None: surf['MATL_ID'] = Mid
        if Col != None: surf['COLOR'] = Col
        if Thi != None: surf['THICKNESS'] = Thi
        if Bac != None: surf['BACKING'] = Bac
        if Geo != None: surf['GEOMETRY'] = Geo
        if Fyi != None: surf['FYI'] = Fyi
        if Len != None: surf['LENGTH'] = Len
        if LeaPat != None: surf['LEAK_PATH'] = LeaPat
        if Hrrpua != None: surf['HRRPUA'] = Hrrpua
        if qramp != None: surf['RAMP_Q'] = qramp
        if Rgb != None: surf['RGB'] = Rgb
        if adiabatic: surf['ADIABATIC'] = True
        if VOLUME_FLOW != None: surf['VOLUME_FLOW'] = VOLUME_FLOW
        if VEL_T != None: surf['VEL_T'] = VEL_T
        self.surfs[ID] = surf
        
        
    def addTIME(self, T_END=0.0, T_BEGIN=0.0):
        """Adds a time key to internal attribute time
        
        Adds a time key to internal attribute time. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        T_END : float, optional
            Time to end the simulation (default None)
        T_BEGIN : float, optional
            Time to begin the simulation (default None)
        """
        
        time = defaultdict(bool)
        time['T_BEGIN'] = T_BEGIN
        time['T_END'] = T_END
        self.time['ID'] = time
        
        
    def addVENT(self, ID, SURF_ID, XB=None, CTRL_ID=None, MB=None,
                IOR=None):
        """Adds a vent key to internal attribute vents
        
        Adds a vent key to internal attribute vents. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            String identifier for the vent
        SURF_ID : str
            String identifier specifying the surface of the vent
        XB : float array(6), optional
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates (default None)
        CTRL_ID : str, optional
            String identifier for control determining if the vent is
            active (default None)
        MB : str, optional
            String specifying short-hand position of axis (default None)
        IOR : int, optional
            Integer specifying the orientation of the vent
            (default None)
        """
        
        vent = defaultdict(bool)
        vent['ID'] = ID
        vent['SURF_ID'] = SURF_ID
        if XB is not None:
            if type(XB) is list: XB = np.array(XB)
            vent['XB'] = XB
        if CTRL_ID != None: vent['CTRL_ID'] = CTRL_ID
        if MB != None: vent['MB'] = MB
        if IOR != None: vent['IOR'] = IOR
        if self.vents[ID]:
            counter = self.vents[ID]['counter']
            counter += 1
            self.vents["%s-%0.0f"%(ID, counter)] = vent
            self.vents[ID]['counter'] = counter
        else:
            vent['counter'] = 0
            self.vents[ID] = vent
            
            
    def addZONE(self, ID, XB, LEAK_AREA=None):
        """Adds a zone key to internal attribute zones
        
        Adds a zone key to internal attribute zones. Optional parameters
        that are specified as None will not be explicitly specified in
        a generated input file. These values at runtime will default to
        current FDS default parameters.
        
        Parameters
        ----------
        ID : str
            String identifier for the zone
        XB : float array(6)
            Six component array containing X_min, X_max, Y_min, Y_max,
            Z_min, Z_max coordinates (default None)
        LEAK_AREA : float array(N), optional
            Leakage area to each pressure zone
        """
        
        zone = defaultdict(bool)
        zone['ID'] = ID
        zone['XB'] = XB
        if LEAK_AREA != None: zone['LEAK_AREA'] = LEAK_AREA
        self.zones[ID] = zone
        
        
    def calculateMeshCells(self):
        """Returns a list of mesh keys and number of cells in each mesh
        
        Returns
        -------
        list
            List of mesh keys
        list
            List of number of cells
        """
        
        meshes = []
        numCells = []
        meshKeys = list(self.meshes.keys())
        try:
            meshKeys.remove('unknownCounter')
        except:
            pass
        for key in meshKeys:
            IJK = self.meshes[key]['IJK']
            numCells.append(IJK[0]*IJK[1]*IJK[2])
            meshes.append(key)
        return meshes, numCells
    
    
    def checkOverlappingMESH(self):
        """Returns True if any meshes are overlapping else False
        
        Returns
        -------
        bool
            True if any meshes are overlapping, else False
        """
        
        def in_hull(p,hull):
            if not isinstance(hull,scsp.Delaunay):
                hull = scsp.Delaunay(hull)
            return hull.find_simplex(p)>=0
        def pointsFromXB(XB,extend=[0.05, -0.05, 0.05, -0.05, 0, 0]):
            pts = [[XB[0]+extend[0],XB[2]+extend[2],XB[4]+extend[4]],
                   [XB[0]+extend[0],XB[2]+extend[2],XB[5]+extend[5]],
                   [XB[0]+extend[0],XB[3]+extend[3],XB[4]+extend[4]],
                   [XB[0]+extend[0],XB[3]+extend[3],XB[5]+extend[5]],
                   [XB[1]+extend[1],XB[2]+extend[2],XB[4]+extend[4]],
                   [XB[1]+extend[1],XB[2]+extend[2],XB[5]+extend[5]],
                   [XB[1]+extend[1],XB[3]+extend[3],XB[4]+extend[4]],
                   [XB[1]+extend[1],XB[3]+extend[3],XB[5]+extend[5]]]
            return pts
        meshHulls = defaultdict(bool)
        for key in list(self.meshes.keys()):
            pts = pointsFromXB(self.meshes[key]['XB'])
            meshHull = scsp.Delaunay(pts)
            meshHulls[key] = meshHull
        overlap = False
        for key1 in list(self.meshes.keys()):
            for key2 in list(self.meshes.keys()):
                if (key1 != key2):
                    extend = [0.05, -0.05, 0.05, -0.05, 0, 0]
                    if ('east' in key2): extend = [0.05, 0.1, 0.05, -0.05, 0, 0]
                    if ('west' in key2): extend = [-0.1, -0.05, 0.05, -0.05, 0, 0]
                    if ('north' in key2): extend = [0.05, -0.05, 0.05, 0.1, 0, 0]
                    if ('south' in key2): extend = [0.05, -0.05, -0.1, -0.05, 0, 0]
                    pts = pointsFromXB(self.meshes[key2]['XB'], extend=extend)
                    for p in pts:
                        if in_hull(p, meshHulls[key1]):
                            overlap = True
        return overlap
    
    
    def dictFromLine(self, line, lineType, types):
        """Returns a dictionary with keys and values from a namelist
        
        Parameters
        ----------
        line : str
            String namelist line
        lineType : str
            String type of namelist
        types : dict
            Dictionary containing dictionaries of namelist key types
        
        Returns
        -------
        defaultdict
            Dictionary containing keys from the namelist line
        """
        
        lineDict = defaultdict(bool)
        keys = self.splitLineIntoKeys(line)
        for key in keys:
            keyID, keyID2, keyType, keyValue = self.interpretKey(key, lineType, types)
            #print(keyID, keyID2, keyType, keyValue)
            if keyType == 'string':
                keyValue = keyValue.split("'")[1]
            elif keyType == 'float':
                keyValue = float(keyValue.replace(' ', '').replace(',','').replace('/',''))
            elif keyType == 'int':
                keyValue = int(keyValue.replace(' ', '').replace(',','').replace('/',''))
            elif keyType == 'bool':
                keyValue = keyValue.split(".")[1]
            elif ('list' in keyType) and ('ind' not in keyType) and ('row' not in keyType):
                vals = []
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                for t in keyValues:
                    if 'string' in keyType: preprocess = t.split("'")[1]
                    if 'float' in keyType: preprocess = float(t.replace(' ', '').replace(',','').replace('/',''))
                    if 'int' in keyType: preprocess = int(t.replace(' ', '').replace(',','').replace('/',''))
                    vals.append(preprocess)
                keyValue = vals
            elif ('list' in keyType) and ('ind' in keyType) and ('row' not in keyType):
                #print(keyID, keyID2, keyType, keyValue)
                regex1 = r"(\(.{0,3}):(.{0,3}\))"
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                if 'string' in keyType: keyValues = [x.split("'")[1] for x in keyValues]
                if 'float' in keyType: keyValues = [float(x) for x in keyValues]
                if 'int' in keyType: keyValues = [int(x) for x in keyValues]
                tmp = re.search(regex1, keyID)
                if tmp is not None:
                    ar1 = [int(x) for x in tmp.groups()[0].replace('(','').split(':')]
                    ar2 = [int(x) for x in tmp.groups()[1].replace(')','').split(':')]
                else:
                    (ar1, ar2) = ([1], [len(keyValues)])
                tmp = np.zeros((np.max([ar1, ar2]), 1), dtype='object')
                for i in range(0, tmp.shape[0]):
                    tmp[i-1, 0] = keyValues[i-1]
                keyValue = tmp

            elif ('list' in keyType) and ('ind' not in keyType) and ('row' in keyType):
                vals = []
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                for t in keyValues:
                    if 'string' in keyType: preprocess = t
                    if 'float' in keyType: preprocess = float(t.replace(' ', '').replace(',','').replace('/',''))
                    if 'int' in keyType: preprocess = int(t.replace(' ', '').replace(',','').replace('/',''))
                    vals.append(preprocess)
                keyValue = vals
            elif ('matrix' in keyType):
                #print(keyID, keyID2, keyType, keyValue)
                regex1 = r"(\(.{0,3});(.{0,3}\))"
                while (keyValue[-1] == ' ') or (keyValue[-1] == ',') or (keyValue[-1] == '/'):
                    keyValue = keyValue[:-1]
                keyValues = keyValue.split(",")
                if 'string' in keyType: keyValues = [x.split("'")[1] for x in keyValues]
                if 'float' in keyType: keyValues = [float(x) for x in keyValues]
                tmp = re.search(regex1, keyID)
                if tmp is not None:
                    ar1 = [int(x) for x in tmp.groups()[0].replace('(','').split(':')]
                    ar2 = [int(x) for x in tmp.groups()[1].replace(')','').split(':')]
                    if len(ar1) == 1: ar1 = [ar1[0], ar1[0]]
                    if len(ar2) == 1: ar2 = [ar2[0], ar2[0]]
                else:
                    (ar1, ar2) = ([1, 1], [1, len(keyValues)])
                tmp = np.zeros((np.max(ar1), np.max(ar2)), dtype='object')
                counter = 0
                if ar1[0] == ar1[1]:
                    ar1 = np.array(np.zeros((len(keyValues),)) + ar1[0], dtype=np.int32)
                else:
                    ar1 = list(range(ar1[0], ar1[1]+1))
                if ar2[0] == ar2[1]:
                    ar2 = np.array(np.zeros((len(keyValues),)) + ar2[0], dtype=np.int32)
                else:
                    ar2 = list(range(ar2[0], ar2[1]+1))
                for counter in range(0, len(keyValues)):
                    i = ar1[counter]
                    j = ar2[counter]
                    tmp[i-1, j-1] = keyValues[counter]
                    counter += 1
                keyValue = tmp
                
            else:
                print(lineType.lower(), keyID, keyID2, keyType)
                print(len(keyID))
                print(line)
                print(keys)
                assert False, "Stopped"
            lineDict[keyID2] = keyValue
        return lineDict
    
    
    def dictMerge(self, template, master, path=None):
        """Merges two dictionaries
        
        This function merges two dictionaries into a single dictionary.
        The template dictionary is used as the baseline, and master is
        merged into template. Entries in master will overwrite entries
        in template. Note, nested dictionaries will be merged using the
        same procedure.
        
        Parameters
        ----------
        template : dict or defaultdict
            Baseline dictionary
        master : dict or defaultdict
            Master dictionary. Entries in template will be overwritten
            by entries in master.
        path : str
            Internal variable storing path to current key.
            Used in recursive calls for nested dictionaries.
        
        Returns
        -------
        dict or defaultdict
            Merged dictionary
        """
        
        if path is None: path = []
        for key in master:
            if key in template:
                tCheck = isinstance(template[key], dict)
                mCheck = isinstance(master[key], dict)
                if tCheck and mCheck:
                    self.dictMerge(template[key], master[key], path + [str(key)])
                elif template[key] == master[key]:
                    pass
                else:
                    template[key] = master[key]
            else:
                template[key] = master[key]
        return template
    
    
    def generateFDStext(self, newlines=None, fields=None):
        """Returns str of input file
        
        This function generates the fds input file based on the stored
        attribute dictionaries. The optional input parameters provide
        customization in how the input file is exported. Providing
        a value of None will produce the default configuration.
        
        Parameters
        ----------
        newlines : defaultdict, optional
            Dictionary containing boolean for each field type. If True,
            each key from the namelist will be placed on a new line.
            If False, each key will be placed on the same line.
            (default None)
        fields : list, optional
            List containing the order namelists will be exported to
            the input file. (default None)
        
        Returns
        -------
        str
            text of input file
        """
        
        date = datetime.date.today()
        (year, month, day) = (date.year, date.month, date.day)
        dateStr = "%04.0f-%02.0f-%02.0f"%(year, month, day)
        intro = "Input file generated with python-fds-tools v1"
        types = fdsLineTypes(version=self.version)
        if newlines is None: newlines = self.getNewlineFromTypes()
        if fields is None: fields = self.getDefaultFields()
        if self.meshOrder is False: self.addMPIprocesses(1)
        text = "%s\n"%("!"*72)
        text = "%s%s %s on %s%s%s\n"%(
                text, "!"*5, intro, dateStr, " "*2, "!"*5)
        text = "%s%s\n"%(text, "!"*72)
        
        for field in fields:
            key = self.keyFromLineType(field)
            keyN = "&%s"%(field)
            keyT = getattr(types, field.lower())
            keyD = getattr(self, key)
            if key == 'meshes':
                txt = self.makeMESH(keyD, keyT, order=self.meshOrder)
            elif key == 'ramps':
                txt = self.makeRAMP(keyD)
            else:
                newline1 = newlines[key]
                newline2 = keyD['newline']
                newline = (newline1 or newline2)
                txt = self.makeLinesFromDict(keyD, keyT, keyN, newline)
            text = "%s%s"%(text, txt)
        
        for line in self.customLines:
            text = "%s%s\n"%(text, line)
        
        return text
    
    
    def getDefaultFields(self):
        """Returns default field order
        
        Returns
        -------
        list
            List of default field order
        """
        
        fields = ["HEAD", "TIME", "MISC", "WIND", "INIT", "DUMP", "ZONE", 
                  "PRES", "MESH", "REAC", "RADI", "MATL", "SURF",
                  "RAMP", "OBST", "HOLE", "VENT", "PART", "DEVC",
                  "CTRL", "BNDF", "SLCF", "PROP", "SPEC", "PROF"]
        return fields
    
    
    def getLineType(self, line):
        """Returns namelist key from str line
        
        This function extracts the namelist key from a string line
        
        Parameters
        ----------
        line : str
            String containing the fortran namelist line
        
        Returns
        -------
        str
            String containing fortran namelist type
        """
        
        lineType = line[:4]
        return lineType
    
    
    def getMeshLimits(self):
        """Returns a dictionary containing the extents of defined meshes
        
        This function returns a dictionary containing a key 'XB' with an
        array of the total extents defined in meshes.
        
        Returns
        -------
        dict
            Nested dictionary containing 'Overall'->'XB'->float array(6)
        """
        
        meshLimits = defaultdict(bool)
        limitingXB = [100000, -100000, 100000, -100000, 100000, -100000]
        for key in list(self.meshes.keys()):
            mesh = self.meshes[key]
            XB = mesh['XB']
            limitingXB[0] = min([limitingXB[0], XB[0]])
            limitingXB[1] = max([limitingXB[1], XB[1]])
            limitingXB[2] = min([limitingXB[2], XB[2]])
            limitingXB[3] = max([limitingXB[3], XB[3]])
            limitingXB[4] = min([limitingXB[4], XB[4]])
            limitingXB[5] = max([limitingXB[5], XB[5]])
            meshLimits[key] = mesh
        meshLimits['Overall'] = defaultdict(bool)
        meshLimits['Overall']['XB'] = limitingXB
        return meshLimits
    
    
    def getNewlineFromTypes(self):
        """Returns a dictionary containing default new line parameters
        
        Returns
        -------
        dict
            Dictionary containing default new line parameters
        """
        
        newlines = defaultdict(bool)
        newlines['HEAD'] = False
        newlines['TIME'] = False
        newlines['MISC'] = False
        newlines['INIT'] = True
        newlines['DUMP'] = False
        newlines['ZONE'] = True
        newlines['PRES'] = True
        newlines['MESH'] = False
        newlines['REAC'] = True
        newlines['RADI'] = True
        newlines['MATL'] = True
        newlines['SURF'] = True
        newlines['RAMP'] = False
        newlines['OBST'] = False
        newlines['HOLE'] = False
        newlines['VENT'] = False
        newlines['PART'] = False
        newlines['DEVC'] = False
        newlines['CTRL'] = False
        newlines['BNDF'] = False
        newlines['SLCF'] = False
        newlines['PROP'] = False
        newlines['SPEC'] = False
        return newlines
    
    
    def getPolygonNamesFromFdsFile(self):
        """Returns alist of polygons defined in the fds input file
        
        This function returns a list of polygons defined in the fds
        input file.
        
        Returns
        -------
        list
            List containing names of all obstructions which have
            boundary data available.
        """
        
        names = []
        obstList = list(self.obsts.keys())
        if 'unknownCounter' in obstList:
            obstList.remove('unknownCounter')
        for key in obstList:
            if self.obsts[key]['BNDF_OBST']:
                names.append(self.obsts[key]["ID"])
        names = list(set(names))
        return names
    
    
    def importFile(self, file=None, text=None, textList=None):
        """Adds keys to each namelist from an input file, text, or list
        
        This function will add keys to each namelist from an input file,
        text, or text list.
        
        Parameters
        ----------
        file : str, optional
            String containing path to input file
        text : str, optional
            String containing imported text from an input file
        text : str, optional
            List of strings containing individual namelist lines
        """
        
        if file != None:
            f = self.zopen(file)
            textFDS = f.read()
            textFDS = textFDS.decode("utf-8")
        elif text != None:
            textFDS = text
        elif textList != None:
            textFDS = '\n'.join(textList)
        lines = self.makeFDSLines(textFDS)
        self.parseFDSLines(lines)
        
        
    def interpretKey(self, key, lineType, types):
        """Processes a key from a namelist key pair
        
        This function processes a key from a namelist key pair to
        return the keyID, keyType, and keyValue.
        
        Parameters
        ----------
        key : str
            String containing namelist key pair
        lineType : str
            String containing namelist type
        types : defaultdict
            Dictionary containing types for each key in a namelist type
        
        Returns
        -------
        str
            raw keyID containing all text left of = sign
        str
            regex keyID searching for matrix values left of = sign
        dict
            dictionary containing key types for namelist
        str
           raw keyValue containing all text right of = sign
        """

        keyID = key.split('=')[0].upper()
        keyValue = '='.join(key.split('=')[1:])
        regex1 = r"\(\s*.*\)"
        regex2 = r""
        try:
            keyID2 = re.sub(regex1, regex2, keyID)
        except:
            keyID2 = keyID
        #keyID = keyID.strip()
        #keyID2 = keyID.strip()
        keyID2 = keyID2.replace("\t","")
        while keyID2[-1] == ' ':
            keyID2 = keyID2[:-1]
        while keyID2[0] == ' ':
            keyID2 = keyID2[1:]
        keyType = getattr(types, lineType.lower())[keyID2]
        return keyID, keyID2, keyType, keyValue
    
    
    def keyAssist(self, text, types, dic,
                  internalKeys=['counter'], newline=False):
        """Returns a namelist text line from dictionary inputs.
        
        This function returns a namelist text line based on an input
        dictionary and type dictionary.
        
        Parameters
        ----------
        text : str
            String to which to append namelist fields
        types : dict
            Dictionary containing types for namelist fields
        dic : dict
            Dictionary containing namelist keys and values
        internalKeys : list, optional
            List containing internal software fields not to be exported
            to the text line
        newline : bool, optional
            Flag specifying whether each key in the namelist is to be
            entered on the same of different lines
            
        Returns
        -------
        str
            Updated text string
        """
        keys = list(dic.keys())
        keys.sort()
        if 'ID' in keys:
            keys.insert(0, keys.pop(keys.index('ID')))
            if dic['ID'] is False: dic['ID'] = 'UNKNOWN'
        for key in internalKeys:
            if key in keys:
                keys.remove(key)
        for key2 in keys:
            #print(key2)
            if 'THICKNESS' in key2:
                decimals = 8
            else:
                decimals = 4
            if (types[key2] == 'ignore'):
                pass
            elif (types[key2] == 'string'):
                if dic[key2] is not False:
                    text = "%s%s='%s', "%(text, key2, dic[key2])
            elif (types[key2] == 'float'):
                #print(key2, dic[key2])
                if dic[key2] is not False:
                    text = "%s%s=%s, "%(text, key2, '{:.{prec}f}'.format(dic[key2], prec=decimals))
            elif (types[key2] == 'int'):
                if dic[key2] is not False:
                    text = "%s%s=%0.0f, "%(text, key2, dic[key2])
            elif (types[key2] == 'bool'):
                boolCheck = False
                if (dic[key2] is True): boolCheck = True
                if (dic[key2] == 'TRUE'): boolCheck = True
                if (dic[key2] == '.TRUE.'): boolCheck = True
                if boolCheck:
                    text = "%s%s=.TRUE., "%(text, key2)
                else:
                    text = "%s%s=.FALSE., "%(text, key2)
            elif ('listind' in types[key2]):
                temp = np.array(dic[key2])
                tempTxt = "%s(%0.0f:%0.0f)="%(
                        key2, 1, temp.shape[0])
                if type(temp[0]) == np.float64: temp = [temp]
                for t in temp:
                    for tt in t:
                        if ('string' in types[key2]):
                            tempTxt = "%s '%s',"%(tempTxt, tt)
                        if ('float' in types[key2]):
                            tempTxt = "%s %s,"%(tempTxt, '{:.{prec}f}'.format(tt, prec=decimals))
                        if ('int' in types[key2]):
                            tempTxt = "%s %0.0f,"%(tempTxt, tt)
                text = "%s%s "%(text, tempTxt)
            elif ('list' in types[key2]):
                temp = dic[key2]
                tempTxt = "%s="%(key2)
                if temp is not False:
                    for t in temp:
                        if ('string' in types[key2]):
                            tempTxt = "%s '%s',"%(tempTxt, t)
                        if ('float' in types[key2]):
                            tempTxt = "%s %s,"%(tempTxt, '{:.{prec}f}'.format(t, prec=decimals))
                        if ('int' in types[key2]):
                            tempTxt = "%s %0.0f,"%(tempTxt, t)
                    text = "%s%s "%(text, tempTxt)
            elif ('matrix' in types[key2]):
                temp = np.array(dic[key2])
                sz = temp.shape
                if len(sz) == 1:
                    temp = np.reshape(temp, (temp.shape[0], 1))
                    sz = temp.shape
                ar1 = "(%0.0f:%0.0f,%0.0f:%0.0f)"%(
                        1, sz[1], 1, sz[0])
                tempTxt = "%s%s="%(key2, ar1)
                for t in temp.flatten():
                    if ('string' in types[key2]):
                        tempTxt = "%s '%s',"%(tempTxt, t)
                    if ('float' in types[key2]):
                        tempTxt = "%s %s,"%(tempTxt, '{:.{prec}f}'.format(t, prec=decimals))
                    if ('int' in types[key2]):
                        tempTxt = "%s %0.0f,"%(tempTxt, float(t))
                text = "%s%s "%(text, tempTxt)
                
            else:
                print(keys)
                print(dic)
                print(key2)
                print(types[key2])
                assert False, "Stopped"
            if newline and (types[key2] != 'ignore'):
                text = "%s\n      "%(text)
        #except:
        #    print(keys)
        #    print(dic)
        #    print(types[key2])
        return text
    
    
    def keyFromLineType(self, lineType):
        """Returns internal attribute name from namelist type
        
        Parameters
        ----------
        lineType : str
            String containing namelist type
        
        Returns
        -------
        str
            String containing internal attribute name
        """
        
        if lineType == 'HEAD': key = 'head'
        if lineType == 'DEVC': key = 'devcs'
        if lineType == 'INIT': key = 'inits'
        if lineType == 'OBST': key = 'obsts'
        if lineType == 'VENT': key = 'vents'
        if lineType == 'SURF': key = 'surfs'
        if lineType == 'RAMP': key = 'ramps'
        if lineType == 'CTRL': key = 'ctrls'
        if lineType == 'MESH': key = 'meshes'
        if lineType == 'SLCF': key = 'slcfs'
        if lineType == 'BNDF': key = 'bndfs'
        if lineType == 'TIME': key = 'time'
        if lineType == 'DUMP': key = 'dump'
        if lineType == 'MISC': key = 'misc'
        if lineType == 'ZONE': key = 'zones'
        if lineType == 'REAC': key = 'reacs'
        if lineType == 'MATL': key = 'matls'
        if lineType == 'RADI': key = 'radis'
        if lineType == 'PRES': key = 'pres'
        if lineType == 'HOLE': key = 'holes'
        if lineType == 'PART': key = 'parts'
        if lineType == 'PROP': key = 'props'
        if lineType == 'SPEC': key = 'specs'
        if lineType == 'PROF': key = 'profs'
        if lineType == 'WIND': key = 'winds'
        return key
    
    
    def makeFDSLines(self, textFDS):
        """Returns a list of namelist lines
        
        This function cleans the input file, removing line breaks, and
        splitting the text into lines based on namelist grouping.
        
        Parameters
        ----------
        textFDS : str
            String containg text from an fds input file
        
        Returns
        -------
        list
            List of strings containing namelist lines
        """
        
        linesFDS = [x for x in textFDS.split("&")[1:]]
        for i in range(0, len(linesFDS)):
            line2 = linesFDS[i]
            line2 = '/'.join(line2.split('/')[:-1])
            line2 = line2.replace('\r', ',')
            line2 = line2.replace('\n', ',')
            
            line2 = "%s,"%(line2) if line2[-1] != ',' else line2
            line2 = '%s /'%(line2)
            
            while ',,' in line2: line2 = line2.replace(',,',',')
            while ' ,' in line2: line2 = line2.replace(' ,',',')
            while '  ' in line2: line2 = line2.replace("  ", " ")
            while ',,' in line2: line2 = line2.replace(',,',',')
            line_tmp = list(line2)
            if line_tmp[4] == ',':
                line_tmp[4] = ' '
                line2 = "".join(line_tmp)
                while '  ' in line2: line2 = line2.replace("  ", " ")
            linesFDS[i] = line2
        lineTypes = [x[:4] for x in linesFDS]
        if 'TAIL' in lineTypes:
            ind = np.argwhere([True if x == 'TAIL' else False for x in lineTypes])[0][0]
            linesFDS = linesFDS[:ind]
        return linesFDS
    
    
    def makeLinesFromDict(self, items, types, prefix, newline=False):
        """Returns a str generated from a namelist dictionary
        
        This function generates a text string from a namelist
        dictionary.
        
        Parameters
        ----------
        items : dict
            Dictionary containing key pairs from a namelist group
        types : dict
            Dictionary containing types from a namelist group
        prefix : str
            String containing the namelist type
        newline : bool, optional
            Flag specifying whether each key in the namelist is to be
            entered on the same of different lines
        
        Returns
        -------
        str
            Text containing name list line
        """
        
        text = ''
        keys = list(items.keys())
        keys.sort()
        if 'unknownCounter' in keys: keys.remove('unknownCounter')
        if 'newline' in keys: keys.remove('newline')
        for key in keys:
            text = "%s%s "%(text, prefix)
            text = self.keyAssist(text, types, items[key], newline=newline)
            text = "%s /\n"%(text)
        return text
    
    
    def makeMESH(self, meshes, meshTypes, order=False):
        """Returns a str generated from a meshes namelist dictionary.
        
        Parameters
        ----------
        meshes : dict
            Dictionary containing mesh definitions
        meshTypes : dict
            Dictionary containing types from mesh namelists
        order : list, optional
            Order to output mehes. If False, meshes are not output in
            any particular order. (default False)
        
        Returns
        -------
        str
            Text line generated from dictionary
        """
        
        text = ''
        meshList = list(meshes.keys())
        if 'unknownCounter' in meshList:
            meshList.remove('unknownCounter')
        if (order is not False): meshList = [meshList[x] for x in order]
        for key in meshList:
            text = "%s&MESH "%(text)
            text = self.keyAssist(text, meshTypes, meshes[key])
            text = "%s /\n"%(text)
        return text
    
    
    def makeRAMP(self, ramps):
        """Returns a str generated from a ramps namelist dictionary.
        
        Parameters
        ----------
        ramps : dict
            Dictionary containing ramp definitions
        
        Returns
        -------
        str
            Text line generated from dictionary
        """
        
        text = ''
        for key in list(ramps.keys()):
            ID = ramps[key]['ID']
            makeControl = True
            for F, T in zip(ramps[key]['F'], ramps[key]['T']):
                if makeControl and ramps[key]['CTRL_ID']:
                    text = "%s&RAMP ID='%s', T = %0.4f, F = %0.4f, CTRL_ID='%s'/\n"%(text, ID, T, F, ramps[key]['CTRL_ID'])
                    makeControl = False
                else:
                    text = "%s&RAMP ID='%s', T = %0.4f, F = %0.4f, /\n"%(text, ID, T, F)
        return text
    
    
    def mergeTypeFromLineType(self, lineType):
        """Returns internal merge type based on namelist type.
        
        Parameters
        ----------
        lineType : str
            String containing namelist type
        
        Returns
        -------
        str
            String containing merge type for namelist type
        """
        
        key = 'unknown'
        if lineType == 'HEAD': key = 'merge'
        if lineType == 'DEVC': key = 'enumerate'
        if lineType == 'INIT': key = 'enumerate'
        if lineType == 'OBST': key = 'enumerate'
        if lineType == 'VENT': key = 'enumerate'
        if lineType == 'SURF': key = 'enumerate'
        if lineType == 'RAMP': key = 'append'
        if lineType == 'CTRL': key = 'enumerate'
        if lineType == 'MESH': key = 'enumerate'
        if lineType == 'SLCF': key = 'enumerate'
        if lineType == 'BNDF': key = 'enumerate'
        if lineType == 'TIME': key = 'merge'
        if lineType == 'DUMP': key = 'merge'
        if lineType == 'MISC': key = 'merge'
        if lineType == 'ZONE': key = 'enumerate'
        if lineType == 'REAC': key = 'enumerate'
        if lineType == 'MATL': key = 'enumerate'
        if lineType == 'RADI': key = 'merge'
        if lineType == 'PRES': key = 'merge'
        if lineType == 'HOLE': key = 'enumerate'
        if lineType == 'PART': key = 'enumerate'
        if lineType == 'PROP': key = 'enumerate'
        if lineType == 'SPEC': key = 'enumerate'
        if lineType == 'PROF': key = 'enumerate'
        if lineType == 'WIND': key = 'merge'
        return key
    
    
    def parseFDSLines(self, lines):
        """Adds each line to internal attribute namelist dictionaries.
        
        Parameters
        ----------
        lines : list
            List containing strings of namelist lines
        """
        
        for line in lines:
            lineType = self.getLineType(line)
            key = self.keyFromLineType(lineType)
            types = fdsLineTypes(version=self.version)
            self.parseLine(line, lineType, types, key)
        devcKeys = list(self.devcs.keys())
        devcKeys.remove('unknownCounter')
        for key in devcKeys:
            if self.devcs[key]['INIT_ID']:
                initXYZ = self.inits[self.devcs[key]['INIT_ID']]['XYZ']
                self.devcs[key]['XYZ'] = initXYZ
            else:
                self.devcs[key].pop('INIT_ID')
        
        
    def parseLine(self, line, lineType, types, key):
        """Adds one line to the internal attribute namelist dictionary.
        
        Parameters
        ----------
        line : str
            String containing namelist line
        lineType : str
            String containing namelist line type
        types : dict
            Dictionary containing key types for namelist pair
        key : str
            String containing internal attribute key for namelist line
            type
        """
        #print(line)
        check = True
        try:
            lineDict = self.dictFromLine(line, lineType, types)
        except:
            print("WARNING: Unknown line in input file.\n")
            print("%s\n"%(line))
            check = False
        if check:
            tmp = getattr(self, key)
            mergeType = self.mergeTypeFromLineType(lineType)
            if mergeType == 'merge':
                if not tmp['ID']: tmp['ID'] = defaultdict(bool)
                tmp['ID'] = self.dictMerge(tmp['ID'], lineDict)
                setattr(self, key, tmp)
            elif mergeType == 'append':
                ID = lineDict['ID']
                if tmp[ID]:
                    for keyID2 in list(lineDict.keys()):
                        keyType = getattr(types, lineType.lower())[keyID2]
                        keyValue = lineDict[keyID2]
                        if (keyType == 'listrowfloat'):
                            for v in keyValue:
                                tmp[ID][keyID2].append(v)
                else:
                    tmp[ID] = lineDict
            elif mergeType == 'enumerate':
                ID = lineDict['ID']
                if ID is False:
                    ID = "ID"
                    lineDict["ID"] = ID
                if tmp[ID]:
                    counter = tmp[ID]['counter']
                    if lineDict['ID'] == False: lineDict["ID"] = "%s-%04.0f"%(ID, counter)
                    tmp["%s-%04.0f"%(ID, counter)] = lineDict
                    tmp[ID]['counter'] += 1
                    pass
                else:
                    tmp[ID] = lineDict
                    tmp[ID]['counter'] = 0
            else:
                assert False, "Stopped"
        
        
    def saveModel(self, mpiProcesses, location,
                  allowMeshSplitting=True, splitMultiplier=1.2,
                  meshSplitAxes=[True, True, False]):
        """Saves an fds input file
        
        Input file is generated based on internal attribute namelist
        dictionaries. This functiona also allows splitting of meshes to
        optimize mpi processes balance.
        
        Parameters
        ----------
        mpiProcesses : int
            The number of mpi processes to define in the input file
        location : str
            The path location to save the input file
        allowMeshSplitting : bool, optional
            Flag to enable mesh splitting for balancing mpi processes
            (default is True)
        splitMultiplier : float, optional
            Tolerance used in mesh splitting (default is 1.2)
        meshSplitAxes : list of booleans, optional
            Specifies along which axes the software is allowed to split
            meshes
        """
        self.addMPIprocesses(
                mpiProcesses, allowMeshSplitting=allowMeshSplitting, 
                splitMultiplier=splitMultiplier,
                meshSplitAxes=meshSplitAxes)
        text = self.generateFDStext()
        with open(location, 'w') as f:
            f.write(text)
        print("Input file written to: %s"%(location))
        
        
    def splitLineIntoKeys(self, line2):
        """Returns namelist key pairs from a line.
        
        Parameters
        ----------
        line2 : str
            String containing namelist line
        
        Returns
        -------
        list
            List containing namelist keys
        """
        
        line = line2.replace('\n', ',').replace('\r', ',')
        while (',,' in line) or ('  ' in line):
            line = line.replace(',,', ',').replace('  ', ' ')    
            
        regex1 = r"(\(.{0,3}),(.{0,3}\))"
        regex2 = r"\1;\2"
        try:
            line = re.sub(regex1, regex2, line)
        except:
            pass
        keys = line.split(',')
        keys[0] = keys[0][4:]
        updatedKeys = []
        txt = ''
        for i in range(0,len(keys)):
            if '=' in keys[i]:
                updatedKeys.append(txt)
                txt = keys[i]
            else:
                txt = ','.join([txt,keys[i]])
        updatedKeys.append(txt)
        while '' in updatedKeys:
            updatedKeys.remove('')
        for i, txt in enumerate(updatedKeys):
            while txt[0] == ' ':
                txt = txt[1:]
            updatedKeys[i] = txt
        for i, txt in enumerate(updatedKeys):
            while txt[-1] == ' ' or txt[-1] == ',' or txt[-1] == '/':
                txt = txt[:-1]
            updatedKeys[i] = txt
        return updatedKeys
    
    
    def splitMESHonce(self, mesh, meshSplitAxes):
        """Splits a mesh along its largest axis.
        
        Parameters
        ----------
        mesh : dict
            Dictionary containing information for a single mesh
        meshSplitAxes : list of booleans
            Specifies along which axes the software is allowed to split
            the mesh.
        """
        
        IJK = np.round(mesh['IJK'])
        XB = mesh['XB']
        dxs = [(XB[1]-XB[0])/float(IJK[0]), (XB[3]-XB[2])/float(IJK[1]), (XB[5]-XB[4])/float(IJK[2])]
        ind = np.argmax(IJK)
        IJK_temp = list(IJK)
        while meshSplitAxes[ind] is False:
            IJK_temp = list(IJK_temp)
            IJK_temp[ind] = 0
            ind = np.argmax(IJK_temp)
            if np.sum(IJK_temp) == 0:
                print("Failed to split mesh.")
                break
        
        IJK2 = list(IJK)
        XB2 = list(XB)
        IJK2[ind] = int(IJK[ind]/2)
        if IJK2[ind] % 2 > 0: IJK2[ind] = IJK2[ind]-1
        XB2[int(2*ind+1)] = XB2[int(2*ind)] + dxs[ind]*float(IJK2[ind])
        
        IJK3 = list(IJK)
        XB3 = list(XB)
        IJK3[ind] = IJK[ind] - IJK2[ind]
        XB3[int(2*ind)] = XB2[int(2*ind+1)]
        
        mesh2 = defaultdict(bool)
        mesh2['ID'] = "%s-00"%(mesh["ID"])
        mesh2['XB'] = XB2
        mesh2['IJK'] = IJK2
        
        mesh3 = defaultdict(bool)
        mesh3['ID'] = "%s-01"%(mesh["ID"])
        mesh3['XB'] = XB3
        mesh3['IJK'] = IJK3
        
        self.meshes.pop(mesh['ID'], False)
        self.meshes[mesh2['ID']] = mesh2
        self.meshes[mesh3['ID']] = mesh3
        
        
    def zopen(self, file):
        """Opens a file or zip archive for reading.
        
        Parameters
        ----------
        file : str
            String containing path to file or zip archive
        
        Returns
        -------
        file
            Open binary file for reading
        """
        
        if '.zip' in file:
            zname = '%s.zip'%(file.split('.zip')[0])
            fname = file.split('.zip%s'%(os.sep))[1]
            zip = zipfile.ZipFile(zname, 'r')
            f = zip.open(fname)
        else:
            f = open(file, 'rb')
        return f
    
    
