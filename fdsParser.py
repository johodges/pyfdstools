# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 08:17:40 2019

@author: JHodges
"""

import numpy as np
from collections import defaultdict
import datetime

def dictMerge(a, b, path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                dictMerge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass
            else:
                a[key] = b[key]
        else:
            a[key] = b[key]
    return a

class fdsFileOperations(object):
    def __init__(self):
        self.head = defaultdict(bool)
        self.devcs = defaultdict(bool)
        self.inits = defaultdict(bool)
        self.obsts = defaultdict(bool)
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
        
        self.devcUnknownCounter = 0
        self.obstUnknownCounter = 0
        self.ventUnknownCounter = 0
        self.meshUnknownCounter = 0
        
        self.slcfCounter = 0
        self.bndfCounter = 0
        
    def importFile(self, file=None, text=None, textList=None):
        if file != None:
            with open(file,'r') as f:
                textFDS = f.read()
        elif text != None:
            textFDS = text
        elif textList != None:
            textFDS = '\n'.join(textList)
        lines = self.makeFDSLines(textFDS)
        self.parseFDSLines(lines)
    
    def makeFDSLines(self, textFDS):
        linesFDS = [x.split('/')[0] for x in textFDS.split("&")[1:]]
        for i in range(0, len(linesFDS)):
            line2 = linesFDS[i]
            line2 = "%s,"%(line2) if line2[-1] != ',' else line2
            linesFDS[i] = line2
        return linesFDS
    
    def addOBST(self, ID, XB, SURF_IDS=None, SURF_ID=None, SURF_ID6=None, BNDF_OBST=True):
        obst = defaultdict(bool)
        obst['XB'] = XB
        obst['ID'] = ID
        if SURF_IDS != None: obst['SURF_IDS'] = SURF_IDS
        if SURF_ID != None: obst['SURF_ID'] = SURF_ID
        if SURF_ID6 != None: obst['SURF_ID6'] = SURF_ID6
        if BNDF_OBST: obst['BNDF_OBST'] = True
        self.obsts[ID] = obst
    
    def addHEAD(self, chid, title=None):
        head = defaultdict(bool)
        head['chid'] = chid
        if title != None:
            head['title'] = title
        else:
            head['title'] = title
        self.head = head
    
    def addTIME(self, T_END=0.0, T_BEGIN=0.0):
        time = defaultdict(bool)
        time['T_BEGIN'] = T_BEGIN
        time['T_END'] = T_END
        self.time = time
        
    def addDUMP(self, RENDER_FILE=None, COLUMN_DUMP_LIMIT=False, WRITE_XYZ=False,
                DT_PL3D=None, DT_SL3D=None, DT_SLCF=None, DT_BNDF=None, DT_DEVC=None, DT_CTRL=None):
        dump = defaultdict(bool)
        if RENDER_FILE != None: dump['RENDER_FILE'] = RENDER_FILE
        if COLUMN_DUMP_LIMIT: dump['COLUMN_DUMP_LIMIT'] = COLUMN_DUMP_LIMIT
        if WRITE_XYZ: dump['WRITE_XYZ'] = WRITE_XYZ
        if DT_PL3D != None: dump['DT_PL3D'] = DT_PL3D
        if DT_SL3D != None: dump['DT_SL3D'] = DT_SL3D
        if DT_SLCF != None: dump['DT_SLCF'] = DT_SLCF
        if DT_BNDF != None: dump['DT_BNDF'] = DT_BNDF
        if DT_DEVC != None: dump['DT_DEVC'] = DT_DEVC
        if DT_CTRL != None: dump['DT_CTRL'] = DT_CTRL
        self.dump = dump
    
    def addMESH(self, ID, IJK, XB):
        mesh = defaultdict(bool)
        mesh['ID'] = ID
        mesh['IJK'] = IJK
        mesh['XB'] = XB
        self.meshes[ID] = mesh
        
    def addREAC(self, ID, FUEL=None, FORMULA=None, AIT=None, SY=None, COY=None, HOC=None):
        reac = defaultdict(bool)
        reac['ID'] = ID
        if FUEL != None: reac['FUEL'] = FUEL
        if FORMULA != None: reac['FORMULA'] = FORMULA
        if AIT != None: reac['AUTO_IGNITION_TEMPERATURE'] = AIT
        if SY != None: reac['SOOT_YIELD'] = SY
        if COY != None: reac['CO_YIELD'] = COY
        if HOC != None: reac['HEAT_OF_COMBUSTION'] = HOC
        self.reacs[ID] = reac       
        
    def addMATL(self, ID, Emi=None, Den=None, Con=None, Spe=None, kramp=None, cpramp=None):
        matl = defaultdict(bool)
        matl['ID'] = ID
        if Emi != None: matl['EMISSIVITY'] = Emi
        if Den != None: matl['DENSITY'] = Den
        if Con != None: matl['CONDUCTIVITY'] = Con
        if Spe != None: matl['SPECIFIC_HEAT'] = Spe
        if kramp != None: matl['CONDUCTIVITY_RAMP'] = kramp
        if cpramp != None: matl['SPECIFIC_HEAT_RAMP'] = cpramp
        self.matls[ID] = matl
    
    def addRAMP(self, ID, T, F):
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
    
    def addSURF(self, ID, Mid=None, Col=None, Thi=None, Bac=None, Geo=None,
                Fyi=None, Len=None, LeaPat=None, Hrrpua=None, qramp=None,
                Rgb=None, adiabatic=False):
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
        self.surfs[ID] = surf
        
    def addVENT(self, ID, SURF_ID, XB=None, CTRL_ID=None, MB=None):
        vent = defaultdict(bool)
        vent['ID'] = ID
        vent['SURF_ID'] = SURF_ID
        if XB != None: vent['XB'] = XB
        if CTRL_ID != None: vent['CTRL_ID'] = CTRL_ID
        if MB != None: vent['MB'] = MB
        self.vents[ID] = vent
    
    def addSLCF(self, Qty, PBX=None, PBY=None, PBZ=None, Vec=False, XB=None):
        slcf = defaultdict(bool)
        slcf['ID'] = "SLCF-%05.0f"%(self.slcfCounter)
        slcf['QUANTITY'] = Qty
        if 'PBX' != None: slcf['PBX'] = PBX
        if 'PBY' != None: slcf['PBY'] = PBY
        if 'PBZ' != None: slcf['PBZ'] = PBZ
        if Vec: slcf['VECTOR'] = 'TRUE'
        if XB != None: slcf['XB'] = XB
        self.slcfCounter += 1
        self.slcfs[slcf['ID']] = slcf
    
    def getLineType(self, line):
        lineType = line[:4]
        return lineType
    
    def parseHEAD(self, line):
        head = parseHEAD(line)
        self.head = dictMerge(self.head, head)
    
    def parseDEVC(self, line):
        devc = parseDEVC(line)
        if self.devcs[devc['ID']]:
            self.devcs['NULL%04.0f'%(self.devcUnknownCounter)] = devc
            self.devcUnknownCounter += 1
        else:
            self.devcs[devc['ID']] = devc
    
    def parseINIT(self, line):
        init = parseINIT(line)
        self.init[init['ID']] = init
    
    def parseOBST(self, line):
        obst = parseOBST(line)
        if not obst['ID']:
            self.obsts['NULL%04.0f'%(self.obstUnknownCounter)] = obst
            self.obstUnknownCounter += 1
        else:
            if not self.obsts[obst['ID']]:
                obst['number'] = 0
                self.obsts[obst['ID']] = obst
            else:
                obstNum = self.obsts[obst['ID']]['number']+1
                self.obsts["%s-%04.0f"%(obst['ID'],obstNum)] = obst
                self.obsts[obst['ID']]['number'] = obstNum
                
    def parseVENT(self, line):
        vent = parseVENT(line)
        if not vent['ID']:
            self.vents['NULL%04.0f'%(self.ventUnknownCounter)] = vent
            self.ventUnknownCounter += 1
        else:
            if not self.vents[vent['ID']]:
                vent['number'] = 0
                self.vents[vent['ID']] = vent
            else:
                ventNum = self.vents[vent['ID']]['number']+1
                self.vents["%s-%04.0f"%(vent['ID'],ventNum)] = vent
                self.vents[vent['ID']]['number'] = ventNum

    def parseSURF(self, line):
        surf = parseSURF(line)
        self.surfs[surf['ID']] = surf

    def parseRAMP(self, line):
        ramp = parseRAMP(line)
        if self.ramps[ramp['ID']]:
            Ts = self.ramps[ramp['ID']]['T']
            Ts.append(ramp['T'])
            Fs = self.ramps[ramp['ID']]['F']
            Fs.append(ramp['F'])
            self.ramps[ramp['ID']]['T'] = Ts
            self.ramps[ramp['ID']]['F'] = Fs
        else:
            ramp['T'] = [ramp['T']]
            ramp['F'] = [ramp['F']]
            self.ramps[ramp['ID']] = ramp

    def parseCTRL(self, line):
        ctrl = parseCTRL(line)
        self.ctrls[ctrl["ID"]] = ctrl
    
    def parseMESH(self, line):
        mesh = parseMESH(line)
        if self.meshes[mesh["ID"]]:
            self.meshes["%s-%04.0f"%(mesh["ID"], self.meshUnknownCounter)] = mesh
            self.meshUnknownCounter += 1
        else:
            self.meshes[mesh["ID"]] = mesh
    
    def parseSLCF(self, line):
        slcf = parseSLCF(line)
        slcf['ID'] = "SLCF-%05.0f"%(self.slcfCounter)
        self.slcfs[slcf['ID']] = slcf
        self.slcfCounter += 1
    
    def parseBNDF(self, line):
        bndf = parseBNDF(line)
        bndf['ID'] = "BNDF-%05.0f"%(self.bndfCounter)
        self.bndfs[bndf['ID']] = bndf
        self.bndfCounter += 1
    
    def parseTIME(self, line):
        time = parseTIME(line)
        self.time = dictMerge(self.time, time)
        
    def parseDUMP(self, line):
        dump = parseDUMP(line)
        self.dump = dictMerge(self.dump, dump)
        
    def parseMISC(self, line):
        misc = parseMISC(line)
        self.misc = dictMerge(self.misc, misc)
    
    def parseZONE(self, line):
        zone = parseZONE(line)
        self.zones[zone['ID']] = zone
    
    def parseREAC(self, line):
        reac = parseREAC(line)
        self.reacs[reac['ID']] = reac
    
    def parseMATL(self, line):
        matl = parseMATL(line)
        self.matls[matl['ID']] = matl
    
    def parseFDSLines(self, lines):
        for line in lines:
            lineType = self.getLineType(line)
            if lineType == 'HEAD': self.parseHEAD(line)
            if lineType == 'DEVC': self.parseDEVC(line)
            if lineType == 'INIT': self.parseINIT(line)
            if lineType == 'OBST': self.parseOBST(line)
            if lineType == 'VENT': self.parseVENT(line)
            if lineType == 'SURF': self.parseSURF(line)
            if lineType == 'RAMP': self.parseRAMP(line)
            if lineType == 'CTRL': self.parseCTRL(line)
            if lineType == 'MESH': self.parseMESH(line)
            if lineType == 'SLCF': self.parseSLCF(line)
            if lineType == 'BNDF': self.parseBNDF(line)
            if lineType == 'TIME': self.parseTIME(line)
            if lineType == 'DUMP': self.parseDUMP(line)
            if lineType == 'MISC': self.parseMISC(line)
            if lineType == 'ZONE': self.parseZONE(line)
            if lineType == 'REAC': self.parseREAC(line)
            if lineType == 'MATL': self.parseMATL(line)

        for key in list(self.devcs.keys()):
            if self.devcs[key]['INIT_ID']:
                self.devcs[key]['XYZ'] = self.inits[self.devcs[key]['INIT_ID']]['XYZ']
    
    def generateFDStext(self):
        date = datetime.date.today()
        text = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        text = "%s!!!!! Input file generated with fdsTools v1 on %04.0f-%02.0f-%02.0f          !!!!!\n"%(text, date.year, date.month, date.day)
        text = "%s!!!!! Copyright JENSEN HUGHES %04.0f All right reserved.              !!!!!\n"%(text, date.year)
        text = "%s!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"%(text)
        
        headTxt = makeHEAD(self.head)
        text = "%s%s"%(text, headTxt)
        text = "%s%s"%(text, makeTIME(self.time))
        text = "%s%s"%(text, makeDUMP(self.dump))
        text = "%s%s"%(text, makeZONE(self.zones))

        meshTxt = makeMESH(self.meshes)
        text = "%s%s"%(text, meshTxt)
        text = "%s%s"%(text, makeREAC(self.reacs))
        text = "%s%s"%(text, makeMATL(self.matls))
        
        surfTxt = makeSURF(self.surfs)
        text = "%s%s"%(text, surfTxt)
        
        obstTxt = makeOBST(self.obsts)
        text = "%s%s"%(text, obstTxt)
        
        ventTxt = makeVENT(self.vents)
        text = "%s%s"%(text, ventTxt)
        
        devcTxt = makeDEVC(self.devcs)
        text = "%s%s"%(text, devcTxt)
        
        ctrlTxt = makeCTRL(self.ctrls)
        text = "%s%s"%(text, ctrlTxt)
        
        initTxt = makeINIT(self.inits)
        text = "%s%s"%(text, initTxt)
        
        rampTxt = makeRAMP(self.ramps)
        text = "%s%s"%(text, rampTxt)
        
        bndfTxt = makeBNDF(self.bndfs)
        text = "%s%s"%(text, bndfTxt)
        
        slcfTxt = makeSLCF(self.slcfs)
        text = "%s%s"%(text, slcfTxt)
        
        return text

def makeMATL(matls):
    text = ''
    for key in list(matls.keys()):
        ID = matls[key]['ID']
        FYI = ''
        SPECIFIC_HEAT = ''
        CONDUCTIVITY = ''
        DENSITY = ''
        EMISSIVITY = ''
        if matls[key]['FYI']: FYI = "FYI = '%s', "%(matls[key]['FYI'])
        if matls[key]['SPECIFIC_HEAT']: SPECIFIC_HEAT = "SPECIFIC_HEAT = %0.4f, "%(matls[key]['SPECIFIC_HEAT'])
        if matls[key]['CONDUCTIVITY']: CONDUCTIVITY = "CONDUCTIVITY = %0.4f, "%(matls[key]['CONDUCTIVITY'])
        if matls[key]['DENSITY']: DENSITY = "DENSITY = %0.4f, "%(matls[key]['DENSITY'])
        if matls[key]['EMISSIVITY']: EMISSIVITY = "EMISSIVITY = %0.4f, "%(matls[key]['EMISSIVITY'])
        
        text = "%s&MATL ID = '%s', %s%s%s%s%s /\n"%(text, ID, FYI, SPECIFIC_HEAT, CONDUCTIVITY, DENSITY, EMISSIVITY)
    return text

def makeREAC(reacs):
    text = ''
    for key in list(reacs.keys()):
        ID = reacs[key]['ID']
        FUEL = reacs[key]['FUEL']
        FORMULA = reacs[key]['FORMULA']
        CO_YIELD = ''
        SOOT_YIELD = ''
        IDEAL = ''
        if reacs[key]['CO_YIELD']: CO_YIELD = "CO_YIELD = %0.8f, "%(reacs[key]['CO_YIELD'])
        if reacs[key]['SOOT_YIELD']: SOOT_YIELD = "SOOT_YIELD = %0.8f, "%(reacs[key]['SOOT_YIELD'])
        if reacs[key]['IDEAL']:
            if reacs[key]['IDEAL'] == 'TRUE':
                IDEAL = "IDEAL = .TRUE., "
            elif reacs[key]['IDEAL'] == 'FALSE':
                IDEAL = "IDEAL = .FALSE., "
        text = "%s&REAC ID='%s', FUEL='%s', FORMULA='%s', %s%s%s /\n"%(text, ID, FUEL, FORMULA, CO_YIELD, SOOT_YIELD, IDEAL)
    return text

def makeTIME(time):
    text = ''
    T_BEGIN = ''
    if time['T_BEGIN']: T_BEGIN = "T_BEGIN = %0.4f, "%(time['T_BEGIN'])
    T_END = ''
    if time['T_END']: T_END = "T_END = %0.4f, "%(time['T_END'])
    text = "%s&TIME %s%s /\n"%(text, T_BEGIN, T_END)
    return text

def makeDUMP(dump):
    text = ''
    RENDER_FILE = ''
    COLUMN_DUMP_LIMIT = ''
    DT_CTRL = ''
    DT_DEVC = ''
    DT_BNDF = ''
    DT_SLCF = ''
    DT_SL3D = ''
    DT_RESTART = ''
    if dump['RENDER_FILE']: RENDER_FILE = "RENDER_FILE = '%s', "%(dump['RENDER_FILE'])
    if dump['COLUMN_DUMP_LIMIT']: COLUMN_DUMP_LIMIT = "COLUMN_DUMP_LIMIT = .%s., "%(dump['COLUMN_DUMP_LIMIT'])
    if dump['DT_CTRL']: DT_CTRL = "DT_CTRL = %s, "%(dump['DT_CTRL'])
    if dump['DT_DEVC']: DT_DEVC = "DT_DEVC = %s, "%(dump['DT_DEVC'])
    if dump['DT_BNDF']: DT_BNDF = "DT_BNDF = %s, "%(dump['DT_BNDF'])
    if dump['DT_SLCF']: DT_SLCF = "DT_SLCF = %s, "%(dump['DT_SLCF'])
    if dump['DT_SL3D']: DT_SL3D = "DT_SL3D = %s, "%(dump['DT_SL3D'])
    if dump['DT_RESTART']: DT_RESTART = "DT_RESTART = %s, "%(dump['DT_RESTART'])
    
    text = "%s&DUMP %s%s%s%s%s%s%s%s/\n"%(text, RENDER_FILE, COLUMN_DUMP_LIMIT, DT_CTRL, DT_DEVC, DT_BNDF, DT_SLCF, DT_SL3D, DT_RESTART)
    if text == "&DUMP /\n":
        return ''
    else:
        return text

def makeMISC(misc):
    text = ''
    TMPA = ''
    HUMIDITY = ''
    SUPPRESSION = ''
    BNDF_DEFAULT = ''
    if misc['TMPA']: TMPA = "TMPA = %s, "%(misc['TMPA'])
    if misc['HUMIDITY']: HUMIDITY = "HUMIDITY = %s, "%(misc['HUMIDITY'])
    if misc['SUPPRESSION']:
        if misc['SUPPRESSION'] == 'TRUE':
            SUPPRESSION = "SUPPRESSION = .TRUE., "
        elif misc['SUPPRESSION'] == 'FALSE':
            SUPPRESSION = "SUPPRESSION = .FALSE., "
    if misc['BNDF_DEFAULT']:
        if misc['BNDF_DEFAULT'] == 'TRUE':
            BNDF_DEFAULT = "BNDF_DEFAULT = .TRUE., "
        elif misc['BNDF_DEFAULT'] == 'FALSE':
            BNDF_DEFAULT = "BNDF_DEFAULT = .FALSE., "
    text = "%s&MISC %s%s%s%s /\n"%(text, TMPA, HUMIDITY, SUPPRESSION, BNDF_DEFAULT)
    if text == "&MISC /\n":
        return ''
    else:
        return text

def makeZONE(zones):
    text = ''
    for key in list(zones.keys()):
        ID = zones[key]['ID']
        XB = zones[key]['XB']
        loc ="XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        text = "%s&ZONE ID='%s', %s /\n"%(text, ID, loc)
    return text

def makeBNDF(bndfs):
    text = ''
    for key in list(bndfs.keys()):
        QUANTITY = bndfs[key]['QUANTITY']
        text = "%s&BNDF QUANTITY = '%s' /\n"%(text, QUANTITY)
    return text

def makeSLCF(slcfs):
    text = ''
    for key in list(slcfs.keys()):
        QUANTITY = slcfs[key]['QUANTITY']
        if slcfs[key]['PBX'] is not False: loc = "PBX = %0.4f, "%(slcfs[key]['PBX'])
        if slcfs[key]['PBY'] is not False: loc = "PBY = %0.4f, "%(slcfs[key]['PBY'])
        if slcfs[key]['PBZ'] is not False: loc = "PBZ = %0.4f, "%(slcfs[key]['PBZ'])
        if slcfs[key]['XB'] is not False:
            XB = slcfs[key]['XB']
            loc ="XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        vec = ''
        if slcfs[key]['VECTOR']:
            if slcfs[key]['VECTOR'] == 'TRUE':
                vec = "VECTOR = .TRUE., "
        text = "%s&SLCF QUANTITY = '%s', %s%s /\n"%(text, QUANTITY, loc, vec)
    return text

def makeMESH(meshes):
    text = ''
    for key in list(meshes.keys()):
        ID = meshes[key]['ID']
        IJK = meshes[key]['IJK']
        IJKstr = "%0.0f,%0.0f,%0.0f"%(IJK[0], IJK[1], IJK[2])
        XB = meshes[key]['XB']
        loc ="XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        text = "%s&MESH ID='%s', IJK=%s, %s /\n"%(text, ID, IJKstr, loc)
    return text
        
def makeCTRL(ctrls):
    text = ''
    for key in list(ctrls.keys()):
        ID = ctrls[key]['ID']
        FUNCTION_TYPE = ctrls[key]['FUNCTION_TYPE']
        INPUT_ID = ctrls[key]['INPUT_ID']
        inputIdTxt = 'INPUT_ID ='
        for inp in INPUT_ID:
            "%s '%s',"%(inputIdTxt, inp)
        text = "%s&CTRL ID='%s', FUNCTION_TYPE='%s', %s /\n"%(text, ID, FUNCTION_TYPE, inputIdTxt)
    return text

def makeRAMP(ramps):
    text = ''
    for key in list(ramps.keys()):
        ID = ramps[key]['ID']
        for F, T in zip(ramps[key]['F'], ramps[key]['T']):
            text = "%s&RAMP ID='%s', T = %0.4f, F = %0.4f, /\n"%(text, ID, T, F)
    return text

def makeSURF(surfs):
    text = ''
    for key in list(surfs.keys()):
        ID = surfs[key]['ID']
        MATL_ID = surfs[key]['MATL_ID']
        if surfs[key]['MATL_ID']:
            matlTxt = 'MATL_ID ='
            for matl in MATL_ID:
                matlTxt = "%s '%s',"%(matlTxt, matl)
        else:
            matlTxt = ''
        if surfs[key]['THICKNESS']:
            THICKNESS = surfs[key]['THICKNESS']
            thicknessTxt = 'THICKNESS ='
            for t in THICKNESS:
                thicknessTxt = "%s %0.4f,"%(thicknessTxt, t)
        else:
            thicknessTxt= ''
        COLOR = ''
        if surfs[key]['COLOR']:
            COLOR = "COLOR = '%s', "%(surfs[key]['COLOR'])
        if surfs[key]['RGB']:
            rgb = surfs[key]['RGB']
            COLOR = "RGB = %0.1f, %0.1f, %0.1f, "%(rgb[0], rgb[1], rgb[2])
        BACKING = ''
        if surfs[key]['BACKING']: BACKING = "BACKING = '%s', "%(surfs[key]['BACKING'])
        if surfs[key]['GEOMETRY']:
            GEOMETRY = "GEOMETRY = '%s', "%(surfs[key]['GEOMETRY'])
        else:
            GEOMETRY = ''
        if surfs[key]['FYI']:
            FYI = "FYI = %s, "%(surfs[key]['FYI'])
        else:
            FYI = ''
        if surfs[key]['LENGTH']:
            LENGTH = "LENGTH = %s, "%(surfs[key]['LENGTH'])
        else:
            LENGTH = ''
        if surfs[key]['LEAK_PATH']:
            leakTxt = 'LEAK_PATH ='
            for l in surfs[key]['LEAK_PATH']:
                leakTxt = "%s %0.0f,"%(leakTxt, l)
        else:
            leakTxt = ''
        if surfs[key]['HRRPUA']:
            HRRPUA = "HRRPUA = %0.4f, "%(surfs[key]['HRRPUA'])
        else:
            HRRPUA = ''
        if surfs[key]['RAMP_Q']:
            RAMP_Q = "RAMP_Q = '%s', "%(surfs[key]['RAMP_Q'])
        else:
            RAMP_Q = ''
        if surfs[key]['ADIABATIC']:
            ADB = "ADIABATIC = .TRUE.,"
        else:
            ADB = ""
        text = "%s&SURF ID='%s', %s%s%s %s%s%s%s%s%s%s%s /\n"%(text, ID, COLOR, matlTxt, thicknessTxt, BACKING, GEOMETRY, LENGTH, FYI, HRRPUA, RAMP_Q, leakTxt, ADB)
    return text

def makeVENT(vents):
    text = ''
    for key in list(vents.keys()):
        SURF_ID = vents[key]['SURF_ID']
        ID = vents[key]['ID']
        XB = vents[key]['XB']
        MB = vents[key]['MB']
        if np.any(vents[key]['XB']):
            loc = "XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        if vents[key]['MB']:
            loc ="MB = '%s', "%(MB)
        if vents['CTRL_ID']: CTRL_ID = "CTRL_ID = %s, "%(vents[key]['CTRL_ID'])
        if not vents[key]['CTRL_ID']: CTRL_ID = ""
        text = "%s&VENT ID='%s', SURF_ID='%s', %s, %s/\n"%(text, ID, SURF_ID, loc, CTRL_ID)
    return text

def makeOBST(obsts):
    text = ''
    for key in list(obsts.keys()):
        SURF_IDtxt = ''
        if obsts[key]['SURF_ID']: SURF_IDtxt = "SURF_ID='%s',"%(obsts[key]['SURF_ID'])
        if obsts[key]['SURF_IDS']: SURF_IDtxt = "SURF_IDS=%s,"%(obsts[key]['SURF_IDS'])
        if obsts[key]['SURF_ID6']: SURF_IDtxt = "SURF_ID6='%s',"%(obsts[key]['SURF_ID6'])
        ID = obsts[key]['ID']
        XB = obsts[key]['XB']
        loc = "XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        if obsts[key]['BNDF_OBST']: BNDF_OBST = "BNDF_OBST = .TRUE., "
        if not obsts[key]['BNDF_OBST']: BNDF_OBST = "BNDF_OBST = .FALSE., "
        if obsts[key]['CTRL_ID']: CTRL_ID = "CTRL_ID = %s, "%(obsts['CTRL_ID'])
        if not obsts[key]['CTRL_ID']: CTRL_ID = ""
        text = "%s&OBST ID='%s', %s, %s%s%s/\n"%(text, ID, loc, SURF_IDtxt, BNDF_OBST, CTRL_ID)
    return text

def makeINIT(inits):
    text = ''
    for key in list(inits.keys()):
        XYZ = inits[key]['XYZ']
        XB = inits[key]['XB']
        if 'XYZ' in list(inits[key].keys()): loc = "XYZ = %0.4f,%0.4f,%0.4f"%(XYZ[0], XYZ[1], XYZ[2])
        if 'XB' in list(inits[key].keys()): loc = "XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        ID = inits[key]['ID']
        PART_ID = inits[key]['PART_ID']
        N_PARTICLES = inits[key]['N_PARTICLES']
        text = "%s&INIT ID='%s', %s, N_PARTICLES=%0.0f, PART_ID=%s, /\n"%(text, ID, loc, PART_ID, N_PARTICLES)
    return text

def makeDEVC(devcs):
    text = ''
    for key in list(devcs.keys()):
        XYZ = devcs[key]['XYZ']
        XB = devcs[key]['XB']
        if np.any(devcs[key]['XYZ']): loc = "XYZ = %0.4f,%0.4f,%0.4f"%(XYZ[0], XYZ[1], XYZ[2])
        if np.any(devcs[key]['XB']): loc = "XB = %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f"%(XB[0], XB[1], XB[2], XB[3], XB[4], XB[5])
        ID = devcs[key]['ID']
        QUANTITY = devcs[key]['QUANTITY']
        if devcs[key]['INITIAL_STATE']:
            INITIAL_STATE = "INITIAL_STATE = %s, "%(devcs[key]['INITIAL_STATE'])
        else:
            INITIAL_STATE = ''
        if devcs[key]['INIT_ID']:
            INIT_ID = "INIT_ID = %s, "%(devcs[key]['INIT_ID'])
        else:
            INIT_ID = ''
        if devcs[key]['SETPOINT']:
            SETPOINT = "SETPOINT = %0.4f, "%(devcs[key]['SETPOINT'])
        else:
            SETPOINT = ''
        if devcs[key]['IOR']:
            IOR = "IOR = %0.0f, "%(devcs[key]['IOR'])
        else:
            IOR = ""
        TIME_AVERAGED = ""
        if devcs[key]['TIME_AVERAGED']:
            if devcs[key]['TIME_AVERAGED'] == 'TRUE':
                TIME_AVERAGED = "TIME_AVERAGED = .TRUE., "
        text = "%s&DEVC ID='%s', QUANTITY='%s', %s, %s%s%s%s%s/\n"%(text, ID, QUANTITY, loc, INITIAL_STATE, INIT_ID, SETPOINT, IOR, TIME_AVERAGED)
    return text

def makeHEAD(head):
    chid = head['chid']
    title = head['title']
    if not title: title = chid
    text = "&HEAD CHID = '%s', TITLE = '%s' /\n"%(chid, title)
    return text
  
def parseHEAD(line):
    head = defaultdict(bool)
    if "CHID" in line: head['chid'] = line.split("CHID")[1].split("'")[1]
    if "TITLE" in line: head['title'] = line.split("TITLE")[1].split("'")[1]
    return head

def parseDEVC(line):
    devc = defaultdict(bool)
    if "XYZ" in line:
        tmp = line.split("XYZ")[1]
        tmp = tmp.replace("=",'').split(',')[:3]
        coords = [float(x) for x in tmp]
        devc['XYZ'] = coords
    if "XB" in line:
        tmp = line.split("XB")[1]
        tmp = tmp.replace("=",'').split(',')[:6]
        coords = [float(x) for x in tmp]
        devc['XB'] = coords
    if "ID" in line:
        tmp = line.split("ID")[1]
        devID = tmp.split("'")[1]
        devc['ID'] = devID
    if "INITIAL_STATE" in line:
        tmp = line.split("INITIAL_STATE")[1]
        initialState = tmp.split(".")[1]
        if initialState == "FALSE": initialState = ".FALSE."
        if initialState == "TRUE": initialState = ".TRUE."
        devc['INITIAL_STATE'] = initialState
    if "QUANTITY" in line:
        tmp = line.split("QUANTITY")[1]
        quantity = tmp.split("'")[1]
        devc['QUANTITY'] = quantity
    if "INIT_ID" in line:
        tmp = line.split("INIT_ID")[1]
        initID = tmp.split("'")[1]
        devc['INIT_ID'] = initID
    if "SETPOINT" in line:
        tmp = float(line.split("SETPOINT")[1].split("=")[1].split(',')[0].replace(" ",""))
        devc["SETPOINT"] = tmp
    if "IOR" in line:
        tmp = float(line.split("IOR")[1].split("=")[1].split(',')[0].replace(" ",""))
        devc["IOR"] = tmp
    if "TIME_AVERAGED" in line:
        devc['TIME_AVERAGED'] = line.split("TIME_AVERAGED")[1].split('.')[1]
    return devc

def parseINIT(line):
    init = defaultdict(bool)
    if "XYZ" in line:
        tmp = line.split("XYZ")[1]
        tmp = tmp.replace("=",'').split(',')[:3]
        coords = [float(x) for x in tmp]
        init['XYZ'] = coords
    if "XB" in line:
        tmp = line.split("XB")[1]
        tmp = tmp.replace("=",'').split(',')[:6]
        coords = [float(x) for x in tmp]
        init['XB'] = coords
    if "ID" in line:
        tmp = line.split("ID")[1]
        devID = tmp.split("'")[1]
        init['ID'] = devID
    if "PART_ID" in line:
        tmp = line.split("PART_ID")[1]
        partID = tmp.split("'")[1]
        init['PART_ID'] = partID
    if "N_PARTICLES" in line:
        tmp = line.split("N_PARTICLES")[1]
        tmp = tmp.replace("=",'')
        nParticles = int(tmp.split(",")[0])
        init['N_PARTICLES'] = nParticles
    return init

def parseOBST(line):
    obst = defaultdict(bool)
    if "SURF_ID" in line:
        obst['SURF_ID'] = line.split("SURF_ID")[1].split("'")[1]
    if "ID" in line.replace("SURF_ID",""):
        obst['ID'] = line.replace("SURF_ID","").split("ID")[1].split("'")[1]
    if "XB" in line:
        tmp = np.array([float(x) for x in line.split("XB")[1].replace("=","").split(',')[:6]])
        obst['XB'] =  tmp
    if "BNDF_OBST" in line:
        txt = line.split('BNDF_OBST=')[1].split('.')[1]
        if "TRUE" in txt:
            obst['BNDF_OBST'] = True
        if "FALSE" in txt:
            obst['BNDF_OBST'] = False
    if "CTRL_ID" in line:
        txt = line.split('CTRL_ID')[1].split("'")[1]
        obst['CTRL_ID'] = txt
    return obst

def parseVENT(line):
    vent = defaultdict(bool)
    if "SURF_ID" in line:
        vent['SURF_ID'] = line.split("SURF_ID")[1].split("'")[1]
    if "ID" in line.replace("SURF_ID",""):
        vent['ID'] = line.replace("SURF_ID","").split("ID")[1].split("'")[1]
    if "XB" in line:
        tmp = np.array([float(x) for x in line.split("XB")[1].replace("=","").split(',')[:6]])
        vent['XB'] =  tmp
    if "CTRL_ID" in line:
        tmp = line.split('CTRL_ID')[1].split('=')[1].split(',')[0].split("'")[1]
        vent['CTRL_ID'] = tmp
    if 'MB' in line:
        vent['MB'] = line.split('MB')[1].split("'")[1]
    return vent

def parseSURF(line2):
    line = line2.replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',')
    line = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
    
    keys = line.split(',')
    keys[0] = keys[0].replace('SURF','')
    updatedKeys = []
    txt = ''
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    surf = defaultdict(bool)
    for key in updatedKeys:
        if ("ID" in key) and ('MATL_ID' not in key) and ('VOID' not in key) and ('WIDTH' not in key):
            surf['ID'] = key.split('ID')[1].split("'")[1]
        if "MATL_ID" in key:
            tmp = key.split("=")[1].split("'")[1:]
            if ',' in tmp: tmp.remove(',')
            if '' in tmp: tmp.remove('')
            surf['MATL_ID'] = tmp
        if "COLOR" in key:
            surf['COLOR'] = key.split('COLOR')[1].split("'")[1]
        if "THICKNESS" in key:
            tmp = key.split('=')[1].replace(' ','').split(',')
            if '' in tmp: tmp.remove('')
            surf['THICKNESS'] = [float(x) for x in tmp]                
        if "BACKING" in key:
            surf['BACKING'] = key.split('BACKING')[1].split("'")[1]
        if "GEOMETRY" in key:
            surf['GEOMETRY'] = key.split('GEOMETRY')[1].split("'")[1]
        if "FYI" in key:
            surf['FYI'] = key.split('FYI')[1].split("'")[1]
        if "LENGTH" in key:
            surf['LENGTH'] = float(key.split('LENGTH')[1].replace(' ','').replace('=','').split(',')[0])
        if "LEAK_PATH" in key:
            tmp = key.split('LEAK_PATH')[1].replace(' ','').replace('=','').split(',')
            if '' in tmp: tmp.remove('')
            tmp = [float(x) for x in tmp]
            surf['LEAK_PATH'] = tmp
        if "HRRPUA" in key:
            tmp = key.split('=')[1].split(',')[0].replace(' ','')
            surf['HRRPUA'] = float(tmp)
        if "RAMP_Q" in key:
            tmp = key.split('=')[1].split("'")[1]
            surf['RAMP_Q'] = tmp
        if "RGB" in key:
            surf['RGB'] = [float(x) for x in key.split('=')[1].split(',')[:2]]
        if ("ADIABATIC" in key) and ("ID" not in key):
            surf['ADIABATIC'] = key.split('.')[1]
            
    return surf

def parseRAMP(line):
    keys = line.split(',')
    updatedKeys = []
    txt = ''
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    ramp = defaultdict(bool)
    for key in updatedKeys:
        if ("ID" in key) and ("DEVC_ID" not in key):
            ramp['ID'] = key.split("'")[1]
        elif ("DEVC_ID" in key):
            ramp['DEVC_ID'] = key.split("'")[1]
        elif "T" in key:
            ramp['T'] = float(key.split("T")[1].split('=')[1].split(',')[0].replace(' ',''))
        elif "F" in key:
            ramp['F'] = float(key.split("F")[1].split('=')[1].split(',')[0].replace(' ',''))
    return ramp

def parseCTRL(line2):
    line = line2.replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',')
    line = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
    
    keys = line.split(',')
    keys[0] = keys[0].replace('CTRL ','')
    updatedKeys = []
    txt = ''
    for i in range(0,len(keys)):
        if '=' in keys[i]:
            updatedKeys.append(txt)
            txt = keys[i]
        else:
            txt = ','.join([txt,keys[i]])
    updatedKeys.append(txt)
    ctrl = defaultdict(bool)
    for key in updatedKeys:
        if ("ID" in key) and ('INPUT_ID' not in key):
            ctrl['ID'] = key.split('ID')[1].split("'")[1]
        if "INPUT_ID" in key:
            tmp = key.split("=")[1].split("'")[1:]
            if ',' in tmp: tmp.remove(',')
            if '' in tmp: tmp.remove('')
            if ' ' in tmp: tmp.remove(' ')
            ctrl['INPUT_ID'] = tmp
        if "FUNCTION_TYPE" in key:
            tmp = key.split('=')[1].split("'")[1]
            ctrl['FUNCTION_TYPE'] = tmp
        if ("DELAY" in key) and ("'TIME_DELAY'" not in key):
            tmp = float(key.split('=')[1].split(',')[0].replace(" ",""))
            ctrl['DELAY'] = tmp
    return ctrl

def splitLineIntoKeys(line2):
    line = line2.replace('\n',',').replace(',,',',').replace(',,',',').replace(',,',',')
    line = line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')
    
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
    return updatedKeys

def parseSLCF(line):
    slcf = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'QUANTITY' in key: slcf['QUANTITY'] = key.split('=')[1].split("'")[1]
        if 'PBX' in key: slcf['PBX'] = float(key.split('=')[1].split(',')[0])
        if 'PBY' in key: slcf['PBY'] = float(key.split('=')[1].split(',')[0])
        if 'PBZ' in key: slcf['PBZ'] = float(key.split('=')[1].split(',')[0])
        if 'VECTOR' in key: slcf['VECTOR'] = key.split('=')[1].split('.')[1]
        if "XB" in key:
            slcf['XB'] = [float(x) for x in key.split('XB')[1].replace('=','').split(',')[:6]]
    return slcf

def parseBNDF(line):
    bndf = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'QUANTITY' in key: bndf['QUANTITY'] = key.split('=')[1].split("'")[1]
    return bndf

def parseTIME(line):
    time = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'T_BEGIN' in key: time['T_BEGIN'] = float(key.split('=')[1].split(',')[0])
        if 'T_END' in key: time['T_END'] = float(key.split('=')[1].split(',')[0])
    return time

def parseDUMP(line):
    dump = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'RENDER_FILE' in key: dump['RENDER_FILE'] = key.split("'")[1]
        if 'COLUMN_DUMP_LIMIT' in key: dump['COLUMN_DUMP_LIMIT'] = key.split('.')[1]
        if 'DT_CTRL' in key: dump['DT_CTRL'] = float(key.split('=')[1].split(',')[0])
        if 'DT_DEVC' in key: dump['DT_DEVC'] = float(key.split('=')[1].split(',')[0])
        if 'DT_BNDF' in key: dump['DT_BNDF'] = float(key.split('=')[1].split(',')[0])
        if 'DT_SLCF' in key: dump['DT_SLCF'] = float(key.split('=')[1].split(',')[0])
        if 'DT_SL3D' in key: dump['DT_SL3D'] = float(key.split('=')[1].split(',')[0])
        if 'DT_RESTART' in key: dump['DT_RESTART'] = float(key.split('=')[1].split(',')[0])
    return dump

def parseMISC(line):
    misc = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'TMPA' in key: misc['TMPA'] = float(key.split('=')[1].split(',')[0])
        if 'HUMIDITY' in key: misc['HUMIDITY'] = float(key.split('=')[1].split(',')[0])
        if 'SUPPRESSION' in key: misc['SUPPRESSION'] = key.split('.')[1]
        if 'BNDF_DEFAULT' in key: misc['BNDF_DEFAULT'] = key.split('.')[1]
    return misc

def parseZONE(line):
    zone = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if 'ID' in key: zone['ID'] = line.split("'")[1]
        if 'XB' in key: zone['XB'] = [float(x) for x in key.split('XB')[1].replace('=','').split(',')[:6]]
    return zone

def parseMESH(line):
    mesh = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if "ID" in key:
            mesh['ID'] = key.split('ID')[1].replace('=','').split("'")[1]
        if 'IJK' in key:
            mesh['IJK'] = [float(x) for x in key.split('IJK')[1].replace('=','').split(',')]
        if 'XB' in key:
            mesh['XB'] = [float(x) for x in key.split('XB')[1].replace('=','').split(',')[:6]]
        if 'MPI_PROCESS' in key:
            mesh['MPI_PROCESS'] = key.split('MPI_PROCESS')[1].replace('=','').replace(',','')
    return mesh

def parseREAC(line):
    reac = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if ("ID" in key) and ("IDEAL" not in key): reac['ID'] = key.split('ID')[1].replace('=','').split("'")[1]
        if "FUEL" in key: reac['FUEL'] = key.split("'")[1]
        if "FORMULA" in key: reac['FORMULA'] = key.split("'")[1]
        if "CO_YIELD" in key: reac['CO_YIELD'] = float(key.split('=')[1].split(',')[0])
        if "SOOT_YIELD" in key: reac['SOOT_YIELD'] = float(key.split('=')[1].split(',')[0])
        if "IDEAL" in key: reac['IDEAL'] = key.split('.')
    return reac

def parseMATL(line):
    matl = defaultdict(bool)
    keys = splitLineIntoKeys(line)
    for key in keys:
        if "ID" in key: matl['ID'] = key.split("'")[1]
        if "FYI" in key: matl['FYI'] = key.split("'")[1]
        if ("SPECIFIC_HEAT_RAMP" in key):
            matl['SPECIFIC_HEAT_RAMP'] = key.split('=')[1].split(',')[0]
        elif ("SPECIFIC_HEAT" in key):
            matl['SPECIFIC_HEAT'] = float(key.split('=')[1].split(',')[0])
        if ("CONDUCTIVITY_RAMP" in key):
            matl['CONDUCTIVITY_RAMP'] = key.split('=')[1].split(',')[0]
        elif ("CONDUCTIVITY" in key):
            matl['CONDUCTIVITY'] = float(key.split('=')[1].split(',')[0])
        if "DENSITY" in key: matl['DENSITY'] = float(key.split('=')[1].split(',')[0])
        if "EMISSIVITY" in key: matl['EMISSIVITY'] = float(key.split('=')[1].split(',')[0])
    return matl
