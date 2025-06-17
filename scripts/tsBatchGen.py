#!/usr/bin/env python3
import argparse
import gc
import pathlib
import os
import xml.etree.ElementTree as ETree
import datetime

class MetaSettings:
    def __init__(self):
        self.probeName = None
        self.targetName = None
        self.targetType = None
        self.speciesName = None
        self.cellType = None
        self.voxelDims = (0, 0, 0)
        
    def fromXmlNode(self, element):
        for child in element:
            if child.tag == 'Species':
                self.speciesName = cleanXmlValueString(child.text)
            elif child.tag == 'CellType':
                self.cellType = cleanXmlValueString(child.text)
            elif child.tag == 'TargetName':
                self.targetName = cleanXmlValueString(child.text)
            elif child.tag == 'ProbeName':
                self.probeName = cleanXmlValueString(child.text)
            elif child.tag == 'TargetType':
                self.targetType = cleanXmlValueString(child.text)  
            elif child.tag == 'VoxelDimsNano':
                vx = 0
                vy = 0
                vz = 0
                if 'X' in child.attrib:
                    vx = int(child.attrib['X'])
                if 'Y' in child.attrib:
                    vy = int(child.attrib['Y'])
                if 'Z' in child.attrib:
                    vz = int(child.attrib['Z'])
                self.voxelDims = (vx, vy, vz)
            else:
                print("Tag \"", child.tag, "\" not recognized for Meta")
                
    def mergeVoxelDims(self, other):
        newdims = [0,0,0]
        if self.voxelDims[0] == 0:
            newdims[0] = other.voxelDims[0]
        else:
            newdims[0] = self.voxelDims[0]
        if self.voxelDims[1] == 0:
            newdims[1] = other.voxelDims[1]
        else:
            newdims[1] = self.voxelDims[1]
        if self.voxelDims[2] == 0:
            newdims[2] = other.voxelDims[2]
        else:
            newdims[2] = self.voxelDims[2]
        self.voxelDims = (newdims[0], newdims[1], newdims[2])
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (self.probeName is None):
            self.probeName = other.probeName
        if overwrite or (self.targetName is None):
            self.targetName = other.targetName
        if overwrite or (self.targetType is None):
            self.targetType = other.targetType
        if overwrite or (self.speciesName is None):
            self.speciesName = other.speciesName
        if overwrite or (self.cellType is None):
            self.cellType = other.cellType
        if overwrite:
            self.voxelDims = other.voxelDims
        else:
            self.mergeVoxelDims(other)
        
class CellsegSettings:
    def __init__(self):
        self.templateName = None
        self.lightZMin = 0
        self.lightZMax = 0
        self.nucZMin = 0
        self.nucZMax = 0
        
    def fromXmlNode(self, element):
        for child in element:
            if child.tag == 'PresetName':
                self.templateName = cleanXmlValueString(child.text)
            elif child.tag == 'TransZTrim':
                if 'Min' in child.attrib:
                    self.lightZMin = int(child.attrib['Min'])
                if 'Max' in child.attrib:
                    self.lightZMax = int(child.attrib['Max'])
            elif child.tag == 'NucZTrim':
                if 'Min' in child.attrib:
                    self.nucZMin = int(child.attrib['Min'])
                if 'Max' in child.attrib:
                    self.nucZMax = int(child.attrib['Max'])            
            else:
                print("Tag \"", child.tag, "\" not recognized for cellseg settings")
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (self.templateName is None):
            self.templateName = other.templateName
        if overwrite or (self.lightZMin <= 0):
            self.lightZMin = other.lightZMin
        if overwrite or (self.lightZMax <= 0):
            self.lightZMax = other.lightZMax
        if overwrite or (self.nucZMin <= 0):
            self.nucZMin = other.nucZMin
        if overwrite or (self.nucZMax <= 0):
            self.nucZMax = other.nucZMax
        
class SpotDetectSettings:
    def __init__(self):
        self.thPreset = 0
        self.gaussRad = 0
        
    def fromXmlNode(self, element):
        if 'GaussRad' in element.attrib:
            self.gaussRad = int(element.attrib['GaussRad'])
        for child in element:
            if child.tag == 'ThresholdSettings':
                if 'Preset' in child.attrib:
                    self.thPreset = int(child.attrib['Preset'])  
            else:
                print("Tag \"", child.tag, "\" not recognized for spot detect settings")        
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (self.thPreset == 0):
            self.thPreset = other.thPreset
        if overwrite or (self.gaussRad <= 0):
            self.gaussRad = other.gaussRad 

class QuantSettings:
    def __init__(self):
        self.qNoClouds = True
        self.qCellZero = True
    
    def fromXmlNode(self, element):
        if 'DoClouds' in element.attrib:
            self.qNoClouds = not bool(element.attrib['DoClouds'])
        if 'CellZero' in element.attrib:
            self.qCellZero = bool(element.attrib['CellZero'])
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or self.qNoClouds:
            self.qNoClouds = other.qNoClouds  
        if overwrite or self.qNoClouds:
            self.qCellZero = other.qCellZero     
        
class JobSettings:
    def __init__(self):
        self.cpuCount = 0
        self.ramGigs = 0
        self.timeString = None
        
    def fromXmlNode(self, element):
        if 'CpuCount' in element.attrib:
            self.cpuCount = int(element.attrib['CpuCount'])
        if 'RamGigs' in element.attrib:
            self.ramGigs = int(element.attrib['RamGigs'])
        if 'Time' in element.attrib:
            self.timeString = element.attrib['Time']
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (self.cpuCount <= 0):
            self.cpuCount = other.cpuCount
        if overwrite or (self.ramGigs <= 0):
            self.ramGigs = other.ramGigs
        if overwrite or (self.timeString is None):
            self.timeString = other.timeString
            
class IDistroSettings:
    def __init__(self):
        self.correctionMtxPath = None
        self.forceNoProbe = False
        self.forceNoDPC = True
        self.zTrimMin = 0
        self.zTrimMax = 0
        
        self.jobSettings = JobSettings()
        self.jobSettings.cpuCount = 4
        self.jobSettings.ramGigs = 128
        self.jobSettings.timeString = "12:00:00"
        
    def fromXmlNode(self, element):
        if 'CorrectionMtx' in element.attrib:
            self.correctionMtxPath = element.attrib['CorrectionMtx']
        if 'ForceNoProbe' in element.attrib:
            self.forceNoProbe = bool(element.attrib['ForceNoProbe'])
        if 'ForceNoDPC' in element.attrib:
            self.forceNoDPC = bool(element.attrib['ForceNoDPC'])  
        if 'TrimZMin' in element.attrib:
            self.zTrimMin = int(element.attrib['TrimZMin'])         
        if 'TrimZMax' in element.attrib:
            self.zTrimMax = int(element.attrib['TrimZMax'])         
        for child in element:
            if child.tag == 'JobSettings':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.jobSettings is not None:
                    self.jobSettings.copyDataFrom(childnode, True)
                else:
                    self.jobSettings = childnode                 
            else:
                print("Tag \"", child.tag, "\" not recognized for idistro settings")        
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (not self.forceNoProbe):
            self.forceNoProbe = other.forceNoProbe
        if overwrite or self.forceNoDPC:
            self.forceNoDPC = other.forceNoDPC
        if overwrite or (self.correctionMtxPath is None):
            self.correctionMtxPath = other.correctionMtxPath
        if overwrite or (self.zTrimMin <= 0):
            self.zTrimMin = other.zTrimMin 
        if overwrite or (self.zTrimMax <= 0):
            self.zTrimMax = other.zTrimMax 
        if self.jobSettings is None:
            self.jobSettings = JobSettings()
            self.jobSettings.copyDataFrom(other.jobSettings, True)
        else:
            self.jobSettings.copyDataFrom(other.jobSettings, overwrite)
        
        
class ImageChannelInfo:
    def __init__(self):
        self.channelNumber = 0
        self.metaData = None
        self.spotsSettings = None
        self.quantSettings = None
        self.spotsJobSettings = None
        self.parentBatch = None
        
        self.dirName = None
        self.scriptPath = None
        self.quantResPath = None
    
    def fromXmlNode(self, element):
        if 'ChannelNumber' in element.attrib:
            self.channelNumber = int(element.attrib['ChannelNumber'])
        for child in element:
            if child.tag == 'Meta':
                metanode = MetaSettings()
                metanode.fromXmlNode(child)
                if self.metaData is not None:
                    self.metaData.copyDataFrom(metanode, True)
                else:
                    self.metaData = metanode
            elif child.tag == 'SpotDetectSettings':
                childnode = SpotDetectSettings()
                childnode.fromXmlNode(child)
                if self.spotsSettings is not None:
                    self.spotsSettings.copyDataFrom(childnode, True)
                else:
                    self.spotsSettings = childnode    
            elif child.tag == 'QuantSettings':
                childnode = QuantSettings()
                childnode.fromXmlNode(child)
                if self.quantSettings is not None:
                    self.quantSettings.copyDataFrom(childnode, True)
                else:
                    self.quantSettings = childnode
            elif child.tag == 'SpotJob':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.spotsJobSettings is not None:
                    self.spotsJobSettings.copyDataFrom(childnode, True)
                else:
                    self.spotsJobSettings = childnode              
            else:
                print("Tag \"", child.tag, "\" not recognized for channel settings")
                
    def copyFromParent(self, parent, overwrite):
        if parent is None:
            return

        if self.metaData is None:
            self.metaData = MetaSettings()
        self.metaData.copyDataFrom(parent.metaData, overwrite)
        if self.spotsSettings is None:
            self.spotsSettings = SpotDetectSettings()        
        self.spotsSettings.copyDataFrom(parent.spotsSettings, overwrite)
        if self.quantSettings is None:
            self.quantSettings = QuantSettings()        
        self.quantSettings.copyDataFrom(parent.quantSettings, overwrite)
        if self.spotsJobSettings is None:
            self.spotsJobSettings = JobSettings()         
        self.spotsJobSettings.copyDataFrom(parent.spotsJobSettings, overwrite)

    def copyDataFrom(self, other, overwrite):
        if other is None:
            return        
        self.channelNumber = other.channelNumber
        copyFromParent(other, overwrite) #Contingent upon field names staying the same             
    
class ImageBatchInfo:
    def __init__(self):
        self.name = None
        self.inputDir = None
        self.outputDir = None
        self.channelCount = 0
        self.lightChannel = 0
        self.nucChannel = 0
        self.metaData = None
        self.cellsegSettings = None
        self.spotsSettings = None
        self.quantSettings = None
        self.cellsegJobSettings = None
        self.spotsJobSettings = None
        self.idistroSettings = None
        
        self.batchScriptPath = None
        self.postResScriptPath = None
        
        self.parentSet = None
        self.channels = list()
        
    def fromXmlNode(self, element):
        if 'Name' in element.attrib:
            self.name = element.attrib['Name']
            
        for child in element:
            if (child.tag == 'CommonMeta') or (child.tag == 'Meta'):
                metanode = MetaSettings()
                metanode.fromXmlNode(child)
                if self.metaData is not None:
                    self.metaData.copyDataFrom(metanode, True)
                else:
                    self.metaData = metanode
            elif child.tag == 'CellSegSettings':
                childnode = CellsegSettings()
                childnode.fromXmlNode(child)
                if self.cellsegSettings is not None:
                    self.cellsegSettings.copyDataFrom(childnode, True)
                else:
                    self.cellsegSettings = childnode             
            elif child.tag == 'SpotDetectSettings':
                childnode = SpotDetectSettings()
                childnode.fromXmlNode(child)
                if self.spotsSettings is not None:
                    self.spotsSettings.copyDataFrom(childnode, True)
                else:
                    self.spotsSettings = childnode    
            elif child.tag == 'QuantSettings':
                childnode = QuantSettings()
                childnode.fromXmlNode(child)
                if self.quantSettings is not None:
                    self.quantSettings.copyDataFrom(childnode, True)
                else:
                    self.quantSettings = childnode
            elif child.tag == 'SpotJob':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.spotsJobSettings is not None:
                    self.spotsJobSettings.copyDataFrom(childnode, True)
                else:
                    self.spotsJobSettings = childnode  
            elif child.tag == 'CellSegJob':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.cellsegJobSettings is not None:
                    self.cellsegJobSettings.copyDataFrom(childnode, True)
                else:
                    self.cellsegJobSettings = childnode
            elif child.tag == 'Paths':
                for gchild in child:
                    if gchild.tag == 'ImageDir':
                        self.inputDir = cleanXmlValueString(gchild.text)
                    elif gchild.tag == 'OutputDir':
                        self.outputDir = cleanXmlValueString(gchild.text)
            elif child.tag == 'ChannelInfo':
                if 'ChannelCount' in child.attrib:
                    self.channelCount = int(child.attrib['ChannelCount'])   
                if 'NucChannel' in child.attrib:
                    self.nucChannel = int(child.attrib['NucChannel']) 
                if 'TransChannel' in child.attrib:
                    self.lightChannel = int(child.attrib['TransChannel'])
                for gchild in child:
                    if gchild.tag == 'ImageChannel':
                        chinfo = ImageChannelInfo()
                        chinfo.copyFromParent(self, False)
                        chinfo.fromXmlNode(gchild)
                        chinfo.parentBatch = self
                        self.channels.append(chinfo)
            elif child.tag == 'IntensityDistroSettings':
                childnode = IDistroSettings()
                childnode.fromXmlNode(child)
                if self.idistroSettings is not None:
                    self.idistroSettings.copyDataFrom(childnode, True)
                else:
                    self.idistroSettings = childnode            
            else:
                print("Tag \"", child.tag, "\" not recognized for batch settings")
                
    def copyFromParent(self, parent, overwrite):
        if parent is None:
            return        
        if self.metaData is None:
            self.metaData = MetaSettings()
        self.metaData.copyDataFrom(parent.metaData, overwrite)
        if self.spotsSettings is None:
            self.spotsSettings = SpotDetectSettings()        
        self.spotsSettings.copyDataFrom(parent.spotsSettings, overwrite)
        if self.quantSettings is None:
            self.quantSettings = QuantSettings()        
        self.quantSettings.copyDataFrom(parent.quantSettings, overwrite)
        if self.cellsegJobSettings is None:
            self.cellsegJobSettings = JobSettings()         
        self.cellsegJobSettings.copyDataFrom(parent.cellsegJobSettings, overwrite)    
        if self.spotsJobSettings is None:
            self.spotsJobSettings = JobSettings()         
        self.spotsJobSettings.copyDataFrom(parent.spotsJobSettings, overwrite)
        if self.cellsegSettings is None:
            self.cellsegSettings = CellsegSettings()         
        self.cellsegSettings.copyDataFrom(parent.cellsegSettings, overwrite)
        if self.idistroSettings is None:
            self.idistroSettings = IDistroSettings()         
        self.idistroSettings.copyDataFrom(parent.idistroSettings, overwrite)    
      
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (self.name is None):
            self.name = other.name
        if overwrite or (self.inputDir is None):
            self.inputDir = other.inputDir
        if overwrite or (self.outputDir is None):
            self.outputDir = other.outputDir
        if overwrite or (self.channelCount <= 0):
            self.channelCount = other.channelCount
        if overwrite or (self.lightChannel <= 0):
            self.lightChannel = other.lightChannel
        if overwrite or (self.nucChannel <= 0):
            self.nucChannel = other.nucChannel
        copyFromParent(other, overwrite)
        
    def updateOverrides(self):
        for chinfo in self.channels:
            chinfo.copyFromParent(self, False)

class ImageBatchSet:
    def __init__(self):
        self.metaData = None
        self.cellsegSettings = None
        self.spotsSettings = None
        self.quantSettings = None
        self.cellsegJobSettings = None
        self.spotsJobSettings = None
        self.idistroSettings = None
        
        self.tsDir = None
        self.moduleName = None
        
        self.batches = list()
        
    def fromXmlNode(self, element):
        for child in element:
            if (child.tag == 'CommonMeta') or (child.tag == 'Meta'):
                metanode = MetaSettings()
                metanode.fromXmlNode(child)
                if self.metaData is not None:
                    self.metaData.copyDataFrom(metanode, True)
                else:
                    self.metaData = metanode
            elif child.tag == 'CellSegSettings':
                childnode = CellsegSettings()
                childnode.fromXmlNode(child)
                if self.cellsegSettings is not None:
                    self.cellsegSettings.copyDataFrom(childnode, True)
                else:
                    self.cellsegSettings = childnode             
            elif child.tag == 'SpotDetectSettings':
                childnode = SpotDetectSettings()
                childnode.fromXmlNode(child)
                if self.spotsSettings is not None:
                    self.spotsSettings.copyDataFrom(childnode, True)
                else:
                    self.spotsSettings = childnode    
            elif child.tag == 'QuantSettings':
                childnode = QuantSettings()
                childnode.fromXmlNode(child)
                if self.quantSettings is not None:
                    self.quantSettings.copyDataFrom(childnode, True)
                else:
                    self.quantSettings = childnode
            elif child.tag == 'JobSettings':
                for gchild in child:
                    if gchild.tag == 'SpotJob':
                        childnode = JobSettings()
                        childnode.fromXmlNode(gchild)
                        if self.spotsJobSettings is not None:
                            self.spotsJobSettings.copyDataFrom(childnode, True)
                        else:
                            self.spotsJobSettings = childnode  
                    elif gchild.tag == 'CellSegJob':
                        childnode = JobSettings()
                        childnode.fromXmlNode(gchild)
                        if self.cellsegJobSettings is not None:
                            self.cellsegJobSettings.copyDataFrom(childnode, True)
                        else:
                            self.cellsegJobSettings = childnode
                    elif gchild.tag == 'TrueSpotDir':
                        self.tsDir = cleanXmlValueString(gchild.text)
                    elif gchild.tag == 'MatlabModuleName':
                        self.moduleName = cleanXmlValueString(gchild.text)   
            elif child.tag == 'IntensityDistroSettings':
                childnode = IDistroSettings()
                childnode.fromXmlNode(child)
                if self.idistroSettings is not None:
                    self.idistroSettings.copyDataFrom(childnode, True)
                else:
                    self.idistroSettings = childnode               
            elif child.tag == 'ImageBatch':
                childnode = ImageBatchInfo()
                childnode.copyFromParent(self, False)
                childnode.fromXmlNode(child)
                childnode.parentSet = self
                self.batches.append(childnode)
            else:
                print("Tag \"", child.tag, "\" not recognized for set settings")
                
    def updateOverrides(self):
        for imgbatch in self.batches:
            imgbatch.copyFromParent(self, False)
            imgbatch.updateOverrides()
   
class SingleImage:
    def __init__(self):
        self.name = None
        self.tifPath = None
        self.resultsDir = None
        self.csScriptPath = None
        self.csResPath = None
        
        self.batchInfo = None
        
def cleanXmlValueString(inputStr):
    #Removes quotes from start and end if present
    if inputStr.startswith("\""):
        inputStr = inputStr[1:]
    if inputStr.endswith("\""):
        inputStr = inputStr[:-1]   
    return inputStr

def getdtstr():
    now = datetime.datetime.now()
    return "[" + str(now) + "]"

def genChannelJob(tifImage, channelInfo, tsDir, moduleName):
    chStr = 'CH' + str(channelInfo.channelNumber)
    channelInfo.dirName = chStr
    chDir = os.path.join(tifImage.resultsDir, chStr)
    os.makedirs(chDir, exist_ok=True)
    
    chName = tifImage.name + '_' + chStr
    outStemName = chName + '_spotCall'
    fullOutStem = os.path.join(chDir, outStemName)
    spotScriptPath = os.path.join(chDir, chName + '_spotsJob.sh')
    channelInfo.scriptPath = spotScriptPath
    
    scriptHandle = open(spotScriptPath, 'w')
    scriptHandle.write("#!/bin/bash\n\n")
    scriptHandle.write("module load " + moduleName + "\n")
    scriptHandle.write("if [ -s \"" + tifImage.csResPath + "\" ]; then\n")
    scriptHandle.write("\tbash \"" + os.path.join(tsDir, "TrueSpot_RNASpots.sh") + "\"")
    scriptHandle.write(" -input \"" + tifImage.tifPath + "\"")
    scriptHandle.write(" -outstem \"" + fullOutStem + "\"")
    scriptHandle.write(" -imgname \"" + chName + "\"")
    scriptHandle.write(" -chtotal " + str(tifImage.batchInfo.channelCount))
    if tifImage.batchInfo.lightChannel > 0:
        scriptHandle.write(" -chtrans " + str(tifImage.batchInfo.lightChannel))
    scriptHandle.write(" -chsamp " + str(channelInfo.channelNumber))
    scriptHandle.write(" -cellseg \"" + tifImage.csResPath + "\"")
    if channelInfo.metaData is not None:
        vxSz = channelInfo.metaData.voxelDims
        scriptHandle.write(" -voxelsize \"(" + str(vxSz[0]) + "," + str(vxSz[1]) + "," + str(vxSz[2]) + ")\"")
        if channelInfo.metaData.probeName is not None:
            scriptHandle.write(" -probetype \"" + channelInfo.metaData.probeName + "\"")
        if channelInfo.metaData.targetName is not None:
            scriptHandle.write(" -target \"" + channelInfo.metaData.targetName + "\"") 
        if channelInfo.metaData.targetType is not None:
            scriptHandle.write(" -targettype \"" + channelInfo.metaData.targetType + "\"")
        if channelInfo.metaData.speciesName is not None:
            scriptHandle.write(" -species \"" + channelInfo.metaData.speciesName + "\"")     
        if channelInfo.metaData.cellType is not None:
            scriptHandle.write(" -celltype \"" + channelInfo.metaData.cellType + "\"")  
    if channelInfo.spotsSettings is not None:
        if channelInfo.spotsSettings.gaussRad > 0:
            scriptHandle.write(" -gaussrad " + str(channelInfo.spotsSettings.gaussRad))
        if channelInfo.spotsSettings.thPreset > 0:
            scriptHandle.write(" -sensitivity " + str(channelInfo.spotsSettings.thPreset))
        elif channelInfo.spotsSettings.thPreset == 0:
            scriptHandle.write(" -sensitivity 0")        
        elif channelInfo.spotsSettings.thPreset < 0:
            scriptHandle.write(" -precision " + str(channelInfo.spotsSettings.thPreset * -1))
    scriptHandle.write(" -autominth -automaxth")
    if channelInfo.spotsJobSettings is not None:
        if channelInfo.spotsJobSettings.cpuCount > 0:
            scriptHandle.write(" -threads " + str(channelInfo.spotsJobSettings.cpuCount))
    scriptHandle.write(" -log \"" + os.path.join(chDir, chName + '_spots_mat.log') + "\"")
    scriptHandle.write("\n")
    scriptHandle.write("else\n")
    scriptHandle.write("\techo -e \"Cellseg results not found! Terminating...\"\n")
    scriptHandle.write("\texit 1\n")
    scriptHandle.write("fi\n\n")
    
    #Quant
    scriptHandle.write("if [ -s \"" + fullOutStem + "_callTable.mat\" ]; then\n")
    scriptHandle.write("\tbash \"" + os.path.join(tsDir, "TrueSpot_RNAQuant.sh") + "\"")
    scriptHandle.write(" -runinfo \"" + fullOutStem + "_rnaspotsrun.mat\"")
    if channelInfo.quantSettings is not None:
        if channelInfo.quantSettings.qNoClouds:
            scriptHandle.write(" -noclouds")
        if channelInfo.quantSettings.qCellZero:
            scriptHandle.write(" -cellzero")        
    scriptHandle.write(" -norefilter") 
    scriptHandle.write(" -log \"" + os.path.join(chDir, chName + '_quant_mat.log') + "\"")
    scriptHandle.write("\n")
    scriptHandle.write("else\n")
    scriptHandle.write("\techo -e \"Call table not found! Can't do quant! Terminating...\"\n")
    scriptHandle.write("\texit 1\n")
    scriptHandle.write("fi\n\n")
    channelInfo.quantResPath = fullOutStem + "_quantData.mat"
    
    scriptHandle.close()
    return channelInfo

def genImageJobs(tifImage, batchInfo):
    os.makedirs(tifImage.resultsDir, exist_ok=True)
    tifImage.csResPath = os.path.join(tifImage.resultsDir, "CellSeg_" + tifImage.name + ".mat")
    tifImage.batchInfo = batchInfo
    
    scriptHandle = open(tifImage.csScriptPath, 'w')
    scriptHandle.write("#!/bin/bash\n\n")
    scriptHandle.write("module load " + batchInfo.parentSet.moduleName + "\n")
    scriptHandle.write("if [ ! -s \"" + tifImage.csResPath + "\" ]; then\n")
    scriptHandle.write("\tbash \"" + os.path.join(batchInfo.parentSet.tsDir, "TrueSpot_CellSeg.sh") + "\"")
    scriptHandle.write(" -input \"" + tifImage.tifPath + "\"")
    scriptHandle.write(" -outpath \"" + tifImage.resultsDir + "\"")
    scriptHandle.write(" -imgname \"" + tifImage.name + "\"")
    scriptHandle.write(" -chtotal " + str(batchInfo.channelCount))
    scriptHandle.write(" -chlight " + str(batchInfo.lightChannel))
    scriptHandle.write(" -chnuc " + str(batchInfo.nucChannel))
    if batchInfo.lightChannel <= 0:
        scriptHandle.write(" -nuconly")
    if batchInfo.cellsegSettings is not None:
        if batchInfo.cellsegSettings.templateName is not None:
            scriptHandle.write(" -template \"" + batchInfo.cellsegSettings.templateName + "\"")
        if batchInfo.cellsegSettings.lightZMin > 0:
            scriptHandle.write(" -lightzmin " + str(batchInfo.cellsegSettings.lightZMin)) 
        if batchInfo.cellsegSettings.lightZMax > 0:
            scriptHandle.write(" -lightzmax " + str(batchInfo.cellsegSettings.lightZMax))     
        if batchInfo.cellsegSettings.nucZMin > 0:
            scriptHandle.write(" -nuczmin " + str(batchInfo.cellsegSettings.nucZMin)) 
        if batchInfo.cellsegSettings.nucZMax > 0:
            scriptHandle.write(" -nuczmax " + str(batchInfo.cellsegSettings.nucZMax))
    scriptHandle.write(" -log \"" + os.path.join(tifImage.resultsDir, tifImage.name + '_cellseg_mat.log') + "\"")
    scriptHandle.write("\n")
    scriptHandle.write("else\n")
    scriptHandle.write("\techo -e \"Cellseg results found! Skipping cell segmentation...\"\n")
    scriptHandle.write("fi\n\n")
    
    #Submit jobs for channels
    for chInfo in batchInfo.channels:
        chInfo = genChannelJob(tifImage, chInfo, batchInfo.parentSet.tsDir, batchInfo.parentSet.moduleName)
        chDir = os.path.join(tifImage.resultsDir, chInfo.dirName)
        scriptHandle.write("if [ -s \"" + chInfo.quantResPath + "\" ]; then\n")
        scriptHandle.write("\techo -e \"Quant data already found! Skipping...\"\n")
        scriptHandle.write("else\n")
        scriptHandle.write("\tchmod 774 \"" + chInfo.scriptPath + "\"\n")
        scriptHandle.write("\tsbatch --job-name=\"TS_" + tifImage.name + "_" + chInfo.dirName + "\"")
        if chInfo.spotsJobSettings is not None:
            scriptHandle.write(" --cpus-per-task=" + str(chInfo.spotsJobSettings.cpuCount))
            scriptHandle.write(" --mem=" + str(chInfo.spotsJobSettings.ramGigs) + "g")
            scriptHandle.write(" --time=" + str(chInfo.spotsJobSettings.timeString))
        scriptHandle.write(" --error=\""  + os.path.join(chDir, tifImage.name + "_" + chInfo.dirName + '_tsSlurm.err') + "\"")
        scriptHandle.write(" --out=\""  + os.path.join(chDir, tifImage.name + "_" + chInfo.dirName + '_tsSlurm.out') + "\"")
        scriptHandle.write(" \"" + chInfo.scriptPath + "\"\n")
        scriptHandle.write("fi\n\n")

    scriptHandle.close()
    
    return tifImage

def genIDistro(idisSettings, trgDir, tsDir, moduleName):
    idistroScriptPath = os.path.join(trgDir, 'idistroRun.sh');
    os.makedirs(trgDir, exist_ok=True)
    
    scriptHandle = open(idistroScriptPath, 'w')
    scriptHandle.write("#!/bin/bash\n\n")
    scriptHandle.write("module load " +  moduleName + "\n")
    
    scriptHandle.write("matlab -nodisplay -nosplash -logfile \"")
    scriptHandle.write(os.path.join(trgDir, 'idistroDump.log') + "\"")
    scriptHandle.write(" -r \"cd '" + tsDir + "/src';")
    scriptHandle.write(" Main_QCIntensityDistro('-input', '" + trgDir + "'")
    if idisSettings.zTrimMin > 0:
        scriptHandle.write(", '-zmin', '" + str(idisSettings.zTrimMin) + "'")
    if idisSettings.zTrimMax > 0:
        scriptHandle.write(", '-zmax', '" + str(idisSettings.zTrimMax) + "'")
    if idisSettings.forceNoProbe:
        scriptHandle.write(", '-noprobe'")
    if idisSettings.forceNoDPC:
        scriptHandle.write(", '-nodpc'")
    if (idisSettings.jobSettings is not None) and (idisSettings.jobSettings.cpuCount > 1):
        scriptHandle.write(", '-workers', '" + str(idisSettings.jobSettings.cpuCount) + "'")
    if idisSettings.correctionMtxPath is not None:
        scriptHandle.write(", '-correctionmtx', '" + idisSettings.correctionMtxPath + "'")
        
    scriptHandle.write("); quit;\"\n")    
    scriptHandle.close()
    
    return idistroScriptPath

def genBatch(myBatch):
    print(getdtstr(), "Working on batch " + myBatch.name)
    
    #Get list of tif images
    #https://realpython.com/get-all-files-in-directory-python/
    tifList = list()
    dirHandle = pathlib.Path(myBatch.inputDir)
    for dirItem in dirHandle.iterdir():
        if dirItem.is_file():
            sfx = dirItem.suffix.lower()
            if (sfx == '.tif') or (sfx == '.tiff'):
                myImg = SingleImage()
                myImg.name = dirItem.stem
                myImg.tifPath = os.path.join(myBatch.inputDir, dirItem.name)
                myImg.resultsDir = os.path.join(myBatch.outputDir, myImg.name)
                myImg.csScriptPath = os.path.join(myImg.resultsDir, myImg.name + '_tsCellSeg.sh')
                tifList.append(myImg)
                
    if len(tifList) < 1:
        print(getdtstr(), "\tNo tif images found in", myBatch.inputDir)
        return
    
    #Create output and base batch script
    os.makedirs(myBatch.outputDir, exist_ok=True)
    myBatch.batchScriptPath = os.path.join(myBatch.outputDir, 'tsSlurmBatch.sh');
    batchScriptHandle = open(myBatch.batchScriptPath, 'w')
    batchScriptHandle.write("#!/bin/bash\n\n")
    
    for tifImage in tifList:
        print(getdtstr(), "\tImage found:", tifImage.name)
        tifImage = genImageJobs(tifImage, myBatch)
        batchScriptHandle.write("chmod 774 \"" + tifImage.csScriptPath + "\"\n")
        batchScriptHandle.write("sbatch --job-name=\"TSCS_" + tifImage.name + "\"")
        if myBatch.cellsegJobSettings is not None:
            batchScriptHandle.write(" --cpus-per-task=" + str(myBatch.cellsegJobSettings.cpuCount))
            batchScriptHandle.write(" --mem=" + str(myBatch.cellsegJobSettings.ramGigs) + "g")
            batchScriptHandle.write(" --time=" + str(myBatch.cellsegJobSettings.timeString))
        batchScriptHandle.write(" --error=\""  + os.path.join(tifImage.resultsDir, tifImage.name + '_tscsSlurm.err') + "\"")
        batchScriptHandle.write(" --out=\""  + os.path.join(tifImage.resultsDir, tifImage.name + '_tscsSlurm.out') + "\"")
        batchScriptHandle.write(" \"" + tifImage.csScriptPath + "\"\n\n")        

    batchScriptHandle.close()
    return myBatch

def genPostJobs(myBatch):
    #XML for count dump
    xmlPath = os.path.join(myBatch.outputDir, 'countInfo.xml')
    outHandle = open(xmlPath, 'w')
    outHandle.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    outHandle.write("<BatchSamples>\n")
    for channelInfo in myBatch.channels:
        outHandle.write("\t<SampleChannel")
        outHandle.write(" Name=\"CH" + str(channelInfo.channelNumber) + "\"")
        outHandle.write(" ChannelIndex=\"" + str(channelInfo.channelNumber) + "\"")
        outHandle.write(">\n")
        outHandle.write("\t\t<QueryParams>\n")
        outHandle.write("\t\t\t<QueryParam Key=\"ChannelNo\" Value=\"" + str(channelInfo.channelNumber) + "\"/>\n")
        outHandle.write("\t\t\t<QueryParam Key=\"INameContains\" Value=\"_CH" + str(channelInfo.channelNumber) + "\"/>\n")
        if channelInfo.metaData is not None:
            if channelInfo.metaData.targetName is not None:
                outHandle.write("\t\t\t<QueryParam Key=\"TargetName\" Value=\"" + channelInfo.metaData.targetName + "\"/>\n")
            if channelInfo.metaData.probeName is not None:
                outHandle.write("\t\t\t<QueryParam Key=\"ProbeName\" Value=\"" + channelInfo.metaData.probeName + "\"/>\n")
        outHandle.write("\t\t</QueryParams>\n")
        outHandle.write("\t\t<CountThresholds ThMid=\"0\" ThLow=\"0\" ThHigh=\"0\"/>\n")
        outHandle.write("\t</SampleChannel>\n")
    outHandle.write("</BatchSamples>\n")
    outHandle.close()
    
    #Standard QC
    stem1 = os.path.join(myBatch.outputDir, 'procResQC')
    scriptPath1 = stem1 + ".sh"
    outHandle = open(scriptPath1, 'w')
    outHandle.write("#!/bin/bash\n\n")
    outHandle.write("module load " + myBatch.parentSet.moduleName + "\n")
    outHandle.write("matlab -nodisplay -nosplash -logfile \"")
    outHandle.write(os.path.join(myBatch.outputDir, 'procResMATQC.log') + "\"")
    outHandle.write(" -r \"cd '" + myBatch.parentSet.tsDir + "/src';")
    outHandle.write(" Main_QCSummary('-input', '" + myBatch.outputDir + "'); quit;\"\n")
    outHandle.close()  
    
    #Threshold assessment
    stem2 = os.path.join(myBatch.outputDir, 'procResThreshAssess')
    scriptPath2 = stem2 + ".sh"
    outHandle = open(scriptPath2, 'w')
    outHandle.write("#!/bin/bash\n\n")
    outHandle.write("module load " + myBatch.parentSet.moduleName + "\n")
    outHandle.write("matlab -nodisplay -nosplash -logfile \"")
    outHandle.write(os.path.join(myBatch.outputDir, 'procResThreshAssess.log') + "\"")
    outHandle.write(" -r \"cd '" + myBatch.parentSet.tsDir + "/src';")
    outHandle.write(" Main_AnalyzeBatchThresholds('-input', '" + myBatch.outputDir + "'); quit;\"\n")
    outHandle.close()      
    
    #Count dump (no XML)
    stem3 = os.path.join(myBatch.outputDir, 'procResQuantDumpAuto')
    scriptPath3 = stem3 + ".sh"
    outHandle = open(scriptPath3, 'w')
    outHandle.write("#!/bin/bash\n\n")
    outHandle.write("module load " + myBatch.parentSet.moduleName + "\n")
    outHandle.write("matlab -nodisplay -nosplash -logfile \"")
    outHandle.write(os.path.join(myBatch.outputDir, 'procResQuantDumpAuto.log') + "\"")
    outHandle.write(" -r \"cd '" + myBatch.parentSet.tsDir + "/src';")
    outHandle.write(" Main_DumpQuantResults('-input', '" + myBatch.outputDir + "'); quit;\"\n")
    outHandle.close()    
    
    #Count dump (with XML, no autorun)
    stem4 = os.path.join(myBatch.outputDir, 'procResQuantDumpXML')
    scriptPath4 = stem4 + ".sh"
    outHandle = open(scriptPath4, 'w')
    outHandle.write("#!/bin/bash\n\n")
    outHandle.write("module load " + myBatch.parentSet.moduleName + "\n")
    outHandle.write("matlab -nodisplay -nosplash -logfile \"")
    outHandle.write(os.path.join(myBatch.outputDir, 'procResQuantDumpXML.log') + "\"")
    outHandle.write(" -r \"cd '" + myBatch.parentSet.tsDir + "/src';")
    outHandle.write(" Main_DumpQuantResults('-input', '" + myBatch.outputDir + "',")
    outHandle.write(" '-chdef', '" + xmlPath + "'); quit;\"\n")
    outHandle.close()
    
    #Intensity distribution analysis
    stem5 = os.path.join(myBatch.outputDir, 'idistroJob')
    if myBatch.idistroSettings is not None:
        idistroScriptPath = genIDistro(myBatch.idistroSettings, myBatch.outputDir, myBatch.parentSet.tsDir, myBatch.parentSet.moduleName)
    else:
        idistroScriptPath = None
    
    #Master script
    scriptPath = os.path.join(myBatch.outputDir, 'doProcRes.sh')
    myBatch.postResScriptPath = scriptPath
    outHandle = open(scriptPath, 'w')
    outHandle.write("#!/bin/bash\n\n")
    outHandle.write("chmod 774 \"" + scriptPath1 + "\"\n")
    outHandle.write("sbatch --job-name=\"TSQC_" + myBatch.name + "\" --cpus-per-task=2 --time=8:00:00 --mem=16g")
    outHandle.write(" --error=\"" + stem1 + ".err\"")
    outHandle.write(" --out=\"" + stem1 + ".out\"")
    outHandle.write(" \"" + scriptPath1 + "\"\n")
    
    outHandle.write("chmod 774 \"" + scriptPath2 + "\"\n")
    outHandle.write("sbatch --job-name=\"TSThA_" + myBatch.name + "\" --cpus-per-task=2 --time=2:00:00 --mem=4g")
    outHandle.write(" --error=\"" + stem2 + ".err\"")
    outHandle.write(" --out=\"" + stem2 + ".out\"")
    outHandle.write(" \"" + scriptPath2 + "\"\n")    
    
    outHandle.write("chmod 774 \"" + scriptPath3 + "\"\n")
    outHandle.write("sbatch --job-name=\"TSQDA_" + myBatch.name + "\" --cpus-per-task=2 --time=2:00:00 --mem=8g")
    outHandle.write(" --error=\"" + stem3 + ".err\"")
    outHandle.write(" --out=\"" + stem3 + ".out\"")
    outHandle.write(" \"" + scriptPath3 + "\"\n")
    
    if idistroScriptPath is not None:
        outHandle.write("chmod 774 \"" + idistroScriptPath + "\"\n")
        outHandle.write("sbatch --job-name=\"TSIDis_" + myBatch.name + "\"")
        if myBatch.idistroSettings.jobSettings is not None:
            outHandle.write(" --cpus-per-task=" + str(myBatch.idistroSettings.jobSettings.cpuCount))
            outHandle.write(" --mem=" + str(myBatch.idistroSettings.jobSettings.ramGigs) + "g")
            outHandle.write(" --time=" + myBatch.idistroSettings.jobSettings.timeString)
        else:
            outHandle.write(" --cpus-per-task=2 --time=2:00:00 --mem=8g") 
        outHandle.write(" --error=\"" + stem5 + ".err\"")
        outHandle.write(" --out=\"" + stem5 + ".out\"")
        outHandle.write(" \"" + idistroScriptPath + "\"\n")    
    
    outHandle.close()        
    
    return myBatch

def readBatchXml(xmlpath):
    xmlDoc = ETree.parse(xmlpath)
    treeRoot = xmlDoc.getroot()
    
    if treeRoot.tag != 'ImageSet':
        print(getdtstr(), "Root tag must be \"ImageSet\"! Exiting...")
        return None
    
    batchSet = ImageBatchSet()
    batchSet.fromXmlNode(treeRoot)
    batchSet.updateOverrides()
    
    del(treeRoot)
    del(xmlDoc)
    gc.collect()
    return batchSet
    
def main(args):
    print("TS Batch Job Generator initiated! Version 25.06.09.02")
    print("Input Specification:", args.xmlpath)
    
    print(getdtstr(), "Reading input xml...")
    batchSet = readBatchXml(args.xmlpath)
      
    for batchInfo in batchSet.batches:
        print(getdtstr(), "Working on", batchInfo.name, "...")
        batchInfo = genBatch(batchInfo)
        batchInfo = genPostJobs(batchInfo)
        print()
        print('chmod 774 \"', batchInfo.batchScriptPath, '\"', sep='')
        print('chmod 774 \"', batchInfo.postResScriptPath, '\"', sep='')
        print('bash \"', batchInfo.batchScriptPath, '\"', sep='', flush=True)
        print()
        
if __name__ == "__main__":
    # Args
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("xmlpath", help="Input batch specification as XML")
    parser.add_argument("--help", "-h", "-?", action="help", help="Show this help message and exit.")
    args = parser.parse_args()
    main(args)