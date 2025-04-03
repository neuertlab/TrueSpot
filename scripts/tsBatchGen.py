#!/usr/bin/env python3
import argparse
import gc
import pathlib
import os
import xml.etree.ElementTree as ETree

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
                self.speciesName = child.text
            elif child.tag == 'CellType':
                self.cellType = child.text
            elif child.tag == 'TargetName':
                self.targetName = child.text
            elif child.tag == 'ProbeName':
                self.probeName = child.text
            elif child.tag == 'TargetType':
                self.targetType = child.text  
            elif child.tag == 'VoxelDimsNano':
                vx = 0
                vy = 0
                vz = 0
                if 'X' in child.attr:
                    vx = int(child.attr['X'])
                if 'Y' in child.attr:
                    vy = int(child.attr['Y'])
                if 'Z' in child.attr:
                    vz = int(child.attr['Z'])
                self.voxelDims = (vx, vy, vz)
            else:
                print("Tag \"", child.tag, "\" not recognized for Meta")
                
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
        if overwrite or (self.voxelDims[0] <= 0):
            self.voxelDims[0] = other.voxelDims[0]
        if overwrite or (self.voxelDims[1] <= 0):
            self.voxelDims[1] = other.voxelDims[1] 
        if overwrite or (self.voxelDims[2] <= 0):
            self.voxelDims[2] = other.voxelDims[2]         
        
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
                self.templateName = child.text
            elif child.tag == 'TransZTrim':
                if 'Min' in child.attr:
                    self.lightZMin = int(child.attr['Min'])
                if 'Max' in child.attr:
                    self.lightZMax = int(child.attr['Max'])
            elif child.tag == 'NucZTrim':
                if 'Min' in child.attr:
                    self.nucZMin = int(child.attr['Min'])
                if 'Max' in child.attr:
                    self.nucZMax = int(child.attr['Max'])            
            else:
                print("Tag \"", child.tag, "\" not recognized for cellseg settings")
                
    def copyDataFrom(self, other, override):
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
        self.gaussRad = int(child.attr['GaussRad'])
        for child in element:
            if child.tag == 'ThresholdSettings':
                if 'Preset' in child.attr:
                    self.thPreset = int(child.attr['Preset'])  
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
        if 'DoClouds' in element.attr:
            self.qNoClouds = not bool(element.attr['DoClouds'])
        if 'CellZero' in element.attr:
            self.qCellZero = bool(element.attr['CellZero'])
                
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
        if 'CpuCount' in element.attr:
            self.cpuCount = int(element.attr['CpuCount'])
        if 'RamGigs' in element.attr:
            self.ramGigs = int(element.attr['RamGigs'])
        if 'Time' in element.attr:
            self.timeString = element.attr['Time']
                
    def copyDataFrom(self, other, overwrite):
        if other is None:
            return
        if overwrite or (self.cpuCount <= 0):
            self.cpuCount = other.cpuCount
        if overwrite or (self.ramGigs <= 0):
            self.ramGigs = other.ramGigs
        if overwrite or (self.templateName is None):
            self.timeString = other.timeString    

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
        if 'ChannelNumber' in element.attr:
            self.channelNumber = int(element.attr['ChannelNumber'])
        for child in element:
            if child.tag == 'Meta':
                metanode = MetaSettings()
                metanode.fromXmlNode(child)
                if self.metaData is not None:
                    self.metaData.copyDataFrom(metanode, true)
                else:
                    self.metaData = metanode
            elif child.tag == 'SpotDetectSettings':
                childnode = SpotDetectSettings()
                childnode.fromXmlNode(child)
                if self.spotsSettings is not None:
                    self.spotsSettings.copyDataFrom(childnode, true)
                else:
                    self.spotsSettings = childnode    
            elif child.tag == 'QuantSettings':
                childnode = QuantSettings()
                childnode.fromXmlNode(child)
                if self.quantSettings is not None:
                    self.quantSettings.copyDataFrom(childnode, true)
                else:
                    self.quantSettings = childnode
            elif child.tag == 'SpotJob':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.spotsJobSettings is not None:
                    self.spotsJobSettings.copyDataFrom(childnode, true)
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
        
        self.batchScriptPath = None
        
        self.parentSet = None
        self.channels = list()
        
    def fromXmlNode(self, element):
        if 'Name' in element.attr:
            self.name = element.attr['ChannelNumber']
            
        for child in element:
            if (child.tag == 'CommonMeta') or (child.tag == 'Meta'):
                metanode = MetaSettings()
                metanode.fromXmlNode(child)
                if self.metaData is not None:
                    self.metaData.copyDataFrom(metanode, true)
                else:
                    self.metaData = metanode
            elif child.tag == 'CellSegSettings':
                childnode = CellsegSettings()
                childnode.fromXmlNode(child)
                if self.cellsegSettings is not None:
                    self.cellsegSettings.copyDataFrom(childnode, true)
                else:
                    self.cellsegSettings = childnode             
            elif child.tag == 'SpotDetectSettings':
                childnode = SpotDetectSettings()
                childnode.fromXmlNode(child)
                if self.spotsSettings is not None:
                    self.spotsSettings.copyDataFrom(childnode, true)
                else:
                    self.spotsSettings = childnode    
            elif child.tag == 'QuantSettings':
                childnode = QuantSettings()
                childnode.fromXmlNode(child)
                if self.quantSettings is not None:
                    self.quantSettings.copyDataFrom(childnode, true)
                else:
                    self.quantSettings = childnode
            elif child.tag == 'SpotJob':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.spotsJobSettings is not None:
                    self.spotsJobSettings.copyDataFrom(childnode, true)
                else:
                    self.spotsJobSettings = childnode  
            elif child.tag == 'CellSegJob':
                childnode = JobSettings()
                childnode.fromXmlNode(child)
                if self.cellsegJobSettings is not None:
                    self.cellsegJobSettings.copyDataFrom(childnode, true)
                else:
                    self.cellsegJobSettings = childnode
            elif child.tag == 'Paths':
                for gchild in child:
                    if gchild.tag == 'ImageDir':
                        self.inputDir = gchild.text
                    elif gchild.tag == 'OutputDir':
                        self.outputDir = gchild.text
            elif child.tag == 'ChannelInfo':
                if 'ChannelCount' in child.attr:
                    self.channelCount = int(child.attr['ChannelCount'])   
                if 'NucChannel' in child.attr:
                    self.nucChannel = int(child.attr['NucChannel']) 
                if 'TransChannel' in child.attr:
                    self.lightChannel = int(child.attr['TransChannel'])
                for gchild in child:
                    if gchild.tag == 'ImageChannel':
                        chinfo = ImageChannelInfo()
                        chinfo.copyFromParent(self, false)
                        chinfo.fromXmlNode(gchild)
                        chinfo.parentBatch = self
                        self.channels.append(chinfo)
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
            chinfo.copyFromParent(self, false)

class ImageBatchSet:
    def __init__(self):
        self.metaData = None
        self.cellsegSettings = None
        self.spotsSettings = None
        self.quantSettings = None
        self.cellsegJobSettings = None
        self.spotsJobSettings = None
        
        self.tsDir = None
        self.moduleName = None
        
        self.batches = list()
        
    def fromXmlNode(self, element):
        for child in element:
            if (child.tag == 'CommonMeta') or (child.tag == 'Meta'):
                metanode = MetaSettings()
                metanode.fromXmlNode(child)
                if self.metaData is not None:
                    self.metaData.copyDataFrom(metanode, true)
                else:
                    self.metaData = metanode
            elif child.tag == 'CellSegSettings':
                childnode = CellsegSettings()
                childnode.fromXmlNode(child)
                if self.cellsegSettings is not None:
                    self.cellsegSettings.copyDataFrom(childnode, true)
                else:
                    self.cellsegSettings = childnode             
            elif child.tag == 'SpotDetectSettings':
                childnode = SpotDetectSettings()
                childnode.fromXmlNode(child)
                if self.spotsSettings is not None:
                    self.spotsSettings.copyDataFrom(childnode, true)
                else:
                    self.spotsSettings = childnode    
            elif child.tag == 'QuantSettings':
                childnode = QuantSettings()
                childnode.fromXmlNode(child)
                if self.quantSettings is not None:
                    self.quantSettings.copyDataFrom(childnode, true)
                else:
                    self.quantSettings = childnode
            elif child.tag == 'JobSettings':
                for gchild in child:
                    if gchild.tag == 'SpotJob':
                        childnode = JobSettings()
                        childnode.fromXmlNode(gchild)
                        if self.spotsJobSettings is not None:
                            self.spotsJobSettings.copyDataFrom(childnode, true)
                        else:
                            self.spotsJobSettings = childnode  
                    elif gchild.tag == 'CellSegJob':
                        childnode = JobSettings()
                        childnode.fromXmlNode(gchild)
                        if self.cellsegJobSettings is not None:
                            self.cellsegJobSettings.copyDataFrom(childnode, true)
                        else:
                            self.cellsegJobSettings = childnode
                    elif gchild.tag == 'TrueSpotDir':
                        self.tsDir = gchild.text
                    elif gchild.tag == 'MatlabModuleName':
                        self.moduleName = gchild.text                    
            elif child.tag == 'ImageBatch':
                childnode = ImageBatchInfo()
                childnode.copyFromParent(self, false)
                childnode.fromXmlNode(child)
                childnode.parentSet = self
                self.batches.append(childnode)
            else:
                print("Tag \"", child.tag, "\" not recognized for set settings")
                
    def updateOverrides(self):
        for imgbatch in self.batches:
            imgbatch.copyFromParent(self, false)
            imgbatch.updateOverrides()
   
class SingleImage:
    def __init__(self):
        self.name = None
        self.tifPath = None
        self.resultsDir = None
        self.csScriptPath = None
        self.csResPath = None
        
        self.batchInfo = None
        
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
        scriptHandle.write(" -voxelsize \"(" + vxSz[0] + "," + vxSz[1] + "," + vxSz[2] + ")\"")
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
        if channelInfo.spotsSettings.thPreset < 0:
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
    scriptHandle.write(" -log \"" + os.path.join(chDir, chName + '_quant_mat.log') + "\"")
    scriptHandle.write("\n")
    scriptHandle.write("else\n")
    scriptHandle.write("\techo -e \"Call table not found! Can't do quant! Terminating...\"\n")
    scriptHandle.write("\texit 1\n")
    scriptHandle.write("fi\n\n")
    channelInfo.quantResPath = fullOutStem + "_quantData.mat"
    
    close(scriptHandle)
    return channelInfo

def genImageJobs(tifImage, batchInfo):
    os.makedirs(tifImage.resultsDir, exist_ok=True)
    tifImage.csResPath = os.path.join(tifImage.resultsDir, "CellSeg_" + tifImage.name + ".mat")
    tifImage.batchInfo = batchInfo
    
    scriptHandle = open(tifImage.csResPath, 'w')
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
        scriptHandle.write("\tchmod 775 \"" + chInfo.scriptPath + "\"\n")
        scriptHandle.write("\tsbatch --job-name=\"TS_" + tifImage.name + "_" + chInfo.dirName + "\"")
        if chInfo.spotsJobSettings is not None:
            scriptHandle.write(" --cpus-per-task=" + str(chInfo.spotsJobSettings.cpuCount))
            scriptHandle.write(" --mem=" + str(chInfo.spotsJobSettings.memGigs) + "g")
            scriptHandle.write(" --time=" + str(chInfo.spotsJobSettings.timeString))
        scriptHandle.write("- -error=\""  + os.path.join(chDir, tifImage.name + "_" + chInfo.dirName + '_tsSlurm.err') + "\"")
        scriptHandle.write(" --out=\""  + os.path.join(chDir, tifImage.name + "_" + chInfo.dirName + '_tsSlurm.out') + "\"")
        scriptHandle.write(" \"" + chInfo.scriptPath + "\"\n")
        scriptHandle.write("fi\n\n")

    close(scriptHandle)
    
    return tifImage
            
def genBatch(myBatch):
    print(getdtstr(), "Working on batch", myBatch.name)
    
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
        batchScriptHandle.write("chmod 775 \"" + tifImage.csScriptPath + "\"\n")
        scriptHandle.write("\tsbatch --job-name=\"TSCS_" + tifImage.name + "\"")
        if myBatch.cellsegJobSettings is not None:
            scriptHandle.write(" --cpus-per-task=" + str(myBatch.cellsegJobSettings.cpuCount))
            scriptHandle.write(" --mem=" + str(myBatch.cellsegJobSettings.memGigs) + "g")
            scriptHandle.write(" --time=" + str(myBatch.cellsegJobSettings.timeString))
        scriptHandle.write("- -error=\""  + os.path.join(tifImage.resultsDir, tifImage.name + '_tscsSlurm.err') + "\"")
        scriptHandle.write(" --out=\""  + os.path.join(tifImage.resultsDir, tifImage.name + '_tscsSlurm.out') + "\"")
        scriptHandle.write(" \"" + tifImage.csScriptPath + "\"\n\n")        

    close(batchScriptHandle)
    
    return myBatch

def readBatchXml(xmlpath):
    xmlDoc = ETree.parse(xmlpath)
    treeRoot = xmlDoc.getroot()
    
    if treeRoot.tag != 'ImageSet':
        print(getdtstr(), "Root tag must be \"ImageSet\"! Exiting...")
        return None
    
    batchSet = ImageBatchSet()
    batchSet.fromXmlNode(treeRoot)
    
    del(treeRoot)
    del(xmlDoc)
    gc.collect()
    return batchSet

def main(args):
    print("TS Batch Job Generated initiated! Version 25.04.02.00")
    print("Input Specification:", args.xmlpath)
    
    print(getdtstr(), "Reading input xml...")
    batchSet = readBatchXml(args.xmlpath)
      
    for batchInfo in batchSet:
        print(getdtstr(), "Working on", batchInfo.name, "...")
        batchInfo = genBatch(batchInfo)
        print('chmod 775 \"', batchInfo.batchScriptPath, '\"', sep='')
        print('bash \"', batchInfo.batchScriptPath, '\"', sep='', flush=True)
        
if __name__ == "__main__":
    # Args
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("xmlpath", help="Input batch specification as XML")
    parser.add_argument("--help", "-h", "-?", action="help", help="Show this help message and exit.")
    args = parser.parse_args()
    main(args)