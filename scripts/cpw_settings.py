# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:14:14 2026

@author: bhospelh
"""

class SegSettings:
    def __init__(self):
        self.avgDia = 0.0
        self.minSize = 0.0
        self.normalize = False
        self.modelName = None
        self.useEnsemble = False
        self.cellThreshold = 0.0
        self.flowThreshold = 0.4
        self.do3d = False
        self.zMin = -1
        self.zMax = -1
           
class ImageSettings:
    def __init__(self):
        self.filePath = None
        self.zaxis = -1
        self.caxis = -1
        self.nucCh = -1
        self.cellCh = -1
        self.zxRatio = 1.0
        self.is3d = True
        self.xMin = -1
        self.xMax = -1
        self.yMin = -1
        self.yMax = -1
        self.voxelSize = None
        
class CellposeRun:
    def __init__(self):
        self.nucOutputPath = None
        self.cellOutputPath = None
        self.noOverwrite = False
        self.runNucSeg = True
        self.runCellSeg = True
        self.inputOneBased = False
        self.imgSettings = ImageSettings()
        self.nucSettings = SegSettings()
        self.cellSettings = SegSettings()
        
        self.nucSettings.modelName = "cpsam_v2"
        self.cellSettings.modelName = "cpsam_v2"