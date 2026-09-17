# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:31:10 2026

@author: bhospelh
"""

import argparse
import datetime
import numpy
import gc
import os.path

from cellpose import models
from cellpose.io import imread
#from cellpose.io import imread_3D

from cpw_settings import CellposeRun

def parseDimArg(argstring):
    argstring = argstring.replace(" ","")
    argstring = argstring.replace("(","")
    argstring = argstring.replace(")","")
    splitstr = argstring.split(",")
    dimi = []
    for i in range(len(splitstr)):
        dimi.append(int(splitstr[i]))
    return tuple(dimi)

def getdtstr():
    now = datetime.datetime.now()
    return "[" + str(now) + "]"

def getStackChannel(stack, chdim, chidx):
    #Copied from the Big-FISH wrapper - I still hate it
    if chdim == 0:
        return stack[chidx,:,:,:]
    elif chdim == 1:
        return stack[:,chidx,:,:]
    elif chdim == 2:
        return stack[:,:,chidx,:]
    elif chdim == 3:
        return stack[:,:,:,chidx]
    else:
        return None
    
def outputMask(maskDat, outpath):
    #If 3D, then just adds planes as more rows. Thus Z is needed to read correctly.
    if len(maskDat.shape) > 2:
        Z = maskDat.shape[0]
        Y = maskDat.shape[1]
        X = maskDat.shape[2]
        Yamal = Z * Y
        #maskDat = numpy.swapaxes(maskDat, 0, 2) #Move z to last dim
        #maskDat = numpy.swapaxes(maskDat, 0, 1)
        maskDat = numpy.reshape(maskDat, (Yamal, X))
    numpy.savetxt(outpath, maskDat, delimiter=",")
    
def readImage3D(path):
    imageRaw = imread(path)
    #Determine c axis
    smallestAxis = -1
    smallestAxisSize = 0x7fffffff
    axisCount = len(imageRaw.shape)
    for i in range(axisCount):
        if imageRaw.shape[i] < smallestAxisSize:
            smallestAxisSize = imageRaw.shape[i]
            smallestAxis = i
            
    if smallestAxis != 3:
        imageRaw = numpy.moveaxis(imageRaw, smallestAxis, 3)
        
    return imageRaw

def runCellpose(runparams):
    print(getdtstr(), "Cellpose 4 Wrapper Initialized! Version 26.09.16.00")
    
    #Print parameters for recordkeeping
    print(getdtstr(),"Input File:", runparams.imgSettings.filePath)
    if runparams.imgSettings.nucCh >= 0:
        print(getdtstr(),"Nuclear marker channel:", runparams.imgSettings.nucCh)
    else:
        runparams.runNucSeg = False
        
    if runparams.imgSettings.cellCh >= 0:
        print(getdtstr(),"Cell boundary marker channel:", runparams.imgSettings.cellCh)
    else:
        runparams.runCellSeg = False
        
    print(getdtstr(),"Run nuclear segmentation:", runparams.runNucSeg)
    print(getdtstr(),"Run cell segmentation:", runparams.runCellSeg)
    print(getdtstr(),"Output overwrite block:", runparams.noOverwrite)
    #print(getdtstr(),"Export mask PNGs:", runparams.exportPng)
    
    if runparams.imgSettings.xMin >= 0:
        print(getdtstr(),"X range:", runparams.imgSettings.xMin, "-", runparams.imgSettings.xMax)
    if runparams.imgSettings.yMin >= 0:
        print(getdtstr(),"Y range:", runparams.imgSettings.yMin, "-", runparams.imgSettings.yMax)
    
    if runparams.runNucSeg:
        print(getdtstr(),"Nuclear segmentation output path:", runparams.nucOutputPath)
        print(getdtstr(),"Nuclear segmentation settings ---")
        print("\tAverage diameter:", runparams.nucSettings.avgDia)
        print("\tMinimum size:", runparams.nucSettings.minSize)
        print("\tNormalize:", runparams.nucSettings.normalize)
        print("\tModel name:", runparams.nucSettings.modelName)
        #print("\tUse ensemble model:", runparams.nucSettings.useEnsemble)
        print("\tCell threshold:", runparams.nucSettings.cellThreshold)
        print("\tFlow threshold:", runparams.nucSettings.flowThreshold)
        print("\tRun 3D:", runparams.nucSettings.do3d)
        
    if runparams.runCellSeg:
        print(getdtstr(),"Cell segmentation output path:", runparams.cellOutputPath)
        print(getdtstr(),"Cell segmentation settings ---")
        print("\tAverage diameter:", runparams.cellSettings.avgDia)
        print("\tMinimum size:", runparams.cellSettings.minSize)
        print("\tNormalize:", runparams.cellSettings.normalize)
        print("\tModel name:", runparams.cellSettings.modelName)
        #print("\tUse ensemble model:", runparams.cellSettings.useEnsemble)
        print("\tCell threshold:", runparams.cellSettings.cellThreshold)
        print("\tFlow threshold:", runparams.cellSettings.flowThreshold)
        print("\tRun 3D:", runparams.cellSettings.do3d)
    
    #Load and pre-process image (print z-related parameters, if applicable)
    
    if runparams.imgSettings.is3d:
        print(getdtstr(),"Z/X anisotropic ratio:", runparams.imgSettings.zxRatio)
        if runparams.imgSettings.nucCh >= 0:
            print(getdtstr(),"Nuclear channel Z range:", runparams.nucSettings.zMin, "-", runparams.nucSettings.zMax)
        if runparams.imgSettings.cellCh >= 0:
            print(getdtstr(),"Cell channel Z range:", runparams.cellSettings.zMin, "-", runparams.cellSettings.zMax)
        
        if (runparams.imgSettings.zaxis < 0) and (runparams.imgSettings.caxis < 0):
            imageRaw = readImage3D(runparams.imgSettings.filePath)
            #imageRaw = imread_3D(runparams.imgSettings.filePath)
            runparams.imgSettings.caxis = len(imageRaw.shape) - 1
        else:
            imageRaw = imread(runparams.imgSettings.filePath, z_axis=runparams.imgSettings.zaxis, channel_axis=runparams.imgSettings.caxis)
            
        C = imageRaw.shape[runparams.imgSettings.caxis]
        Z = imageRaw.shape[0]
        Y = imageRaw.shape[1]
        X = imageRaw.shape[2]
        print(getdtstr(), C, "channel", X, "x", Y, "x", Z, "stack loaded!")
    else:
        if runparams.imgSettings.caxis < 0:
            imageRaw = imread(runparams.imgSettings.filePath)
            runparams.imgSettings.caxis = len(imageRaw.shape) - 1
        else:
            imageRaw = imread(runparams.imgSettings.filePath, channel_axis=runparams.imgSettings.caxis)
            
        C = imageRaw.shape[runparams.imgSettings.caxis]
        Z = 1
        Y = imageRaw.shape[0]
        X = imageRaw.shape[1]
        print(getdtstr(), C, "channel", X, "x", Y, "image loaded!")
        
    print(getdtstr(), "Trimming image...")
    if runparams.imgSettings.xMin < 0:
        runparams.imgSettings.xMin = 0
    if runparams.imgSettings.yMin < 0:
        runparams.imgSettings.yMin = 0
    if runparams.imgSettings.xMax <= 0:
        runparams.imgSettings.xMax = X - 1
    if runparams.imgSettings.xMax < runparams.imgSettings.xMin:
        runparams.imgSettings.xMax = runparams.imgSettings.xMin 
    if runparams.imgSettings.yMax <= 0:
        runparams.imgSettings.yMax = Y - 1
    if runparams.imgSettings.yMax < runparams.imgSettings.yMin:
        runparams.imgSettings.yMax = runparams.imgSettings.yMin
        
    x0 = runparams.imgSettings.xMin
    x1 = runparams.imgSettings.xMax + 1
    y0 = runparams.imgSettings.yMin
    y1 = runparams.imgSettings.yMax + 1
    
    ZCT = 0
    ZNT = 0
        
    celldat_3 = None
    nucdat_3 = None
    celldat = None
    nucdat = None
    if runparams.imgSettings.is3d:
        if runparams.nucSettings.zMin < 0:
            runparams.nucSettings.zMin = 0
        if runparams.nucSettings.zMax <= 0:
            runparams.nucSettings.zMax = Z - 1
        if runparams.nucSettings.zMax < runparams.nucSettings.zMin:
            runparams.nucSettings.zMax = runparams.nucSettings.zMin
            
        zn0 = runparams.nucSettings.zMin
        zn1 = runparams.nucSettings.zMax + 1
            
        if runparams.cellSettings.zMin < 0:
            runparams.cellSettings.zMin = 0
        if runparams.cellSettings.zMax <= 0:
            runparams.cellSettings.zMax = Z - 1
        if runparams.cellSettings.zMax < runparams.cellSettings.zMin:
            runparams.cellSettings.zMax = runparams.cellSettings.zMin
            
        zc0 = runparams.cellSettings.zMin
        zc1 = runparams.cellSettings.zMax + 1
            
        if runparams.imgSettings.cellCh >= 0:
            celldat_3 = imageRaw[zc0:zc1, y0:y1, x0:x1, runparams.imgSettings.cellCh]
            if len(celldat_3.shape) > 2:
                celldat = numpy.nanmax(celldat_3, axis=0)
            else:
                celldat = celldat_3
                celldat_3 = None
            celldat = numpy.expand_dims(celldat, axis=0)
            ZCT = celldat_3.shape[0]
            
        if runparams.imgSettings.nucCh >= 0:
            nucdat_3 = imageRaw[zn0:zn1, y0:y1, x0:x1, runparams.imgSettings.nucCh]
            nucdat = numpy.nanmax(nucdat_3, axis=0)
            if len(nucdat_3.shape) > 2:
                nucdat = numpy.nanmax(nucdat_3, axis=0)
            else:
                nucdat = nucdat_3
                nucdat_3 = None
            nucdat = numpy.expand_dims(nucdat, axis=0)
            ZNT = nucdat_3.shape[0]
    else:
        if runparams.imgSettings.cellCh >= 0:
            celldat = imageRaw[y0:y1, x0:x1, runparams.imgSettings.cellCh]
            celldat = numpy.expand_dims(celldat, axis=0)
            ZCT = 1
            
        if runparams.imgSettings.nucCh >= 0:
            nucdat = imageRaw[y0:y1, x0:x1, runparams.imgSettings.nucCh]
            nucdat = numpy.expand_dims(nucdat, axis=0)
            ZNT = 1
            
    XT = (x1 - x0)
    YT = (y1 - y0)

    del(imageRaw)
    gc.collect()

    #Run nuclear segmentation, if requested
    nucSegResult = None
    if runparams.runNucSeg:
        print(getdtstr(), "Now running nuclear segmentation...")
        nmodel = models.CellposeModel(gpu=True, pretrained_model=runparams.nucSettings.modelName, diam_mean=runparams.nucSettings.avgDia)
        if runparams.nucSettings.do3d and (nucdat_3 is not None):
            print(getdtstr(), "Attempting 3D mode")
            dummyChannel = numpy.zeros((ZNT, YT, XT))
            runStack = numpy.stack([nucdat_3, dummyChannel, dummyChannel])
            runStack = numpy.swapaxes(runStack, 0, 1)
            del(dummyChannel)
            nucSegResult = nmodel.eval(runStack, normalize=runparams.nucSettings.normalize, flow_threshold=runparams.nucSettings.flowThreshold, cellprob_threshold=runparams.nucSettings.cellThreshold, min_size=runparams.nucSettings.minSize, anisotropy=runparams.imgSettings.zxRatio, do_3D=True, z_axis=0, channel_axis=1)
            del(nmodel)
            del(runStack)
        else:
            dummyChannel = numpy.zeros((YT, XT))
            ndatFlat = nucdat[0];
            runStack = numpy.stack([ndatFlat, dummyChannel, dummyChannel])
            #runStack = numpy.swapaxes(runStack, 0, 1)
            del(dummyChannel)
            del(ndatFlat)
            nucSegResult = nmodel.eval(runStack, normalize=runparams.nucSettings.normalize, flow_threshold=runparams.nucSettings.flowThreshold, cellprob_threshold=runparams.nucSettings.cellThreshold, min_size=runparams.nucSettings.minSize, channel_axis=0)
            del(nmodel)
            del(runStack)
    gc.collect()
    
    #Run cell segmentation, if requested
    cellSegResult = None
    if runparams.runCellSeg:
        print(getdtstr(), "Now running cell segmentation...")
        cmodel = models.CellposeModel(gpu=True, pretrained_model=runparams.cellSettings.modelName, diam_mean=runparams.cellSettings.avgDia)
        if runparams.nucSettings.do3d and (nucdat_3 is not None):
            print(getdtstr(), "Attempting 3D mode")
            dummyChannel = numpy.zeros((ZCT, YT, XT))
            runStack = numpy.stack([celldat_3, dummyChannel, dummyChannel])
            runStack = numpy.swapaxes(runStack, 0, 1)
            del(dummyChannel)
            cellSegResult = cmodel.eval(runStack, normalize=runparams.cellSettings.normalize, flow_threshold=runparams.cellSettings.flowThreshold, cellprob_threshold=runparams.cellSettings.cellThreshold, min_size=runparams.cellSettings.minSize, anisotropy=runparams.imgSettings.zxRatio, do_3D=True, z_axis=0, channel_axis=1)
            del(cmodel)
            del(runStack)
        else:
            dummyChannel = numpy.zeros((YT, XT))
            if nucdat is not None:
                runStack = numpy.stack([celldat[0], nucdat[0], dummyChannel])
            else:
                runStack = numpy.stack([celldat[0], dummyChannel, dummyChannel])    
            #runStack = numpy.swapaxes(runStack, 0, 1)
            del(dummyChannel)
            cellSegResult = cmodel.eval(runStack, normalize=runparams.cellSettings.normalize, flow_threshold=runparams.cellSettings.flowThreshold, cellprob_threshold=runparams.cellSettings.cellThreshold, min_size=runparams.cellSettings.minSize, channel_axis=0)
            del(cmodel)
            del(runStack)
    gc.collect()
            
    #Output results
    print(getdtstr(), "Outputting results...")
    if nucSegResult is not None:
        maskRaw = nucSegResult[0]
        #Expand to original input size
        if len(maskRaw.shape) > 2:
            nucMask = numpy.zeros((Z, Y, X))
            nucMask[zn0:zn1, y0:y1, x0:x1] = maskRaw[:,:,:]
        else:
            nucMask = numpy.zeros((Y, X))
            nucMask[y0:y1, x0:x1] = maskRaw[:,:]
        outputMask(nucMask, runparams.nucOutputPath)
        
    if cellSegResult is not None:
        maskRaw = cellSegResult[0]
        #Expand to original input size
        if len(maskRaw.shape) > 2:
            cellMask = numpy.zeros((Z, Y, X))
            cellMask[zn0:zn1, y0:y1, x0:x1] = maskRaw[:,:,:]
        else:
            cellMask = numpy.zeros((Y, X))
            cellMask[y0:y1, x0:x1] = maskRaw[:,:]
        outputMask(cellMask, runparams.cellOutputPath)
        
def main(args):
    runparams = CellposeRun()
    runparams.imgSettings.filePath = args.inpath
    
    if args.nuc_out:
        runparams.nucOutputPath = args.nuc_out
    if args.cell_out:
        runparams.cellOutputPath = args.cell_out
    if args.ch_nuc:
        runparams.imgSettings.nucCh = args.ch_nuc
    if args.ch_cell:
        runparams.imgSettings.cellCh = args.ch_cell
    if args.voxelsz:
        runparams.imgSettings.voxelSize = parseDimArg(args.voxelsz)
        runparams.imgSettings.zxRatio = runparams.imgSettings.voxelSize[0] / runparams.imgSettings.voxelSize[2]
        runparams.imgSettings.is3d = True
    if args.pixelsz:
        runparams.imgSettings.voxelSize = parseDimArg(args.pixelsz)
        runparams.imgSettings.is3d = False
    if args.in3d:
        runparams.imgSettings.is3d = True
    if args.in2d:
        runparams.imgSettings.is3d = False
    if args.nonucseg:
        runparams.runNucSeg = False
    if args.nocellseg:
        runparams.runCellSeg = False
    if args.onecoords:
        runparams.inputOneBased = True
    if args.caxis:
        runparams.imgSettings.caxis = args.caxis
    if args.zaxis:
        runparams.imgSettings.zaxis = args.zaxis
    if args.xmin:
        runparams.imgSettings.xMin = args.xmin
    if args.xmax:
        runparams.imgSettings.xMax = args.xmax
    if args.ymin:
        runparams.imgSettings.yMin = args.ymin
    if args.ymax:
        runparams.imgSettings.yMax = args.ymax
    if args.navgdia:
        runparams.nucSettings.avgDia = args.navgdia
    if args.cavgdia:
        runparams.cellSettings.avgDia = args.cavgdia
    if args.nminsz:
        runparams.nucSettings.minSize = args.nminsz
    if args.cminsz:
        runparams.cellSettings.minSize = args.cminsz
    if args.nnorm:
        runparams.nucSettings.normalize = True
    if args.cnorm:
        runparams.cellSettings.normalize = True
    if args.norm:
        runparams.nucSettings.normalize = True
        runparams.cellSettings.normalize = True
    if args.nmodel:
        runparams.nucSettings.modelName = args.nmodel
    if args.cmodel:
        runparams.cellSettings.modelName = args.cmodel
    if args.ncth:
        runparams.nucSettings.cellThreshold = args.ncth
    if args.ccth:
        runparams.cellSettings.cellThreshold = args.ccth
    if args.nfth:
        runparams.nucSettings.flowThreshold = args.nfth
    if args.cfth:
        runparams.cellSettings.flowThreshold = args.cfth
    if args.n3d:
        runparams.nucSettings.do3d = True
    if args.c3d:
        runparams.cellSettings.do3d = True
    if args.nzmin:
        runparams.nucSettings.zMin = args.nzmin
    if args.nzmax:
        runparams.nucSettings.zMax = args.nzmax
    if args.czmin:
        runparams.cellSettings.zMin = args.czmin
    if args.czmax:
        runparams.cellSettings.zMax = args.czmax
        
        
    if not runparams.imgSettings.filePath:
        print("ERROR: Input path is required!")
        return
    
    if not os.path.isfile(runparams.imgSettings.filePath):
        print("ERROR: Input \"", runparams.imgSettings.filePath, "\" does not exist!")
        return
    
    #Adjust for one-based coordinates
    if runparams.inputOneBased:
        if runparams.imgSettings.nucCh > 0:
            runparams.imgSettings.nucCh -= 1
        if runparams.imgSettings.cellCh > 0:
            runparams.imgSettings.cellCh -= 1
        if runparams.imgSettings.caxis > 0:
            runparams.imgSettings.caxis -= 1
        if runparams.imgSettings.zaxis > 0:
            runparams.imgSettings.zaxis -= 1
        if runparams.imgSettings.xMin > 0:
            runparams.imgSettings.xMin -= 1
        if runparams.imgSettings.xMax > 0:
            runparams.imgSettings.xMax -= 1
        if runparams.imgSettings.yMin > 0:
            runparams.imgSettings.yMin -= 1
        if runparams.imgSettings.yMax > 0:
            runparams.imgSettings.yMax -= 1
        if runparams.nucSettings.zMin > 0:
            runparams.nucSettings.zMin -= 1
        if runparams.nucSettings.zMax > 0:
            runparams.nucSettings.zMax -= 1
        if runparams.cellSettings.zMin > 0:
            runparams.cellSettings.zMin -= 1
        if runparams.cellSettings.zMax > 0:
            runparams.cellSettings.zMax -= 1
            
    #Update any end coordiates to be exclusive
    # if runparams.imgSettings.xMax >= 0:
    #     runparams.imgSettings.xMax += 1
    # if runparams.imgSettings.yMax >= 0:
    #     runparams.imgSettings.yMax += 1
    # if runparams.nucSettings.zMax >= 0:
    #     runparams.nucSettings.zMax += 1
    # if runparams.cellSettings.zMax >= 0:
    #     runparams.cellSettings.zMax += 1
    
    #Generate output paths from input if not provided
    indir = os.path.dirname(runparams.imgSettings.filePath)
    inname = os.path.basename(runparams.imgSettings.filePath)
    lastdot = inname.rfind('.')
    if lastdot >= 0:
        inname = inname[:lastdot]
    if runparams.runNucSeg and not runparams.nucOutputPath:
        runparams.nucOutputPath = os.path.join(indir, inname + "_nucmask.csv")
        print("WARNING: Nuclear segmentation requested, but output path was not provided.")
        print("\tSet to:", runparams.nucOutputPath)
        
    if runparams.runCellSeg and not runparams.cellOutputPath:
        runparams.cellOutputPath = os.path.join(indir, inname + "_cellmask.csv")
        print("WARNING: Cell segmentation requested, but output path was not provided.")
        print("\tSet to:", runparams.cellOutputPath)
        
    
    #Attempt to determine if input is 2D or 3D without loading
    if runparams.nucSettings.zMax > 1:
        runparams.imgSettings.is3d = True
    if runparams.cellSettings.zMax > 1:
        runparams.imgSettings.is3d = True
        
    runCellpose(runparams)
        
    
if __name__ == "__main__":
#if True:
    # Args
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("inpath", help="Input tif image path.")
    parser.add_argument("--nuc_out", help="Output path for nuclear segmentation mask (csv)")
    parser.add_argument("--cell_out", help="Output path for cell segmentation mask (csv)")
    parser.add_argument("--ch_nuc", type=int, help="Index of nuclear marker channel.")
    parser.add_argument("--ch_cell", type=int, help="Index of cell marker channel.")
    parser.add_argument("--voxelsz", help="Voxel size (for 3D stack) in nm formatted 'z,y,x' (example: 300,65,65)")
    parser.add_argument("--pixelsz", help="Pixel size (for 2D image) in nm formatted 'y,x' (example: 65,65)")
    parser.add_argument("--in3d", action="store_true", help="Input is 3D image stack (implied by use of --voxelsz)")
    parser.add_argument("--in2d", action="store_true", help="Input is 2D image (implied by use of --pizelsz)")
    parser.add_argument("--nonucseg", action="store_true", help="Skip nuclear segmentation")
    parser.add_argument("--nocellseg", action="store_true", help="Skip cell segmentation")
    parser.add_argument("--onecoords", action="store_true", help="All input coordinates/indices (axes, image positions) are one-based (Default is zero-based)")
    parser.add_argument("--caxis", type=int, help="Channel axis for multi-dim input image/stack")
    parser.add_argument("--zaxis", type=int, help="Z axis for multi-dim input image/stack")
    parser.add_argument("--xmin", type=int, help="Starting x coordinate for input region to evaluate")
    parser.add_argument("--xmax", type=int, help="Ending x coordinate (inclusive) for input region to evaluate")
    parser.add_argument("--ymin", type=int, help="Starting y coordinate for input region to evaluate")
    parser.add_argument("--ymax", type=int, help="Ending y coordinate (inclusive) for input region to evaluate")
    parser.add_argument("--navgdia", type=float, help="Average nuclear diameter (pixels)")
    parser.add_argument("--cavgdia", type=float, help="Average cell diameter (pixels)")
    parser.add_argument("--nminsz", type=float, help="Minimum nucleus area (pixels)")
    parser.add_argument("--cminsz", type=float, help="Minimum cell area (pixels)")
    parser.add_argument("--nnorm", action="store_true", help="Normalize input for nuclear segmentation")
    parser.add_argument("--cnorm", action="store_true", help="Normalize input for cell segmentation")
    parser.add_argument("--norm", action="store_true", help="Normalize input for cell and nuclear segmentation")
    parser.add_argument("--nmodel", help="Name/Path of nuclear model to use (Default: cpsam_v2)")
    parser.add_argument("--cmodel", help="Name/Path of cell model to use (Default: cpsam_v2)")
    parser.add_argument("--ncth", type=float, help="Set cell probability threshold for nuclear segmentation (Default: 0)")
    parser.add_argument("--ccth", type=float, help="Set cell probability threshold for cell segmentation (Default: 0)")
    parser.add_argument("--nfth", type=float, help="Set flow threshold for nuclear segmentation (Default: 0.4)")
    parser.add_argument("--cfth", type=float, help="Set flow threshold for cell segmentation (Default: 0.4)")
    parser.add_argument("--n3d", action="store_true", help="Do nuclear segmentation in 3D")
    parser.add_argument("--c3d", action="store_true", help="Do cell segmentation in 3D")
    parser.add_argument("--nzmin", type=int, help="Starting z slice (inclusive) for nuclear segmentation")
    parser.add_argument("--nzmax", type=int, help="Ending z slice (inclusive) for nuclear segmentation")
    parser.add_argument("--czmin", type=int, help="Starting z slice (inclusive) for cell segmentation")
    parser.add_argument("--czmax", type=int, help="Ending z slice (inclusive) for cell segmentation")
    parser.add_argument("--help", "-h", "-?", action="help", help="Show this help message and exit.")
    args = parser.parse_args()
    main(args)
