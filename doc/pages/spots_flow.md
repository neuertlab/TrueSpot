# Spot Detection Module - Internal Flow
This page provides a brief overview of the modules within the spot calling module. This is to make it slightly easier for anyone looking into the code to find functions of interest.

[Back](../dochome.md)

![General diagram of spot detection module flow](./main_rnaspots_modules_240305.png)

## Node Summaries
Below is a short description of each node in the general flow chart shown above. Classes are generally used as namespaces rather than instantiatable object definitions, with some exceptions.

**Main_RNASpots** - A MATLAB function file that deals with parsing arguments and storing them in a data structure that the pipeline core will accept. It also takes care of some "sanity checks".

**Adapter_RNASpots** - A MATLAB function file that primarily serves to funnel inputs from different user interfaces into the spot detect pipeline. It provides a few more sanity checks that may or may not be handled by the input interface. It was created in case a GUI was ever added in addition to the command line interface.

**RNASpotsRun** - An instantiatable class for storing parameters describing a spot detection run and some of its results. The definition used to be more strict with every parameter being assigned to its own instance variable, but besides being an issue for clutter, it made it harder to store and load an object's contents in a mat file in such a way that they could be easily read without the `RNASpotsRun.m` class definition.

**RNA_Pipeline_Core** - The core wrapper function that runs the actual spot detection and automated thresholding. It takes an instance of `RNASpotsRun` as its main argument. It can also accept input images that have been preloaded by another MATLAB script, as would be useful in a GUI.

**RNA_Threshold_SpotDetector** - A class containing static wrapper functions for the spot detection process. Mostly left in for compatibility purposes. 

**RNADetection** - A class containing static functions for the initial maxima detection and threshold scanning. If you are interested in how spot calls are made, this is where to look.

**RNAThreshold** - A class containing static functions for automated thresholding, preset management, wrapping of automated thresholding, and debug plots. If you are interested in the automated thresholding process, this is where to look. `estimate_threshold()` is the function that makes the threshold suggestion. It accepts a struct containing all parameters and input data (ie. spot count tables). A default version of the input struct can be generated with `genEmptyThresholdParamStruct()`. The output is another struct that contains the threshold pool, selected threshold, intermediate data, and all parameters used.

**Seglr2** - A class containing static functions for fitting a 2-part linear spline to a curve. There is a lot of experimental stuff in there.

**RNA_Threshold_Common** - The containing class for a lot of general use static functions for spot detection and thresholding, though many of these have since been moved to more specific modules. This file still contains the code for LoG filter generation and application as well as dead pixel detection. It also contains some experimental clustering stuff and old versions of the maxima detection functions. Don't use those. They are functional, but the output is a mess.

**RNAUtils** - A class containing short static functions for common utility use.
