# TrueSpot Batch XMLs
This page provides a comprehensive overview of the XML structure and fields recognized by TrueSpot Lite and `tsBatchGen.py`. Parameters are organized hierarchically and most options correlate directly to one or more command line interface parameters.

[Back](../dochome.md)

An example input XML for TrueSpot Lite can be found [here](../sampleXmls/TSLite_Test.xml).

Parameters are organized into blocks. The outermost block has the tag `ImageSet`. An `ImageSet` can contain any number of batches (`ImageBatch`) and blocks specifying parameters common to all batches. TrueSpot Lite will only read the first `ImageSet` in the file and it will process ALL batches within that set. Having multiple `ImageSet`s in one file is not recommended anyway since XML parsers tend to balk at files that don't have a single defined root element. 

## ImageBatch
An image batch is a group of images (specified by input table or directory) to run together. It is assumed that all images in a batch are of the same target and have the same channel settings.

The `ImageBatch` block has two mandatory child blocks: `Paths` and `ChannelInfo`. The other blocks are optional common blocks that apply to the whole set unless overridden by channel specifications.

```
<ImageBatch Name="{BATCH_NAME}">
	<CommonMeta>...</CommonMeta>
	<CellSegSettings>...</CellSegSettings>
	<SpotDetectSettings>...</SpotDetectSettings>
	<QuantSettings>...</QuantSettings>
	<Paths>...</Paths>
	<ChannelInfo>...</ChannelInfo>
</ImageBatch>
```

### ImageBatch: Paths
`Paths` has no attributes, two mandatory children (`Input` or `ImageDir` - these are synonymous, and `OutputDir`), and one optional child (`ControlPath`).

```
<Paths>
	<Input>"{INPUT FILE/DIR PATH}"</Input>
	<OutputDir>"{OUTPUT DIR PATH}"</OutputDir>
	<ControlPath>"{CONTROL IMAGE PATH}"</ControlPath>
</Paths>
```

`Input` should be the path to a directory of TIF files. These should be the TIF files that contain the sample channels to analyze. TrueSpot Lite will attempt to process all TIF files in the target directory (NON recursively). Therefore, any images that do not follow the specifications in `ChannelInfo` or other blocks should be moved to a different directory.

`Output` specifies the output directory. Individual inner directory and file names are generated from the image file names, so only a place to put this output needs to be specified. If the specified output folder does not exist at run time, TrueSpot Lite will create it.

At this time, `ControlPath` only interprets a single path to a single control image.

### ImageBatch: ChannelInfo
The `ChannelInfo` block contains information about the channels in the input images. It is also the block containing sample channel specific parameters in the form of `ImageChannel` blocks.

`ChannelInfo` has only one mandatory attribute (`ChannelCount`) and five optional ones. It does not have any children besides its `ImageChannel` blocks.

```
<ChannelInfo ChannelCount="{INT}" NucChannel="{INT}" TransChannel="{INT}" ControlChannelCount="{INT}" ControlNucChannel="{INT}" ControlTransChannel="{INT}">
	<ImageChannel>...{Sample channel A)...</ImageChannel>
	<ImageChannel>...{Sample channel B)...</ImageChannel>
	...
</ChannelInfo>
```

`ChannelCount` is the number of channels in every image stack in the input directory. The third party TIF parser we use does not read any data inside a TIF that indicates how many channels there are. Therefore in order to split the stacks constituting each channel, the total number of channels needs to be input manually.

`NucChannel` is the 1-based channel index for the channel containing nuclear marker data such as DAPI. This is only required for cell segmentation and some rendering.

`TransChannel` is the 1-based channel index for the channel containing TRANS light data. This is used for cell segmentation, background extraction (for Spot Detect threshold flooring) and some rendering.

`ControlChannelCount`, `ControlNucChannel`, and `ControlTransChannel` are the same as above, but for the control image(s) if provided.

## ImageChannel
An `ImageChannel` block contains parameters specific to a single sample channel. `ImageChannel` has no mandatory children, though it has one mandatory attribute (`ChannelNumber`, which specifies its 1-based index). An `ImageChannel` block must be present for TrueSpot Lite to run a channel as a sample, even if there are no channel specific parameters.

```
<ImageChannel ChannelNumber="{INT}">
	<Meta>...</Meta>
	<CellSegSettings>...</CellSegSettings>
	<SpotDetectSettings>...</SpotDetectSettings>
	<QuantSettings>...</QuantSettings>
</ImageChannel>
```

Note that at the channel level, the metadata block is just called `Meta`.

## Common Blocks
These can be included at the set, batch, or channel level. Any blocks appearing above the channel level will be assumed to be common, overridden by any parameter setting specified down the branch.

### Metadata
The metadata block (with tag `CommonMeta` at the set or batch level, and `Meta` at the channel level) can be used to specify metadata tags to save with the run data. These can be helpful for record keeping and results output. The metadata block has no top level attributes and seven optional children.

```
<(Common)Meta>
	<Species>"{SPECIES NAME}"</Species>
	<CellType>"{CELL TYPE/LINE}"</CellType>
	<TargetName>"{SAMPLE TARGET NAME}"</TargetName>
	<ProbeName>"{PROBE NAME}"</ProbeName>
	<TargetType>"{SAMPLE TARGET MOLECULE TYPE}"</TargetType>
	<VoxelDimsNano X="{FLOAT}" Y="{FLOAT}" Z="{FLOAT}"/>
	<PointDimsNano X="{FLOAT}" Y="{FLOAT}" Z="{FLOAT}"/>
</(Common)Meta>
```

`PixelDimsNano` can be substituted for `VoxelDimsNano` with only `X` and `Y` attributes in the case of 2D images. `VoxelDimsNano`/`PixelDimsNano` and `PointDimsNano` are specified in nanometers. These values are not used by TrueSpot for processing, but are handy for scaling, rendering, and recording keeping.

### CellSegSettings
The `CellSegSettings` block specifies parameters for cell segmentation, if it is to be run. The [parameters](./cellseg_allargs.md) are equivalent to what is used for the command line interface. The full arguments page can be referenced for descriptions. This whole block is optional, and all children are optional for the block.

```
<CellSegSettings>
	<PresetName>"{PRESET NAME}"</PresetName>
	<TransZTrim Min="{INT}" Max="{INT}"/>
	<NucZTrim Min="{INT}" Max="{INT}"/>
	<CellSize Min="{INT}" Max="{INT}"/>
	<NucSize Min="{INT}" Max="{INT}"/>
	<XTrim>"{INT}"</XTrim>
	<YTrim>"{INT}"</YTrim>
	<NucZRange>"{INT}"</NucZRange>
	<NucThSample>"{INT}"</NucThSample>
	<NucCutoff>"{FLOAT}"</NucCutoff>
	<NucDXY>"{FLOAT}"</NucDXY>
	<Options ExportCellMaskToFormat="{png | tif}" ExportNucMaskToFormat="{png | tif}" Overwrite="{BOOL}" DumpSettingsToText="{BOOL}"/>
</CellSegSettings>
```

### SpotDetectSettings
`SpotDetectSettings` contains parameters for spot detection and automated threshold selection. As with `CellSegSettings`, the block is not required and it does not have any required children or attributes. However, we recommend using `GaussRad` to set a Gaussian filter radius smaller than 7 for images with voxel sizes significantly larger than 65nm. We also recommend using `Workers` to attempt multithreading to speed up the threshold scan.

See the [full command line argument list](./spots_allargs.md) for descriptions of parameters.

```
<SpotDetectSettings GaussRad="{FLOAT}" Workers="{INT}" NoDPC="{BOOL}">
	<ZTrim Min="{INT}" Max="{INT}"/>
	<Options DumpRunParamsToText="{BOOL}" Overwrite="{BOOL}" CsvDump="{All | ThresholdRange | SelectedThreshold}" CsvZeroBased="{BOOL}"/>
	<ThresholdSettings Preset="{INT}" ScanMin="{INT}" ScanMax="{INT}">
		<WindowSettings Min="{FLOAT}" Max="{FLOAT}" Increment="{FLOAT}"/>
		<MADFactor Min="{FLOAT}" Max="{FLOAT}"/>
		<Weights FitRightIntersect="{FLOAT}" MedMad="{FLOAT}" Fit="{FLOAT}"/>
		<MiscOptions IncludeRawCurve="{BOOL}" IncludeDiffCurve="{BOOL}" StDevFactor="{FLOAT}" LogMode="{All | None | FitOnly}"/>
	<ThresholdSettings/>
</SpotDetectSettings>
```

The `Preset` attribute of `ThresholdSettings` takes an integer from -5 to +5 specifying a preset level. 0 is default. Positive numbers represent increased sensitivity and negative numbers represent increased precision.

### QuantSettings
A `QuantSettings` block specifies quantification/fit parameters. It has no children, instead all parameters are specified as attributes.

See the [full command line argument list](./quant_allargs.md) for descriptions of parameters.

```
<QuantSettings DoClouds="{BOOL}" DoRefilter="{BOOL}" CellZero="{BOOL}" ZtoXY="{FLOAT}" FitRadXY="{INT}" FitRadZ="{INT}" ManualThreshold="{INT}"/>
```

Note that if the voxel size is specified in metadata, the Z to XY ratio can be automatically calculated.
The `ManualThreshold` option is intended as an override. We recommend letting quant run without a threshold specified in most cases however as this allows it to fit spots over a wider range of possible thresholds meaning quant does not have to be run again should the desired threshold be lowered.
