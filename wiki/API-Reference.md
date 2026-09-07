# API Reference

Everything below is reachable as `pyfdstools.<name>`. Each routine
carries a numpydoc docstring, so `help(fds.query2dAxisValue)` gives the
full signature and argument list.

```python
import pyfdstools as fds
help(fds.query2dAxisValue)
```

## Finding files and metadata

| Routine | Purpose |
| --- | --- |
| `getFileList(dir, chid, ext)` | list a case's files with an extension |
| `getSmvFile(dir, chid)` | the case's smokeview file, raising if absent |
| `getEndianness(dir, chid)` | byte order the output was written in |
| `getDatatypeByEndianness(dtype, endianness)` | apply that byte order to a numpy dtype |
| `zopen(path)`, `zreadlines(path)` | open or read a file that may be inside an archive |
| `parseSMVFile(smvFile)` | parse a smokeview file into grids, obstructions and file records |
| `readSLCFquantities(chid, dir)` | what slices the case wrote |
| `readBoundaryQuantities(dir, chid)` | what boundary files the case wrote |
| `getBndeQuantities(smvFile)` | what `&GEOM` boundary element files the case wrote |
| `buildQTYstring(chid, dir, qty)` | the file name suffix FDS used for a quantity |

## 2-D slices

| Routine | Purpose |
| --- | --- |
| `query2dAxisValue(dir, chid, qty, axis, value, ...)` | a plane across every mesh |
| `read2dSliceFile(file, chid, ...)` | one slice file with its coordinates |
| `readSingleSlcfFile(file, ...)` | one slice file, raw |
| `readSLCFtimes(file)` | the timestamps in a slice file |
| `readSLCFheader(f, endianness)` | quantity, units and extents from a slice header |
| `getAxisAndValueFromXB(XB, grid, cen)` | which plane a slice lies in |
| `getAxisFromLims(lims)` | the same, from cell extents alone |
| `renderSliceCsvs(data, chid, outdir)` | write each frame to csv |
| `getBlockedCellsInPlane(dir, chid, axis, value)` | mask of cells blocked by obstructions |

## 3-D slices and plot3D

| Routine | Purpose |
| --- | --- |
| `readSLCF3Ddata(chid, dir, qty, ...)` | a 3-D slice across every mesh |
| `readPlot3Ddata(chid, dir, time)` | a plot3D snapshot across every mesh |
| `findSliceLocation(grid, data, axis, value)` | take a plane out of a volume |
| `extractPoint(point, grid, data)` | the history at the nearest grid point |
| `readP3Dfile(file)`, `writeP3Dfile(file, data)` | one plot3D file |
| `readXYZfile(file)`, `writeXYZfile(file, grid)` | one grid file |
| `rearrangeGrid(grid)` | flat grid to meshgrid form |
| `getAbsoluteGrid(grids, makeUniform=False)` | build a grid spanning several meshes |
| `getLimsFromGrid(grid)` | the bounding box of an absolute grid |
| `printExtents(grid, data)` | print the coordinate and value ranges |

## Boundary data

| Routine | Purpose |
| --- | --- |
| `queryBndf(dir, chid, fdsPath, qtys, axis, value)` | a boundary plane across every mesh |
| `importBoundaryFile(file, ...)` | one boundary file as patches |
| `readBoundaryHeader(file)` | quantity, units and patch count |
| `readBoundaryFile(file)` | one boundary file, raw |
| `extractMaxBndfValues(...)` | peak value per named polygon over time |
| `buildAbsPatch(patches, ...)` | assemble patches onto one grid |
| `fdspatch` | one patch: `data`, `lims`, `orientation`, `buildSpace()`, `average()` |
| `readGcfFile(file)` | `&GEOM` vertices, faces and surfaces |
| `readBeFile(file)`, `writeBeFile(...)` | `&GEOM` boundary element data |

## Particles and smoke3D

| Routine | Purpose |
| --- | --- |
| `importParticles(dir, chid)` | particles from every mesh |
| `importParticle(file)` | one particle file |
| `importParticle_meta(file)` | particle classes without the coordinates |
| `extractS3dValues(dir, chid, decode=True)` | smoke3D fields in physical units |
| `readS3dFile(file)` | one smoke3D file, raw bytes |
| `writeS3dFile(file, times, data, qty, dx)` | write a smoke3D file |

## CSV output

| Routine | Purpose |
| --- | --- |
| `load_csv(dir, chid, suffix='_devc')` | any of the case's csv outputs |
| `load_hrr(file)` | an `_hrr.csv` by path |

## Time averaging and filtering

| Routine | Purpose |
| --- | --- |
| `timeAverage(data, times, window, ...)` | boxcar average over a moving window |
| `kalmanFilter(z, Q, R)` | Kalman filter a 1-D series |
| `slcfTimeAverage(file, dt, ...)` | average one slice file |
| `slcfsTimeAverage(dir, chid, qty, dt, ...)` | average every slice of a quantity |
| `bndfTimeAverage(file, dt, ...)` | average one boundary file |
| `bndfsTimeAverage(dir, chid, qty, dt, ...)` | average every boundary file of a quantity |

## Plotting

| Routine | Purpose |
| --- | --- |
| `plotSlice(x, z, data, axis, ...)` | render a 2-D field |
| `visualizePlot3D(x, z, T, U, V, W, HRR, ...)` | plot the temperature of a plot3D slice |
| `smvVisual(obstructions, surfaces, namespace, ...)` | render obstructions in 3-D |
| `buildSMVcolormap(percentile=None, width=None)` | the smokeview colormap |
| `getVTcolors()`, `getJHcolors()` | categorical color sequences |
| `getPlotColors(n)` | generate n distinct colors |
| `maxValuePlot(...)`, `maxValueCSV(...)` | plot or write per-group maxima |

## Writing output files

| Routine | Purpose |
| --- | --- |
| `writeSlice(...)` | write a derived slice, optionally registering it in a smokeview file |
| `writeSliceToSmv(...)` | register a slice file in a smokeview file |
| `writeSLCFheader(...)`, `writeSLCFTime(...)` | the pieces of a slice file |
| `writeBndfHeader(...)`, `writeBndfPatchInfo(...)` | the pieces of a boundary file |
| `buildBndfSmvLine(...)` | the smokeview record for a boundary file |

## Input files

| Routine | Purpose |
| --- | --- |
| `fdsFileOperations()` | read, edit, generate and write `.fds` files |
| `fdsLineTypes(version)` | datatype of every namelist parameter |
| `checkDevices(model, smvFile)` | devices buried inside obstructions |

`fdsFileOperations` methods worth knowing:

| Method | Purpose |
| --- | --- |
| `importFile(file=None, text=None, textList=None)` | read a model |
| `saveModel(mpiProcesses, location)` | write it out |
| `generateFDStext(precision=15)` | the text without writing a file |
| `addHEAD`, `addTIME`, `addMESH`, `addSURF`, `addMATL`, `addOBST`, `addVENT`, `addDEVC`, `addSLCF`, `addBNDF`, `addREAC` | add a namelist |
| `calculateMeshCells()` | mesh names and cell counts |
| `checkOverlappingMESH()` | whether any two meshes overlap |
| `sortDEVCs()` | order the devices as FDS requires |

## Model uncertainty

| Routine | Purpose |
| --- | --- |
| `readErrorTable(fdsVersion)` | the published bias and standard deviation table |
| `getQuantities(fdsVersion)` | the quantities in that table |
| `calculatePercentile(values, qty, percentile, ...)` | apply the uncertainty |
| `plotPercentile(value, qty, ...)` | plot the resulting distribution |

## Fire engineering helpers

| Routine | Purpose |
| --- | --- |
| `astFromGhf(ghf, h, e, Tgauge=20)` | adiabatic surface temperature from gauge heat flux |
| `getTwoZone(z, val)` | reduce a vertical profile to a two-layer model |
| `pointsFromXB(XB, extend)` | the eight corners of an `XB` |
| `in_hull(p, hull)` | whether points fall inside a convex hull |
| `pts2polygons(groups)` | convex hulls from groups of points |

## ParaView export

Requires `pip install "pyfdstools[paraview]"`.

| Routine | Purpose |
| --- | --- |
| `exportSl3dDataToVtk`, `exportSl2dDataToVtk` | slices |
| `exportBndfDataToVtk`, `exportBndeDataToVtk` | boundary data |
| `exportPrt5DataToVtk` | particles |
| `exportS3dDataToVtk` | smoke3D |
| `obstToStl(dir, chid, outDir)` | obstructions as an STL |

## VTKHDF

Requires `pip install "pyfdstools[vtk]"`.

| Routine | Purpose |
| --- | --- |
| `query2dAxisValue_vtkhdf(...)` | a 2-D slice from FDS VTKHDF output |

## Examples

| Routine | Purpose |
| --- | --- |
| `runExamples(scripts=None, raiseOnError=False)` | run the bundled example scripts |
| `getExamplesDirectory()` | where they live |
| `EXAMPLE_SCRIPTS` | their names |

## Deprecated

These emit a `DeprecationWarning` and will be removed.

| Deprecated | Use instead |
| --- | --- |
| `query2dAxisValueXYZ` | `query2dAxisValue` |
| `readSLCF3DdataXYZ` | `readSLCF3Ddata` |
| `readSLCF2Ddata` | `query2dAxisValue` |
| `timeAverage2` | `timeAverage` |
| `getFileListFromResultDir` | `getFileList` (it is now an alias) |
| `writeBeFile_old` | `writeBeFile` |
