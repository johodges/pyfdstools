# Boundary Data

`&BNDF` output records a quantity on the faces of every solid surface,
written to one `.bf` file per mesh per quantity. `&GEOM` surfaces write
`.be` boundary element files instead; those are covered at the end.

## What did the case write?

```python
import os
import pyfdstools as fds

dataDir = os.path.join(os.path.dirname(fds.__file__), 'examples', 'data')
workingDir = os.path.join(dataDir, 'case001.zip')
chid = 'case001'

quantities = fds.readBoundaryQuantities(workingDir, chid)
print(list(quantities.keys()))   # ['WALL TEMPERATURE']
```

The values are the `.bf` files holding each quantity.

## Reading a boundary plane

`queryBndf` reads every patch facing a chosen direction at a chosen
coordinate and assembles them onto one grid.

```python
fdsFilePath = os.path.join(dataDir, 'case001.fds')

datas, times = fds.queryBndf(
    workingDir, chid, fdsFilePath,
    ['WALL TEMPERATURE'],   # quantities, as a list
    -2,                     # axis: faces whose normal points in -y
    4.4)                    # coordinate of the surface

entry = datas['WALL TEMPERATURE']
print(entry['DATA'].shape)   # (13, 21, 600) -> (N, M, NT)
print(entry['UNITS'])        # C
```

| Key | Contents |
| --- | --- |
| `X` | first in-plane coordinate, `array(N, M)` |
| `Z` | second in-plane coordinate, `array(N, M)` |
| `DATA` | values, `array(N, M, NT)` |
| `UNITS` | units of the quantity |

`times` is shared across the quantities requested.

### The signed axis

Unlike slice queries, boundary queries take a **signed** axis. A surface
has two sides and they carry different data — the wall of a compartment
has one face inside and one outside.

| `axis` | Reads |
| --- | --- |
| `1` | faces whose outward normal points in +x |
| `-1` | faces whose outward normal points in −x |
| `2`, `-2` | the same for y |
| `3`, `-3` | the same for z |

Passing an unsigned axis where the geometry needs a signed one gives you
the faces on the wrong side, which usually shows up as an unexpectedly
cold field.

Passing an unsupported value raises:

```
ValueError: axis must be one of -3, -2, -1, 1, 2, or 3; received 7
```

Note that `queryBndf` takes the quantities as a **list** but the axis and
value as **scalars**.

### Why the input file path?

`queryBndf` reads the `.fds` file to recover the obstruction names and
extents, which it uses to work out which patches belong to which
surface. Pass the path to the input file itself, not to its directory.

## Plotting

The returned arrays go straight into `plotSlice`:

```python
import numpy as np
import matplotlib.pyplot as plt

entry = datas['WALL TEMPERATURE']
frame = entry['DATA'][:, :, -1]

fig, ax = fds.plotSlice(
    entry['X'], entry['Z'], frame, axis=-2,
    clabel='Wall temperature (%s)' % (entry['UNITS']),
    qnty_mn=20, qnty_mx=1000,
    cbarticks=[20, 200, 400, 600, 800, 1000])
fig.savefig('wall_temperature.png')
```

Faces not covered by any patch are `NaN`, so reduce with `np.nanmax`
rather than `np.max`.

## Peak values over time

`extractMaxBndfValues` returns the peak of a boundary quantity over the
whole surface at each output time — the usual input to a "did anything
get hot enough to ignite" check.

```python
times, maxValues, names = fds.extractMaxBndfValues(
    fdsFilePath, smvFilePath, workingDir, chid,
    ['WALL TEMPERATURE'],
    tStart=0, tEnd=120, tInt=1, tBand=3)
```

`pyfdstools/examples/extract_boundary_max.py` is a worked example that
also plots the result and writes it to csv.

## Reading one boundary file

For the raw per-mesh patches:

```python
bndfFile = fds.getFileList(workingDir, chid, 'bf')[0]
smvFile = fds.getSmvFile(workingDir, chid)
smvData = fds.parseSMVFile(smvFile)

times, patches, units = fds.importBoundaryFile(
    bndfFile, smvFile=smvFile, gridNum=0, smvData=smvData)

for patch in patches:
    print(patch.orientation, patch.lims, patch.data.shape)
```

Each patch is an `fdspatch`:

| Attribute or method | Contents |
| --- | --- |
| `data` | `array(NX, NY, NT)` of values |
| `lims` | `[xmin, xmax, ymin, ymax, zmin, zmax]` |
| `orientation` | signed axis the patch faces |
| `buildSpace()` | fills `x`, `y`, `z` node coordinate grids |
| `average(inds)` | mean over the given time indices |
| `extractPoints()` | flattens to coordinates, values and orientations |

To read just the header without the data:

```python
quantity, shortName, units, npatch = fds.readBoundaryHeader(bndfFile)
```

## Assembling patches yourself

`buildAbsPatch` is what `queryBndf` uses internally, and is available if
you have selected patches by some other rule:

```python
xGrid, zGrid, values = fds.buildAbsPatch(
    patches, xmin, xmax, ymin, ymax, zmin, zmax,
    dx, dz, axis=-2)
```

It validates its inputs and raises `ValueError` with the offending patch
named if a patch does not sit on the assembled grid, or if
cell-centered and node-centered patches are mixed.

## Time averaging

`bndfsTimeAverage` averages every boundary file of a quantity and writes
new `.bf` files plus a smokeview file registering them, so the averaged
field can be opened in smokeview:

```python
outFiles, outQty, refFiles, newSmvFile = fds.bndfsTimeAverage(
    workingDir, chid, 'WALL TEMPERATURE', dt=30, outDir='./averaged')
```

See [Time Averaging](Time-Averaging).

## &GEOM boundary elements

Unstructured `&GEOM` surfaces write `.be` files against a `.gcf`
geometry file rather than `.bf` files.

```python
smvFile = fds.getSmvFile(workingDir, chid)
available = fds.getBndeQuantities(smvFile)

for name, info in available.items():
    print(name, info['quantity'], info['gridfile'])

vertices, faces, surfaces = fds.readGcfFile(gcfFile)
times, values, header = fds.readBeFile(beFile)
```

`values` is `array(NV, NT)`, one column per output time, indexed by
geometry vertex. `case002` in the bundled data has boundary element
output if you want something to try this on.

Writing a derived `.be` file is supported by `writeBeFile`, but
registering it in the smokeview file is not — `appendNewBeFileToSMV` is
an unimplemented stub. Add the `BNDE` record by hand, or follow what
`bndfsTimeAverage` does for rectangular boundary files.

## Exporting to ParaView

```python
fds.exportBndfDataToVtk(chid, workingDir, outDir='./vtk')
fds.exportBndeDataToVtk(chid, workingDir, outDir='./vtk')
```

Requires the `paraview` extra.
