# 3-D Slice Data

An `&SLCF` whose `XB` spans a volume rather than a plane produces a 3-D
slice. `readSLCF3Ddata` reads every mesh's contribution and interpolates
it onto one grid.

## Reading a volume

```python
import os
import numpy as np
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'

grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'TEMPERATURE')

print(grid.shape)     # (25, 11, 25, 3) -> (NX, NY, NZ, xyz)
print(values.shape)   # (25, 11, 25, 6) -> (NX, NY, NZ, NT)
print(units)          # C
```

Note the argument order: `readSLCF3Ddata` takes `chid` before
`workingDir`, unlike `query2dAxisValue`.

| Return value | Contents |
| --- | --- |
| `grid` | `array(NX, NY, NZ, 3)` of absolute coordinates |
| `values` | `array(NX, NY, NZ, NT)` of slice values |
| `times` | timestamp of each frame |
| `units` | units of the quantity |

Cells outside every mesh are `NaN`.

## When no 3-D slice exists

If the case wrote no 3-D slice of the requested quantity,
`readSLCF3Ddata` prints what the case does contain and returns
`(False, False, False, False)`:

```python
grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'NOT A REAL QUANTITY')

if grid is False:
    raise SystemExit('quantity not present in this case')
```

## Extracting a plane

Giving `axis` and `value` extracts just that plane. This is much cheaper
than reading the whole volume and slicing it yourself, because only the
meshes intersecting the plane are read:

```python
grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'TEMPERATURE', axis=3, value=1.0)

print(grid.shape)     # (25, 11, 2) -> (N, M, two in-plane coordinates)
print(values.shape)   # (25, 11, 6) -> (N, M, NT)
```

Meshes that do not span the requested plane are skipped, and `verbose=True`
reports which:

```python
grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'TEMPERATURE', axis=3, value=1.0, verbose=True)
```

## Slicing a volume you already read

If you already have the whole volume in memory, `findSliceLocation`
takes a plane out of it without going back to disk:

```python
grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'TEMPERATURE')

x, z, frame = fds.findSliceLocation(grid, values[:, :, :, -1], axis=3, value=1.0)

fig, ax = fds.plotSlice(x, z, frame, axis=3, clabel='Temperature (C)')
```

It picks the nearest grid plane rather than interpolating, so the plane
you get is the one closest to `value`.

## Extracting a point

```python
history = fds.extractPoint([2.0, 1.5, 1.0], grid, values)
```

Returns the time history at the nearest grid point, and warns if that
point is more than 0.25 m away from the one you asked for.

## Dumping the volume to csv

`pyfdstools/examples/dump_3d_slice_to_csv.py` writes a 3-D slice out as
one csv per z-plane per frame:

```bash
python dump_3d_slice_to_csv.py \
    --chid case001 --quantity TEMPERATURE \
    --working_dir data/case001.zip
```

## Speed and memory

A 3-D slice over a long run is large: `NX × NY × NZ × NT` float32 values
per mesh, plus the assembled array. Two things help:

* **Restrict the time range.** Passing `time` and `dt` reads only the
  frames inside the window.
* **Extract the plane at read time.** Passing `axis` and `value` avoids
  materialising the volume at all.

```python
grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'TEMPERATURE',
    time=60, dt=20, axis=3, value=1.0)
```

`saveTimesFile=True` caches each slice file's timestamps in a csv next
to it, which speeds up repeated reads of a case on disk.

## Exporting to ParaView

For interactive 3-D work, export to VTK and open the result in ParaView:

```python
fds.exportSl3dDataToVtk(chid, workingDir, outDir='./vtk')
```

Requires the `paraview` extra. See
[Exporting to ParaView](Exporting-to-ParaView).

## Deprecated

`readSLCF3DdataXYZ` does the same job using the `.xyz` files FDS writes
when `WRITE_XYZ=.TRUE.` is set on `&DUMP`. It emits a
`DeprecationWarning`; use `readSLCF3Ddata`, which reads the mesh grids
from the smokeview file and so works on any case.
