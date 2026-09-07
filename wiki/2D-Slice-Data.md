# 2-D Slice Data

FDS writes each `&SLCF` to a `.sf` file per mesh. A slice with one
degenerate dimension is 2-D; this page covers those. For slices covering
a volume see [3-D Slice Data](3D-Slice-Data).

## Finding out what is available

```python
import os
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'

quantities, files, dimensions, meshes, centers, units = \
    fds.readSLCFquantities(chid, workingDir)
```

| Return value | Contents |
| --- | --- |
| `quantities` | FDS quantity in each file |
| `files` | path of each `.sf` file |
| `dimensions` | `[iX, eX, iY, eY, iZ, eZ]` extent in cell indices |
| `meshes` | mesh string parsed from the file name |
| `centers` | `True` if the slice is cell-centered |
| `units` | units of the quantity |

To list the planes as coordinates rather than cell indices, ask for a
plane you know does not exist — the error path prints them:

```python
fds.query2dAxisValue(workingDir, chid, 'TEMPERATURE', 3, 1e9)
```

```
Warning did not find a 2-D slice of TEMPERATURE on axis 3 at value 1000000000.0000
Available slices of qty TEMPERATURE:
	Axis	Value
	1	2.5500
	-1	-1.0000
Note an axis of -1 indicates a 3-D slice which is not currently supported in this function.
```

An axis of `-1` means that file holds a 3-D slice.

## Reading a plane across every mesh

`query2dAxisValue` is the routine you want almost always. It finds every
`.sf` file holding the requested quantity in the requested plane,
assembles them onto one grid, and returns the result.

```python
data, units = fds.query2dAxisValue(
    workingDir, chid, 'TEMPERATURE', axis=1, value=2.55, time=30, dt=60)

x = data['x']          # array(N, M) first in-plane coordinate
z = data['z']          # array(N, M) second in-plane coordinate
values = data['datas'] # array(N, M, NT)
times = data['times']  # timestamp of each frame
```

### Arguments

| Argument | Meaning |
| --- | --- |
| `axis` | axis the plane is normal to: 1 = x, 2 = y, 3 = z |
| `value` | coordinate of the plane along that axis |
| `time`, `dt` | frame selection, see [Concepts](Concepts) |
| `atol` | tolerance when matching `value` against the planes present (default `1e-8`) |
| `printInfo` | print each file as it is read |
| `verbose` | print progress per mesh |

### Matching the coordinate

`value` is matched against the slice coordinates with `atol`, which
defaults to `1e-8` — effectively exact. If your slice was written at
`2.5499999` and you ask for `2.55`, loosen it:

```python
data, units = fds.query2dAxisValue(
    workingDir, chid, 'TEMPERATURE', 1, 2.55, atol=1e-3)
```

### Cells outside every mesh

Anywhere the absolute grid is not covered by a mesh comes back as `NaN`.
Use the `nan`-aware reductions:

```python
import numpy as np
print('peak', np.nanmax(values))
print('mean', np.nanmean(values[:, :, -1]))
```

## Reading a single slice file

When you want one mesh's raw data on its own grid:

```python
slcfFile = fds.getFileList(workingDir, chid, 'sf')[0]
x, z, values, times, coords = fds.read2dSliceFile(slcfFile, chid)
```

`coords` is `[xmin, xmax, ymin, ymax, zmin, zmax]` of the slice; exactly
one axis has equal bounds, which is the plane it lies in.
`read2dSliceFile` returns `None` if the file holds a 3-D slice.

Pass `cen=True` for a cell-centered slice to have the coordinates
shifted to cell centers, and `grid=` to supply a mesh grid you already
have rather than letting it read the mesh's `.xyz` file.

For the rawest form, which keeps the singleton axis and does no
coordinate work at all:

```python
lims, values, times = fds.readSingleSlcfFile(slcfFile)
# values.shape == (NX+1, NY+1, NZ+1, NT), one axis of length 1
```

## Running averages

Passing `dt` without `time` applies a running mean of that width to
every frame, which is the usual way to smooth a noisy field before
plotting:

```python
x, z, smoothed, times, coords = fds.read2dSliceFile(
    slcfFile, chid, dt=20.0)
```

Each output frame is the mean of the frames whose timestamps fall within
`t ± dt/2`. See [Time Averaging](Time-Averaging) for writing the
averaged field back out as a slice file smokeview can open.

## Extracting a plane from a 3-D slice

If the case wrote a 3-D slice rather than a 2-D one, take a plane out of
it instead:

```python
grid, values, times, units = fds.readSLCF3Ddata(
    chid, workingDir, 'TEMPERATURE', axis=3, value=1.0)
```

Giving `axis` and `value` makes `readSLCF3Ddata` extract only that plane,
which avoids assembling the whole volume. See
[3-D Slice Data](3D-Slice-Data).

## Blocked cells

Obstructions occupy part of a slice plane. To mask them out:

```python
blocked = fds.getBlockedCellsInPlane(workingDir, chid, axis=3, value=0.0)
values[blocked[..., None].astype(bool).repeat(values.shape[2], axis=2)] = np.nan
```

`blocked` is `1` where an obstruction blocks the plane and `0`
elsewhere, sized on the same absolute grid.

## Writing csv

```python
fds.renderSliceCsvs(data, chid, outdir='.')
```

One file per frame, named `<chid>_<time>.csv`. Rows are the second
in-plane coordinate, columns the first.

To do it yourself for a single frame:

```python
import pandas as pd

frame = pd.DataFrame(data['datas'][:, :, -1].T,
                     index=data['z'][0, :],
                     columns=data['x'][:, 0])
frame.to_csv('final_frame.csv')
```

## Writing a derived slice back out

You can compute a new field and write it as a slice file that smokeview
will display alongside the originals. See
[Time Averaging](Time-Averaging) for the mechanics, and
`pyfdstools/examples/make_slice_from_devc.py` for a worked example that
builds a radiative heat flux slice from a grid of devices.

## The command line example

`pyfdstools/examples/dump_2d_slice_to_csv.py` does all of the above and
takes its parameters as arguments:

```bash
python dump_2d_slice_to_csv.py \
    --chid case002 --quantity TEMPERATURE \
    --axis 3 --value 7.2 --time 30 --dt -1 \
    --working_dir data/case002.zip
```

## Deprecated

`query2dAxisValueXYZ` and `readSLCF2Ddata` do the same job using the
`.xyz` files FDS writes when `WRITE_XYZ=.TRUE.` is set on `&DUMP`. They
emit a `DeprecationWarning`; use `query2dAxisValue`, which reads the mesh
grids from the smokeview file and so works on any case.
