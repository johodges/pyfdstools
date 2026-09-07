# Quick Start

This page reads a 2-D temperature slice out of the FDS case bundled with
pyfdstools, plots it, and writes it to csv. Every snippet runs as-is
after `pip install pyfdstools` — no FDS installation and no simulation
of your own is needed.

## The bundled cases

pyfdstools ships several completed FDS cases so that the examples and
the test suite have something real to read. They live in
`pyfdstools/examples/data` and are read straight out of their zip
archives.

| Case | Meshes | Contains |
| --- | --- | --- |
| `case001` | 1 | slices, boundary, plot3D, particles, smoke3D, csv |
| `case002` | 4 | slices, `&GEOM` boundary elements, VTK output |
| `hfg_slice` | 4 | heat flux gauge devices and slices |
| `stretched_mesh_example` | 14 | stretched meshes (`&TRNZ`) |
| `visibility_adjustment` | 1 | a visibility slice and smoke3D |

Locate one like this:

```python
import os
import pyfdstools as fds

dataDir = os.path.join(os.path.dirname(fds.__file__), 'examples', 'data')
workingDir = os.path.join(dataDir, 'case001.zip')
chid = 'case001'
```

## What is in the case?

Before reading anything, ask the case what it wrote.

```python
quantities, files, dimensions, meshes, centers, units = \
    fds.readSLCFquantities(chid, workingDir)

for quantity, unit in sorted(set(zip(quantities, units))):
    print(quantity, unit)
```

```
TEMPERATURE C
U-VELOCITY m/s
V-VELOCITY m/s
W-VELOCITY m/s
```

## Read a 2-D slice

`query2dAxisValue` selects a plane by the axis it is normal to
(1 = x, 2 = y, 3 = z) and the coordinate along that axis, then assembles
every mesh that contributes to it onto one grid.

```python
data, units = fds.query2dAxisValue(
    workingDir, chid, 'TEMPERATURE',
    axis=1, value=2.55,     # the plane x = 2.55 m
    time=30, dt=60)         # mean over 0 s to 60 s

print(data['datas'].shape)  # (11, 25, 1) -> (n_y, n_z, n_times)
print(units)                # C
```

The returned dictionary holds four entries:

| Key | Contents |
| --- | --- |
| `x` | first in-plane coordinate, `array(N, M)` |
| `z` | second in-plane coordinate, `array(N, M)` |
| `datas` | values, `array(N, M, NT)` |
| `times` | timestamp of each frame |

Leave `time` and `dt` out to get every frame the case wrote:

```python
data, units = fds.query2dAxisValue(workingDir, chid, 'TEMPERATURE', 1, 2.55)
print(data['datas'].shape)   # (11, 25, 121)
```

## Plot it

```python
import matplotlib.pyplot as plt

fig, ax = fds.plotSlice(
    data['x'], data['z'], data['datas'][:, :, -1], axis=1,
    clabel='Temperature (%s)' % (units),
    qnty_mn=0, qnty_mx=1000, cbarnumticks=11)

fig.savefig('case001_temperature.png')
plt.show()
```

`plotSlice` defaults to the smokeview colormap and picks a figure aspect
ratio from the extent of the slice. See [Plotting](Plotting) for the
full set of options.

## Write it to csv

```python
fds.renderSliceCsvs(data, chid, outdir='.')
```

One csv per frame, named `<chid>_<time>.csv`, with the second in-plane
coordinate as the row index and the first as the column headers.

## Put it together

```python
import os
import matplotlib.pyplot as plt
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'
quantity = 'TEMPERATURE'
axis, value = 1, 2.55

data, units = fds.query2dAxisValue(
    workingDir, chid, quantity, axis, value, time=30, dt=60)

fds.renderSliceCsvs(data, chid, os.getcwd())

fig, ax = fds.plotSlice(
    data['x'], data['z'], data['datas'][:, :, -1], axis,
    clabel='%s (%s)' % (quantity, units),
    qnty_mn=0, qnty_mx=1000, cbarnumticks=11)
fig.savefig('%s_%s.png' % (chid, quantity))
plt.show()
```

## Point it at your own case

Replace `workingDir` with the directory holding your results and `chid`
with your `CHID`. A directory and a zip archive of that directory are
interchangeable everywhere in pyfdstools:

```python
workingDir = '/path/to/my/results'          # a directory
workingDir = '/path/to/my/results.zip'      # or the same thing zipped
```

If the plane you asked for does not exist, `query2dAxisValue` prints the
planes that do and returns `(None, None)` — see
[Troubleshooting](Troubleshooting).

## Where next

* [2-D Slice Data](2D-Slice-Data) — cell-centered slices, multi-mesh
  cases, reading one file at a time
* [Concepts](Concepts) — what the absolute grid is and why it matters
* [Boundary Data](Boundary-Data) — surface temperatures and heat fluxes
* [Input Files](Input-Files) — generating and editing `.fds` files
