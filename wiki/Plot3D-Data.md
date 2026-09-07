# Plot3D Data

Plot3D output is a snapshot of five gas phase quantities over the whole
mesh at a set of output times. FDS writes one `.q` file per mesh per
time, plus one `.xyz` file per mesh holding the grid.

## Enabling it in FDS

```
&DUMP DT_PL3D=30., WRITE_XYZ=.TRUE. /
```

`WRITE_XYZ` is required: the `.q` files hold no coordinates, so the
`.xyz` files are the only source of the grid.

## Reading a snapshot

```python
import os
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'

grid, values = fds.readPlot3Ddata(chid, workingDir, time=30.0)

print(grid.shape)     # (25, 11, 25, 3)
print(values.shape)   # (25, 11, 25, 5)
```

Note the argument order: `chid` comes before `workingDir`.

The output time nearest the one you ask for is used. For a multi-mesh
case each mesh's own `.q` file for that time is read.

## The five quantities

The last axis of `values` is fixed by the plot3D format:

| Index | Quantity | Units |
| --- | --- | --- |
| 0 | temperature | C |
| 1 | u-velocity | m/s |
| 2 | v-velocity | m/s |
| 3 | w-velocity | m/s |
| 4 | HRRPUV | kW/m3 |

```python
temperature = values[:, :, :, 0]
hrrpuv = values[:, :, :, 4]
```

By default FDS writes those five; `QUANTITIES` on `&DUMP` can change
what occupies each slot, in which case the table above no longer
applies.

## Taking a plane

```python
x, z, frame = fds.findSliceLocation(grid, values[:, :, :, 0], axis=2, value=4.4)

fig, ax = fds.plotSlice(x, z, frame, axis=2,
                        clabel='Temperature (C)', qnty_mn=20, qnty_mx=200)
fig.savefig('pl3d_temperature.png')
```

`plot3d=True` returns all five fields at once instead:

```python
x, z, T, U, V, W, HRR = fds.findSliceLocation(
    grid, values, axis=2, value=4.4, plot3d=True)
```

Which is what `visualizePlot3D` takes:

```python
fds.visualizePlot3D(x, z, T, U, V, W, HRR, qnty_mn=20, qnty_mx=200)
```

`plotSlice` gives far more control; `visualizePlot3D` is a convenience
wrapper that plots the temperature field only.

## Printing the extents

```python
fds.readPlot3Ddata(chid, workingDir, 30.0, verbose=True)
```

prints the coordinate and value ranges of each mesh as it is read, which
is a quick sanity check that the case is what you think it is.

## Reading one file

```python
grid, gridHeader = fds.readXYZfile(xyzFile)
values, dataHeader = fds.readP3Dfile(qFile)
```

Both return flat arrays with the header separately; `rearrangeGrid`
turns the grid into meshgrid form:

```python
xGrid, yGrid, zGrid = fds.rearrangeGrid(grid)
```

## Writing plot3D files

```python
fds.writeP3Dfile(path, values)     # values.shape == (NX, NY, NZ, 5)
fds.writeXYZfile(path, grid)
```

Useful for feeding a derived field to smokeview, which will display any
`.q` file named to match the case's `.xyz` files.

## Errors

`readPlot3Ddata` raises `FileNotFoundError` if the case has no `.xyz`
files, with a message reminding you to set `WRITE_XYZ=.TRUE.`, or if a
mesh has no `.q` file at any time.

## The command line example

```bash
python extract_2d_slice_from_pl3d.py \
    --chid case001 --quantity TEMPERATURE \
    --axis 2 --value 4.4 --time 30 \
    --working_dir data/case001.zip
```
