# Concepts

Five ideas explain most of the pyfdstools API. Reading this page once
will save you guessing at argument names later.

## CHID and the working directory

Every pyfdstools reader takes the same two arguments first:

```python
fds.query2dAxisValue(workingDir, chid, ...)
fds.readSLCF3Ddata(chid, workingDir, ...)   # note: this one is reversed
fds.load_csv(workingDir, chid, ...)
```

* **`chid`** is the `CHID` from the `&HEAD` line of the input file. FDS
  names every output file after it.
* **`workingDir`** is the directory holding those output files.

The argument order is not consistent across the package —
`readSLCF3Ddata` and `readPlot3Ddata` take `chid` first, the rest take
the directory first. Check the signature if a call misbehaves.

### Zip archives

Anywhere a directory is accepted, a zip archive of that directory works
too:

```python
workingDir = '/path/to/results'
workingDir = '/path/to/results.zip'
```

pyfdstools reads members out of the archive without extracting it, which
is convenient for keeping a run's output compressed. Internally, files
inside an archive are named
`<archive>.zip<separator><name inside the archive>`; if you see a path
like that in an error message, that is what it means.

Reading from an archive is slower for large slice files, because a
compressed member has to be decompressed from the start rather than
seeked into. For repeated queries over a big case, extract it first.

### RESULTS_DIR

If the input file sets `RESULTS_DIR` on `&DUMP`, FDS writes its output
into that subdirectory rather than next to the input file. Point
pyfdstools at the directory holding the **input file**; it reads the
input file and follows `RESULTS_DIR` itself. This only works for
directories, not archives.

## Meshes and the absolute grid

An FDS case is divided into meshes, and each mesh writes its own output
files. A slice through a multi-mesh case is therefore spread across
several files with different grids.

The high-level readers hide this. `query2dAxisValue`,
`readSLCF3Ddata`, `readPlot3Ddata` and `queryBndf` each:

1. find every file that contributes to what you asked for,
2. build an **absolute grid** spanning all of them,
3. interpolate each mesh's data onto it,
4. return one array.

Cells not covered by any mesh come back as `NaN`, so use `np.nanmax`,
`np.nanmean` and friends rather than their plain counterparts:

```python
peak = np.nanmax(data['datas'])
```

The absolute grid is built from the union of the mesh coordinates, so a
case with meshes of different resolutions gets a non-uniform absolute
grid. Pass `makeUniform=True` to `getAbsoluteGrid` if you need a
rectangular one at the finest resolution.

### Reading one mesh at a time

When you want the raw per-mesh data instead, the lower-level readers
give it to you:

```python
for slcfFile in fds.getFileList(workingDir, chid, 'sf'):
    lims, data, times = fds.readSingleSlcfFile(slcfFile)
```

## Axis and value

A plane is named by the axis it is **normal to** and its coordinate
along that axis:

| `axis` | Plane | In-plane coordinates returned as |
| --- | --- | --- |
| 1 | x = value | `x` is y, `z` is z |
| 2 | y = value | `x` is x, `z` is z |
| 3 | z = value | `x` is x, `z` is y |

The returned coordinate arrays are always called `x` and `z` whichever
plane you asked for; for a horizontal slice (`axis=3`) the array named
`z` holds y coordinates. `plotSlice` takes the same `axis` argument and
labels the figure accordingly.

Boundary queries also accept a **signed** axis, where the sign selects
which side of a surface to read: `axis=-2` reads the faces whose outward
normal points in −y. See [Boundary Data](Boundary-Data).

## Time selection

Every reader takes the same pair of optional arguments:

| `time` | `dt` | Result |
| --- | --- | --- |
| omitted | omitted | every frame in the file |
| given | omitted | the single frame nearest `time` |
| given | given | the mean of the frames within `time ± dt/2` |
| omitted | given | a running mean of width `dt` at every frame |

The last form is only supported by `read2dSliceFile` and
`query2dAxisValue`.

Frames are selected by timestamp, not by counting records, so a case
whose output interval changed part way through is handled correctly.

```python
# every frame
data, units = fds.query2dAxisValue(workingDir, chid, qty, 1, 2.55)

# the frame nearest t = 30 s
data, units = fds.query2dAxisValue(workingDir, chid, qty, 1, 2.55, time=30)

# the mean over 0 s to 60 s
data, units = fds.query2dAxisValue(workingDir, chid, qty, 1, 2.55,
                                   time=30, dt=60)
```

## Cell-centered data

FDS writes slices either at cell nodes or at cell centers, depending on
whether `CELL_CENTERED=.TRUE.` was set on the `&SLCF` line. Smokeview
records which in its own file, and pyfdstools reads that flag rather
than guessing:

```python
quantities, files, dims, meshes, centers, units = \
    fds.readSLCFquantities(chid, workingDir)

for f, c in zip(files, centers):
    print(os.path.basename(f), 'cell-centered' if c else 'node-centered')
```

A cell-centered slice has one fewer value per axis than the node grid,
so the coordinates are shifted by half a cell before the data are mapped
onto the absolute grid.

## Endianness

FDS writes its binary output in the byte order of the machine it ran on
and records that order in a `.end` file. pyfdstools reads it
automatically via `getEndianness`, so results produced on a different
architecture read back correctly. Little-endian is assumed when no
`.end` file is present, which is right for every mainstream platform.
