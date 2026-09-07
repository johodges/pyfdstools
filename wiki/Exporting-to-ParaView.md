# Exporting to ParaView

`pyfdstools.paraview` converts every FDS output type to VTK, so the
results can be explored interactively in
[ParaView](https://www.paraview.org/).

## Installing

The export routines need `pyevtk` and `numpy-stl`, which are not core
dependencies:

```bash
python -m pip install "pyfdstools[paraview]"
```

Until they are installed, reaching for one of these names raises an
ImportError naming the extra:

```
ImportError: pyfdstools.paraview requires optional dependencies which are
not installed (No module named 'evtk'). Install them with
'pip install pyfdstools[paraview]'.
```

Nothing else about pyfdstools changes; importing the package works
either way.

## Exporting

Each routine takes the same arguments and writes one file per timestep
plus a `.pvd` index. Open the `.pvd` in ParaView and it loads the whole
time series with the timeline populated.

```python
import os
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'
outDir = './vtk'

fds.exportSl3dDataToVtk(chid, workingDir, outDir=outDir)   # 3-D slices
fds.exportSl2dDataToVtk(chid, workingDir, outDir=outDir)   # 2-D slices
fds.exportBndfDataToVtk(chid, workingDir, outDir=outDir)   # boundary files
fds.exportBndeDataToVtk(chid, workingDir, outDir=outDir)   # &GEOM boundary
fds.exportPrt5DataToVtk(chid, workingDir, outDir=outDir)   # particles
fds.exportS3dDataToVtk(chid, workingDir)                   # smoke3D
```

| Argument | Effect |
| --- | --- |
| `outtimes` | times to export; every output time when omitted |
| `outDir` | where to write; alongside the results when omitted |
| `binary` | appended binary data rather than ASCII, default `True` |
| `dtype` | VTK datatype name, chosen from the data when omitted |
| `ftype` | `'ImageData'` or `'RectilinearGrid'`, slices only |

### Uniform or stretched meshes

`ftype='ImageData'` (the default) stores only an origin and a spacing,
which loads faster and uses less disk. It is only correct for uniformly
spaced meshes. For a case using `&TRNX`, `&TRNY` or `&TRNZ`, use
`'RectilinearGrid'`, which stores the coordinates explicitly:

```python
fds.exportSl3dDataToVtk(chid, workingDir, outDir=outDir,
                        ftype='RectilinearGrid')
```

### Exporting a subset of times

A long run produces a lot of files. Restrict them:

```python
import numpy as np

fds.exportSl3dDataToVtk(chid, workingDir, outDir=outDir,
                        outtimes=np.arange(0, 601, 30))
```

### Smoke3D

Smoke3D data are stored as single bytes and are exported that way by
default, which keeps the files small. Pass `decode=True` for physical
units:

```python
fds.exportS3dDataToVtk(chid, workingDir, decode=True, dtype='Float32')
```

## Geometry

ParaView cannot read FDS obstructions, so export them as an STL to show
alongside the field data:

```python
fds.obstToStl(workingDir, chid, outDir=outDir)
```

This writes `<chid>.stl` covering every rectangular obstruction in the
case. `&GEOM` surfaces are already triangulated; export them with
`exportBndeDataToVtk`.

## In ParaView

1. Open the `.pvd` file. ParaView loads every timestep and enables the
   timeline.
2. Add the `.stl` with **File > Open** and set it to a solid color to
   show the geometry.
3. For a 3-D slice, **Volume** representation gives a smoke-like render;
   **Slice** and **Contour** filters cut through it.
4. For particles, apply a **Glyph** filter and scale by a particle
   quantity.
5. Colour by a quantity from the toolbar, then **Rescale to Custom Data
   Range** and fix the range so it does not rescale on every timestep.

## VTKHDF

Recent FDS versions can write VTKHDF directly. To read that back into
the same dictionary layout `query2dAxisValue` returns:

```bash
python -m pip install "pyfdstools[vtk]"
```

```python
data, units = fds.query2dAxisValue_vtkhdf(
    workingDir, chid, 'TEMPERATURE', axis=1, value=2.55, time=30)
```

The result has the same `x`, `z`, `datas` and `times` keys, so
downstream code does not need to know which source it came from.

## The standalone script

`pyfdstools/example_paraview_export.py` predates the `paraview` module
and duplicates part of it. It is kept as a worked example of the VTK
file format; use `pyfdstools.paraview` for real work.
