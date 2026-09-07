# Particle and Smoke3D Data

## Particle output

`&PART` with `&DUMP DT_PART` writes Lagrangian particle positions and
quantities to one `.prt5` file per mesh.

### Reading every mesh

```python
import os
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'

particles, byTime = fds.importParticles(workingDir, chid)
```

| Return value | Contents |
| --- | --- |
| `particles['tags']` | per particle history, keyed by particle tag |
| `particles['classes']` | per particle class metadata and counts |
| `particles['times']` | every output time across all meshes |
| `byTime` | which particle tags exist at each time |

A particle is tracked by its tag, so its history survives crossing a
mesh boundary.

### Reading one file

```python
prt5File = fds.getFileList(workingDir, chid, 'prt5')[0]
particles, times = fds.importParticle(prt5File)
```

### What classes and quantities are present

```python
particles, times = fds.importParticle(prt5File)

for pid, info in particles['classes'].items():
    if pid == 'times':
        continue
    print(pid, info['qtyNames'], info['qtyUnits'])
```

### Cross-checking against the .bnd file

`importParticle_meta` reads the same data but takes its timestamps from
the companion `.prt5.bnd` metadata file and checks each one against the
timestamp in the binary record, raising `ValueError` if they disagree:

```python
particles, times = fds.importParticle_meta(prt5File)
```

Use it when a run was interrupted and the binary may be inconsistent
with its metadata. Otherwise use `importParticle`, which does not need
the `.bnd` file to be present.

### Exporting to ParaView

Particles are far easier to inspect interactively than as arrays:

```python
fds.exportPrt5DataToVtk(chid, workingDir, outDir='./vtk')
```

Requires the `paraview` extra. See
[Exporting to ParaView](Exporting-to-ParaView).

## Smoke3D output

`&DUMP` writes smoke3D (`.s3d`) files for soot density, HRRPUV and
temperature. These are the fields smokeview renders as smoke and fire.

Smoke3D data are stored as single bytes, run-length encoded, scaled to a
fixed range per quantity. `extractS3dValues` decodes them back to
physical units.

```python
values, times = fds.extractS3dValues(workingDir, chid)

for meshNum in values:
    for quantity, field in values[meshNum].items():
        print(meshNum, quantity, field.shape)
```

`values` is keyed by mesh number, then by quantity. `times` holds the
output times. If the case wrote no smoke3D output, both come back as
`None`.

### Raw bytes

Pass `decode=False` to get the stored bytes without scaling, which is
what you want if you intend to write them back out unchanged:

```python
values, times = fds.extractS3dValues(workingDir, chid, decode=False)
```

### How the scaling works

| Quantity | Stored range maps to |
| --- | --- |
| `SOOT DENSITY` | derived from the mass extinction coefficient and the cell size |
| `HRRPUV` | 0 to 1200 kW/m3 |
| `TEMPERATURE` | ambient to 2000 C |

Because the range is fixed, smoke3D data are lossy: a field above the
range is clipped, and the resolution is the range divided by 254. Use
slice output when you need the actual values.

### Writing smoke3D files

```python
fds.writeS3dFile(path, times, field, quantity='TEMPERATURE', dx=0.1)
```

`dx` is a representative cell size, needed for the soot density scaling.
`writeS3dFile_debug` writes the same file while comparing each record
against a reference file, which is how a format mismatch gets located.

`pyfdstools.examples.examplePostProcessVisibility` is a worked example:
it reads a visibility slice, rescales it for a different mass extinction
factor, and writes the result back out.

### Exporting to ParaView

```python
fds.exportS3dDataToVtk(chid, workingDir)
```

Exports the raw bytes by default, which keeps the files small; pass
`decode=True` for physical units.
