# Troubleshooting

## `ImportError: pyfdstools.paraview requires optional dependencies...`

```
ImportError: pyfdstools.paraview requires optional dependencies which are
not installed (No module named 'evtk'). Install them with
'pip install pyfdstools[paraview]'.
```

The ParaView export and VTKHDF submodules wrap large third party
libraries which are not core dependencies. Install the extra the message
names:

```bash
python -m pip install "pyfdstools[paraview]"
python -m pip install "pyfdstools[vtk]"
python -m pip install "pyfdstools[all]"
```

Importing pyfdstools itself never requires them, so this only appears
when you reach for a name they provide.

## `AttributeError: module pyfdstools has no attribute ...`

The name does not exist. If the message goes on to say that an optional
submodule could not be imported, the name may be one of theirs — install
the extra as above. Otherwise check the spelling against the
[API Reference](API-Reference).

## `Warning did not find a 2-D slice of ... on axis ... at value ...`

`query2dAxisValue` returned `(None, None)` because the case has no slice
in that plane. It prints the planes that do exist:

```
Available slices of qty TEMPERATURE:
	Axis	Value
	1	2.5500
	-1	-1.0000
```

Three common causes:

* **Wrong plane.** Pick one from the list.
* **A coordinate that is close but not equal.** The match tolerance
  defaults to `1e-8`. Loosen it: `atol=1e-3`.
* **The slice is 3-D.** An axis of `-1` means that file covers a volume.
  Use `readSLCF3Ddata(chid, workingDir, qty, axis=..., value=...)`
  instead, which extracts a plane from it.

## `No 3D slice data found for qty ...`

`readSLCF3Ddata` returned `(False, False, False, False)`. It prints
every quantity the case wrote; check the spelling — FDS quantity names
are exact, including spaces (`'WALL TEMPERATURE'`, not
`'Wall Temperature'`).

Guard for it:

```python
grid, values, times, units = fds.readSLCF3Ddata(chid, workingDir, qty)
if grid is False:
    raise SystemExit('%s was not written by this case' % (qty))
```

## `ValueError: axis must be one of -3, -2, -1, 1, 2, or 3`

`queryBndf` takes the axis as a **scalar**, not a list, and it may be
signed:

```python
fds.queryBndf(workingDir, chid, fdsFilePath, ['WALL TEMPERATURE'], -2, 4.4)
#                                            ^ list              ^ scalar
```

where `fdsFilePath` is the path to the case's `.fds` input file.

The sign selects which side of the surface to read. See
[Boundary Data](Boundary-Data).

## `FileNotFoundError: No smokeview (.smv) file for chid ... was found`

Either the `chid` is wrong, or the working directory is. Check what is
actually there:

```python
import os
import zipfile

if workingDir.endswith('.zip'):
    print(zipfile.ZipFile(workingDir).namelist())
else:
    print(os.listdir(workingDir))
```

If the input file sets `RESULTS_DIR` on `&DUMP`, point pyfdstools at the
directory holding the **input file**, not the output subdirectory; it
follows `RESULTS_DIR` itself.

## `FileNotFoundError: No xyz file for chid ... Add WRITE_XYZ=.TRUE.`

Plot3D output carries no coordinates of its own, so `readPlot3Ddata`
needs the `.xyz` files. Re-run with:

```
&DUMP DT_PL3D=30., WRITE_XYZ=.TRUE. /
```

Slice and boundary readers do **not** need this — they take the grid
from the smokeview file.

## `FileNotFoundError: No csv file matching '...' for chid ...`

The case did not write that output. `_devc.csv` only exists if the case
had `&DEVC` lines; `_ctrl.csv` only if it had `&CTRL`.

## Everything comes back as NaN

Cells outside every mesh are `NaN` by design. If the whole array is
`NaN`:

* The plane may lie outside the domain. Check the mesh extents:
  ```python
  smv = fds.parseSMVFile(fds.getSmvFile(workingDir, chid))
  for trnx, trny, trnz in smv['grids']:
      print(trnx[[0, -1], 1], trny[[0, -1], 1], trnz[[0, -1], 1])
  ```
* The run may have produced no output at the requested time.

Use `np.nanmax`, `np.nanmean` and friends rather than their plain
counterparts throughout.

## The values look wrong by a constant

If you are comparing against numbers from an earlier pyfdstools release,
several calculations were corrected in v0.0.24 and now give different
answers:

| Routine | Was |
| --- | --- |
| slice time averaging with `time` and `dt` | summed N frames and divided by N−1 |
| `astFromGhf` | mixed Celsius and Kelvin; a zero flux did not return the gauge temperature |
| `getTwoZone` | reversed an ascending profile, reporting an upper layer cooler than the lower |
| `readPlot3Ddata` | read one mesh's `.q` file for every mesh |

Re-run the calculation rather than reconciling against the old numbers.

## `ValueError: setting an array element with a sequence` while parsing

Fixed in v0.0.24. Stretched meshes (`&TRNX`, `&TRNY`, `&TRNZ`) insert
extra lines into the smokeview grid records, which earlier releases
could not parse. Upgrade.

## `AttributeError: module 'numpy' has no attribute 'trapz'`

Fixed in v0.0.24. `np.trapz` was removed in numpy 2.0. Upgrade
pyfdstools, or pin `numpy<2`.

## Reading is slow

* **Reading from a zip archive.** Compressed members cannot be seeked
  into, so a large slice file is decompressed from the start on every
  read. Extract the case first if you are querying it repeatedly.
* **Reading every frame when you need one.** Pass `time` and `dt`.
* **Assembling a volume to take a plane out of it.** Pass `axis` and
  `value` to `readSLCF3Ddata`.
* **Re-scanning timestamps.** Pass `saveTimesFile=True` to cache them in
  a csv beside each slice file.

## No figure appears

On a machine with no display, matplotlib needs a non-interactive
backend:

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
```

or set `MPLBACKEND=Agg`. `fig.savefig` still works; `plt.show()` becomes
a no-op.

## A plot's colors jump between frames

Fix the color scale. Without `qnty_mn` and `qnty_mx`, each frame is
scaled to its own range:

```python
fds.plotSlice(data['x'], data['z'], data['datas'][:, :, i], axis,
              qnty_mn=20, qnty_mx=1000)
```

## Smokeview does not show a file I wrote

Smokeview only displays files a record in the smokeview file names. The
`slcfsTimeAverage` and `bndfsTimeAverage` routines write a new smokeview
file for you and return its path — open that one, not the original. If
you wrote the output file yourself, use `writeSliceToSmv` or
`buildBndfSmvLine` to add the record.

`&GEOM` boundary element files are the exception:
`appendNewBeFileToSMV` is an unimplemented stub, so the `BNDE` record has
to be added by hand.

## An input file round-trips with quotes around a number

The parameter is not in `fdsLineTypes`, so it is treated as a string.
Add it to the relevant `getXXXXtypes` method in
`pyfdstools/fdsTypes.py`:

```python
surfTypes['MY_PARAMETER'] = 'float'
```

and please open a pull request.

## Iterating a namelist group hits an integer

Multi-entry groups carry a bookkeeping `'unknownCounter'` entry:

```python
for key in model.meshes:
    if key == 'unknownCounter':
        continue
    ...
```

## Still stuck

Open an issue at
https://github.com/johodges/pyfdstools/issues with the pyfdstools
version (`fds.__version__`), the FDS version that produced the results,
and the full traceback.
