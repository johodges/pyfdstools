# Contributing

## Getting set up

```bash
git clone https://github.com/johodges/pyfdstools
cd pyfdstools
python -m venv .venv
source .venv/bin/activate          # .venv\Scripts\activate on Windows
python -m pip install -e ".[all,test]"
```

The editable install means your changes take effect without
reinstalling.

## Running the tests

```bash
python -m pytest tests -m "not slow"    # fast: ~20 s
python -m pytest tests                  # everything, including the examples
```

The suite runs against the FDS cases bundled in
`pyfdstools/examples/data`, so it needs no FDS installation and no
network access.

| Marker | Contents |
| --- | --- |
| default | unit tests, roughly 20 s |
| `slow` | runs each bundled example end to end, roughly 45 s |

Run one file, or one test:

```bash
python -m pytest tests/test_slices.py -v
python -m pytest tests/test_slices.py::test_query2dAxisValue -v
```

## Layout

```
pyfdstools/
  smokeviewParser.py       parses .smv files
  extractPlot3Ddata.py     slice and plot3D readers, plotSlice
  extractBoundaryData.py   .bf readers, the fdspatch class
  extractGeomData.py       &GEOM geometry and boundary elements
  extractParticleData.py   .prt5 readers
  extractS3D.py            smoke3D readers and writers
  extractCSVdata.py        FDS csv output
  fdsFileOperations.py     reading and writing .fds input files
  fdsTypes.py              namelist parameter datatypes
  fdsErrorCalculation.py   model uncertainty
  utilities.py             shared helpers
  colorSchemes.py          colormaps and color sequences
  paraview.py              VTK export          (optional: pyevtk, numpy-stl)
  vtkhdf_plugin.py         VTKHDF reader       (optional: vtk)
  examples/                runnable example scripts and the bundled cases
tests/                     the pytest suite
wiki/                      the source of this wiki
.github/workflows/         CI
```

## Optional dependencies

`paraview` and `vtkhdf_plugin` are loaded on demand through a module
level `__getattr__` in `pyfdstools/__init__.py`, so `import pyfdstools`
works without `vtk`, `pyevtk` or `numpy-stl`.

**Do not add an unconditional import of an optional library to any other
module.** Doing so breaks every user who installed only the core
requirements. The `minimal-install` CI job exists to catch this: it
installs with no extras and asserts none of them are importable while
the package still works.

## Continuous integration

Four workflows run on every push and pull request.

| Workflow | Jobs |
| --- | --- |
| `tests.yml` | unit tests on Python 3.9–3.13 and on Linux, Windows and macOS; a minimal-install job; a job pinning both the oldest supported dependencies and numpy 2 |
| `examples.yml` | the bundled examples end to end on Linux and Windows, and the optional extras; also weekly on a schedule |
| `lint.yml` | pyflakes, `compileall` on the oldest and newest Python, and a distribution build checked with twine |
| `publish.yml` | builds and publishes to PyPI on a `v*` tag |

`lint.yml` **fails** on undefined names, because those raise `NameError`
the first time the branch containing them runs. Unused imports and
unused locals are reported but do not fail the build.

Reproduce the lint job locally:

```bash
python -m pip install pyflakes
python -m pyflakes pyfdstools tests | grep "undefined name '"
```

## Writing tests

Fixtures in `tests/conftest.py` give you the bundled cases:

| Fixture | Contents |
| --- | --- |
| `case001_zip` | single mesh, slices, boundary, plot3D, particles, smoke3D |
| `case002_zip` | four meshes, `&GEOM` output |
| `hfg_zip` | four meshes, heat flux gauges |
| `stretched_zip` | fourteen stretched meshes |
| `case001_dir` | `case001` extracted to a directory |
| `outdir` | an empty temporary directory |
| `data_dir` | the bundled data directory |

`case001_zip` and `case001_dir` are the same case reached two ways;
readers take different code paths for each, so both are worth covering.

Assert on what the numbers mean, not just on shapes:

```python
def test_readSingleSlcfFile_time_average_matches_mean(case001_zip):
    """The averaged frame must be the mean of the frames in the window."""
    sf = sorted(fds.getFileList(case001_zip, 'case001', 'sf'))[0]
    _, allData, times = fds.readSingleSlcfFile(sf)
    times = np.asarray(times)

    _, averaged, _ = fds.readSingleSlcfFile(sf, time=30.0, dt=20.0)

    inds = np.where((times >= 20.0) & (times <= 40.0))[0]
    expected = allData[:, :, :, inds].mean(axis=3)
    assert np.allclose(averaged[:, :, :, 0], expected, rtol=1e-5, atol=1e-4)
```

When you fix a bug, add the test that would have caught it, and say in
its docstring what the old behaviour was.

## Style

The codebase is not uniformly formatted and there is no formatter in CI.
Match the surrounding code:

* four spaces, no tabs
* `camelCase` for functions and variables, which is what the package uses
  throughout
* numpydoc docstrings on every public function, with `Parameters`,
  `Returns` and `Raises` where they apply
* `is None` rather than `== None`
* raise an exception with an actionable message rather than printing and
  returning `None`

## Deprecating rather than removing

Scripts written against this package are often years old. When a routine
is superseded, keep it and emit a warning:

```python
warnings.warn(
    "oldRoutine is deprecated and will be removed in a future "
    "release; use newRoutine instead.",
    DeprecationWarning, stacklevel=2)
```

and note the replacement in the docstring with a `.. deprecated::`
block.

## Releasing

1. Bump `__version__` in `pyfdstools/_version.py`.
2. Commit and push.
3. Tag it: `git tag v0.0.24 && git push origin v0.0.24`.

`publish.yml` verifies the tag matches `__version__`, builds the sdist
and wheel, runs twine check and publishes to PyPI through trusted
publishing. Configure the project on PyPI to trust the workflow before
the first release.

## Editing this wiki

The pages live in `wiki/` in the main repository, so they are reviewed
alongside code changes. To publish them:

```bash
git clone https://github.com/johodges/pyfdstools.wiki.git
cp wiki/*.md pyfdstools.wiki/
cd pyfdstools.wiki && git add -A && git commit -m "Update wiki" && git push
```

Run any snippet you add before committing it — every example on this
wiki was executed against the bundled cases.
