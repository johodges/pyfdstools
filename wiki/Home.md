# pyfdstools

A Python library for building and post-processing [Fire Dynamics
Simulator](https://github.com/firemodels/fds) (FDS) models.

pyfdstools reads every binary output format FDS writes — slice files,
boundary files, plot3D, particles, smoke3D, geometry — assembles them
across meshes onto a single grid, and hands you numpy arrays. It also
reads and writes FDS input files as Python objects, so a model can be
generated, edited and re-saved from a script.

It reads results straight out of a zip archive, so simulation output can
stay compressed.

## New here?

* **[Installation](Installation)** — pip, source, and the optional extras
* **[Quick Start](Quick-Start)** — read a slice, plot it, dump it to csv
* **[Concepts](Concepts)** — CHID, working directories, meshes and the
  absolute grid

## Reading results

| Page | FDS output |
| --- | --- |
| [2-D Slice Data](2D-Slice-Data) | `.sf` files from `&SLCF` with a fixed coordinate |
| [3-D Slice Data](3D-Slice-Data) | `.sf` files from `&SLCF` covering a volume |
| [Boundary Data](Boundary-Data) | `.bf` files from `&BNDF`, and `.be` from `&GEOM` |
| [Plot3D Data](Plot3D-Data) | `.q` and `.xyz` files from `&DUMP` |
| [Device and CSV Data](Device-and-CSV-Data) | `_devc.csv`, `_hrr.csv` and friends |
| [Particle and Smoke3D Data](Particle-and-Smoke3D-Data) | `.prt5` and `.s3d` files |

## Working with results

* **[Plotting](Plotting)** — `plotSlice`, colormaps, highlighting a threshold
* **[Time Averaging](Time-Averaging)** — averaging windows, and writing
  averaged fields back out for smokeview
* **[Exporting to ParaView](Exporting-to-ParaView)** — VTK export of every
  output type
* **[Model Uncertainty](Model-Uncertainty)** — applying the FDS validation
  guide bias and standard deviation

## Building models

* **[Input Files](Input-Files)** — reading, editing, generating and writing
  `.fds` files

## Reference

* **[API Reference](API-Reference)** — every public routine, grouped by task
* **[Troubleshooting](Troubleshooting)** — what the common errors mean
* **[Contributing](Contributing)** — running the tests, the CI, the layout

## Citing

> Hodges, J. L., pyFDStools: A Python Package to Assist in Developing and
> Post-Processing Data Produced Through the Computational Fluid Dynamics
> Software Fire Dynamics Simulator, (2020), GitHub repository,
> https://github.com/johodges/pyfdstools.
