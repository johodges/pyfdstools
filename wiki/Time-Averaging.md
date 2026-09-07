# Time Averaging

Instantaneous fields from an LES simulation are noisy. pyfdstools
averages at three levels: in the reader, over an array in memory, and
over whole output files that are written back out for smokeview.

## While reading

Every reader takes `time` and `dt`. Passing both returns a single frame
averaged over `time ± dt/2`:

```python
data, units = fds.query2dAxisValue(
    workingDir, chid, 'TEMPERATURE', 1, 2.55, time=60, dt=30)
# one frame, the mean of every frame between t = 45 s and t = 75 s
```

Passing `dt` **without** `time` applies a running mean of that width to
every frame:

```python
x, z, smoothed, times, coords = fds.read2dSliceFile(
    slcfFile, chid, dt=20.0)
```

Frames are selected by timestamp, so a case whose output interval
changed part way through still averages over the right window.

> **A note on older results.** Releases before v0.0.24 summed the frames
> in the window and divided by one fewer than the count, which biased
> every averaged slice — by up to 9 °C on the bundled `case001`. If you
> have numbers produced by an earlier release, re-run the averaging.

## Over an array

`timeAverage` applies a boxcar average to an array whose last axis is
time:

```python
import numpy as np

data, units = fds.query2dAxisValue(
    workingDir, chid, 'TEMPERATURE', 1, 2.55)   # every frame

values = data['datas']                          # (N, M, NT)
times = np.asarray(data['times'])

averaged, outTimes = fds.timeAverage(values, times, window=10.0)
```

The data are first interpolated onto a uniform time base built from the
smallest positive timestep present, then convolved with a rectangular
window. That means the input does not need to be evenly spaced.

| Argument | Effect |
| --- | --- |
| `window` | width of the averaging window in seconds |
| `outdt` | timestep of the returned series; the internal uniform base is used when `<= 0` |
| `smoothEnds` | average the first and last half-windows over the partial window available |
| `queryTime` | average over a single window centred here and return one frame |

### The ends

With `smoothEnds=False` (the default), the first and last half-window of
the output are copied from the interpolated series **unaveraged**. This
is deliberate — it keeps the output the same length as the input — but
it means the ends are noisier than the middle. Either set
`smoothEnds=True`, or trim them:

```python
halfWindow = int(round(window / (times[1] - times[0])))
interior = averaged[:, :, halfWindow:-halfWindow]
```

If `window` exceeds the length of the series, `timeAverage` prints a
warning and returns the input unchanged.

### Resampling at the same time

```python
averaged, outTimes = fds.timeAverage(values, times, window=10.0, outdt=5.0)
```

### A single window

```python
frame, t = fds.timeAverage(values, times, window=10.0, queryTime=60.0)
```

### Other filters

`kalmanFilter` smooths a 1-D series without the phase lag a boxcar
introduces, which suits a device trace better than a field:

```python
smoothed = fds.kalmanFilter(devices['"TC-1"'].values, Q=1e-5, R=0.25)
```

Raise `Q` to follow the measurements more closely, raise `R` to smooth
harder.

`timeAverage2` is a deprecated exponentially weighted variant; it emits a
`DeprecationWarning`. Use `timeAverage`.

## Averaging whole output files

The results of the above live in memory. To get an averaged field into
smokeview, the averaged data have to be written back out as output files
**and registered in a smokeview file** — smokeview will not display a
file it does not know about.

Both routines do this for you and return the path of the new smokeview
file to open.

### Slices

```python
outFiles, outQty, refFiles, newSmvFile = fds.slcfsTimeAverage(
    workingDir, chid, 'TEMPERATURE', dt=30, outDir='./averaged')

print('open this in smokeview:', newSmvFile)
```

Every `.sf` file holding the quantity is averaged, across all meshes.
`outdt` resamples the output to a coarser interval, which shrinks the
files considerably:

```python
outFiles, outQty, refFiles, newSmvFile = fds.slcfsTimeAverage(
    workingDir, chid, 'TEMPERATURE', dt=30, outdt=10, outDir='./averaged')
```

For a single file:

```python
outFile, outQty = fds.slcfTimeAverage(slcfFile, dt=30, outFile='avg.sf')
```

### Boundary files

```python
outFiles, outQty, refFiles, newSmvFile = fds.bndfsTimeAverage(
    workingDir, chid, 'WALL TEMPERATURE', dt=30, outDir='./averaged')
```

and for a single file:

```python
outFile, outQty, shortName, units = fds.bndfTimeAverage(
    bndfFile, dt=30, outFile='avg.bf')
```

The averaged quantity is named after the original plus the window, so
that it is distinguishable in smokeview. Pass `outQty` to choose the
name yourself.

`pyfdstools.examples.exampleBndfTimeAverage` is a runnable version of
this against the bundled case.

## Choosing a window

There is no universal answer, but:

* Long enough to average over the largest turbulent structures — for a
  compartment fire, tens of seconds rather than a few.
* Short compared with the timescale of what you are trying to see. A
  30 s window will flatten a flashover that develops in 20 s.
* When comparing against experimental data, match the averaging the
  instrument applied. A shielded thermocouple has a time constant of its
  own.

The FDS User's Guide discusses this in the context of comparing
predictions with measurements.
