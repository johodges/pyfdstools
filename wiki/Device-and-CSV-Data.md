# Device and CSV Data

FDS writes several csv files per case: `_devc.csv` for `&DEVC` output,
`_hrr.csv` for the global heat release rate budget, `_ctrl.csv` for
`&CTRL` states, and so on. `load_csv` reads any of them into a pandas
data frame.

## Reading

```python
import os
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')
chid = 'case001'

devices = fds.load_csv(workingDir, chid, '_devc')
print(devices.columns.tolist()[:4])
# ['Time', '"U"', '"burn"', '"con"']
```

The `suffix` argument selects which file:

| Suffix | File |
| --- | --- |
| `'_devc'` | device output (the default) |
| `'_hrr'` | heat release rate budget |
| `'_ctrl'` | control function states |
| `'_steps'` | timestep history |
| `'_cpu'` | per-routine timings |

```python
hrr = fds.load_csv(workingDir, chid, '_hrr')
print(hrr.columns.tolist()[:5])
# ['Time', 'HRR', 'HRR_OX', 'Q_RADI', 'Q_CONV']
```

Column names come straight from FDS, so device columns keep the quotes
FDS writes around them. Strip them if you prefer:

```python
devices.columns = [c.strip('"') for c in devices.columns]
```

## Why not pandas.read_csv?

FDS csv files have two header rows — units then labels — and they
contain non-ASCII characters in the unit row (`°C`, `kW/m²`). Some
outputs also have a trailing comma on every line. `load_csv` detects
the header length, repairs the encoding and drops the trailing field,
which `pd.read_csv` on its own does not.

If a file has an unusual header, pin it explicitly:

```python
devices = fds.load_csv(workingDir, chid, '_devc', labelRow=1)
```

`skipcols` drops columns by index, which is useful when a file has a
malformed one:

```python
devices = fds.load_csv(workingDir, chid, '_devc', skipcols=[5])
```

## Plotting a time history

```python
import matplotlib.pyplot as plt

hrr = fds.load_csv(workingDir, chid, '_hrr')

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(hrr['Time'], hrr['HRR'], linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('HRR (kW)')
fig.tight_layout()
fig.savefig('hrr.png')
```

## Selecting a group of devices

Devices are usually named systematically, so selecting a group is a
string match on the column names:

```python
thermocouples = [c for c in devices.columns if 'TC' in c]
peak = devices[thermocouples].max().max()
```

`maxValuePlot` and `maxValueCSV` plot and write a set of such groups:

```python
import numpy as np

times = devices['Time'].values
values = devices[thermocouples].values
names = [c.strip('"') for c in thermocouples]

fds.maxValuePlot(times, values, names, 'thermocouples.png',
                 vName='Temperature (C)')
fds.maxValueCSV(times, values, names, 'thermocouples')
```

## Smoothing a noisy trace

Device output from a turbulent simulation is noisy. Two options:

```python
# Kalman filter: no phase lag, tune with Q (process) and R (measurement)
smoothed = fds.kalmanFilter(devices['"TC-1"'].values, Q=1e-5, R=0.25)

# Boxcar average over a window, via timeAverage
values = devices[thermocouples].values.T[:, None, :]   # (N, 1, NT)
averaged, outTimes = fds.timeAverage(values, times, window=10.0)
```

See [Time Averaging](Time-Averaging).

## Two-zone reduction

A vertical array of thermocouples can be reduced to the equivalent
two-layer model from the FDS verification guide:

```python
import numpy as np

elevations = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
profile = devices[['"TC-%0.1f"' % z for z in elevations]].values[-1, :]

lower, upper, interface = fds.getTwoZone(elevations, profile)
print('lower layer %.1f C, upper layer %.1f C, interface at %.2f m'
      % (lower, upper, interface))
```

The elevations may be given in either direction; the ordering is
detected from the data.

## Building a slice from a grid of devices

If you laid out devices on a regular grid, they can be interpolated onto
a slice plane and written out as a slice file smokeview will display.
`pyfdstools/examples/make_slice_from_devc.py` does exactly this for a
radiometer array, and
`pyfdstools.examples.exampleAddGasPhaseHeatFluxSlice` is the same code as
a callable function.

## Adiabatic surface temperature

Heat flux gauge output can be converted to adiabatic surface
temperature:

```python
ghf = devices['"gauge-1"'].values          # kW/m2
ast = fds.astFromGhf(ghf, h=0.01, e=0.9, Tgauge=20.0)
```

`h` is the convective heat transfer coefficient in kW/m2-K, `e` the
gauge emissivity, and `Tgauge` the gauge temperature in Celsius. The
result is in Celsius, and equals `Tgauge` when the gauge flux is zero.
