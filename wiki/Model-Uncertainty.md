# Model Uncertainty

The FDS Validation Guide publishes, for each validated quantity, a model
**bias** and a relative standard deviation derived from comparing FDS
against a body of experiments. pyfdstools ships those tables and applies
them, which turns a single predicted value into a distribution.

This is the basis of the probabilistic approach described in the guide:
rather than asking "does the prediction exceed the criterion", you ask
"what value should I use so that the criterion is met with 95 %
confidence".

## What is tabulated

```python
import pyfdstools as fds

print(fds.getQuantities())
```

```
['HGL Temperature, Forced Ventilation',
 'HGL Temperature, Natural Ventilation',
 'HGL Temperature, No Ventilation',
 ...]
```

29 quantities are tabulated for FDS 6.7.1. Three versions ship:

```python
fds.getQuantities(fdsVersion='6.2.0')
fds.getQuantities(fdsVersion='6.7.1')   # the default
fds.getQuantities(fdsVersion='6.7.4')
```

Asking for a version that is not shipped tells you which are:

```
FileNotFoundError: No FDS error table for version 9.9.9.
Available versions: 6.2.0, 6.7.1, 6.7.4
```

The tables themselves are csv files under `pyfdstools/fdsErrorTables`;
read one directly if you want the bias and standard deviation:

```python
table, quantities = fds.readErrorTable(fdsVersion='6.7.1')
print(table.loc['HGL Temperature, Natural Ventilation'][['Bias', 'sigmaM']])
```

## Applying the uncertainty

```python
predicted = [100.0, 250.0, 400.0]

adjusted = fds.calculatePercentile(
    predicted,
    'HGL Temperature, Natural Ventilation',
    percentile=0.95)

print(adjusted)   # [112.47, 281.17, 449.88]
```

Each predicted value is divided by the bias to correct the systematic
offset, then treated as the mean of a normal distribution whose standard
deviation is the mean times the published relative standard deviation.
The value at the requested percentile of that distribution is returned.

Asking for an untabulated quantity raises with the list of valid names:

```python
fds.calculatePercentile([100.0], 'NOT A QUANTITY', 0.95)
```

```
ValueError: Quantity 'NOT A QUANTITY' is not in the FDS 6.7.1 error table.
Known quantities: HGL Temperature, Forced Ventilation, ...
```

## Applying it to a field

`calculatePercentile` takes an array, so a whole slice can be adjusted
by flattening and reshaping:

```python
import numpy as np

data, units = fds.query2dAxisValue(
    workingDir, chid, 'TEMPERATURE', 1, 2.55, time=60, dt=30)

frame = data['datas'][:, :, 0]
mask = np.isfinite(frame)

adjusted = np.full_like(frame, np.nan)
adjusted[mask] = fds.calculatePercentile(
    frame[mask], 'HGL Temperature, Natural Ventilation', 0.95)

fig, ax = fds.plotSlice(data['x'], data['z'], adjusted, 1,
                        clabel='95th percentile temperature (C)')
```

Whether that is meaningful is a judgement call: the tabulated
uncertainty was derived for a specific measurement — hot gas layer
temperature, in this case — and applying it point by point to a field
assumes the same uncertainty holds everywhere. Read the Validation
Guide's discussion before relying on it.

## Visualising the distribution

```python
import matplotlib.pyplot as plt

fig, ax = fds.plotPercentile(
    250.0, 'HGL Temperature, Natural Ventilation')
fig.savefig('uncertainty.png', dpi=300)
```

Plots the probability density and the cumulative distribution around the
prediction, with the predicted value marked. `colors=` overrides the two
curve colors.

## Worked example

`pyfdstools/examples/error_calculation.py` runs the whole thing from the
command line:

```bash
python error_calculation.py \
    --values 100 200 300 \
    --quantity "HGL Temperature, Natural Ventilation" \
    --percentile 0.95
```

## Reference

Overholt, K. J., and McGrattan, K. B., *Fire Dynamics Simulator
Technical Reference Guide, Volume 3: Validation*, NIST Special
Publication 1018-3. The bias and relative standard deviation tables are
in the chapter on model uncertainty.
