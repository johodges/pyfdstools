# Plotting

`plotSlice` renders any 2-D field on a coordinate grid — a slice, a
boundary plane, a plane taken out of a plot3D snapshot. It returns the
matplotlib figure and axes, so anything it does not offer you can do
afterwards.

## The basics

```python
import os
import matplotlib.pyplot as plt
import pyfdstools as fds

workingDir = os.path.join(
    os.path.dirname(fds.__file__), 'examples', 'data', 'case001.zip')

data, units = fds.query2dAxisValue(
    workingDir, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)

fig, ax = fds.plotSlice(
    data['x'], data['z'], data['datas'][:, :, -1], axis=1,
    clabel='Temperature (%s)' % (units))

fig.savefig('slice.png', dpi=300)
plt.show()
```

The `axis` argument is only used to label the figure: for `axis=1` the
horizontal axis is labelled `y (m)` and the vertical `z (m)`. Pass
`xlabel` and `zlabel` to override.

## Fixing the color scale

Leaving the limits out scales each frame to its own range, which makes
an animation flicker. Fix them when comparing frames or cases:

```python
fig, ax = fds.plotSlice(
    data['x'], data['z'], data['datas'][:, :, -1], 1,
    qnty_mn=20, qnty_mx=1000,
    cbarnumticks=11,
    clabel='Temperature (C)')
```

| Argument | Effect |
| --- | --- |
| `qnty_mn`, `qnty_mx` | limits of the color scale |
| `cbarnumticks` | number of evenly spaced colorbar ticks |
| `cbarticks` | explicit tick locations, overrides `cbarnumticks` |
| `tickDecimals` | decimal places on the tick labels |
| `levels` | number of contour levels, or explicit level values |

```python
fig, ax = fds.plotSlice(
    data['x'], data['z'], data['datas'][:, :, -1], 1,
    qnty_mn=20, qnty_mx=1000,
    cbarticks=[20, 100, 200, 400, 600, 800, 1000],
    tickDecimals=0)
```

## Extending the scale

`extend` says which end of the scale continues past the limits, and is
drawn as an arrow on the colorbar.

| Value | Meaning |
| --- | --- |
| `'both'` | both ends (the default) |
| `'below'` or `'min'` | the low end only |
| `'above'` or `'max'` | the high end only |
| `'neither'` | neither |

`'below'` and `'above'` are pyfdstools spellings of matplotlib's `'min'`
and `'max'`; both work.

## Highlighting a threshold

To make one value stand out — a tenability criterion, an ignition
temperature — `highlightValue` draws a black band across the colormap at
that value:

```python
fig, ax = fds.plotSlice(
    data['x'], data['z'], data['datas'][:, :, -1], 1,
    qnty_mn=20, qnty_mx=1000,
    highlightValue=200, highlightWidth=3,
    clabel='Temperature (C)')
```

A colorbar tick is added at the highlighted value and any tick too close
to it is dropped. `highlightWidth` is the half-width of the band in
colormap entries and defaults to 2.

## Colormaps

The default is the smokeview blue-cyan-green-yellow-red ramp, so that
figures match what smokeview shows. Any matplotlib colormap works too:

```python
x, z, frame = data['x'], data['z'], data['datas'][:, :, -1]

fds.plotSlice(x, z, frame, 1, cmap='viridis')
fds.plotSlice(x, z, frame, 1, cmap='SMV')       # the default, named
fds.plotSlice(x, z, frame, 1, cmap=fds.buildSMVcolormap())
```

For line plots, two categorical sequences are provided:

```python
colors = fds.getVTcolors()    # 13 colors
colors = fds.getJHcolors()    # 16 colors
colors = fds.getPlotColors(25)  # generate as many as you need
```

## Contours or an image

`plotSlice` draws filled contours by default. Two alternatives:

```python
fds.plotSlice(x, z, frame, 1, linecontour=True)   # line contours
fds.plotSlice(x, z, frame, 1, contour=False)      # imshow, no interpolation
```

The image form shows the cells as they are, which is what you want when
checking mesh resolution; the contour form is smoother for presentation.
With `contour=False` and `extend` set to one end, values past the other
end are masked out.

## Layout

| Argument | Effect |
| --- | --- |
| `figsize` | figure size in inches; derived from the slice aspect ratio when omitted |
| `figsizeMult` | length of the shorter figure axis, default 4 |
| `fs` | font size, default 16 |
| `title` | axes title |
| `reverseXY` | swap the horizontal and vertical axes |
| `xmn`, `xmx`, `zmn`, `zmx` | axis limits; data extent when omitted |
| `fixXLims`, `fixZLims` | whether those limits are applied |
| `addCbar` | draw the colorbar at all |

## Drawing several slices on one figure

Pass an existing figure and axes:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

for ax, time in zip(axes, [30, 60, 90]):
    data, units = fds.query2dAxisValue(
        workingDir, 'case001', 'TEMPERATURE', 1, 2.55, time=time, dt=10)
    fds.plotSlice(data['x'], data['z'], data['datas'][:, :, 0], 1,
                  fig=fig, ax=ax, qnty_mn=20, qnty_mx=1000,
                  title='t = %d s' % (time),
                  addCbar=(ax is axes[-1]),
                  clabel='Temperature (C)')

fig.savefig('slices.png', dpi=300)
```

## Getting the mappable

For a colorbar you place yourself, or an animation that updates the
image in place:

```python
fig, ax, im = fds.plotSlice(x, z, frame, 1, returnIm=True, addCbar=False)
cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
```

## Animations

```python
import matplotlib.animation as animation

data, units = fds.query2dAxisValue(workingDir, 'case001', 'TEMPERATURE', 1, 2.55)

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
frames = []
for i in range(0, data['datas'].shape[2]):
    fds.plotSlice(data['x'], data['z'], data['datas'][:, :, i], 1,
                  fig=fig, ax=ax, qnty_mn=20, qnty_mx=1000,
                  addCbar=(i == 0),
                  title='t = %.0f s' % (data['times'][i]))
    fig.savefig('frame_%04d.png' % (i), dpi=150)
    ax.clear()
```

Keep `qnty_mn` and `qnty_mx` fixed across frames, or the colors will
jump from frame to frame.

## Running headless

On a machine with no display — CI, a compute node — select the Agg
backend before importing pyplot:

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
```

or set `MPLBACKEND=Agg` in the environment. `plt.show()` then does
nothing and `fig.savefig` still works.

## 3-D geometry

`smvVisual` renders a case's obstructions as a 3-D figure:

```python
surfaces, obstructions = fds.buildSMVgeometry(smvFile)
fig, ax = fds.smvVisual(obstructions, surfaces, 'mycase',
                        limits=[0, 15, 0, 8, 0, 5])
```

Faces are colored from the surface definitions in the smokeview file.
For anything interactive, export to ParaView instead.
