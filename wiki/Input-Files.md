# Input Files

`fdsFileOperations` reads an FDS input file into Python objects, lets you
edit them, and writes the model back out. That makes parametric studies,
bulk edits and generated geometry a matter of a few lines rather than a
text-templating exercise.

## Reading

```python
import os
import pyfdstools as fds

path = os.path.join(os.path.dirname(fds.__file__),
                    'examples', 'data', 'case001.fds')

model = fds.fdsFileOperations()
model.importFile(path)

print(model.head['ID']['CHID'])                 # case001
print([k for k in model.meshes if k != 'unknownCounter'])   # ['MESH-001']
print(list(model.surfs.keys()))                 # ['UPHOLSTERY', 'WALL', 'ignitor']
```

You can also parse text you already have:

```python
model.importFile(text=fdsText)
model.importFile(textList=listOfLines)
```

## How the model is stored

Each namelist group is an attribute holding a dictionary keyed by the
entry's `ID`:

| Attribute | Namelist |
| --- | --- |
| `head`, `time`, `misc`, `dump`, `pres`, `comb`, `clip` | single-entry groups, keyed `'ID'` |
| `meshes`, `surfs`, `matls`, `obsts`, `holes`, `vents`, `devcs`, `ctrls`, `props`, `parts`, `specs`, `reacs`, `ramps`, `radis`, `winds`, `slcfs`, `bndfs`, `isof`, `profs`, `inits`, `zones`, `geom`, `hvac`, `mult`, `move`, `sm3d`, `trnx`, `trny`, `trnz`, `tabl`, `catf` | multi-entry groups |

The spelling is not uniform — some of these are plural and some are not.
`sorted(vars(model))` lists them all for the version you have.

```python
mesh = model.meshes['MESH-001']
print(mesh['IJK'])   # [24, 12, 20]
print(mesh['XB'])    # [0.0, 4.8, 2.0, 4.4, 0.0, 2.4]
```

### The bookkeeping keys

Multi-entry groups carry an extra `'unknownCounter'` entry, used to name
namelists that have no `ID` of their own. **Skip it when iterating:**

```python
for key in model.meshes:
    if key == 'unknownCounter':
        continue
    print(key, model.meshes[key]['XB'])
```

Device keys additionally get a `DEVICE-nnnnnn-` prefix when the model is
written out, which imposes the ordering FDS requires. The prefix is
internal; the `ID` written to the file comes from the entry's own `ID`.

## Editing

```python
# retarget a heat release rate
model.surfs['UPHOLSTERY']['HRRPUA'] = 750.0

# move an obstruction
for key in model.obsts:
    if key == 'unknownCounter':
        continue
    model.obsts[key]['XB'][4] += 0.1

# extend the run
model.time['ID']['T_END'] = 600.0
```

## Building from scratch

```python
model = fds.fdsFileOperations()
model.addHEAD('compartment', title='Single compartment test')
model.addTIME(T_END=600.0)
model.addMESH('MESH-1', [40, 40, 40], [0.0, 4.0, 0.0, 4.0, 0.0, 4.0])
model.addSURF('BURNER', Hrrpua=1000.0)
model.addVENT('FIRE', 'BURNER', XB=[1.5, 2.5, 1.5, 2.5, 0.0, 0.0])
model.addOBST('TABLE', [1.0, 3.0, 1.0, 3.0, 0.7, 0.8], SURF_ID='INERT')
model.addDEVC('TC-1', 'TEMPERATURE', XYZ=[2.0, 2.0, 2.0])
model.addSLCF('TEMPERATURE', PBY=2.0)
model.addBNDF('WALL TEMPERATURE')

model.saveModel(1, 'compartment.fds')
```

The `addXXXX` methods cover the common namelists; anything they do not
cover can be added by writing into the dictionary directly, or through
`customLines`:

```python
model.customLines.append("&PROP ID='sprinkler', QUANTITY='SPRINKLER LINK TEMPERATURE' /")
```

## Writing

```python
model.saveModel(mpiProcesses=1, location='out.fds')
```

`mpiProcesses` distributes the meshes across MPI processes by setting
`MPI_PROCESS` on each `&MESH`.

To get the text without writing a file:

```python
text = model.generateFDStext(precision=6)
```

`precision` is the number of decimal places used for real values.
Generating text does not modify the model, and generating it twice gives
the same result.

## Inspecting a model

```python
meshes, cells = model.calculateMeshCells()
print(sum(cells), 'cells across', len(meshes), 'meshes')

if model.checkOverlappingMESH():
    print('warning: meshes overlap')
```

`checkOverlappingMESH` shrinks each mesh slightly before testing, so
meshes that merely share a face are not reported. A mesh whose name
contains `east`, `west`, `north` or `south` is instead grown on that
side, so a mesh which should abut its neighbour there but does not is
reported.

## Parametric studies

The pattern the class exists for:

```python
import copy
import os

base = fds.fdsFileOperations()
base.importFile('base_case.fds')

for hrrpua in [250.0, 500.0, 750.0, 1000.0]:
    model = copy.deepcopy(base)
    chid = 'case_hrrpua_%04d' % (hrrpua)
    model.head['ID']['CHID'] = chid
    model.surfs['BURNER']['HRRPUA'] = hrrpua
    model.saveModel(4, os.path.join('cases', '%s.fds' % (chid)))
```

Copy the base model rather than re-importing it each time only if the
import is slow; `deepcopy` is safe here because the model is plain
dictionaries and lists.

## Adding devices from data

`pyfdstools.examples.exampleAddOccupantFedDevices` reads occupant
positions from a csv and adds a fractional effective dose device at each
one:

```python
for i, (x, y, z) in enumerate(positions):
    model.addDEVC('FED-%04d' % (i), 'FED', XYZ=[x, y, z])
```

The same pattern covers thermocouple trees, radiometer arrays and
anything else laid out from data.

## Namelist datatypes

`fdsLineTypes` declares the datatype of every parameter of every
namelist, which is how the parser knows that `IJK` is a list of integers
and `SURF_ID` a quoted string:

```python
types = fds.fdsLineTypes(version='6.7.4')
print(types.surf['HRRPUA'])    # float
print(types.mesh['IJK'])       # listint
```

A parameter absent from these tables is still read, but is treated as a
string and so written back out quoted. If a numeric parameter of yours
round-trips with quotes around it, add it to the relevant
`getXXXXtypes` method in `pyfdstools/fdsTypes.py`.

Pass `version=` to `fdsFileOperations` for a case written for an older
FDS release:

```python
model = fds.fdsFileOperations(version='6.7.1')
```

## Verifying an input file

`checkDevices` reports any `&DEVC` whose `XYZ` falls inside an
obstruction, where FDS would record the solid rather than the gas:

```python
model = fds.fdsFileOperations()
model.importFile('mycase.fds')

buried = fds.checkDevices(model, smvFile='mycase.smv')
for name in buried:
    print('device inside an obstruction:', name)
```

It needs a smokeview file from a completed run, because FDS snaps
obstructions to the mesh and it is the snapped extents that matter.

## The command line example

```bash
python read_and_write_input_files.py \
    --path data/case001.fds --outdir generated --precision 4
```
