"""Tests for reading, editing and writing FDS input files."""

import os
import warnings

import numpy as np
import pytest

import pyfdstools as fds


@pytest.fixture()
def case001_fds(data_dir):
    return os.path.join(data_dir, 'case001.fds')


def test_import_file(case001_fds):
    model = fds.fdsFileOperations()
    model.importFile(case001_fds)
    assert model.head['ID']['CHID'] == 'case001'
    assert len(model.meshes) > 0
    assert len(model.surfs) > 0


def test_round_trip_preserves_key_namelists(case001_fds, outdir):
    """Writing a model back out and re-reading it must not lose data."""
    original = fds.fdsFileOperations()
    original.importFile(case001_fds)

    out = os.path.join(outdir, 'roundtrip.fds')
    original.saveModel(1, out)
    assert os.path.exists(out)

    reread = fds.fdsFileOperations()
    reread.importFile(out)

    assert reread.head['ID']['CHID'] == original.head['ID']['CHID']
    for group in ('meshes', 'surfs', 'obsts', 'vents', 'devcs', 'reacs'):
        assert (len(getattr(reread, group))
                == len(getattr(original, group))), group


def test_round_trip_preserves_mesh_geometry(case001_fds, outdir):
    original = fds.fdsFileOperations()
    original.importFile(case001_fds)
    out = os.path.join(outdir, 'roundtrip.fds')
    original.saveModel(1, out)
    reread = fds.fdsFileOperations()
    reread.importFile(out)

    for key in original.meshes:
        if key == 'unknownCounter':
            continue
        assert np.allclose(np.array(original.meshes[key]['XB'], dtype=float),
                           np.array(reread.meshes[key]['XB'], dtype=float))
        assert (list(original.meshes[key]['IJK'])
                == list(reread.meshes[key]['IJK']))


def test_generateFDStext_is_valid_input(case001_fds):
    model = fds.fdsFileOperations()
    model.importFile(case001_fds)
    text = model.generateFDStext(precision=6)
    assert '&HEAD' in text
    assert '&MESH' in text
    assert '&SURF' in text
    # Every namelist opened must be closed.
    assert text.count('&') == text.count('/')


def test_build_a_model_from_scratch(outdir):
    model = fds.fdsFileOperations()
    model.addHEAD('synthetic', title='built by the test suite')
    model.addTIME(T_END=10.0)
    model.addMESH('MESH-1', [10, 10, 10], [0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
    model.addSURF('BURNER', Hrrpua=500.0)
    model.addVENT('BURNER-VENT', 'BURNER',
                  XB=[0.4, 0.6, 0.4, 0.6, 0.0, 0.0])
    model.addDEVC('TC-1', 'TEMPERATURE', XYZ=[0.5, 0.5, 0.5])

    out = os.path.join(outdir, 'synthetic.fds')
    model.saveModel(1, out)

    reread = fds.fdsFileOperations()
    reread.importFile(out)
    assert reread.head['ID']['CHID'] == 'synthetic'
    assert 'MESH-1' in reread.meshes
    assert 'BURNER' in reread.surfs
    assert reread.time['ID']['T_END'] == 10.0


def test_addMESH_records_geometry():
    model = fds.fdsFileOperations()
    model.addMESH('M1', [4, 5, 6], [0.0, 1.0, 0.0, 2.0, 0.0, 3.0])
    assert list(model.meshes['M1']['IJK']) == [4, 5, 6]
    assert list(model.meshes['M1']['XB']) == [0.0, 1.0, 0.0, 2.0, 0.0, 3.0]


def test_calculateMeshCells():
    model = fds.fdsFileOperations()
    model.addMESH('M1', [10, 10, 10], [0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
    model.addMESH('M2', [5, 5, 5], [1.0, 2.0, 0.0, 1.0, 0.0, 1.0])
    meshes, cells = model.calculateMeshCells()
    assert sorted(cells) == [125, 1000]


def test_import_from_text(case001_fds):
    with open(case001_fds, 'r') as f:
        text = f.read()
    fromText = fds.fdsFileOperations()
    fromText.importFile(text=text)
    fromPath = fds.fdsFileOperations()
    fromPath.importFile(case001_fds)
    assert fromText.head['ID']['CHID'] == fromPath.head['ID']['CHID']
    assert len(fromText.meshes) == len(fromPath.meshes)


def test_generateFDStext_does_not_mutate_the_model(case001_fds):
    """Generating text must leave the model it reads unchanged.

    generateFDStext used to subscript each namelist defaultdict for a
    'newline' key, which inserted a phantom entry into every collection
    it touched.
    """

    model = fds.fdsFileOperations()
    model.importFile(case001_fds)
    before = {group: set(getattr(model, group).keys())
              for group in ('surfs', 'obsts', 'vents', 'matls')}
    model.generateFDStext(precision=6)
    after = {group: set(getattr(model, group).keys())
             for group in before}
    assert before == after


def test_generateFDStext_is_idempotent(case001_fds):
    """Repeated calls must produce the same text and the same model.

    sortDEVCs renames the device keys to impose the output ordering FDS
    requires. It used to append a fresh prefix on every call, so the
    keys grew without bound and the model drifted each time text was
    generated.
    """

    model = fds.fdsFileOperations()
    model.importFile(case001_fds)

    first = model.generateFDStext(precision=6)
    keysAfterFirst = list(model.devcs.keys())
    second = model.generateFDStext(precision=6)
    keysAfterSecond = list(model.devcs.keys())

    assert first == second
    assert keysAfterFirst == keysAfterSecond
    assert max(len(k) for k in keysAfterSecond) < 80


def test_sortDEVCs_puts_aspiration_devices_last():
    """FDS requires an aspiration device after the devices it samples."""
    model = fds.fdsFileOperations()
    model.addDEVC('SAMPLE-1', 'TEMPERATURE', XYZ=[0.5, 0.5, 0.5])
    model.addDEVC('ASPIRATOR', 'ASPIRATION', XYZ=[0.5, 0.5, 0.5])
    model.addDEVC('SAMPLE-2', 'TEMPERATURE', XYZ=[0.5, 0.5, 0.6])
    model.sortDEVCs()

    ordered = [k for k in sorted(model.devcs.keys())
               if k.startswith('DEVICE-')]
    quantities = [model.devcs[k]['QUANTITY'] for k in ordered]
    assert quantities.index('ASPIRATION') == len(quantities) - 1


# ---------------------------------------------------------------------
# Parameters newer than the fdsTypes tables
# ---------------------------------------------------------------------

def test_unknown_parameter_does_not_drop_the_namelist_line(outdir):
    """An undeclared parameter must not take its whole line with it.

    interpretKey returns False as the type, which made dictFromLine
    raise TypeError on "'list' in False"; the bare except in parseLine
    swallowed that and discarded the entire namelist line, leaving only
    "WARNING: Unknown line in input file".
    """

    src = os.path.join(outdir, 'unknown.fds')
    with open(src, 'w') as f:
        f.write("&HEAD CHID='unknown' /\n")
        f.write("&MESH ID='M1', IJK=10,10,10, XB=0.0,1.0,0.0,1.0,0.0,1.0 /\n")
        f.write("&PRES CHECK_POISSON=.TRUE., NOT_YET_IN_FDSTYPES=.TRUE. /\n")
        f.write("&TIME T_END=1.0 /\n")

    model = fds.fdsFileOperations()
    with pytest.warns(UserWarning, match='NOT_YET_IN_FDSTYPES'):
        model.importFile(src)

    # The declared parameter on the same line survives.
    assert model.pres['ID']['CHECK_POISSON'] in (True, 'TRUE', '.TRUE.')


@pytest.mark.parametrize('parameter,value', [
    ('NEW_LOGICAL', '.TRUE.'),
    ('NEW_FLOAT', '1.25'),
    ('NEW_STRING', "'abc'"),
])
def test_unknown_parameter_round_trips_verbatim(parameter, value, outdir):
    """Whatever was written is what comes back out."""
    src = os.path.join(outdir, 'roundtrip.fds')
    with open(src, 'w') as f:
        f.write("&HEAD CHID='rt' /\n")
        f.write("&MESH ID='M1', IJK=10,10,10, XB=0.0,1.0,0.0,1.0,0.0,1.0 /\n")
        f.write("&SURF ID='S1', HRRPUA=500.0, %s=%s /\n" % (parameter, value))
        f.write("&TIME T_END=1.0 /\n")

    model = fds.fdsFileOperations()
    with pytest.warns(UserWarning):
        model.importFile(src)

    out = os.path.join(outdir, 'out.fds')
    model.saveModel(1, out)

    surfLine = [l for l in open(out) if l.startswith('&SURF')]
    assert len(surfLine) == 1
    assert '%s=%s' % (parameter, value) in surfLine[0]
    # The declared parameter on the same line is unaffected.
    assert 'HRRPUA=500' in surfLine[0]


def test_unknown_parameter_is_reported_once_per_namelist(outdir):
    """One warning per parameter, not one per line that uses it."""
    src = os.path.join(outdir, 'repeated.fds')
    with open(src, 'w') as f:
        f.write("&HEAD CHID='repeated' /\n")
        f.write("&MESH ID='M1', IJK=10,10,10, XB=0.0,1.0,0.0,1.0,0.0,1.0 /\n")
        for i in range(0, 5):
            f.write("&SURF ID='S%d', FUTURE_PARAM=1.0 /\n" % (i))
        f.write("&TIME T_END=1.0 /\n")

    model = fds.fdsFileOperations()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.importFile(src)
    matching = [w for w in caught if 'FUTURE_PARAM' in str(w.message)]
    assert len(matching) == 1
    assert ('SURF', 'FUTURE_PARAM') in model.unknownParameters


def test_importFile_closes_the_file(case001_fds):
    """importFile leaked a handle per model imported."""
    model = fds.fdsFileOperations()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.importFile(case001_fds)
    assert [w for w in caught if 'unclosed' in str(w.message)] == []
