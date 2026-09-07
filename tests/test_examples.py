"""Smoke tests for the example scripts bundled with pyfdstools.

These run the examples end to end against the bundled cases. They are
slow relative to the rest of the suite, so they are marked and can be
deselected with ``-m 'not slow'``.
"""

import os
import subprocess
import sys

import pytest

import pyfdstools as fds


@pytest.mark.slow
@pytest.mark.parametrize('script', fds.EXAMPLE_SCRIPTS)
def test_example_runs(script):
    examplesDir = fds.getExamplesDirectory()
    assert os.path.exists(os.path.join(examplesDir, script))

    env = os.environ.copy()
    env['MPLBACKEND'] = 'Agg'
    # Run against the working tree rather than any installed copy.
    repoRoot = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env['PYTHONPATH'] = repoRoot + os.pathsep + env.get('PYTHONPATH', '')

    # Not every example accepts --outdir, so let them write to their
    # default pyfdstools/examples/generated directory, which is ignored
    # by git.
    result = subprocess.run(
        [sys.executable, script],
        cwd=examplesDir, capture_output=True, text=True, env=env,
        timeout=600)

    assert result.returncode == 0, (
        "%s failed\nstdout:\n%s\nstderr:\n%s"
        % (script, result.stdout, result.stderr))
    assert 'Traceback' not in result.stderr


@pytest.mark.slow
def test_runExamples_reports_failures():
    """runExamples must surface a failing example rather than swallow it."""
    with pytest.raises(RuntimeError, match='failed'):
        fds.runExamples(scripts=['this_example_does_not_exist.py'],
                        raiseOnError=True)
