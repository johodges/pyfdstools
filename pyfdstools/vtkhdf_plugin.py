import vtk
from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import os

from .fdsFileOperations import fdsFileOperations
from .utilities import getDatatypeByEndianness, getEndianness
from .utilities import getFileListFromZip, getFileList, zopen, zreadlines
from .smokeviewParser import parseSMVFile

def round_if_needed(arr, tol):
    if tol is None or tol <= 0:
        return arr
    return np.round(arr / tol) * tol


def get_point_array(part, array_name):
    arr = part.GetPointData().GetArray(array_name)
    if arr is None:
        raise ValueError(f"Point-data array '{array_name}' not found in partition.")
    return vtk_to_numpy(arr)


def partition_to_structured_block(part, array_name, tol=1e-10, plane=None):
    """
    Convert one partition into a structured 2D block.

    Returns
    -------
    block : dict
        {
            "plane": "XY" | "XZ" | "YZ",
            "coord_names": ("x","y") or similar,
            "normal_name": "z" or similar,
            "x": 2D grid of first in-plane coord,
            "z": 2D grid of second in-plane coord,
            "datas": 2D scalar values,
            "normal_value": scalar approximately constant for this partition,
        }
    """
    pts = vtk_to_numpy(part.GetPoints().GetData())
    vals = get_point_array(part, array_name)

    if plane is None:
        plane, a0, a1, an = detect_surface_plane_from_points(pts, tol=tol)
    else:
        if plane == 3:
            a0, a1, an = 0, 1, 2
        elif plane == 2:
            a0, a1, an = 0, 2, 1
        elif plane == 1:
            a0, a1, an = 1, 2, 0
        else:
            raise ValueError(f"Unsupported plane '{plane}'")

    c0 = pts[:, a0]
    c1 = pts[:, a1]
    cn = pts[:, an]

    c0r = round_if_needed(c0, tol)
    c1r = round_if_needed(c1, tol)

    u0 = np.unique(c0r)
    u1 = np.unique(c1r)

    n0 = len(u0)
    n1 = len(u1)

    if n0 * n1 != len(vals):
        raise ValueError(
            f"Partition is not a complete structured block for plane {plane}: "
            f"n0={n0}, n1={n1}, n0*n1={n0*n1}, npts={len(vals)}"
        )

    i0 = {v: i for i, v in enumerate(u0)}
    i1 = {v: i for i, v in enumerate(u1)}

    V = np.full((n1, n0), np.nan, dtype=vals.dtype)

    for cc0, cc1, vv in zip(c0r, c1r, vals):
        j = i1[cc1]
        i = i0[cc0]
        V[j, i] = vv

    X, Y = np.meshgrid(u0, u1, indexing="xy")

    coord_labels = {
        3: (("x", "y"), "z"),
        2: (("x", "z"), "y"),
        1: (("y", "z"), "x"),
    }

    coord_names, normal_name = coord_labels[plane]

    return {
        "plane": plane,
        "coord_names": coord_names,
        "normal_name": normal_name,
        "x": X,
        "z": Y,
        "datas": V,
        "normal_value": float(np.mean(cn)),
    }

def cached_surface_to_structured_blocks(surface, array_name, tol=1e-10, plane=None):
    """
    Convert every partition in a cached surface into structured 2D blocks.
    If plane is None, it is auto-detected from the first non-empty partition.
    """
    blocks = []
    detected_plane = plane

    for p in range(surface.GetNumberOfPartitions()):
        part = surface.GetPartition(p)
        if part is None:
            continue

        pts = vtk_to_numpy(part.GetPoints().GetData())
        if pts.shape[0] == 0:
            continue

        if detected_plane is None:
            detected_plane, *_ = detect_surface_plane_from_points(pts, tol=tol)

        block = partition_to_structured_block(
            part,
            array_name=array_name,
            tol=tol,
            plane=detected_plane,
        )
        block["partition"] = p
        blocks.append(block)

    if not blocks:
        raise ValueError("No non-empty partitions found.")

    return blocks


def stitch_structured_blocks(blocks, tol=1e-10):
    """
    Stitch structured partition blocks into one global 2D mesh.

    Missing regions are filled with NaN.
    """
    if not blocks:
        raise ValueError("No blocks to stitch.")

    plane = blocks[0]["plane"]
    coord_names = blocks[0]["coord_names"]
    normal_name = blocks[0]["normal_name"]

    for b in blocks[1:]:
        if b["plane"] != plane:
            raise ValueError("Blocks have inconsistent detected planes.")

    all_c0 = np.concatenate([b["x"][0, :] for b in blocks])
    all_c1 = np.concatenate([b["z"][:, 0] for b in blocks])

    all_c0 = round_if_needed(all_c0, tol)
    all_c1 = round_if_needed(all_c1, tol)

    g0 = np.unique(all_c0)
    g1 = np.unique(all_c1)

    n0 = len(g0)
    n1 = len(g1)

    map0 = {v: i for i, v in enumerate(g0)}
    map1 = {v: i for i, v in enumerate(g1)}

    dtype = np.result_type(*[b["datas"].dtype for b in blocks])
    Vg = np.full((n1, n0), np.nan, dtype=dtype)

    for b in blocks:
        X = round_if_needed(b["x"], tol)
        Y = round_if_needed(b["z"], tol)
        V = b["datas"]

        for j in range(V.shape[0]):
            for i in range(V.shape[1]):
                jj = map1[Y[j, i]]
                ii = map0[X[j, i]]
                Vg[jj, ii] = V[j, i]

    Xg, Yg = np.meshgrid(g0, g1, indexing="xy")

    return {
        "plane": plane,
        "coord_names": coord_names,
        "normal_name": normal_name,
        "x": Xg,
        "z": Yg,
        "datas": Vg,
        "normal_value": float(np.mean([b["normal_value"] for b in blocks])),
    }


def cached_surface_to_structured_mesh(surface, array_name, tol=1e-10, plane=None):
    """
    Convenience wrapper: convert a cached surface into one stitched
    structured mesh for contourf.
    """
    blocks = cached_surface_to_structured_blocks(
        surface,
        array_name=array_name,
        tol=tol,
        plane=plane,
    )
    return stitch_structured_blocks(blocks, tol=tol)


def get_step_times(reader):
    info = reader.GetOutputInformation(0)

    key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()

    if info.Has(key):
        n = info.Length(key)
        return np.array([info.Get(key, i) for i in range(n)])

    return None

def query2dAxisValue_vtkhdf(workingDir, chid, quantity, axis, value, time=None, dt=None):
    endianness = getEndianness(workingDir, chid)
    datatype = getDatatypeByEndianness(np.float32, endianness)
    
    if '.zip' in workingDir:
        smvFile = getFileListFromZip(workingDir, chid, 'smv')[0]
    else:
        smvFile = getFileList(workingDir, chid, 'smv')[0]
    smvOutputs = parseSMVFile(smvFile)
    
    smv_grids = smvOutputs['grids']
    smv_slcf = smvOutputs['files']['SLICES']
    
    resultDir = workingDir
    if '.zip' not in workingDir:
        fdsFileName = getFileList(workingDir, chid, 'fds')[0]
        fdsFile = fdsFileOperations()
        fdsFile.importFile(fdsFileName)
        if fdsFile.dump['ID'] is not False:
            if fdsFile.dump['ID']['RESULTS_DIR'] is not False:
                resultDir = workingDir + os.sep + fdsFile.dump['ID']['RESULTS_DIR'] + os.sep
    
    fn = resultDir + os.sep + chid + "_Z_pos_%08d.vtkhdf"%(value*100)
    
    #colors = vtkNamedColors()
    
    # Set up reader
    reader = vtkHDFReader()
    reader.SetFileName(fn)
    reader.UseCacheOn()
    reader.Update()
    
    # Get queried time
    timesteps = get_step_times(reader)
    if time is not None:
        step = np.argmin(abs(timesteps-time))
    else:
        step = -1
    reader.SetStep(step)
    reader.Modified()
    
    # Set up surface
    surface = vtkDataSetSurfaceFilter()
    surface.SetInputConnection(reader.GetOutputPort())
    surface.SetNonlinearSubdivisionLevel(0)
    #surface.FastModeOn()
    surface.Update()
    output = surface.GetOutputDataObject(0)
    #cached_surface = clone_data_object(output)
    
    # Convert to numpy array
    mesh = cached_surface_to_structured_mesh(
        output,
        array_name=quantity,
        tol=1e-10,
        plane=axis,
    )
    return mesh
