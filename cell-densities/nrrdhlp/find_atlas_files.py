import os
import re

from voxcell import VoxelData

from .mtypes import parse_mtype

fn_hierarchy = "hierarchy.json"
fn_regions = "brain_regions.nrrd"


def load_all_densities(root, density_kw="densities"):
    fns = [os.path.join(root, x) for x in os.listdir(root) if density_kw in x]

    mtypes = [parse_mtype(os.path.split(x)[1]) for x in fns]
    data = [(mtype, VoxelData.load_nrrd(fn)) for mtype, fn in zip(mtypes, fns)
            if mtype is not None]
    return dict(data)


def find_hierarchy(root, format="path"):
    expected_fn_hierarchy = os.path.join(root, fn_hierarchy)
    if not os.path.isfile(expected_fn_hierarchy):
        return None
    if format == "path":
        return expected_fn_hierarchy
    if format == "json":
        with open(expected_fn_hierarchy, 'r') as fid:
            import json
            hier_json = json.load(fid)
        if 'msg' in hier_json:
            hier_json = hier_json['msg'][0]
        return hier_json
    if format == "voxcell":
        from voxcell.nexus.voxelbrain import RegionMap
        return RegionMap.load_json(expected_fn_hierarchy)
    raise ValueError("Unknown format spec: {0}".format(format))


def find_regions(root, format="path"):
    expected_fn_regions = os.path.join(root, fn_regions)
    if not os.path.isfile(expected_fn_regions):
        return None
    if format == "path":
        return expected_fn_regions
    elif format == "voxcell":
        from voxcell import VoxelData
        return VoxelData.load_nrrd(expected_fn_regions)
    raise ValueError("Unknown format spec: {0}".format(format))
