import os
import numpy


lst_broken_regions = ['VISrll', 'VISlla', 'VISmma', 'VISmmp', 'VISm']  # REMOVE THIS ONCE HIERARCHY IS FIXED!


def find_root(j, root_acronym):
    if j['acronym'] == root_acronym:
        return j
    if 'children' in j:
        for child in j['children']:
            cand = find_root(child, root_acronym)
            if cand is not None:
                return cand
    return None


def find_almost_leaves(j, target_lvl, conjunction_func=numpy.max, partition=False):
    if 'children' not in j or len(j['children']) == 0:
        lvl = 0
        lst_ids = []
    else:
        res_lvl, lst_ids = zip(*[find_almost_leaves(child, target_lvl, conjunction_func=conjunction_func,
                                                    partition=partition)
                                 for child in j['children'] if child['acronym'] not in lst_broken_regions])  # AND THIS
        lst_ids = numpy.hstack(list(lst_ids)).tolist()
        lvl = conjunction_func(list(res_lvl)) + 1
    if lvl == target_lvl:
        if partition:
            lst_ids = []
        lst_ids.append(j['acronym'])
    return lvl, numpy.unique(lst_ids)


def list_cortex_regions(hierarchy, root_acronym="Isocortex"):
    if isinstance(hierarchy, str) or isinstance(hierarchy, os.PathLike):
        if os.path.isdir(hierarchy):
            from .find_atlas_files import find_hierarchy
            json_data = find_hierarchy(hierarchy, format="json")
            assert json_data is not None, "No hierarchy found under {0}".format(hierarchy)
        else:
            with open(hierarchy, 'r') as fid:
                import json
                json_data = json.load(fid)
            if 'msg' in json_data:
                json_data = json_data['msg'][0]
    else:
        json_data = hierarchy
    root = find_root(json_data, root_acronym)
    assert root is not None, "No region {0} found".format(root_acronym)
    return list(find_almost_leaves(root, 1, conjunction_func=numpy.max)[1])
