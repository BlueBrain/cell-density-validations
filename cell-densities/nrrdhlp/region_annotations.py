import copy
import json

import numpy
import pandas
import voxcell
from scipy.spatial import KDTree

class AnnotationWrapper(object):
    def __init__(self, forge, method_dict={}, use_hier_file=None, use_anno_file=None):

        if use_anno_file:
            self._ann = voxcell.voxel_data.VoxelData.load_nrrd(use_anno_file)
            pass
        else:
            from .atlas_from_forge import annotation_from_forge
            self._ann = annotation_from_forge(forge)

        if use_hier_file:
            self.hier_file = use_hier_file
            self.hier = voxcell.nexus.voxelbrain.RegionMap.load_json(use_hier_file)
        else:
            from .atlas_from_forge import hierarchy_from_forge
            self.hier_file, self.hier = hierarchy_from_forge(forge)

        self.voxel_counts = pandas.Series(self._ann.raw.flatten()).value_counts()
        self._input_voxel_counts = self.voxel_counts.copy()

        Y, X, Z = numpy.meshgrid(range(self._ann.shape[1]),
                                 range(self._ann.shape[0]),
                                 range(self._ann.shape[2]))
        self.ann_df = pandas.DataFrame({
            "x": X.flatten().astype(float),
            "y": Y.flatten().astype(float),
            "z": Z.flatten().astype(float),
            "d": self._ann.raw.flatten()
        }).set_index("d")
        self.methods = method_dict.copy()

    @staticmethod
    def __nearest_n_neighbor_extrapolator__(child_coords, reg_coords, n=6):
        kd = KDTree(child_coords.values)
        distance_arr, idx_arr = kd.query(reg_coords.values, n)

        neighbor_reg_ids = child_coords.index.values[idx_arr]
        local_indices = numpy.repeat(numpy.arange(len(idx_arr)).reshape((-1, 1)), idx_arr.shape[1], axis=1)
        df = pandas.DataFrame(
            {
                "lcl_idx": local_indices.flatten(),
                "neighbor_ids": neighbor_reg_ids.flatten(),
                "inverse_distance": 1.0 / distance_arr.flatten()
            }
        )
        evaluated = df.groupby(["lcl_idx", "neighbor_ids"])["inverse_distance"].sum().unstack("neighbor_ids")
        evaluated = evaluated.idxmax(axis="columns").sort_index()
        return evaluated
    
    def lookup_method(self, reg_id):
        reg_acronym = self.hier.get(reg_id, "acronym")
        if reg_acronym in self.methods:
            return self.methods[reg_acronym]
        parent = self.hier._parent[reg_id]
        if parent is None:
            return "create"
        return self.lookup_method(parent)
    
    def make_new_region_id(self):
        return int(numpy.max(list(self.hier._data.keys())) + 1)

    def make_new_region(self, parent_id, **kwargs):
        if not hasattr(self, "replacement_ids"): self.replacement_ids = {}
        parent_name = self.hier.get(parent_id, "name")
        parent_acronym = self.hier.get(parent_id, "acronym")
        reg_name = parent_name + ": Other"
        reg_acronym = parent_acronym + "_O"
        reg_id = self.make_new_region_id()

        self.hier._parent[reg_id] = parent_id
        self.hier._children[parent_id].append(reg_id)
        self.hier._children[reg_id] = []

        # ADD ATTRIBUTES HERE
        entry_dict = {
            "id": reg_id,
            "acronym": reg_acronym,
            "name": reg_name,
            "parent_structure_id": parent_id,
            "color_hex_triplet": self.hier.get(parent_id, "color_hex_triplet"),
            "ontology_id": 1,  # fixed value
            "hemisphere_id": 3,  # fixed value
            "children": []
        }
        entry_dict.update(kwargs)
        self.hier._data[reg_id] = entry_dict
        self.replacement_ids[parent_id] = reg_id
        return reg_id

    @property
    def ann(self):
        ann_out = numpy.zeros_like(self._ann.raw)
        _ann_df = self.ann_df.astype(int)
        ann_out[_ann_df["x"], _ann_df["y"], _ann_df["z"]] = self.ann_df.index.values
        ann_out = voxcell.VoxelData(ann_out, self._ann.voxel_dimensions, offset=self._ann.offset)
        return ann_out
    
    def all_children_of(self, reg_id):
        child_ids = list(self.hier.find(reg_id, "id", with_descendants=True))
        child_ids.remove(reg_id)
        return child_ids
    
    def child_region_indices(self, reg_id):
        child_ids = self.all_children_of(reg_id)
        child_coords = self.ann_df.loc[numpy.isin(self.ann_df.index, child_ids)]
        return child_coords

    def extrapolate_to_leaf_regions(self, reg_id, n=6):
        child_coords = self.child_region_indices(reg_id)
        reg_coords = self.ann_df.loc[reg_id]

        evaluated = self.__nearest_n_neighbor_extrapolator__(child_coords, reg_coords, n=n)

        self.ann_df.index.values[self.ann_df.index == reg_id] = evaluated.values
        self.voxel_counts = self.voxel_counts.add(evaluated.value_counts(), fill_value=0)
        self.voxel_counts[reg_id] = 0
    
    def create_new_leaf_region(self, reg_id):
        new_reg_id = self.make_new_region(reg_id)
        L = numpy.sum(self.ann_df.index == reg_id)
        self.ann_df.index.values[self.ann_df.index == reg_id] = new_reg_id
        self.voxel_counts[new_reg_id] = L
        self.voxel_counts[reg_id] = 0
        self._input_voxel_counts[new_reg_id] = 0
    
    def __recursive_check__(self, reg_id, method=None, **kwargs):
        # Seems to be number of voxels belonging to reg_id
        count_direct = self.voxel_counts.get(reg_id, 0)
        child_ids = self.hier._children[reg_id]
        # Seems to be number of voxels belonging to children of reg_id
        count_indirect = 0
        log = ""
        for _id in child_ids:
            _count, _log = self.__recursive_check__(_id, method=method, **kwargs)
            count_indirect += _count
            log += _log

        method = method or self.lookup_method(reg_id)
        if method == "extrapolate":
            if (count_direct > 0) and (count_indirect > 0):
                log += """Region {0} has {1} voxels, {2} of which ({3}%) are directly associated\n
                """.format(self.hier.get(reg_id, "name"),
                           count_direct + count_indirect,
                           count_direct,
                           (100.0 * count_direct) / (count_direct + count_indirect))
                self.extrapolate_to_leaf_regions(reg_id, **kwargs)
        elif method == "create":
            # For creation, we don't care about voxels of children, only whether there are children
            # or not.
            if (count_direct > 0) and (len(child_ids) > 0):
                self.create_new_leaf_region(reg_id)
        else:
            raise ValueError("Unknown method: {0}".format(method))
        return count_direct + count_indirect, log
    
    def launch_fix(self, **kwargs):
        root_ids = [_id for _id, v in self.hier._parent.items() if v is None]
        log = ""
        for _id in root_ids:
            _, _log = self.__recursive_check__(_id, **kwargs)
            log += _log

        return log

    def update_json_new_hier(self, orig_hier) -> dict:
        out_hier = copy.deepcopy(orig_hier)

        def _append_new_regions(node):
          if isinstance(node, dict):
              [_append_new_regions(ii) for ii in node["children"]]
              if node["id"] in self.replacement_ids.keys():
                  node["children"].append(self.hier._data[self.replacement_ids[node["id"]]])
          return node

        return _append_new_regions(out_hier)

    def fix(self, fn_ann_out, fn_hier_out, fn_log, use_mba_hierarchy=False, **kwargs):

        log_str = self.launch_fix(**kwargs)
        self.ann.save_nrrd(fn_ann_out)

        # Only create a new hierarchy file if new regions were actually created
        if hasattr(self, "replacement_ids"):
            # Note that the "mba_hierarchy_v3l23split.json" file is the default
            # one that Nexus downloads as the "hierarchy" file, which is
            # produced as an OUTPUT of the Atlas-pipeline. Instead, we want to
            # work with the older "hierarchy_l23split.json", which has a
            # slightly different structure and is an INPUT to the
            # Atlas-pipeline. In particular, the latter has a "header" while
            # the former does not, and the latter uses a JSON indent of 1 space
            # while the former uses 2 spaces.
            #
            # Also note that we must re-load the hierarchy files, since the the
            # `self.hier` RegionMap object does not have the precise structure
            # we want.
            if use_mba_hierarchy:
                # This option is for using the just-downloaded "mba_..." hierarchy version
                with open('mba_hierarchy_v3l23split.json', 'r') as f:
                    original_hierarchy = json.load(f)
                out_hierarchy = self.update_json_new_hier(original_hierarchy)
                with open(fn_hier_out, "w") as fid:
                    json.dump(out_hierarchy, fid, indent=2)
            else:
                # This option is for using the original "hierarchy_l23split.json" hierarchy version
                with open(self.hier_file, 'r') as f:
                    original_hierarchy = json.load(f)
                # Remove the "header"
                original_hierarchy = copy.deepcopy(original_hierarchy['msg'][0])
                out_hierarchy = self.update_json_new_hier(original_hierarchy)
                # Re-add the "header" so it is as close as possible to "hierarchy_l23split.json"
                out_hierarchy = {
                    "success": True,
                    "id": 0,
                    "start_row": 0,
                    "num_rows": 1,
                    "total_rows": 1,
                    "msg": [ out_hierarchy ]}
                with open(fn_hier_out, "w") as fid:
                    json.dump(out_hierarchy, fid, indent=1)

        log_str += "Changes in leaf volumes:\n\n"
        voxel_count_changes = self.voxel_counts.divide(self._input_voxel_counts)
        leaf_ids = [_id for _id in voxel_count_changes.index
                    if len(self.hier._children.get(_id, [1])) == 0]
        voxel_count_changes = voxel_count_changes[leaf_ids]
        voxel_count_changes = voxel_count_changes[voxel_count_changes != 1.0]
        for reg_id, value in voxel_count_changes.items():
            log_str += """{0} from {1} to {2} ({3}%)
            """.format(self.hier.get(reg_id, "name"),
                       self._input_voxel_counts[reg_id],
                       self.voxel_counts[reg_id],
                       100 * (value - 1.0))
        with open(fn_log, "w") as fid:
            fid.write(log_str)
