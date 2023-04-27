import numpy
import pandas
import voxcell
from scipy.spatial import KDTree

from .atlas_from_forge import hierarchy_from_forge, annotation_from_forge

class AnnotationWrapper(object):
    def __init__(self, forge, method_dict={}):
        self._ann = annotation_from_forge(forge)
        self.hier = hierarchy_from_forge(forge)
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
        return int(numpy.max(list(self.hier._data.keys())) + 1 + numpy.random.randint(10))
    
    def make_new_region(self, parent_id, **kwargs):
        parent_name = self.hier.get(parent_id, "name")
        parent_acronym = self.hier.get(parent_id, "acronym")
        reg_name = parent_name + ": Other"
        reg_acronym = parent_acronym + "_O"
        reg_id = self.make_new_region_id()

        self.hier._parent[reg_id] = parent_id
        self.hier._children[parent_id].append(reg_id)
        self.hier._children[reg_id] = []

        entry_dict = {
            "id": reg_id,
            "acronym": reg_acronym,
            "name": reg_name,
            "parent_structure_id": parent_id
        }
        entry_dict.update(kwargs)
        self.hier._data[reg_id] = entry_dict
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
        count_direct = self.voxel_counts.get(reg_id, 0)
        child_ids = self.hier._children[reg_id]
        count_indirect = 0
        log = ""
        for _id in child_ids:
            _count, _log = self.__recursive_check__(_id, method=method, **kwargs)
            count_indirect += _count
            log += _log

        if (count_direct > 0) and (count_indirect > 0):
            log += """Region {0} has {1} voxels, {2} of which ({3}%) are directly associated\n
            """.format(self.hier.get(reg_id, "name"), 
                       count_direct + count_indirect,
                       count_direct,
                       (100.0 * count_direct) / (count_direct + count_indirect))
            method = method or self.lookup_method(reg_id)
            if method == "extrapolate":
                self.extrapolate_to_leaf_regions(reg_id, **kwargs)
            elif method == "create":
                self.create_new_leaf_region(reg_id)
            else:
                raise ValueError("Unknow method: {0}".format(method))
        return count_direct + count_indirect, log
    
    def launch_extrapolation(self, **kwargs):
        root_ids = [_id for _id, v in self.hier._parent.items() if v is None]
        log = ""
        for _id in root_ids:
            _, _log = self.__recursive_check__(_id, **kwargs)
            log += _log

        return log
    
    def fix(self, fn_ann_out, fn_hier_out, fn_log, **kwargs):
        log_str = self.launch_extrapolation(**kwargs)
        self.ann.save_nrrd(fn_ann_out)
        
        with open(fn_hier_out, "w") as fid:
            fid.write("This is where the updated hierarchy would go, but I don't know how to write it")

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
