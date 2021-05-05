import numpy

from voxcell import VoxelData
from .masker import Masker


def absmax(arr_in):
    return arr_in[numpy.argmax(numpy.abs(arr_in))]


class Tester(object):

    def __init__(self, hier_fn, annotation_fn, nrrd_dict, circuits_dict=None):
        self.masker = Masker(hier_fn, annotation_fn)
        self.nrrd = nrrd_dict
        self.check_inputs()
        self.file_test_count = {}
        self.file_passed_count = {}
        self.__cached_vols__ = {}
        self.circuits = {}
        if circuits_dict is not None:
            import bluepy.v2 as bluepy
            for k, v in circuits_dict.items():
                self.circuits[k] = bluepy.Circuit(v)

    def check_inputs(self):
        """
        Performs an initial test on the provided inputs
        """
        for nrrdspec in self.nrrd.values():
            assert isinstance(nrrdspec, str) or (isinstance(nrrdspec, list) and len(nrrdspec) > 0)
            # Should also test that the annotation volume is compatible with the nrrd files (voxel size etc)

    @staticmethod
    def test_against_tolerance(A, B, tolerance_spec):
        """
        Tests whether the error between A and B is within the specified tolerance
        :param A: 1-dimensional numpy.array
        :param B: 1-dimensional numpy.array of the same length
        :param tolerance_spec: A dict specifying the maximal tolerated relative or absolute error
        :return: err1 (float): Maximal element-wise absolute error,
                 err2 (float): Maximal element-wise relative error,
                 passed (bool): True if error does not exceed any tolerance
        """
        err1 = absmax(A - B)
        err2 = absmax((A - B) / (A + B + 1E-100))
        if len(tolerance_spec) == 0:
            return err1, err2, err1 == 0
        if "absolute" in tolerance_spec:
            if numpy.abs(err1) > tolerance_spec["absolute"]:
                return err1, err2, False
        if "relative" in tolerance_spec:
            if numpy.abs(err2) > tolerance_spec["relative"]:
                return err1, err2, False
        return err1, err2, True

    @staticmethod
    def neuron_dens_voxel_data(xyz, reference):
        bins = [numpy.arange(l + 1) for l in reference.shape]
        xyz = reference.positions_to_indices(xyz)
        H = numpy.histogramdd([xyz["x"], xyz["y"], xyz["z"]], bins=bins)[0]
        return VoxelData(1E9 * H / reference.voxel_volume,
                         reference.voxel_dimensions,
                         offset=reference.offset)

    def neuron_density_volume_from_circuit(self, specifications):
        circ_name = specifications["circuit"]
        if circ_name not in self.circuits:
            raise ValueError("Circuit {0} unknown!".format(circ_name))
        circ = self.circuits[circ_name]
        grp = specifications.get("group", None)
        fltrs = specifications.get("filters", {})
        props = circ.cells.get(group=grp, properties=["x", "y", "z"] + list(fltrs.keys()))

        for k, v in fltrs.items():
            if isinstance(v, list):
                valid = numpy.in1d(props[k].values, v)
            else:
                valid = (props[k] == v).values
            props = props.iloc[valid]

        return self.neuron_dens_voxel_data(props[["x", "y", "z"]], self.masker.annotation)

    @staticmethod
    def mask_and_sum(A, B, mask, multiply_with_volume):
        """
        :param A: List of VoxelData objects or floats (or both mixed)
        :param B: Same as A, but can be different length
        :param mask: bitwise mask, same shape as VoxelData in A and B
        :param multiply_with_volume: if true, the VoxelData in A and B are multiplied with their voxel volumes
        :return:
         Asum: 1-dimensional numpy.array. If all objects in A and B are VoxelData, then the length of Asum will be
         equal to the number of nonzero entries in mask. Entries are the sum over elements in A for masked voxels.
         If any float or int is provided in A or B, the length of Asum will be 1 and it will contain the sum also over
         voxels.
         Bsum: See Asum, but for B
        """
        perform_idv_voxel_test = all([isinstance(_x, VoxelData) for _x in A + B])

        def voxels_to_1darray(vol):
            if mask is None:
                out = vol.raw.astype(float).flatten()
            else:
                out = vol.raw[mask].astype(float)
            if not perform_idv_voxel_test:
                out = numpy.sum(out, keepdims=True)
            if multiply_with_volume:
                out = out * vol.voxel_volume / 1E9
            return out

        A = [voxels_to_1darray(_x) if isinstance(_x, VoxelData) else numpy.array([_x]) for _x in A]
        B = [voxels_to_1darray(_x) if isinstance(_x, VoxelData) else numpy.array([_x]) for _x in B]

        return numpy.vstack(A).sum(axis=0), numpy.vstack(B).sum(axis=0)

    def evaluate(self, specs):
        """
        :param specs: Evaluate the specified test. Both on a per-voxel basis and for the sum over voxels.
        :return:
          Asum (float), the sum of masked voxels in "left"
          Bsum (float), the sum of masked voxels in "right"
          errspec (dict): keys specify the type of error, values the error value
          passed (bool): True if all error is within tolerance
        """
        if "mask" in specs:
            mask = self.masker.mask(specs["mask"])
        else:
            mask = None

        modality = specs.get("modality", "density")
        is_cell_count = (modality == "count")

        A = [self.load_nrrd(_x) for _x in specs["left"]]
        B = [self.load_nrrd(_x) for _x in specs["right"]]

        A, B = self.mask_and_sum(A, B, mask, is_cell_count)
        abs_err_sum, rel_err_sum, passed = self.test_against_tolerance(A.sum(keepdims=True),
                                                                       B.sum(keepdims=True),
                                                                       specs.get("tolerance", {}))
        err_spec = {"Absolute sum": abs_err_sum, "Relative sum": rel_err_sum}
        if len(A) > 1:
            abs_err, rel_err, passed_vxl = self.test_against_tolerance(A, B, specs.get("tolerance", {}))
            err_spec["Absolute max over voxels"] = abs_err
            err_spec["Relative max over voxels"] = rel_err
            passed = passed & passed_vxl

        for fl in specs["left"] + specs["right"]:
            if isinstance(fl, str):
                self.file_test_count[fl] = self.file_test_count.setdefault(fl, 0) + 1
                self.file_passed_count[fl] = self.file_passed_count.setdefault(fl, 0) + int(passed)
            elif isinstance(fl, dict):
                self.file_test_count[fl["circuit"]] = self.file_test_count.setdefault(fl["circuit"], 0) + 1
                self.file_passed_count[fl["circuit"]] = self.file_passed_count.setdefault(fl["circuit"], 0) + int(passed)

        return A.sum(), B.sum(), err_spec, passed

    def load_nrrd(self, to_load):
        """
        Load an nrrd file or the sum of several
        :param to_load: list of nrrd file names or a single nrrd file name
        :return: VoxelData object. If a single nrrd file was specified that one is loaded. If a list was specified,
        the sum of all of them is returned.
        """
        if isinstance(to_load, int) or isinstance(to_load, float):
            return to_load
        to_load_str = str(to_load)
        if to_load_str not in self.__cached_vols__:
            if isinstance(to_load, dict):
                self.__cached_vols__[to_load_str] = self.neuron_density_volume_from_circuit(to_load)
            else:
                if to_load not in self.nrrd:
                    raise ValueError("No entry for {0} found".format(to_load))

                filespec = self.nrrd[to_load]
                if isinstance(filespec, str):
                    self.__cached_vols__[to_load_str] = VoxelData.load_nrrd(filespec)
                elif isinstance(filespec, list):
                    loaded = VoxelData.load_nrrd(filespec[0])  # Checked for len == 0 before!
                    raw = loaded.raw.copy()
                    for _fn in filespec[1:]:
                        raw += VoxelData.load_nrrd(_fn).raw
                    self.__cached_vols__[to_load_str] = VoxelData(raw, loaded.voxel_dimensions, offset=loaded.offset)

        return self.__cached_vols__[to_load_str]

    def test_count_str(self):
        out = ""
        for fl in self.nrrd.keys():
            out = out + "{0}: Passed {1} out of {2} tests\n".format(fl, self.file_passed_count.get(fl, 0),
                                                                    self.file_test_count.get(fl, 0))
        return out
