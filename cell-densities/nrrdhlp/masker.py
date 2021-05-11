import numpy

from voxcell import VoxelData
from voxcell.nexus.voxelbrain import RegionMap

layer_match_pat = "@.*[lL]ayer {0}[ab]?$"


class Masker(object):
    """
    A helper class that generates bitwise masks for brain volumes according to custom specifications
    """

    def __init__(self, hier_fn, annotation_fn):
        """
        :param hier_fn: path to hierarchy.json file
        :param annotation_fn: path to brain region annotation volume file
        """
        self.rmap = RegionMap.load_json(hier_fn)
        self.annotation = VoxelData.load_nrrd(annotation_fn)

    def mask(self, specs):
        """
        :param specs: Specification of a brain region mask
        :return: bitwise mask according to specifications
        """
        lst_ids = self.target_ids(specs)
        return numpy.in1d(self.annotation.raw.flat, lst_ids).reshape(self.annotation.raw.shape)

    def target_ids(self, specs):
        if isinstance(specs, str):
            if specs.startswith("Layer"):  # TODO: This can be done better. Maybe a dict?
                specs = specs[5:]
                return list(self.rmap.find(layer_match_pat.format(specs), "name",
                                           with_descendants=True))
            return list(self.rmap.find(specs, "acronym", with_descendants=True))
        elif isinstance(specs, list):
            conjunction, rest = specs
            lst_res = [self.target_ids(_rest) for _rest in rest]
            if len(lst_res) == 0:
                return []
            base_result = lst_res[0]
            if conjunction == "and":
                for _res in lst_res:
                    base_result = numpy.intersect1d(base_result, _res)
            elif conjunction == "or":
                for _res in lst_res:
                    base_result = numpy.union1d(base_result, _res)
            else:
                raise ValueError("Mask specification error: {0}".format(str(specs)))
            return base_result
        else:
            raise ValueError("Mask specification error: {0}".format(str(specs)))
