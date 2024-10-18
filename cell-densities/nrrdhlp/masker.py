# SPDX-License-Identifier: Apache-2.0
from importlib.metadata import distribution
import numpy

from voxcell import VoxelData
from voxcell.nexus.voxelbrain import RegionMap

# layer_match_pat = "@.*[lL]ayer {0}[ab]?$"
layer_match_pat = "@.*{0}[ab]?$"
layer_match_prop = "acronym"

from . import atlas_from_forge
from .atlas_from_forge import download_and_read, annotation_from_forge, hierarchy_from_forge

default_endpoint = atlas_from_forge.default_endpoint


class Masker(object):
    """
    A helper class that generates bitwise masks for brain volumes according to custom specifications
    """

    def __init__(self, hier, annotation):
        """
        :param hier: path to hierarchy.json file, or a RegionMap
        :param annotation: path to brain region annotation volume file, or the annotation volume
        """
        if isinstance(hier, RegionMap):
            self.rmap = hier
        else:
            self.rmap = RegionMap.load_json(hier)
        if isinstance(annotation, VoxelData):
            self.annotation = annotation
        else:
            self.annotation = VoxelData.load_nrrd(annotation)

    def mask(self, specs):
        """
        :param specs: Specification of a brain region mask
        :return: bitwise mask according to specifications
        """
        lst_ids = self.target_ids(specs)
        return numpy.isin(self.annotation.raw, lst_ids)

    def target_ids(self, specs):
        if isinstance(specs, str):
            if specs.startswith("Layer"):  # TODO: This can be done better. Maybe a dict?
                specs = specs[5:]
                return list(self.rmap.find(layer_match_pat.format(specs), layer_match_prop,
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
    
    @classmethod
    def from_knowledge_graph(cls, yml_config_or_forge, endpoint=None, token=None):
        import getpass
        from kgforge.core import KnowledgeGraphForge
        
        if isinstance(yml_config_or_forge, KnowledgeGraphForge):
            forge = yml_config_or_forge
        else:
            print("Setting up kgforge...")
            if token is None:
                token = getpass.getpass()
            if endpoint is None:
                endpoint = default_endpoint
            forge = KnowledgeGraphForge(yml_config_or_forge,
                                token=token,
                                endpoint=endpoint, 
                                bucket="bbp/atlas")
        
        ann = annotation_from_forge(forge)
        hier = hierarchy_from_forge(forge)
        return cls(hier, ann)
        
