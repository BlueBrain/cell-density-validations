from importlib.metadata import distribution
import numpy

from voxcell import VoxelData
from voxcell.nexus.voxelbrain import RegionMap

# layer_match_pat = "@.*[lL]ayer {0}[ab]?$"
layer_match_pat = "@.*{0}[ab]?$"
layer_match_prop = "acronym"

default_endpoint = "https://bbp.epfl.ch/nexus/v1"
production_bbp_atlas = "https://bbp.epfl.ch/neurosciencegraph/data/4906ab85-694f-469d-962f-c0174e901885"


def download_and_read(forge, thing, reader, distribution_id=None):
    if distribution_id is None: distribution_id = 0
    if not isinstance(thing.distribution, list):
        distribution = [thing.distribution]
    else:
        distribution = thing.distribution
    print("Downloading: {0}".format(distribution[distribution_id].name))
    forge.download(thing, "distribution.contentUrl", ".", overwrite=True)

    return reader(distribution[distribution_id].name)


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

        if token is None:
            token = getpass.getpass()
        if endpoint is None:
            endpoint = default_endpoint
        
        if isinstance(yml_config_or_forge, KnowledgeGraphForge):
            forge = yml_config_or_forge
        else:
            print("Setting up kgforge...")
            forge = KnowledgeGraphForge(yml_config_or_forge,
                                token=token,
                                endpoint=endpoint, 
                                bucket="bbp/atlas")
        
        atlas_release = forge.retrieve(production_bbp_atlas)
        atlas_release._store_metadata["_rev"]
        parcellation_ontology = forge.retrieve(atlas_release.parcellationOntology.id, cross_bucket=True)
        parcellation_volume = forge.retrieve(atlas_release.parcellationVolume.id)
        print("...done!")
        print("Starting downloads...")
        ann = download_and_read(forge, parcellation_volume, VoxelData.load_nrrd)
        hier = download_and_read(forge, parcellation_ontology, RegionMap.load_json)
        print("...done!")
        return cls(hier, ann)
        
