import pandas
import voxcell
import getpass

from kgforge.core import KnowledgeGraphForge

from .masker import Masker
from .atlas_from_forge import download_and_read, default_endpoint, query_for_cell_densities

def intersection_size(rmap, tgt_ids):
    def func_to_apply(region_label):
        file_ids = rmap.find(region_label, "acronym", with_descendants=True)
        return len(file_ids.intersection(tgt_ids)) > 0
    return func_to_apply

def add_info(forge, out_df, row):
    # NOTE: E/I info and additional classification is supposed to be added here
    # by looking up info on the m/e-types in the ontology! For now hard coded, 
    # but that is to be replaced!
    out_df["mtype"] = row["mtype"][0]
    out_df["etype"] = row["etype"][0]
    
    if "PC" in row["mtype"][0]:
        out_df["synapse_class"] = "EXC"
    elif row["mtype"] == "Excitatory":
        out_df["synapse_class"] = "EXC"
    elif row["mtype"] == "L4_SSC":
        out_df["synapse_class"] = "EXC"
    else:
        out_df["synapse_class"] = "INH"

def build_from_row(forge, masker, msk, sub_regions, resources):
    def func_to_return(row):
        i = row.name
        V = download_and_read(forge, resources[i], voxcell.VoxelData.load_nrrd)
        v = V.raw[msk]

        _df = pandas.DataFrame({"sub-region": sub_regions, "count": v * V.voxel_volume / (1000 ** 3)})
        _df = _df.groupby("sub-region").agg("sum").reset_index()
        _df["sub-region"] = [masker.rmap.get(_x, "acronym") for _x in _df["sub-region"]]
        add_info(forge, _df, row)
        
        return _df
    return func_to_return

def cell_counts_in_region(region_acronym, yml_config, token=None, endpoint=None):
    if endpoint is None:
        endpoint = default_endpoint
    if token is None:
        token = getpass.getpass()
    
    forge = KnowledgeGraphForge(yml_config,
                                token=token,
                                endpoint=endpoint, 
                                bucket="bbp/atlas")
    
    df, resources = query_for_cell_densities(forge)
    masker = Masker.from_knowledge_graph(forge, token=token)

    msk = masker.mask(region_acronym)
    tgt_ids = masker.target_ids(region_acronym)
    sub_regions = masker.annotation.raw[msk]
    rel_df = df.loc[df["brainLocation_brainRegion_label"].apply(intersection_size(masker.rmap, tgt_ids))]

    func = build_from_row(forge, masker, msk, sub_regions, resources)

    composition = pandas.concat([func(row[1]) for row in rel_df.iterrows()], axis=0)
    cols = list(composition.columns)
    cols.remove("count")
    composition = composition.set_index(cols)

    return composition
