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

    return distribution[distribution_id].name, reader(distribution[distribution_id].name)


def query_for_cell_densities(forge):
    resources = forge.search({"type":"METypeDensity","atlasRelease":{"@id":production_bbp_atlas}}, limit=1000, debug=False)
    reshaped_resources = forge.reshape(resources, keep=["id","type", "annotation.hasBody.id", "annotation.hasBody.label",
                     "brainLocation.brainRegion.id", "brainLocation.brainRegion.label", "distribution.atLocation.location"])
    df = forge.as_dataframe(reshaped_resources, nesting="_")

    df["mtype"] = df.apply(lambda row: (row.annotation[0]["hasBody"]["label"],row.annotation[0]["hasBody"]["id"]) , axis=1) 
    df["etype"] = df.apply(lambda row: (row.annotation[1]["hasBody"]["label"],row.annotation[1]["hasBody"]["id"]) , axis=1) 
    type_column = df.pop('type')
    mtype_column = df.pop('mtype')
    etype_column = df.pop('etype')

    df.insert(0, 'type', type_column)
    df.insert(1, 'mtype', mtype_column)
    df.insert(2, 'etype', etype_column)
    df.drop(columns="annotation")

    return df, resources


def annotation_from_forge(forge):
    from voxcell import VoxelData
    atlas_release = forge.retrieve(production_bbp_atlas)
    atlas_release._store_metadata["_rev"]
    parcellation_volume = forge.retrieve(atlas_release.parcellationVolume.id)
    print("Starting download of annotation file...")
    _, ann = download_and_read(forge, parcellation_volume, VoxelData.load_nrrd)
    print("...done!")
    return ann


def hierarchy_from_forge(forge):
    from voxcell.nexus.voxelbrain import RegionMap
    atlas_release = forge.retrieve(production_bbp_atlas)
    atlas_release._store_metadata["_rev"]
    parcellation_ontology = forge.retrieve(atlas_release.parcellationOntology.id, cross_bucket=True)
    print("Starting download of hierarchy file...")
    hier_file, hier = download_and_read(forge, parcellation_ontology, RegionMap.load_json)
    print("...done!")
    return hier_file, hier
