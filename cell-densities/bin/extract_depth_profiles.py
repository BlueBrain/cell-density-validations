import voxcell
import pandas
import numpy
import glob
import os
import json

from nrrdhlp.masker import Masker
from nrrdhlp.hierarchyhlp import *
from nrrdhlp.find_atlas_files import *


def make_region_profiles(masker, depths, regions, density_fn, d_bins):
    print("Loading and analyzing {0}".format(density_fn))
    name = os.path.split(os.path.splitext(density_fn)[0])[1].replace("[cell_density]", "")
    density = voxcell.VoxelData.load_nrrd(density_fn)
    res = []
    idx = []
    for region in regions:
        msk = masker.mask(region)
        d_idx = numpy.digitize(depths.raw[:, :, :, 1][msk], d_bins) - 1
        d = density.raw[msk]
        dens_smpls = [d[d_idx == i] for i in range(len(d_bins) - 1)]
        res.append(list(map(numpy.mean, dens_smpls)))
        idx.append((name, region))
    midx = pandas.MultiIndex.from_tuples(idx, names=["type", "region"])
    return pandas.DataFrame(numpy.vstack(res).transpose(), columns=midx, index=d_bins[:-1])


def configure():
    import sys
    if os.path.isfile(sys.argv[1]):
        with open(sys.argv[1], "r") as fid:
            cfg = json.load(fid)
        root = cfg.get("atlas",
                       os.path.split(sys.argv[1])[0])
    else:
        cfg = {}
        root = sys.argv[1]
    density_files = cfg.get("densities",
                            glob.glob(os.path.join(root, "[cell_density*"))
                            )
    hier_fn = cfg.get("hierarchy", find_hierarchy(root))
    with open(hier_fn, "r") as fid:
        hier = json.load(fid)["msg"][0]
    regions = cfg.get("regions", list_cortex_regions(hier))
    annotation_fn = cfg.get("annotations", find_regions(root))
    depth_fn = cfg.get("depths", os.path.join(root, "[PH]1.nrrd"))
    depths = voxcell.VoxelData.load_nrrd(depth_fn)
    mx_depths = cfg.get("max_depth", numpy.nanmax(depths.raw))
    nbins = cfg.get("nbins", 200)

    if len(sys.argv) > 2:
        out_fn = sys.argv[2]
    else:
        out_fn = cfg.get("out_fn", "layer_profiles.pkl")

    masker = Masker(hier_fn, annotation_fn)
    return density_files, hier, masker, depths, regions, out_fn, mx_depths, nbins


def main():
    density_files, hier, masker, depths, regions, out_fn, mx_depths, nbins = configure()
    d_bins = numpy.linspace(0, mx_depths + 1E-9, nbins + 1)

    res = []
    for fn in density_files:
        res.append(make_region_profiles(masker, depths, regions, fn, d_bins))

    pandas.concat(res, axis=1).to_pickl(out_fn)


if __name__ == "__main__":
    main()

