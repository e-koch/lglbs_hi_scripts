#!/usr/bin/env python
'''
Split the staged M31 line MS into the 4 strip-tile MSes for separate
joint-mosaic imaging (see m31_tile_definitions.py and
casa_on_fir_imaging_guide.md sec 8.1).

Pipeline placement: runs AFTER the standard phangsPipeline staging of the
full m31 target and BEFORE per-tile imaging. Output MSes are named
{tile}_{config}_{product}.ms in the same imaging directory, which is where
the imaging handlers look for them given dir_key.txt (m31_N -> m31).

Selecting only the M31LARGE_* tile fields also drops the NGC185/NGC205
fields carried along in the 14A tracks (the purpose of the original
version of this script).

Usage (via casa -c):
    casa ... -c split_m31_hilores_fields.py <tile|all> [config] [product]

    tile     : m31_1 .. m31_4, or 'all' to loop over every tile
    config   : interferometer config string (default: A+B+C+D)
    product  : line product (default: hilores)

Notes:
- keepflags=False drops fully-flagged rows: fewer rows to grid (typically
  10-30% for VLA data) and a smaller MS to stage to node-local disk.
- datacolumn='all' propagates whichever of DATA/CORRECTED exists in the
  staged MS.
'''

import os
import sys

from casatasks import split, listobs, flagdata

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m31_tile_definitions import M31_TILES, field_selection_string

# ----------------------------------------------------------------------
# Arguments and paths
# ----------------------------------------------------------------------

# Args arrive after the CASA-internal ones; read from the end as in the
# other fir_imaging scripts. Allow 1-3 trailing args.
args = sys.argv[-3:]
this_config = "A+B+C+D"
this_product = "hilores"

if args[-1] in list(M31_TILES) + ['all']:
    tile_arg = args[-1]
elif args[-2] in list(M31_TILES) + ['all']:
    tile_arg, this_product = args[-2], args[-1]
else:
    tile_arg, this_config, this_product = args[0], args[1], args[2]

tiles = list(M31_TILES) if tile_arg == 'all' else [tile_arg]

imaging_path = "/home/ekoch/scratch/VLAXL_imaging/imaging/m31"
ms_path = os.path.join(imaging_path,
                       "m31_{0}_{1}.ms".format(this_config, this_product))

if not os.path.exists(ms_path):
    print("ERROR: {0} not found".format(ms_path))
    sys.exit(1)

# ----------------------------------------------------------------------
# Verify the expected field names exist in the MS
# ----------------------------------------------------------------------

try:
    from casatools import msmetadata
    msmd = msmetadata()
    msmd.open(ms_path)
    ms_fields = set(msmd.fieldnames())
    msmd.close()
except ImportError:
    ms_fields = None
    print("WARNING: casatools msmetadata unavailable; skipping field check")

if ms_fields is not None:
    for tile in tiles:
        missing = [f for f in field_selection_string(tile).split(",")
                   if f not in ms_fields]
        if missing:
            print("ERROR: {0}: fields not in MS FIELD table: {1}".format(
                tile, missing))
            print("MS fields are: {0}".format(sorted(ms_fields)))
            sys.exit(1)

# ----------------------------------------------------------------------
# Split each tile
# ----------------------------------------------------------------------

for tile in tiles:

    fields = field_selection_string(tile)
    n_fields = len(fields.split(","))

    output_ms = os.path.join(
        imaging_path, "{0}_{1}_{2}.ms".format(tile, this_config, this_product))

    if os.path.exists(output_ms):
        print("Found existing {0}. Removing and re-splitting.".format(output_ms))
        os.system("rm -r {0}".format(output_ms))

    print("Splitting {0}: {1} fields from {2}".format(tile, n_fields, ms_path))
    print("Output: {0}".format(output_ms))

    split(
        vis=ms_path,
        outputvis=output_ms,
        field=fields,
        datacolumn="all",
        keepflags=False,
    )

    _ = listobs(output_ms,
                listfile="{0}.listobs.txt".format(output_ms),
                overwrite=True)

    # Record the remaining flag fraction for the imaging bookkeeping
    # summary = flagdata(vis=output_ms, mode='summary')
    # if summary is not None and 'flagged' in summary:
    #     frac = summary['flagged'] / summary['total']
    #     print("{0}: {1:.1f}% of remaining rows flagged".format(
    #         tile, 100 * frac))

print("M31 tile split complete: {0}".format(", ".join(tiles)))
