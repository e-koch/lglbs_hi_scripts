'''
Canonical M31 strip-tile definitions for A+B+C+D imaging.

The LGLBS M31 mosaic is a 49-pointing hex strip: 14 columns alternating 4/3
pointings, 17.4 arcmin spacing, running ~197' x 53' at PA ~37 deg E of N
(the disk major axis). Full-frame joint imaging (12800^2 at 1.5") is not
schedulable on fir (PSF gridding alone >65 hr/channel; see
casa_on_fir_imaging_guide.md). We therefore image 4 column-group tiles
separately and recombine with a linear mosaic in postprocessing.

Tile design (guide sec 8.1):
- Cuts placed between hex columns at u = -45', -0.5', +47' along the strip
  axis (u = position along PA 53 deg from the RA axis).
- 'core' fields: pointings whose column lies inside the tile cut.
- 'boundary' fields: adjacent-column pointings within 21' (the EVLA L-band
  PB radius to pblimit=0.25; PB FWHM 31.7' at 1.42 GHz) of a cut. These are
  gridded by BOTH neighboring tiles so the seam regions are jointly
  deconvolved. Total double-gridded rows ~18%.
- Because tiles share boundary-field data, seam noise between adjacent tile
  images is ~fully correlated: the linear mosaic is a blend, NOT an
  inverse-variance stack (no sqrt(2) gain exists in overlaps).
- Total double-gridded overhead: 70 field-instances vs 49 pointings = ~43%
  extra gridded rows across the 4 tiles (each cut shares the full adjacent
  column, since column spacing 15.1' < 21' PB radius). This is the price of
  jointly-deconvolved seams; per-tile jobs stay at 17-18 fields.
- imsize: core-pointing bounding box + 21' PB + 5' guard per side at 1.5"
  cell, rounded up to FFT-friendly sizes. These are reference values; if
  the pipeline's dynamic sizing chooses differently, prefer these.

Keys already wired up: target_definitions.txt (m31_1..4 phase centers),
dir_key.txt (all tiles -> m31 directory), linearmosaic_definitions.txt
(m31 -> m31_1..4), distance_key.txt.

Field naming assumes the MS FIELD table uses the project names
M31LARGE_0 .. M31LARGE_48 (as in "VLA X-Large Sources.pst"). The split
script verifies this against the MS before selecting.
'''


def _f(indices):
    return ["M31LARGE_{0}".format(i) for i in indices]


M31_TILES = {
    'm31_1': {
        # SW end of the strip (columns u = -97' to -52')
        'phasecenter': 'J2000 00h38m03.59s +40d06m58.5s',
        'imsize': [4608, 4608],
        'core_fields': _f([35, 21, 7, 0, 36, 22, 8, 37, 23, 1, 9, 38, 24, 10]),
        'boundary_fields': _f([39, 25, 11, 2]),
    },
    'm31_2': {
        # columns u = -37' to -6.8'
        'phasecenter': 'J2000 00h41m07.71s +40d51m57.2s',
        'imsize': [4608, 4320],
        'core_fields': _f([39, 2, 25, 11, 40, 12, 26, 3, 41, 13, 27]),
        'boundary_fields': _f([38, 24, 10, 14, 28, 42]),
    },
    'm31_3': {
        # columns u = +8.4' to +38.8'
        'phasecenter': 'J2000 00h43m33.02s +41d28m14.9s',
        'imsize': [4096, 4096],
        'core_fields': _f([14, 42, 28, 4, 15, 43, 29, 16, 30, 44]),
        'boundary_fields': _f([41, 3, 13, 27, 5, 17, 31, 45]),
    },
    'm31_4': {
        # NE end of the strip (columns u = +54' to +100')
        'phasecenter': 'J2000 00h46m07.25s +42d07m57.7s',
        'imsize': [4800, 4608],
        'core_fields': _f([5, 17, 31, 45, 18, 32, 46, 6, 19, 33, 47, 20, 34, 48]),
        'boundary_fields': _f([16, 30, 44]),
    },
}


def tile_fields(tile):
    '''All fields (core + boundary) to be split into and gridded by a tile.'''
    d = M31_TILES[tile]
    return d['core_fields'] + d['boundary_fields']


def field_selection_string(tile):
    '''Comma-joined CASA field selection for the tile.'''
    return ",".join(tile_fields(tile))


if __name__ == "__main__":
    # Sanity summary + checks
    all_core = []
    for name in sorted(M31_TILES):
        d = M31_TILES[name]
        print("{0}: {1} core + {2} boundary = {3} fields, imsize={4}, "
              "phasecenter={5}".format(
                  name, len(d['core_fields']), len(d['boundary_fields']),
                  len(tile_fields(name)), d['imsize'], d['phasecenter']))
        all_core += d['core_fields']
    assert len(all_core) == 49, "Core fields must partition all 49 pointings"
    assert len(set(all_core)) == 49, "A pointing appears in two tiles' cores"
    n_inst = sum(len(tile_fields(t)) for t in M31_TILES)
    print("OK: 49 pointings partitioned; {0} field-instances "
          "(~{1:.0f}% double-gridding overhead)".format(
              n_inst, 100.0 * (n_inst - 49) / 49))
