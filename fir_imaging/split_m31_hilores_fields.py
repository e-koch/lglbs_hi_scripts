#!/usr/bin/env python
import os
import sys
from casatasks import split, listobs

# Initial line staging did not exclude the NGC185/205 fields in the 14A tracks.

ms_path = "/home/ekoch/scratch/VLAXL_imaging/imaging/m31/m31_A+B+C+D_hilores.ms"
output_ms = "/home/ekoch/scratch/VLAXL_imaging/imaging/m31/m31_hilores_M31fields.ms"

if not os.path.exists(ms_path):
    print(f"ERROR: {ms_path} not found")
    sys.exit(1)

print(f"Splitting M31* fields from {ms_path}")
print(f"Output: {output_ms}")

split(
    vis=ms_path,
    outputvis=output_ms,
    field="M31*",
    datacolumn="all"
)

_ = listobs(output_ms, listfile=f"{output_ms}.listobs.txt", overwrite=True)
