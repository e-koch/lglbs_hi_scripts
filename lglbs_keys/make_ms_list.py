
from pathlib import Path

targets = ['m31']

configs = ['A', 'B', 'C', 'D']


# ms_path = Path("/reduction/erickoch/LGLBS/line_imaging/MeasurementSets/")
ms_path = Path("/home/ekoch/nearline/rrg-eros-ab/ekoch/calibrated/")

# e.g. M31_C_20A-346.sb38095502.eb38174408.58988.69745849537.speclines.ms.split.tar

for target in targets:
    for config in configs:

        order_dict = {}
        for this_ms_name in (ms_path / target.upper()).glob(f"*_{config}_*specline*.tar"):

            sdm_name = this_ms_name.name.split(".speclines")[0].split("_")[-1]

            project, sbname, ebname = sdm_name.split(".")[:3]
            mjd = float(".".join(sdm_name.split(".")[3:]))

            track_name = f"{sdm_name}.speclines.ms.split"

            order_dict[mjd] = {
                "project": project,
                "track_name": track_name
            }

        print("")
        print("")
        print(f"{target} {config}")
        # Print string in order:
        for idx, mjd in enumerate(sorted(order_dict.keys())):
            print(f"{target} {order_dict[mjd]['project']} all {config} {idx+1} {target}/{order_dict[mjd]['track_name']}")