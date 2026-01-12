import subprocess
import os

WBT_EXE = r"E:\pythonProject\spu\whitebox_tools.exe"
ESRI_PNTR_PATH = r"E:\pythonProject\spu\d8_esri_pointer.tif"
UAA_PATH = r"E:\pythonProject\spu\upstream_area_m2.tif"

os.makedirs(os.path.dirname(UAA_PATH), exist_ok=True)

subprocess.run(
    [
        WBT_EXE,
        '--run=D8FlowAccumulation',
        f'--input={ESRI_PNTR_PATH}',
        f'--output={UAA_PATH}',
        '--out_type=catchment area',   # stays ONE argument here
        '--pntr',
        '--esri_pntr',
        '-v'
    ],
    check=True
)

print("Wrote:", UAA_PATH)
