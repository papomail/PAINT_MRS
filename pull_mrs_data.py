from pathlib import Path
import shutil



data_folder = Path('/Volumes/NMR/LWP_data/LWP_raw_data/LWP_data_3T/INSPIRE')
basis_folder = Path('/Users/papo/Sync/Projects/PAINT_MRS_CSI/3_0T_basis_threonine_no_MM')

sdat_files = [f for f in sorted(data_folder.rglob('*')) if ("act.sdat" in f.name.lower() and "csi" in f.name.lower())]
for f in sdat_files:
    target = Path.cwd() / Path('/'.join(f.parts[3:]))
    origin = f.resolve()
    print(f'{origin=}\n{target=}\n\n')
    shutil.copytree(origin, target)
    
    