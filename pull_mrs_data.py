from pathlib import Path
import shutil



# data_folder = Path('/Volumes/NMR/LWP_data/LWP_raw_data/LWP_data_3T/INSPIRE')
# data_folder = Path('/Users/papo/Sync/MRdata/IoN_Piglet/MRS_RAY')
# basis_folder = Path('/Users/papo/Sync/Projects/PAINT_MRS_CSI/3_0T_basis_threonine_no_MM')
data_folder = Path('/Users/patxi/Desktop/LWP_Local')
# mrs_files = [f for f in sorted(data_folder.rglob('*')) if ("act.sdat" in f.name.lower() and "csi" in f.name.lower())]
sdat_files = [file for file in sorted(data_folder.rglob('*')) if ('.sdat' in str(file.resolve).lower() )]
# spar_files = [file.rename(file.with_suffix('.SPAIR')) for file in sdat_files]
spar_files = [file.with_suffix('.SPAR') for file in sdat_files]

mrs_files = [j for i in zip(sdat_files, spar_files) for j in i]

print(mrs_files)
# mrs_files = [f for f in sorted(data_folder.rglob('*')) if ("act." in f.name.lower() and "wip_31p_isis" in f.name.lower())]

for f in mrs_files:
    # print(f.parts[-3:])
    # target = Path.cwd() / 'INSPIRE_MRS_DATA' / Path('/'.join(f.parts[-3:]))
    target = Path.cwd() / 'INSPIRE_MRS_DATA' / f.parts[-3]

    target.mkdir(parents=True, exist_ok=True)
    origin = f.resolve()
    print(f'{origin=}\n{target=}\n\n')
    try:
        shutil.copy(origin, target)
    except:
        continue
    

