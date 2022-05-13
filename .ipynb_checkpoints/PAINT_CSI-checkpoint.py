# Import modules
from pathlib import Path
import pwd
import suspect
import numpy as np
import datetime 
import pandas as pd
from tqdm import tqdm 
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from fpdf import FPDF
from pdf2image import convert_from_path
import warnings 
import pickle
from pprint import pprint


# Definitions
from suspect import MRSData, rotation_matrix, transformation_matrix

import os
import numpy

spar_types = {
    "floats": ["ap_size", "lr_size", "cc_size", "ap_off_center", "lr_off_center",
               "cc_off_center", "ap_angulation", "lr_angulation", "cc_angulation",
               "image_plane_slice_thickness", "slice_distance", "spec_col_lower_val",
               "spec_col_upper_val", "spec_row_lower_val", "spec_row_upper_val",
               "spectrum_echo_time", "echo_time"],
    "integers": ["samples", "rows", "synthesizer_frequency", "offset_frequency",
                 "sample_frequency", "echo_nr", "mix_number", "t0_mul_direction",
                 "repetition_time", "averages", "volumes",
                 "volume_selection_method", "nr_of_slices_for_multislice",
                 "spec_num_col", "spec_num_row", "num_dimensions", "TSI_factor",
                 "spectrum_inversion_time", "image_chemical_shift",
                 "t0_mu1_direction"],
    "strings": ["scan_id", "scan_date", "patient_name", "patient_birth_date",
                "patient_position", "patient_orientation", "nucleus",
                "volume_selection_enable", "phase_encoding_enable", "t1_measurement_enable",
                "t2_measurement_enable", "time_series_enable", "Spec.image in plane transf",
                "spec_data_type", "spec_sample_extension", "spec_col_extension",
                "spec_row_extension", "echo_acquisition", "resp_motion_comp_technique",
                "de_coupling", "equipment_sw_verions", "examination_name"],
}


def load_sdat(sdat_filename, spar_filename=None, spar_encoding=None):
    # if the spar filename is not supplied, assume it is in the same folder as
    # the sdat and only differs in the extension
    if spar_filename is None:
        path, ext = os.path.splitext(sdat_filename)
        # match the capitalisation of the sdat extension
        if ext == ".SDAT":
            spar_filename = path + ".SPAR"
        elif ext == ".sdat":
            spar_filename = path + ".spar"

    with open(spar_filename, 'r', encoding=spar_encoding) as fin:
        parameter_dict = {}
        for line in fin:
            # ignore empty lines and comments starting with !
            if line != "\n" and not line.startswith("!"):
                key, value = map(str.strip, line.split(":", 1))
                if key in spar_types["floats"]:
                    parameter_dict[key] = float(value)
                elif key in spar_types["integers"]:
                    parameter_dict[key] = int(value)
                elif key in spar_types["strings"]:
                    parameter_dict[key] = value
                else:
                    pass
                    #print("{} : {}".format(key, value))

    dt = 1 / parameter_dict["sample_frequency"]

    with open(sdat_filename, 'rb') as fin:
        raw_bytes = fin.read()

    floats = _vax_to_ieee_single_float(raw_bytes)
    data_iter = iter(floats)
    complex_iter = (complex(r, -i) for r, i in zip(data_iter, data_iter))
    raw_data = numpy.fromiter(complex_iter, "complex64")
    raw_data = numpy.reshape(raw_data, (parameter_dict["rows"], parameter_dict["samples"])).squeeze()

    # calculate transformation matrix
    voxel_size = numpy.array([parameter_dict["lr_size"],
                              parameter_dict["ap_size"],
                              parameter_dict["cc_size"]])
    position_vector = numpy.array([parameter_dict["lr_off_center"],
                                   parameter_dict["ap_off_center"],
                                   parameter_dict["cc_off_center"]])

    A = numpy.eye(3)
    for a,ang in enumerate(["lr_angulation", "ap_angulation", "cc_angulation"]):
        axis = numpy.zeros(3)
        axis[a] = 1
        A = A @ rotation_matrix(parameter_dict[ang]/180*numpy.pi,axis)
    e1 = A[:,0]
    e1 = e1 / numpy.linalg.norm(e1)
    e2 = A[:,1]
    e2 = e2 / numpy.linalg.norm(e2)
    
    transform = transformation_matrix(e1,
                                    e2,
                                    position_vector,
                                    voxel_size)


    return MRSData(raw_data,
                   dt,
                   parameter_dict["synthesizer_frequency"] * 1e-6,
                   te=parameter_dict["echo_time"],
                   tr=parameter_dict["repetition_time"],
                   transform=transform), parameter_dict


def _vax_to_ieee_single_float(data):
    """Converts a float in Vax format to IEEE format.
    Data should be a single string of chars that have been read in from
    a binary file. These will be processed 4 at a time into float values.
    Thus the total number of byte/chars in the string should be divisible
    by 4.
    Notes
    -----
    Based on VAX data organization in a byte file, we need to do a bunch of
    bitwise operations to separate out the numbers that correspond to the
    sign, the exponent and the fraction portions of this floating point
    number
    role :      S        EEEEEEEE      FFFFFFF      FFFFFFFF      FFFFFFFF
    bits :      1        2      9      10                               32
    bytes :     byte2           byte1               byte4         byte3
    Returns
    -------
    f : array
        Contains floats in IEEE format
    """
    f = []
    nfloat = int(len(data) / 4)
    for i in range(nfloat):

        byte2 = data[0 + i*4]
        byte1 = data[1 + i*4]
        byte4 = data[2 + i*4]
        byte3 = data[3 + i*4]

        # hex 0x80 = binary mask 10000000
        # hex 0x7f = binary mask 01111111

        sign = (byte1 & 0x80) >> 7
        expon = ((byte1 & 0x7f) << 1) + ((byte2 & 0x80) >> 7)
        fract = ((byte2 & 0x7f) << 16) + (byte3 << 8) + byte4

        if sign == 0:
            sign_mult = 1.0
        else:
            sign_mult = -1.0

        if 0 < expon:
            # note 16777216.0 == 2^24
            val = sign_mult * (0.5 + (fract/16777216.0)) * pow(2.0, expon - 128.0)
            f.append(val)
        elif expon == 0 and sign == 0:
            f.append(0)
        else:
            f.append(0)
            # may want to raise an exception here ...

    return f


def csi_ref(series):
    for serie in series:
        num_objects = len(serie.sp_objects)
        if num_objects > 127:     
            s1,s2 = serie.sp_objects[0].csi_mat[0], serie.sp_objects[0].csi_mat[1]
            print(f'csi matrix size {[s1,s2]}')
            mymat = np.zeros([s1,s2])
            for ii in tqdm(range(int(num_objects/2)),f'checking CSI data'):
                col = ii//s1
                row = ii%s1
#                 print(f'current voxel {[row,col]}')
#                 dat = serie.sp_objects[ii].Kspace[0]
                ref = serie.sp_objects[ii+int(num_objects/2)].Kspace[0]
                metric = np.sum(np.abs(ref.flatten()))
                mymat[row,col] = metric  
                
            mymat = mymat/np.max(mymat)
            
            plt.figure()    
            c = plt.imshow(mymat, cmap ='Greens', interpolation ='nearest', origin ='lower') 
            plt.colorbar(c) 
            plt.title('CSI structure',fontweight ="bold") 
#             plt.show() 
            plt.savefig('csi_int_new.png')

            return mymat









act_filename = Path('/Users/patxi/Sync/MRdata/IoN Piglet/MRS RAY/LWP733_D4/CSI/C5985446_WIP_CSI_PB_auto_TR2S_SENSE_13_2_raw_act.SDAT')
ref_filename = Path('/Users/patxi/Sync/MRdata/IoN Piglet/MRS RAY/LWP733_D4/CSI/C5985446_WIP_CSI_PB_auto_TR2S_SENSE_13_2_raw_ref.SDAT')
# xml_filename = Path('/Users/patxi/Sync/MRdata/IoN Piglet/MRS RAY/LWP733_D4/CSI/C5985446_WIP CSI_PB_auto_TR2S_SENSE_13_2.xml') 

act, act_info = load_sdat(act_filename)
ref, ref_info = load_sdat(ref_filename)

plot = px.line(ref[0])


csi_mat = (8,8)




out_folder = Path.cwd()/'Tarquin_results'
out_folder.mkdir(parents=True, exist_ok=True)

fit_params = {
                'output_pdf': out_folder/'tarquin_SV_PRESS.pdf', 
                'basis_csv': Path.cwd()/'3_0T_basis_threonine_no_MM'
                }
 

from matplotlib import pyplot as plt
plt.plot(act[40])
plt.show()





fit_results = suspect.io.tarquin.process(act[40], ref[40], fit_params)

# fit_results = suspect.io.tarquin.process(act[50])



# print(f'{fit_params=}')
# print(f'{fit_results=}')

