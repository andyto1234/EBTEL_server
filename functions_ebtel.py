import numpy as np
from functions_pixels import *
from tqdm import tqdm
from functions_pickle import *
from scipy.io import readsav
import glob
import re
import sunpy
import multiprocessing as mp
from functools import partial


def get_var(date_string):
    with open(date_string+'length.txt', "r") as f:
        length = f.readlines()
        length = [float(l.replace('\n','')) for l in length]

    with open(date_string+'mean_field.txt', "r") as f:
        gauss = f.readlines()
        gauss = [float(l.replace('\n','')) for l in gauss]

    with open(date_string+'heating.txt', "r") as f:
        heating = f.readlines()
        heating = [float(l.replace('\n','')) for l in heating]
    flines = restore(date_string+'flines')
    print('Order of variables: length, gauss, heating')
    print(f'List length of loop length: {len(length)}; Magnetic field: {len(gauss)}; Heating: {len(heating)}')
    
    return length, gauss, heating, flines

def check_file(python_index, file_list):
    file_heating = [float('.'.join(('0',re.findall(r'\d+', i)[-2]))) for i in file_list]
    file_length = [int(re.findall(r'\d+', i)[-1]) for i in file_list]
    synth_file_from_idl = np.array(file_length)+(np.array(file_heating)*1e4)
    synth_file_real_data = np.array(length)+(np.array(heating)*1e4)

    idl_index = list(synth_file_from_idl).index(min(synth_file_from_idl, key=lambda x:abs(x-synth_file_real_data[python_index])))
    # print(f'real heating: {heating_a[python_index]}, idl heating:{file_heating[idl_index]}')
    # print(f'real length: {length_a[python_index]}, idl length:{file_length[idl_index]}')
    return idl_index


def fline_multi(length, flines, aia_submap, file_list, i):
    if length[i] > 5e6:
        idl_index = check_file(i, file_list)
        pix_x, pix_y = filter_pix(flines[i].coords, aia_submap)
        intensity = readsav(file_list[idl_index])['int'][445]
        print(i)
        return pix_x, pix_y, intensity
    else:
        pass


if __name__ == "__main__":
    date = '20110415/'
    aia_submap = restore(date+'aia_submap')
    eis_fixed = restore(date+'eis_map')
    aia = restore(date+'aia_map')
    pfss_in = restore(date+'pfss_in')

    blank_data = np.zeros(len(eis_fixed.data[0:])*len(eis_fixed.data[0])).reshape(len(eis_fixed.data[0:]),len(eis_fixed.data[0]))
    failed_list = []
    length, gauss, heating, flines = get_var(date)
    file_list_multi = glob.glob(date+"simulation_results/*.sav")

    with mp.Pool(processes = 50) as p:
        fline_partial = partial(fline_multi, length, flines, aia_submap, file_list_multi)
        results = p.map(fline_partial, range(len(flines)))
    for result in results:
        blank_data[result[1], result[0]] += result[2]

    # for i in tqdm(range(len(flines))):
    #     if length[i] > 5e6:
    #         idl_index = check_file(i, file_list)
    #         try:
    #             pix_x, pix_y = filter_pix(flines[i].coords, eis_fixed)
    #             blank_data[pix_y, pix_x] += readsav(file_list[idl_index])['int'][400]
    #         except:
    #             failed_list.append(i)

    synth_map_multi = sunpy.map.Map(blank_data, eis_fixed.meta)
    save(date+'synth_eis', synth_map_multi)
    save(date+'failed_list', file_list_multi)
