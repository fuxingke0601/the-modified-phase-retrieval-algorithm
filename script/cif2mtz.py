import sys
from iotbx import reflection_file_reader, mtz

filename = sys.argv[1].split('.')[0]
file = reflection_file_reader.any_reflection_file(file_name=filename+'.hkl')
arrays = file.as_miller_arrays()
array1 = arrays[0].average_bijvoet_mates()
array2 = arrays[1].average_bijvoet_mates()
dataset = array1.as_mtz_dataset(column_root_label='Ical', column_types='I')
dataset.add_miller_array(miller_array=array1.f_sq_as_f(),column_root_label='Fcal', column_types='F')
dataset.add_miller_array(miller_array=array2.discard_sigmas(),column_root_label='Iobs', column_types='I')
array3 = array2.french_wilson().discard_sigmas()
dataset.add_miller_array(miller_array=array3,column_root_label='Fobs', column_types='F')
dataset.mtz_object().write(filename+'.mtz')
