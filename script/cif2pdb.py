import sys
from iotbx import cif

filename = sys.argv[1].split('.')[0]
xray_structure = cif.reader(file_path=filename+'.cif').build_crystal_structures()
xray_structure = xray_structure[list(xray_structure)[0]]
#xray_structure.show_summary().show_scatterers()
out_string = xray_structure.as_pdb_file()

with open(filename+'.pdb','w') as fp:
    fp.write(out_string)
fp.close()
