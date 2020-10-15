import numpy as np
from hycom.info import read_field_names,read_field_grid_names

NAN_TH = 2**99  # Nan threshold

def subset_hycom_field(input_file: str, output_file:str, fields: list, layers=[]):
    """
    This function will create a new set of .a and .b files, with a subset number of fields and layers.
    input_file: str
        Complete path to the .a or .b hycom output file. It can also be the name of the file without the extension
    output_file: str
        Name of the file to use as output.
    fields: list
        List of strings with the fields names to keep.
    layers: list
        List of integers that represent the z-index of each of the layers to keep. If empty, all layers are kept.
    """
    # Selecting the proper name
    if input_file.endswith('.a') or input_file.endswith('.b'):
        input_file = input_file[:-2]
    if output_file.endswith('.a') or output_file.endswith('.b'):
        output_file = output_file[:-2]

    a_input_file = input_file+'.a'
    b_input_file = input_file+'.b'
    b_file = open(b_input_file, 'r')

    new_b_input_file = output_file+'.b'

    # Validate that the field names requested are available in the hycom file
    all_fields = read_field_names(b_input_file)
    if len(fields) == 0:
        fields = all_fields
        # print(F"Reading all the fields in the file: {fields}")
    if not(np.all([field in all_fields for field in fields])):
        print(F"Warning!!!!! Fields {[field for field in fields if not(field in all_fields)]} are not"
              F" in the hycom file {input_file}, removing them from the list.")
        fields = [field for field in fields if field in all_fields ]

    # Reading the header file (first 4 lines, just general info)
    b_file_lines = b_file.readlines()

    lon_size = int(b_file_lines[7].strip().split()[0])
    lat_size = int(b_file_lines[8].strip().split()[0])
    layer_size = lon_size*lat_size
    # size of each layer (it seems all the layers have the same size)
    npad = 4096-np.mod(layer_size, 4096)

    # Data of the new bfile
    new_b_file_lines = b_file_lines[:10]

    # Opening input and output files
    a_file = open(a_input_file, 'rb')
    a_output_file = open(F"{output_file}.a", 'wb')

    # Looking for the starting locations for each layer and each field
    field_loc = {field: [] for field in fields}
    for line_idx, cur_line in enumerate(b_file_lines[9:]):
        field = cur_line.split()[0].strip()
        if field in field_loc:
            if np.any(len(field_loc[field]) == np.array(layers)):
                new_b_file_lines.append(cur_line)
                offset = (line_idx-1) * (layer_size + npad) * 4
                a_file.seek(offset)
                a_output_file.write(a_file.read((layer_size + npad) * 4))
            field_loc[field].append(line_idx)

    file_output_b = open(F"{output_file}.b","w")
    file_output_b.write("".join(new_b_file_lines))

    # Closing both files
    a_file.close()
    b_file.close()
    a_output_file.close()

def read_hycom_fields(file_name: str, fields: list, layers=[], replace_to_nan=True):
    """
    Reads hycom files (.a and .b) and returns the desired fields in a dictionary.
        file_name: str
            Complete path to the .a or .b hycom output file. It can also be the name of the file without the extension
        fields: list
            List of strings with the fields names to read. If empty, the function reads all the fields
        layers: list
            List of integers that represent the z-index of the layers to read. If empty, all the layers will be read.
        replace_to_nan: boolean
            Indicates if the nan values should be replaced with numpy nan
    """
    # Selecting the proper name
    if file_name.endswith('.a') or file_name.endswith('.b'):
        file_name = file_name[:-2]

    a_file_name = file_name+'.a'
    b_file_name = file_name+'.b'
    b_file = open(b_file_name, 'r')

    # Validate that the field names requested are available in the hycom file
    all_fields = read_field_names(b_file_name)
    if len(fields) == 0:
        fields = all_fields
        # print(F"Reading all the fields in the file: {fields}")
    if not(np.all([field in all_fields for field in fields])):
        print(F"Warning!!!!! Fields {[field for field in fields if not(field in all_fields)]} are not"
              F" in the hycom file {file_name}, removing them from the list.")
        fields = [field for field in fields if field in all_fields ]

    # Reading the header file (first 4 lines, just general info)
    b_file_lines = b_file.readlines()

    hycom_ver = b_file_lines[4].strip().split()[0]
    exp_num = b_file_lines[5].strip().split()[0]
    lon_size = int(b_file_lines[7].strip().split()[0])
    lat_size = int(b_file_lines[8].strip().split()[0])
    layer_size = lon_size*lat_size
    # size of each layer (it seems all the layers have the same size)
    npad = 4096-np.mod(layer_size, 4096)

    # Printing basic information
    print(F"Hycom version: {hycom_ver}, Experiment: {exp_num}")
    for cur_line in range(3):
        print(b_file_lines[cur_line].strip())

    # Looking for the starting locations for each layer and each field
    field_loc = {field: [] for field in fields}
    for line_idx, cur_line in enumerate(b_file_lines[9:]):
        field = cur_line.split()[0].strip()
        if field in field_loc:
            field_loc[field].append(line_idx)

    # Counting the number of layers for each field.
    num_layers = {field: len(field_loc[field]) for field in fields}
    print(F"Dims lon: {lon_size}, lat: {lat_size}")

    # Read layers for each field
    a_file = open(a_file_name, 'rb')

    # Define the layers that are going to be retrieved for each field
    if len(layers) != 0:
        layers_per_field = {field: [layer for layer in layers if layer in range(num_layers[field])] for field in fields}
    else:
        layers_per_field = {field: range(num_layers[field]) for field in fields}

    # Create the dictionary that will contain the np arrays with the fields information
    np_fields = {field: np.zeros((len(layers_per_field[field]), lat_size, lon_size)) for field in fields}

    for field in fields:
        print(F"\tReading layers {layers_per_field[field]} for field {field}. Total layers: {num_layers[field]}")

    # For each field read the proper section for each layer, from the binary file
    for field in fields:
        for cur_layer_idx, cur_layer in enumerate(layers_per_field[field]):
            offset = (field_loc[field][cur_layer]-1) * (layer_size+npad)*4
            a_file.seek(offset)
            cur_layer_data = np.fromfile(file_name+'.a', dtype='>f', count=layer_size, offset=offset)
            if replace_to_nan:
                cur_layer_data[cur_layer_data > NAN_TH] = np.nan
            np_fields[field][cur_layer_idx, :, :] = np.reshape(cur_layer_data, (lat_size, lon_size))

    # Closing both files
    a_file.close()
    b_file.close()
    return np_fields


def read_hycom_grid(file_name: str, fields: list,  replace_to_nan=True):
    """
    Reads hycom grid files (.a and .b) and returns the desired fields in a dictionary.
        file_name: str
            Complete path to the .a or .b hycom output file. It can also be the name of the file without the extension
        fields: list
            List of strings with the fields names to read. If empty, the function reads all the fields
        replace_to_nan: boolean
            Indicates if the nan values should be replaced with numpy nan
    """
    # Selecting the proper name
    if file_name.endswith('.a') or file_name.endswith('.b'):
        file_name = file_name[:-2]

    a_file_name = file_name+'.a'
    b_file_name = file_name+'.b'
    b_file = open(b_file_name, 'r')

    # Validate that the field names requested are available in the hycom file
    all_fields = read_field_grid_names(b_file_name)
    if len(fields) == 0:
        fields = all_fields
        # print(F"Reading all the fields in the file: {fields}")
    if not(np.all([field in all_fields for field in fields])):
        print(F"Warning!!!!! Fields {[field for field in fields if not(field in all_fields)]} are not"
              F" in the hycom file {file_name}, removing them from the list.")
        fields = [field for field in fields if field in all_fields ]

    # Reading the header file (first 4 lines, just general info)
    b_file_lines = b_file.readlines()

    lon_size = int(b_file_lines[0].strip().split()[0])
    lat_size = int(b_file_lines[1].strip().split()[0])
    layer_size = lon_size*lat_size
    # size of each layer (it seems all the layers have the same size)
    npad = 4096-np.mod(layer_size, 4096)


    # Looking for the starting locations for each layer and each field
    field_loc = {field: [] for field in fields}
    for line_idx, cur_line in enumerate(b_file_lines[3:]):
        field = cur_line.split()[0].strip()
        field = field[:-1]
        #print(field)
        if field in field_loc:
            field_loc[field].append(line_idx)

    # Counting the number of layers for each field.
    num_layers = {field: len(field_loc[field]) for field in fields}
    print(F"Dims lon: {lon_size}, lat: {lat_size},num_layers: {num_layers}")

    # Read layers for each field
    a_file = open(a_file_name, 'rb')

    # Define the layers that are going to be retrieved for each field
    layers=[0]
    if len(layers) != 0:
        layers_per_field = {field: [layer for layer in layers if layer in range(num_layers[field])] for field in fields}
    else:
        layers_per_field = {field: range(num_layers[field]) for field in fields}

    # Create the dictionary that will contain the np arrays with the fields information
    np_fields = {field: np.zeros((lat_size, lon_size)) for field in fields}

    for field in fields:
        print(F"\tReading layers {layers_per_field[field]} for field {field}. Total layers: {num_layers[field]}")

    # For each field read the proper section for each layer, from the binary file
    for field in fields:
        for cur_layer_idx, cur_layer in enumerate(layers_per_field[field]):
            offset = (field_loc[field][cur_layer]) * (layer_size+npad)*4
            a_file.seek(offset)
            cur_layer_data = np.fromfile(file_name+'.a', dtype='>f', count=layer_size, offset=offset)
            if replace_to_nan:
                cur_layer_data[cur_layer_data > NAN_TH] = np.nan
            np_fields[field][:, :] = np.reshape(cur_layer_data, (lat_size, lon_size))

    # Closing both files
    a_file.close()
    b_file.close()

    return np_fields

def read_hycom_depth(file_name:str, idm, jdm, replace_to_nan=True):
   """
   Read depth file from HYCOM (2D field)
   file_name: string
         name of the file with path
   idm, jdm: int
         size of the 2D field
   replace_to_nan: bool
         if True replace huge values by NaN
   """
   layer_size=idm*jdm
   depth=np.fromfile(file_name,dtype='>f',count=layer_size)
   if replace_to_nan:
      depth[depth > NAN_TH] = np.nan
   bathy=np.reshape(depth,(jdm,idm))
   print(F"Bathymetry Min value:{np.nanmin(bathy)}, Max value:{np.nanmax(bathy)}")
   return bathy

def write_hycom_depth(field,file_name:str,description:str):
   """
   write HYCOM depth file
   field: 2 D array of bathymetry
   file_name: str name of the file
   description: str first line of b file
   """
   # size of each layer
   idm=field.shape[1] ; jdm=field.shape[0]
   layer_size=idm*jdm
   npad = 4096-np.mod(layer_size, 4096)

   # Selecting the proper name
   if file_name.endswith('.a') or file_name.endswith('.b'):
      file_name = file_name[:-2]

   a_file_name = file_name+'.a'
   b_file_name = file_name+'.b'

   ## change precision of field
   field32=np.array(field,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout=open(a_file_name,'wb')
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   fout.close()

   # writing .b file
   field32[field32 > 1e5]=np.nan
   b_file=[description,F"idm/jdm= {idm} {jdm}","   ","   ","   ",\
           F"min, max depth = {np.nanmin(field32)} {np.nanmax(field32)}"]

   fbout=open(b_file_name,'w')
   fbout.writelines("%s\n" % i for i in b_file)
   fbout.close()

def write_hycom_grid(file_name:str,plon, plat, qlon,  qlat,  ulon, ulat, \
		    vlon, vlat, pang,pscx, pscy, qscx,  qscy, uscx, uscy,  vscx, vscy, cori, pasp,mapflg):
   """
   write HYCOM grid file
   file_name: str name of the file
   fields to provide
   """
   # size of each layer
   idm=plon.shape[1] ; jdm=plon.shape[0]
   layer_size=idm*jdm
   npad = 4096-np.mod(layer_size, 4096)

   # Selecting the proper name
   if file_name.endswith('.a') or file_name.endswith('.b'):
      file_name = file_name[:-2]

   a_file_name = file_name+'.a'
   b_file_name = file_name+'.b'

   ## open file
   fout=open(a_file_name,'wb')

   ## plon
   field32=np.array(plon,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## plat
   field32=np.array(plat,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## qlon
   field32=np.array(qlon,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## qlat
   field32=np.array(qlat,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## ulon
   field32=np.array(ulon,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## ulat
   field32=np.array(ulat,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## vlon
   field32=np.array(vlon,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## vlat
   field32=np.array(vlat,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## pang
   field32=np.array(pang,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## pscx
   field32=np.array(pscx,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## pscy
   field32=np.array(pscy,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## qscx
   field32=np.array(qscx,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## qscy
   field32=np.array(qscy,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## uscx
   field32=np.array(uscx,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## uscy
   field32=np.array(uscy,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## vscx
   field32=np.array(vscx,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## vscy
   field32=np.array(vscy,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))


   ## cori
   field32=np.array(cori,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   ## pasp
   field32=np.array(pasp,dtype='>f4') ## float32 big_endian
   ## writing .a file
   fout.write(np.reshape(field32,layer_size))
   if (npad != 4096):
      fout.write(np.zeros([npad],dtype='>f4'))

   fout.close()


   # writing .b file
   b_file=[F" {idm} 'idm   ' = longitudinal array size",\
           F" {jdm} 'jdm   ' = latitudinal array size",\
           F"  {mapflg}  'mapflg' = map flag (0=mercator,10=panam,12=ulon-panam)",\
           "plon:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(plon),np.nanmax(plon)),\
           "plat:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(plat),np.nanmax(plat)),\
           "qlon:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(qlon),np.nanmax(qlon)),\
           "qlat:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(qlat),np.nanmax(qlat)),\
           "ulon:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(ulon),np.nanmax(ulon)),\
           "ulat:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(ulat),np.nanmax(ulat)),\
           "vlon:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(vlon),np.nanmax(vlon)),\
           "vlat:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(vlat),np.nanmax(vlat)),\
           "pang:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(pang),np.nanmax(pang)),\
           "pscx:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(pscx),np.nanmax(pscx)),\
           "pscy:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(pscy),np.nanmax(pscy)),\
           "qscx:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(qscx),np.nanmax(qscx)),\
           "qscy:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(qscy),np.nanmax(qscy)),\
           "uscx:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(uscx),np.nanmax(uscx)),\
           "uscy:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(uscy),np.nanmax(uscy)),\
           "vscx:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(vscx),np.nanmax(vscx)),\
           "vscy:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(vscy),np.nanmax(vscy)),\
           "cori:  min,max = {0:15.5E} {1:15.5E}".format(np.nanmin(cori),np.nanmax(cori)),\
           "pasp:  min,max = {0:15.5f} {1:15.5f}".format(np.nanmin(pasp),np.nanmax(pasp))]

   fbout=open(b_file_name,'w')
   fbout.writelines("%s\n" % i for i in b_file)
   fbout.close()

def write_hycom_relax(filename:str,fieldname:str,density,field):
   """
   subroutine to write HYCOM relax files
   filename:string of  absolute path+filename
   fieldname: string 'intf' or 'temp' or 'saln'
   density: 1-D float array density layer values
   field: 2-D float array
   """
   # size of each layer
   idm=field.shape[1] ; jdm=field.shape[0]
   kdm=field.shape[2] ; tdm=field.shape[3]
   layer_size=idm*jdm
   npad = 4096-np.mod(layer_size, 4096)

   # Selecting the proper name
   if filename.endswith('.a') or filename.endswith('.b'):
      filename = filename[:-2]

   a_file_name = filename+'.a'
   b_file_name = filename+'.b'

   if (fieldname == 'intf'):
      long_name='Interface Depths' ; short_name = ' int '

   if (fieldname == 'temp'):
      long_name = 'Potential Temperature' ; short_name = ' tem '

   if (fieldname == 'saln'):
      long_name = 'Salinity' ;  short_name = ' sal '

   ## open file
   fout=open(a_file_name,'wb')

   ## relax .a file
   for t in np.arange(tdm):
      for k in np.arange(kdm):
         field32=np.array(field[:,:,k,t],dtype='>f4') ## float32 big_endian
         ## writing .a file
         fout.write(np.reshape(field32,layer_size))
         if (npad != 4096):
            fout.write(np.zeros([npad],dtype='>f4'))


   fout.close()

   # writing .b file
   field[field > 1e20]=np.nan
   fbout=open(b_file_name,'w')
   fbout.write(long_name+"\n")
   fbout.write("  \n")
   fbout.write("  \n")
   fbout.write("  \n")
   fbout.write(F"i/jdm = {idm},{jdm}  \n")
   for t in np.arange(tdm):
      for k in np.arange(kdm):
          fbout.write(short_name+":  month,layer,dens,range = {0:02d} {1:02d} {2:6.3f} {3:12.5E} {4:12.5E} \n".\
             format(t+1,k+1,density[k],np.nanmin(field[:,:,k,t]),np.nanmax(field[:,:,k,t])))


   fbout.close()
