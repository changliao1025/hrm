from pye3sm.mosart.mesh.mosart_extract_contributing_cells_by_outlet_id import mosart_extract_contributing_cells_by_outlet_id
from pye3sm.mosart.mesh.mosart_extract_by_cellid import mosart_extract_by_cellid
from pyearth.toolbox.reader.text_reader_string import text_reader_string

lCellID = 16060

#debug
#lCellID = 18345

sFilenamae_mosart_in='/compyfs/inputdata/rof/mosart/MOSART_NLDAS_8th_20160426.nc'
lCellID_outlet_in=lCellID
sFilename_mosart_out='/qfs/people/liao313/workspace/python/hrm/data/output/mosart/mosart_mississippi.nc'
sFilename_cellid_out='/qfs/people/liao313/workspace/python/hrm/data/output/mosart/mosart_mississippi_cellid.txt'
#save the id to a file
iFlag_save_to_file = 0
if iFlag_save_to_file == 1:
    aCellID = mosart_extract_contributing_cells_by_outlet_id(sFilenamae_mosart_in, lCellID_outlet_in)
    #save to a file
    with open(sFilename_cellid_out, 'w') as f:
        for item in aCellID:
            f.write("%s\n" % item)

else:
    #read the id from a file
    dummy = text_reader_string(sFilename_cellid_out)
    aCellID = dummy[:,0].astype(float).astype(int)

mosart_extract_by_cellid(sFilenamae_mosart_in, sFilename_mosart_out, aCellID)