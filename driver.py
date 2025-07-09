import h5py
import netCDF4 as nc
import geopandas
import os
import numpy as np
import sys
import json
import scipy.sparse as sparse

def create_database_mp(grp,ID,X,Y,md,cid_mapping):
 #open access to Duke HB database for macroscale polygon
 fpduke = nc.Dataset(os.path.join(md['rdir'],f'experiments/simulations/baseline/{ID}/input_file.nc')) # open the HB database
 #create macroscale polygon group
 mpgrp = grp.create_group("tile:1,is:%d,js:%d" % (X,Y))
 #mpgrp = grp.create_group("%d" % (ID,))
 #metadata
 mtdgrp = mpgrp.create_group("metadata")
 mtdgrp['ilat'] = X
 mtdgrp['ilon'] = Y
 mtdgrp['latitude'] = fpduke['metadata'].latitude
 mtdgrp['longitude'] = fpduke['metadata'].longitude
 mtdgrp['frac'] = (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np.float64)
 #mtdgrp['frac'][:] = np.ones(fpduke['parameters']['hru'][:].size)/fpduke['parameters']['hru'][:].size
 #TMP
 mtdgrp['tid'] = (fpduke['parameters']['hru'][:]+1).astype(np.int32) #temporary
 print((fpduke['parameters']['hru'][:]+1).astype(np.int32))
 mtdgrp['tile'] = (fpduke['parameters']['hru'][:]+1).astype(np.int32)
 mtdgrp['type'] = 3*np.ones(fpduke['parameters']['hru'][:].size).astype(np.int32) #temporary
 #soil
 sgrp = mpgrp.create_group("soil")
 #bl
 sgrp['bl'] = 0.05*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #br
 sgrp['br'] = 0.05*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #bsw
 sgrp['bsw'] = 0.05*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #bwood
 sgrp['bwood'] = 0.05*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_awc_lm2
 sgrp['dat_awc_lm2'] = 0.1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_chb
 sgrp['dat_chb'] = (fpduke['parameters']['BB'][:,0]).astype(np.float64) #units???
 #dat_emis_dry
 sgrp['dat_emis_dry'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_emis_sat
 sgrp['dat_emis_sat'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_heat_capacity_dry
 sgrp['dat_heat_capacity_dry'] = 1.1e+06*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_k_sat_ref
 sgrp['dat_k_sat_ref'] = (fpduke['parameters']['SATDK'][:,0]).astype(np.float64) #units???
 #dat_psi_sat_ref
 sgrp['dat_psi_sat_ref'] = (fpduke['parameters']['SATPSI'][:,0]).astype(np.float64) #units??
 #dat_refl_dry_dif
 sgrp['dat_refl_dry_dif'] = 0.333*np.ones((2,fpduke['parameters']['hru'][:].size)).astype(np.float64) ###
 #dat_refl_dry_dir
 sgrp['dat_refl_dry_dir'] = 0.333*np.ones((2,fpduke['parameters']['hru'][:].size)).astype(np.float64) ###
 #dat_refl_sat_dif
 sgrp['dat_refl_sat_dif'] = 0.333*np.ones((2,fpduke['parameters']['hru'][:].size)).astype(np.float64) ###
 #dat_refl_sat_dir
 sgrp['dat_refl_sat_dir'] = 0.333*np.ones((2,fpduke['parameters']['hru'][:].size)).astype(np.float64) ###
 #dat_tf_depr
 sgrp['dat_tf_depr'] = 2.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_thermal_cond_dry
 sgrp['dat_thermal_cond_dry'] = 0.21*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_thermal_cond_exp
 sgrp['dat_thermal_cond_exp'] = 5.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_thermal_cond_sat
 sgrp['dat_thermal_cond_sat'] = 1.5*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_thermal_cond_scale
 sgrp['dat_thermal_cond_scale'] = 0.5*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_thermal_cond_weight
 sgrp['dat_thermal_cond_weight'] = 0.7*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_w_sat
 sgrp['dat_w_sat'] = (fpduke['parameters']['MAXSMC'][:,0]).astype(np.float64)
 #dat_z0_momentum
 sgrp['dat_z0_momentum'] = 0.01*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #depth_to_bedrock
 sgrp['depth_to_bedrock'] = (fpduke['parameters']['m'][:]).astype(np.float64)
 #frac
 sgrp['frac'] = (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np.float64) ###
 (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np    .float64)
 #gw_hillslope_length
 sgrp['gw_hillslope_length'] = 1000*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_hillslope_relief
 sgrp['gw_hillslope_relief'] = 300*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_hillslope_zeta_bar
 sgrp['gw_hillslope_zeta_bar'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_perm
 sgrp['gw_perm'] = 4.77204e-16*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_res_time
 sgrp['gw_res_time'] = 5.184e+06*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_scale_length
 sgrp['gw_scale_length'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_scale_perm
 sgrp['gw_scale_perm'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_scale_relief
 sgrp['gw_scale_relief'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_scale_soil_depth
 sgrp['gw_scale_soil_depth'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #gw_soil_e_depth
 sgrp['gw_soil_e_depth'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###!
 #sgrp['gw_soil_e_depth'][:] = 14.0
 #hidx_j
 sgrp['hidx_j'] = (fpduke['parameters']['hband'][:]).astype(np.int32)+1
 #sgrp['hidx_j'][:] = 1
 nhband = int(np.max(fpduke['parameters']['hband'][:])+1)
 #hidx_k
 sgrp['hidx_k'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.int32) #temporary
 #irrigation
 sgrp['irrigation'] = 0.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #ksat_0cm
 sgrp['ksat_0cm'] = (fpduke['parameters']['SATDK'][:,0]).astype(np.float64) #units???
 #sgrp['ksat_0cm'][:] = 0.0109779
 #ksat_200cm
 sgrp['ksat_200cm'] = (fpduke['parameters']['SATDK'][:,-1]).astype(np.float64) #units???
 #sgrp['ksat_200cm'][:] = 0.000429401
 #landuse
 sgrp['landuse'] = 3*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #microtopo
 sgrp['microtopo'] = 1e+36*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #ntile
 ntile = len(fpduke['parameters']['hru'][:])
 sgrp['ntile'] = np.float64(ntile) ###!
 #pann
 sgrp['pann'] = 1000.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #rsa_exp_global
 sgrp['rsa_exp_global'] = 1.5*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64)
 #soil_depth
 sgrp['soil_depth'] = 2.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #tann
 sgrp['tann'] = 286.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #tile
 sgrp['tile'] = (fpduke['parameters']['hru'][:]).astype(np.int32) + 1#temporary
 #sgrp['tile'] = (fpduke['parameters']['hru'][:]).astype(np.int32) #temporary
 #tile_elevation
 sgrp['tile_elevation'] = (fpduke['parameters']['dem'][:]).astype(np.float64)
 #tile_hlsp_elev
 sgrp['tile_hlsp_elev'] = (fpduke['parameters']['hand'][:]).astype(np.float64)
 #tile_hlsp_frac
 sgrp['tile_hlsp_frac'] = (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np.float64) ####!
 #sgrp['tile_hlsp_frac'][:] = np.ones(fpduke['parameters']['hru'][:].size)/fpduke['parameters']['hru'][:].size
 #####Convert watershed properties
 width = sparse.csr_matrix((fpduke['wmatrix_Basin1']['data'][:],
                            fpduke['wmatrix_Basin1']['indices'][:],
                            fpduke['wmatrix_Basin1']['indptr'][:]),
                            shape=(nhband,nhband),dtype=np.float64)
 width = width.todense() #the definition needs some work... currently assuming the hband connectivity defines the width
 hbands = (fpduke['parameters']['hband'][:]).astype(np.int32)
 tile_area = fpduke['parameters']['area'][:]
 hband_width = np.copy(np.diagonal(width, offset=1))
 hband_width = np.concatenate((hband_width[0][np.newaxis],hband_width,hband_width[-1][np.newaxis]))
 hband_width = (hband_width[1:] + hband_width[0:-1])/2
 hband_area = []
 for hband in range(len(hband_width)):
     m = hbands == hband
     hband_area.append(np.sum(tile_area[m]))
 hband_area = np.array(hband_area)
 hband_length = hband_area/hband_width
 tile_length = np.zeros(tile_area.size)
 tile_width = np.zeros(tile_area.size)
 for hband in range(len(hband_width)):
     m = hbands == hband
     tile_length[m] = hband_length[hband]
     tile_width[m] = tile_area[m]/tile_length[m]
 tile_hpos = np.zeros(tile_area.size)
 for hband in range(len(hband_width)):
     m = hbands == hband
     if hband > 0:tile_hpos[m] = (np.sum(hband_length[0:hband]) + hband_length[hband]/2)/np.sum(hband_length)
     else:tile_hpos[m] = (hband_length[hband]/2)/np.sum(hband_length)
 
 #tile_hlsp_hpos
 sgrp['tile_hlsp_hpos'] = (tile_hpos).astype(np.float64)
 #tile_hlsp_length
 sgrp['tile_hlsp_length'] = (tile_length).astype(np.float64)
 sgrp['tile_hlsp_length'][:] = 100 #TMP
 #tile_hlsp_slope
 sgrp['tile_hlsp_slope'] = (fpduke['parameters']['slope'][:]).astype(np.float64)
 sgrp['tile_hlsp_slope'][:] = 0.01 #TMP
 #tile_hlsp_width
 sgrp['tile_hlsp_width'] = (tile_width/tile_width[0]).astype(np.float64)
 sgrp['tile_hlsp_width'][:] = 1 #TMP
 #vegn
 sgrp['vegn'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #wtd
 sgrp['wtd'] = 10.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #lake
 #lgrp = mpgrp.create_group("lake")
 #glacier
 #ggrp = mpgrp.create_group("glacier")
 #close netcdf file
 #river network
 rgrp = mpgrp.create_group("stream_network")
 rgrp['nc'] = fpduke['stream_network']['downstream_channels'][:].shape[0]
 rgrp['ninlets'] = fpduke['stream_network']['inlets'][:].shape[0]
 rgrp['noutlets'] = fpduke['stream_network']['outlets'][:].shape[0]
 for var in fpduke['stream_network'].variables:
    rgrp[var] = fpduke['stream_network'][var][:].T
 #Learn mapping of cid to gfdl cell structure
 rgrp['gfdl2cid'] = cid_mapping['gfdl2cid'][:].T
 rgrp['cid2gfdl'] = cid_mapping['cid2gfdl'][:].T
 #Learn mapping of ucids that the river can flow to downstream in a given time step
 ucids_downstream = np.unique(fpduke['stream_network']['downstream_channels'][:,1,:])
 ucids_downstream = ucids_downstream[ucids_downstream!=-1]
 ucids_downstream = ucids_downstream[ucids_downstream!=-9999]
 tmp = np.zeros(rgrp['cid2gfdl'][:].T.shape[0]).astype(np.int32)
 tmp[:] = -9999
 for i in range(ucids_downstream.size):
     tmp[ucids_downstream[i]-1] = i+1
 rgrp['ucids_downstream'] = tmp[:].T
 rgrp['nmps'] = cid_mapping['cid2gfdl'][:].shape[0]
 rgrp['nucids_d'] = ucids_downstream.size
 rgrp['cid'] = ID
 #Determine the maximum number of channels for all the cids that it empties out into
 maxnc = 0
 for ucid in ucids_downstream:
   cdir = '%s/%s' % (md['rdir'],'experiments/simulations/baseline/%d' % ucid)
   fp = nc.Dataset('%s/input_file.nc' % cdir)
   c_length = fp['stream_network']['length'][:]
   fp.close()
   if c_length.size > maxnc:maxnc = c_length.size
 downstream_c_length = np.zeros((ucids_downstream.size,maxnc))
 downstream_c_length[:] = -9999
 for ucid in ucids_downstream:
   cdir = '%s/%s' % (md['rdir'],'experiments/simulations/baseline/%d' % ucid)
   fp = nc.Dataset('%s/input_file.nc' % cdir)
   c_length = fp['stream_network']['length'][:]
   downstream_c_length[tmp[ucid-1]-1,:c_length.size] = c_length[:]
   fp.close()
 rgrp['downstream_c_length'] = downstream_c_length[:].T
 rgrp['nc_d'] = maxnc

 fpduke.close()

 return

mdfile = sys.argv[1] # metadata file
metadata = json.load(open(mdfile,'r')) # read in metadata

#create output file
os.system(f'rm {os.path.join(metadata["rdir"],r"ptiles.h5")}') # remove file if it exists
fp = h5py.File(os.path.join(metadata["rdir"],r"ptiles.h5"), 'w') # create file
grp = fp.create_group("grid_data") # create group

#iterate through the different macroscale polygons
df = geopandas.read_file(os.path.join(metadata["rdir"],r'data/shp/domain.shp')) # read in the domain shapefile
print(df)
xs = np.arange(np.max(df.X.values)).astype(np.int32)+1
ys = np.arange(np.max(df.Y.values)).astype(np.int32)+1
faces = np.arange(np.max(df.TILE.values)).astype(np.int32)+1
#gfdl cell structure 2 cid
gfdl2cid = np.zeros((np.max(faces),np.max(xs),np.max(ys))).astype(np.int32)
gfdl2cid[:] = -9999
for i in range(len(df)):
    gfdl2cid[df.TILE.values[i]-1,df.X.values[i]-1,df.Y.values[i]-1] = df.ID.values[i]
#cid to gfdl cell structure
cid2gfdl = np.zeros((df.ID.values.size,3))
cid2gfdl[:] = -9999
for i in range(len(df)):
    cid2gfdl[i][:] = [df.TILE.values[i],df.X.values[i],df.Y.values[i]]
cid_mapping = {'gfdl2cid':gfdl2cid,'cid2gfdl':cid2gfdl}
nmp = len(df['ID']) # number of macroscale polygons
for imp in range(nmp): # iterate through the macroscale polygons
    ID = df['ID'][imp] # macroscale polygon ID
    X = df['X'][imp] # macroscale polygon X index
    Y = df['Y'][imp] # macroscale polygon Y index
    print(ID,X,Y)
    create_database_mp(grp,ID,X,Y,metadata,cid_mapping) # create the macroscale polygon database

#Close file
fp.close() # close file
    

