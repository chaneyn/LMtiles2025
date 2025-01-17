import h5py
import netCDF4 as nc
import geopandas
import os
import numpy as np
import scipy.sparse as sparse

def create_database_mp(grp,ID,X,Y):
 #open access to Duke HB database for macroscale polygon
 fpduke = nc.Dataset('/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/GFDL_TEST/experiments/simulations/baseline/%d/input_file.nc' % ID)
 #create macroscale polygon group
 #mpgrp = grp.create_group("tile:1,is:%d,js:%d" % (X+1,Y+1))
 mpgrp = grp.create_group("%d" % (ID,))
 #metadata
 mtdgrp = mpgrp.create_group("metadata")
 mtdgrp['ilat'] = X
 mtdgrp['ilon'] = Y
 mtdgrp['latitude'] = fpduke['metadata'].latitude
 mtdgrp['longitude'] = fpduke['metadata'].longitude
 mtdgrp['frac'] = (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np.float64)
 mtdgrp['tid'] = (fpduke['parameters']['hru'][:]+1).astype(np.int32) #temporary
 mtdgrp['tile'] = (fpduke['parameters']['hru'][:]).astype(np.int32)
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
 sgrp['dat_emis_sat'] = 1.1e+06*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_heat_capacity_dry
 sgrp['dat_heat_capacity_dry'] = 0.05*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_k_sat_ref
 sgrp['dat_k_sat_ref'] = (fpduke['parameters']['SATDK'][:,0]).astype(np.float64) #units???
 #dat_psi_sat_ref
 sgrp['dat_psi_sat_ref'] = (fpduke['parameters']['SATPSI'][:,0]).astype(np.float64) #units??
 #dat_refl_dry_dif
 sgrp['dat_refl_dry_dif'] = 0.333*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_refl_dry_dir
 sgrp['dat_refl_dry_dir'] = 0.333*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_refl_sat_dif
 sgrp['dat_refl_sat_dif'] = 0.333*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #dat_refl_sat_dir
 sgrp['dat_refl_sat_dir'] = 0.333*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
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
 sgrp['dat_z0_momentum'] = 0.7*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) ###
 #depth_to_bedrock
 sgrp['depth_to_bedrock'] = (fpduke['parameters']['m'][:]).astype(np.float64)
 #frac
 sgrp['frac'] = (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np.float64) ###
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
 #hidx_j
 sgrp['hidx_j'] = (fpduke['parameters']['hband'][:]).astype(np.int32)+1
 nhband = int(np.max(fpduke['parameters']['hband'][:])+1)
 #hidx_k
 sgrp['hidx_k'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.int32) #temporary
 #irrigation
 sgrp['irrigation'] = 0.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #ksat_0cm
 sgrp['ksat_0cm'] = (fpduke['parameters']['SATDK'][:,0]).astype(np.float64) #units???
 #ksat_200cm
 sgrp['ksat_200cm'] = (fpduke['parameters']['SATDK'][:,-1]).astype(np.float64) #units???
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
 sgrp['tile'] = (fpduke['parameters']['hru'][:]).astype(np.int32) #temporary
 #tile_elevation
 sgrp['tile_elevation'] = (fpduke['parameters']['dem'][:]).astype(np.float64)
 #tile_hlsp_elev
 sgrp['tile_hlsp_elev'] = (fpduke['parameters']['hand'][:]).astype(np.float64)
 #tile_hlsp_frac
 sgrp['tile_hlsp_frac'] = (fpduke['parameters']['area'][:]/np.sum(fpduke['parameters']['area'][:])).astype(np.float64) ####!
 #####Convert watershed properties
 width = sparse.csr_matrix((fpduke['wmatrix_Basin1']['data'][:],
                            fpduke['wmatrix_Basin1']['indices'][:],
                            fpduke['wmatrix_Basin1']['indptr'][:]),
                            shape=(nhband,nhband),dtype=np.float64)
 width = width.todense() #the definition needs some work... currently assuming the hband connectivity defines the width
 print(width)
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
 #tile_hlsp_slope
 sgrp['tile_hlsp_slope'] = (fpduke['parameters']['slope'][:]).astype(np.float64)
 #tile_hlsp_width
 sgrp['tile_hlsp_width'] = (tile_width/tile_width[0]).astype(np.float64)
 #vegn
 sgrp['vegn'] = 1*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #wtd
 sgrp['wtd'] = 10.0*np.ones(fpduke['parameters']['hru'][:].size).astype(np.float64) #temporary
 #lake
 #lgrp = mpgrp.create_group("lake")
 #glacier
 #ggrp = mpgrp.create_group("glacier")
 #close netcdf file
 fpduke.close()

 return

#create output file
os.system('rm test.h5')
fp = h5py.File('test.h5', 'w')

#iterate through the different macroscale polygons
df = geopandas.read_file('/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/GFDL_TEST/data/shp/domain.shp')
nmp = len(df['ID'])
for imp in range(nmp):
    ID = df['ID'][imp]
    X = df['X'][imp]
    Y = df['Y'][imp]
    print(ID,X,Y)
    create_database_mp(fp,ID,X,Y)

#Close file
fp.close()
    

