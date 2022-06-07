#!/usr/bin/env python3

'''Test changing subX file order of dimensions to see if I can use cdo operators
to make a grid text file to re-grid.'''

# #Open a random subX file. Get grid description though CDO
# aa = a.transpose('lead','S', 'model','Y','X')
# bb = aa.rename_dims({'lead':'time'})

# bb.drop('S')

# cc = bb.ETo[:,:,:,:,:]



# cc = bb.drop_dims('S')
# dd = cc.drop_dims('model')



# #Convert to an xarray object
# var_OUT = xr.Dataset(
#     data_vars = dict(
#         ETo = (['time','Y','X'], bb.ETo[:,0,0,:,:]),
#     ),
#     coords = dict(
#         X = bb.X.values,
#         Y = bb.Y.values,
#         time = bb.lead.values),
#     attrs = dict(
#         Description = 'Reference crop evapotranspiration (mm/day)'),
# )                    

# #Save as a netcdf for later processing
# var_OUT.to_netcdf(path = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/test_/test_change_dim.nc4', mode ='w')
# print(f'Saved file into {home_dir}.')
