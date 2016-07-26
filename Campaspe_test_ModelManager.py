import Groundwater
from ModelInterface.flopyInterface import flopyInterface
# MM is short for model manager

testMM = Groundwater.GWModelManager()

testMM.load_GW_model(r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\Campaspe_packaged.pkl")

print "************************************************************************"
print " Updating recharge boundary "

old_val = testMM.GW_build['Campaspe'].parameters.param['magic_rain']
new_val = 0.05

rch = {}
rch[0] = testMM.GW_build['Campaspe'].boundaries.bc['Rain_reduced']['bc_array'][0] / old_val * new_val

testMM.GW_build['Campaspe'].boundaries.create_model_boundary_condition('Rain_reduced', 'recharge', bc_static=True)
testMM.GW_build['Campaspe'].boundaries.assign_boundary_array('Rain_reduced', rch)

print " Include irrigation in the recharge array"




print " Updating river conditions "



print "************************************************************************"
print " Build and run MODFLOW model "


modflow_model = flopyInterface.ModflowModel(testMM.GW_build['Campaspe'])

modflow_model.runMODFLOW()

print " Return the stream-aquifer exchange for reaches as list "

print " Return the average depth to the GW table in areas as list "


#modflow_model.viewHeads()
#modflow_model.viewHeads2()