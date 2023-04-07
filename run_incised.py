import ccr
# import ccr.incised as inc

c_riv, c_well, flows = ccr.incised_vertical.model_output(base_dir='test1_newBCs_Sce0')

print('River concentration: ' + str(c_riv))
print('Well concentration: ' + str(c_well))
print(flows)
