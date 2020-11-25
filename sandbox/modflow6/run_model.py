import os
import matplotlib.pyplot as plt

from xmipy import XmiWrapper



mf6_dll = "d:\checkouts\Athesis\sandbox\modflow6\lib\libmf6.dll"
mf6_config_file = "d:\checkouts\Athesis\sandbox\modflow6\data\mfsim.nam"

model_dir = os.path.dirname(mf6_config_file)

mf6 = XmiWrapper(lib_path=mf6_dll, working_directory=model_dir)
mf6.set_int("ISTDOUTTOFILE", 0)



mf6.initialize(mf6_config_file)
mf6.update()

# get head and reshape top layer for plotting
head_tag = mf6.get_var_address("X", "HOOGHOUDT")
head = mf6.get_value_ptr(head_tag)

head_3d = head.reshape(1,258,258)
plt.imshow(head_3d[0,:,:])
plt.colorbar()
plt.show()

mf6.finalize()