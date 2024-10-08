import numpy as np
from datetime import datetime as dt



now = dt.now()
_ID = now.strftime(f"%d%m%Y%H%M%S")

#----------------SETTINGS----------------#
binary_path = f'set_files/BINARY_NAME.dat'  # path to simulated binary file
np_array_path = f"set_files/NUMPY_NAME{_ID}.npy" # the target path to the folder were the resulting NumPy file should be saved
histogram_import_number = 100 # number of histograms to save from the simulated binary file
#----------------SETTINGS----------------#
hist_dim = (60,60) # dimensions of 2D Dwell Time Histograms as they have been simulated with the 2D-Fit. Should not be changed


print("converting . . . ")

with open(binary_path,'r',encoding="ISO-8859-1") as f, open(np_array_path, 'ab') as f2:
    for i in range(histogram_import_number):
        name = f.read(1000)
        name = name.split('.')[0]
        name = name.split('/')[-1]  
        data = np.fromfile(f, dtype=np.float64, count=int(hist_dim[0] * hist_dim[1]))
        hist = np.transpose(data.reshape(hist_dim))
        arr = np.array([name, hist], dtype=object)
        np.save(f2, arr)
f2.close
f.close

print(">----------------COMPLETED----------------<")


