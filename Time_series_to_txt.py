import numpy as np

ts_length = 10000000
binary_path = f'set_files/BINARY_NAME.dat' # path to the binary file
txt_path = f'set_files/TEXT_NAME.txt' # path to the save location of the text file 

with open(binary_path,'r',encoding="ISO-8859-1") as f, open(txt_path, 'ab') as f1:
    name = f.read(1000) # skip header containg label
    data = np.fromfile(f, dtype=np.uint16,count = ts_length)
    np.savetxt(f1, data)

f.close()
f1.close()