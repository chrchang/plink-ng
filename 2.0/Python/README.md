This provides a basic Python API for pgenlib  (See [python_api.txt](python_api.txt) for details.)


### Build instructions
To build this library you will first need to clone the repository:

```
# clone repo
git clone https://github.com/chrchang/plink-ng
# go to python folder
cd plink-ng/2.0/Python 
```


Then install Cython and NumPy:
```
pip3 install "cython>=0.29.21" "numpy>=1.19.0"
```

and then build and install the package
```
python3 setup.py build_ext
python3 setup.py install
```


##### Example usage:
```
#write a 2 sample file
import numpy as np
import pgenlib as pg

with pg.PgenWriter("test.pgen".encode("utf-8"), 2, 3, False) as writer:
	writer.append_alleles(np.array([0,1,1,1],dtype=np.int32))
	writer.append_alleles(np.array([0,1,0,0],dtype=np.int32))
	writer.append_alleles(np.array([0,0,0,0],dtype=np.int32))

```
