This provides a basic Python API for pgenlib; see [python_api.txt](python_api.txt) for details.


### Build instructions
PyPI: `pip install Pgenlib`

To build from GitHub instead, clone the repository:

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
python3 setup.py build_clib build_ext -i
python3 setup.py install
```

You can test the package with `pytest`.

##### Example usage:
```
#write a 2 sample file
import numpy as np
import pgenlib as pg

with pg.PgenWriter("test.pgen".encode("utf-8"), 2, variant_ct=3, nonref_flags=False) as writer:
	writer.append_alleles(np.array([0,1,1,1],dtype=np.int32))
	writer.append_alleles(np.array([0,1,0,0],dtype=np.int32))
	writer.append_alleles(np.array([0,0,0,0],dtype=np.int32))

```

See tests/test_pgenlib.py for more sophisticated examples.
