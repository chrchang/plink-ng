This provides a basic Python API for pgenlib  (See [python_api.txt](python_api.txt) for details.)


##### Build this with this.
Cython and NumPy must be installed.
```
python3 setup.py build_ext
[sudo] python3 setup.py install
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
