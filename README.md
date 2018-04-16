WHAT
====
Implementations of fast exact k-means algorithms as described in http://arxiv.org/abs/1602.02514 and implementations of turbo-charged mini-batch k-means as described in http://arxiv.org/pdf/1602.02934

for interfaces
- (LIB) Shared library with accompanying C++ header file
- (EX) Command-line exectuble
- (PY) Python library


REQUIREMENTS
============
Minimal installation requirements:
- C++ compiler supporting C++11
- Linux operating system

Optional but recommended:
- BLAS implementation, we recommend this one : http://www.openblas.net/

Specific to Python library:
- Python and Cython


CONFIGURATION
=============
In `Makefile`, set `USEBLAS` to either `NO` or `YES`
if `USEBLAS = YES`, then set `LIBBLASDIR`, `INCBLASDIR` (unless blas paths will be found automatically)


BUILDING
========
- For (LIB) and (EX) and (PY) : `make all`
- For (EX) : `make main`
- For (LIB) : `make lib`

USING
=====
(EX) If succesfully installed, you should find an executable in directory bin
Run the executable with -h flag to see the options

(LIB) You need to add lib directory to your LD_LIBRARY_PATH : put the following line in your ``~/.bashrc` file:
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/kmeans/lib
```
(PY) If successfully installed, you should be able to `import kmeans` when in directory lib.
To use from a different directory,

(a) as per (LIB), and

(b) add the path to lib to your python path, either by:
```
export PYTHONPATH=${PYTHONPATH}:/path/to/kmeans/lib
```
or directly in your python script :
```
import sys
sys.path.insert(0,'/path/to/kmeans/lib')
```
Example use is found in `examples/examples.py`



DOESN'T WORK?
=============
Please contact me at jnewling@idiap.ch
