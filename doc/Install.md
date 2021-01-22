# Installation Guide

## Prerequisites:

`liqa` has been tested on **python** 3.5, 3.6, 3.7 and **R** 3.5.2. If you don’t know the version of python you can check it by:
```bash
$ python --version
# Python 3.7.3
$ R --version
# R version 3.5.2
```

Following packages are required for running `liqa`: 

  - Python:
    * Pysam
    * Numpy
    * Scipy
    * Lifelines
  - R:
    * gcmr
    * betareg

Among the required packages, `math`, `sys`, `os`, `re` and `time` are included in the python standard library, meaning that they should be already installed with python.

## PyPI  
The recommended way to install `liqa` is using [pip](https://pip.pypa.io/en/stable/):

```bash
$ pip install liqa
```
This will pull and install the latest stable release from [PyPi](https://pypi.org/).

If you do not have permission (when you get a permission denied error), you should install `liqa` by 

```bash
$ pip install --user liqa
```

**Note**: you need to make sure that the `pip` is for python3，or we should install `liqa` by
```bash 
python3 -m pip install liqa 
#or
pip3 install liqa
```


## Github  
Download the package from [Github](https://github.com/WGLab/LIQA) and install it locally:

```bash
git clone https://github.com/WGLab/LIQA.git
cd LIQA
pip install .
```

## Test
You can check whether the installation is complete by:
```bash
$ liqa

#Please specify task (liqa -task <task>):
#
#        refgene:   preprocess reference file
#
#        quantify:   quantify isoform expression
#
#        diff:   detect differential splicing gene/isoform

```
