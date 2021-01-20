# Installation Guide

## Prerequisites:

`liqa` has been tested on python 3.7.x. If you don’t know the version of python you can check it by:
```python
>>>import platform
>>>platform.python_version()
#3.7.3
```
### The required packages for running liqa are listed below:
	* python packages:
		+ pysam
		+ numpy
		+ scipy
    	+ lifelines
    	+ copulas


## PyPI  
The recommended way to install **liqa** is using [pip](https://pip.pypa.io/en/stable/):

```bash
$ pip install liqa
```
This will pull and install the latest stable release from [PyPi](https://pypi.org/).

If you do not have permission (when you get a permission denied error), you should install liqa by 

```bash
$ pip install --user liqa
```

## Github  
Download the package from [Github](https://github.com/WGLab/LIQA) and install it locally:

```bash
git clone https://github.com/WGLab/LIQA.git
cd LIQA
pip install .
```

## Test
You can check whether the installation in complete or not by:
```bash
$ liqa

#Please specify task (liqa -task <task>):

#        refgene:   preprocess reference file

#        quantify:   quantify isoform expression

#        diff:   detect differential splicing gene/isoform

```
