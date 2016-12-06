# RADis
Pipeline for performing in-silico RADseq experiments

## Installation
RADis is composed of Python and Shell scripts, and should therefore work across most systems. The Shell scripts rely on the Z shell (zsh), so that must be available on the system. All Python code should be compatible with either Python 2 or Python 3.

[Heng Li's version of GNU sort](https://github.com/lh3/foreign/tree/master/sort), with support for alphanumerical sorting of genome coordinate files, is distributed with this software. I've modified the Makefile so that the executable name does not conflict with the existing version of sort, in case you elect to move it to your system `$PATH`. You must compile this software for certain RADis software to work.

```shell
make
```
