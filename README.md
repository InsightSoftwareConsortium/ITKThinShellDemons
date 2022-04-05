![Build Status](https://github.com/InsightSoftwareConsortium/ITKThinShellDemons/workflows/Build,%20test,%20package/badge.svg)
![PyPi Version](https://img.shields.io/pypi/v/itk-thinshelldemons.svg)
![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)

# ITKThinShellDemons

This module implements the Thin Shell Demons regularization proposed in

> Thin Shell Demons  
> Zhao Q, Price T, Pizer S, Niethammer M, Alterovitz R, Rosenman J  
> MIUA 2015  



> :warning: **This module requires to be compiled against an ITK version with**  
> - [PointSetToPointSetMetricWithIndexv4](https://github.com/InsightSoftwareConsortium/ITK/pull/2385)   
> 

<p align="center">
<img src="https://user-images.githubusercontent.com/1044135/158479969-7313ed94-c5fb-4803-ae8d-2c0631893664.png" width="400" height="300">
</p>

## V4 version (TSHDv4)

The v4 version is implemented in ThinShellDemonsMetricv4. With this approach different
transformation can be combined with the thin shell regularization. For C++ examples see:
1. [Affine](./test/itkThinShellDemonsTestv4_Affine.cxx)
2. [DisplacementField](./test/itkThinShellDemonsTestv4_Displacement.cxx)
3. [SyN Diffeomorphism](./test/itkThinShellDemonsTestv4_SyN.cxx)

For Python samples please see:
1. [Affine](./examples/test_tsd_affine.py)
2. [DisplacementField](./examples/test_tsd_displacement.py)
3. [SyN Diffeomorphism](./examples/test_tsd_syn.py)
4. [BSpline](./examples/test_tsd_bspline.ipynb)

Please be aware that the regularization is implemented on the fixed mesh. This is due to the 
nomenclature in of the moving transform being computed from fixed image to moving image. For 
re-sampling the moving image is pulled back to the fixed image domain. For point sets
it's easier to compute the push forward than the pull back (for some transformation the 
inverse might be difficult compute or not available). Hence, the nomenclature matches
the ITK approach for images and computes a transform from fixed to moving domain. However to
register a point set it is the fixed point set that is transformed to the moving point set 
domain and it is computationally more efficient to regularize on the fixed mesh.


## Requirements

Compile against ITK.

## Authors
Qingyu Zhao (original version)  
Samuel Gerber (v4 version and updates to original)  
Pranjal Sahu
