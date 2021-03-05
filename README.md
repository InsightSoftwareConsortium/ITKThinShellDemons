![Build Status](https://github.com/InsightSoftwareConsortium/ITKThinShellDemons/workflows/Build,%20test,%20package/badge.svg)
![PyPi Version](https://img.shields.io/pypi/v/itk-thinshelldemons.svg)
![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)

# ITKThinShellDemons

This module implements the Thin Shell Demons regularization proposed in

> Thin Shell Demons  
> Zhao Q, Price T, Pizer S, Niethammer M, Alterovitz R, Rosenman J  
> MIUA 2015  


There are two implementations of this approach available:
1. (THDv4 ) An implementation taking advantage of the pointset to pointset registration
   capabilities of the ITK v4 registration framework.
2. (THD) An implementation started by Qingyu Zhao and updated to fix some remaiing task
   and run with ITK 5


## V4 version (THDv4)

The v4 version is implemnted in ThinShellDemonsMetricv4. With this approach different
transformation can be combined with the thin shell regularization. For examples see:
1. [Affine](./test/itkThinShellDemonsTestv4_Affine.cxx)
2. [DisplacementField](./test/itkThinShellDemonsTestv4_Displacement.cxx)
3. [SyN Diffeomorphism](./test/itkThinShellDemonsTestv4_SyN.cxx)

Please be aware that the regularization is implented on the fixed mesh. This is due to the 
nomenclature in of the moving transform being computed from fixed image to moving image. For 
resampling the moving image is pulled back to the fixed image domain. For point sets
it's easier to compute the oush forward than the pull back (for some transformation the 
inverse might be difficult compute or not available). Hence, the nomenclature matches
the ITK apprach for images and computes a transform from fixed to moving domain. Howver to
register a point set it is the fixed point set that is transformed to the moving point set 
domain and it is computationally more efficent to regularize on the fixed mesh.

**The implementation currently supports the thin shell regularization. Missing features are a 
confidence matching weighting and geomtric feature matching (currently only point distances 
are used).**


## Original version (THD)

This is an implementation that fits within the pre v4 ITK registration workflow.
The implementation is not very generic at this point and miss details to be able
to interact with other registration components in teh ITK registration ecosystem.
However, different optimizers can be used within this method.

**The implementation supports the thin shell regularization, geometric feature matching
(currently updated at each iteration during optimtization). Missing features are a
confidence matching weighting.**

For an example see [ThinShellDemonsTest](./test/itkThinShellDemonsTest.cxx)

The implementation has the followings components:

1. Transformation (itkMeshDisplacementTransform)

The class "MeshDisplacementTransformation" defines a finite dimensional vector
space on mesh vertices. Its private member m_VectorField is a 1D parameter
array in the form of [x_1,y_1,z_1,x_2,y_2,z_2,...], where the subscript denote
the index of the vertex.  A mesh has to be initially associated with a transformation
object to serve as a template. The template essentially designates the number of
vertices, so that m_VectorField can be initialized and allocated with a correct
size (# of vertices * 3)

2. Metric (itkMeshToMeshMetric -> itkThinShellDemonsMetric)

MeshToMeshMetric: This class is templated over the type of PointsetToPointsetMetric.
This class serves as the basis for all kinds of mesh-to-mesh metrics (in some sense
computing the similarity between two meshes). It expects a mesh-to-mesh transformation
to be plugged in. This class computes an objective function value (also with its
derivative w.r.t. the transformation parameters) that measures a registration
metric between the fixed mesh and the moving mesh.

ThinShellDemonsMetric: This Class inherits the basic MeshToMeshMetric. It expects a
mesh-to-mesh transformaton to be plugged in. This class computes a metric value, which
is a combination of geometric feature matching quality and the Thin Shell deformation
Energy. This metric computation part (objective function) is the core of the Thin Shell
Demons algorithm. When initializing a metric object of this class with two meshes,
the metric object first pre-computes geometric feature matching between the two meshes.
~~The matching results stay the same during the optimization process~~. *This was changed
to be updated during optimization and curvature matching was implemnted.*

3. Optimizer

Different thin shell energy approximation leads to different objective function
formulations, thereby requiring different optimizers. The current objective
function adopts a quadratic form. Therefore, Conjugate Gradient or LBFGS is a
preferable optimizer.

4. Registration Method (itkMeshToMeshRegistrationMethod)

This class is templated over the pointset-to-pointset registration method. Users
will create an object of this class to perform Thin Shell Demons. See the test
example for usage.


## Requirements

Compile against ITK with ITKVtkGlue module built.

## Authors
Qingyu Zhao (original version)  
Samuel Gerber (v4 version and updateds to original)
