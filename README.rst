ITKThinShellDemons
=================================

.. image:: https://github.com/InsightSoftwareConsortium/ITKThinShellDemons/workflows/Build,%20test,%20package/badge.svg
    :alt:    Build Status

.. image:: https://img.shields.io/pypi/v/itk-thinshelldemons.svg
    :target: https://pypi.python.org/pypi/itk-thinshelldemons
    :alt: PyPI Version

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://github.com/InsightSoftwareConsortium/ITKThinShellDemons/blob/master/LICENSE
    :alt: License

Overview
--------

This is a module for [ITK](http://itk.org): The Insight Toolkit for Segmentation and
Registration. It is designed to work with the ITKv4
modular system and  to be places it ITK/Modules/External or uses as a
Remote module.


Getting Started
===============

The official ITKv4 Wiki documentation on adding an external module is here:
http://www.itk.org/Wiki/ITK_Release_4/Modularization/Add_an_external_module_(external_module)
http://www.itk.org/Wiki/ITK/Policy_and_Procedures_for_Adding_Remote_Modules


Thin Shell Demons External Module
=================================

Thin Shell Demons External Module is a mesh-to-mesh registration external module.
Reference: "Thin Shell Demons: Zhao Q, Price T, Pizer S, Niethammer M, Alterovitz R, Rosenman J, MIUA 2015

Following the standard ITK registration workflow, this module has the followings components (files)
1. Transformation (itkMeshDisplacementTransform)

	The class "MeshDisplacementTransformation" defines a finite dimensional vector space on mesh vertices. Its private member m_VectorField is a 1D parameter array in the form of [x_1,y_1,z_1,x_2,y_2,z_2,...], where the subscript denote the index of the vertex.

	A mesh has to be initially associated with a transformation object to serve as a template. The template essentially designates the number of vertices, so that m_VectorField can be initialized and allocated with a correct size (# of vertices * 3)
)

2. Metric (itkMeshToMeshMetric -> itkThinShellDemonsMetric)

	MeshToMeshMetric: This class is templated over the type of PointsetToPointsetMetric. This class serves as the basis for all kinds of mesh-to-mesh metrics (in some sense computing the similarity between two meshes). It expects a mesh-to-mesh transformation to be plugged in. This class computes an objective function value (also with its derivative w.r.t. the transformation parameters) that measures a registration metric between the fixed mesh and the moving mesh.

	ThinShellDemonsMetric: This Class inherits the basic MeshToMeshMetric. It expects a mesh-to-mesh transformaton to be plugged in. This class computes a metric value, which is a combination of geometric feature matching quality and the Thin Shell deformation Energy. This metric computation part (objective function) is the core of the Thin Shell Demons algorithm. When initializing a metric object of this class with two meshes, the metric object first pre-computes geometric feature matching between the two meshes. The matching results stay the same during the optimization process.

3. Optimizer

	Different thin shell energy approximation leads to different objective function formulations, thereby requiring different optimizers. The current objective function adopts a quadratic form. Therefore, Conjugate Gradient is a preferable optimizer.

4. Registration Method (itkMeshToMeshRegistrationMethod)

	This class is templated over the pointset-to-pointset registration method. Users will create an object of this class to perform Thin Shell Demons. See the test example for usage.


License
=======

This software is distributed under the Apache License. Please see
LICENSE for details.


Author
======


Qingyu Zhao
