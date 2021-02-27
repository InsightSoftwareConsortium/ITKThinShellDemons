/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkMeshToMeshMetric_h
#define itkMeshToMeshMetric_h

#include "itkImageBase.h"
#include "itkTransform.h"
#include "itkSingleValuedCostFunction.h"
#include "itkMacro.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

namespace itk
{
/** \class MeshToMeshMetric
 *
 * \brief This class is templated over the type of PointsetToPointsetMetric. This class serves as the basis for all kinds of mesh-to-mesh metrics (in some sense computing the similarity between two meshes). It expects a mesh-to-mesh transformation to be plugged in. This class computes an objective function value (also with its derivative w.r.t. the transformation parameters) that measures a registration metric between the fixed mesh and the moving mesh.
 *
 */

template< typename TFixedMesh,  typename TMovingMesh >
class ITK_TEMPLATE_EXPORT MeshToMeshMetric:public SingleValuedCostFunction
{
public:

  /** Standard class typedefs. */
  typedef MeshToMeshMetric   Self;
  typedef SingleValuedCostFunction Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Type used for representing point components  */
  typedef Superclass::ParametersValueType CoordinateRepresentationType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshToMeshMetric, SingleValuedCostFunction);

  /**  Type of the moving Mesh. */
  typedef TMovingMesh                           MovingMeshType;
  typedef typename TMovingMesh::PixelType       MovingMeshPixelType;
  typedef typename MovingMeshType::ConstPointer MovingMeshConstPointer;

  /**  Type of the fixed Mesh. */
  typedef TFixedMesh                           FixedMeshType;
  typedef typename FixedMeshType::ConstPointer FixedMeshConstPointer;

  /** Constants for the Mesh dimensions */
  itkStaticConstMacro(MovingMeshDimension, unsigned int,
                      TMovingMesh::PointDimension);
  itkStaticConstMacro(FixedMeshDimension, unsigned int,
                      TFixedMesh::PointDimension);

  typedef typename FixedMeshType::PointsContainer::ConstIterator     FixedPointIterator;
  typedef typename FixedMeshType::PointDataContainer::ConstIterator  FixedPointDataIterator;

  typedef typename MovingMeshType::PointsContainer::ConstIterator    MovingPointIterator;
  typedef typename MovingMeshType::PointDataContainer::ConstIterator MovingPointDataIterator;

  /**  Type of the Transform Base class */
  typedef Transform< CoordinateRepresentationType,
                     itkGetStaticConstMacro(MovingMeshDimension),
                     itkGetStaticConstMacro(FixedMeshDimension) > TransformType;

  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;

  /**  Type of the measure. */
  typedef Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef Superclass::ParametersType ParametersType;

  /** Get/Set the Fixed Mesh.  */
  itkSetConstObjectMacro(FixedMesh, FixedMeshType);
  itkGetConstObjectMacro(FixedMesh, FixedMeshType);

  /** Get/Set the Moving Mesh.  */
  itkSetConstObjectMacro(MovingMesh, MovingMeshType);
  itkGetConstObjectMacro(MovingMesh, MovingMeshType);

  /** Connect the Transform. */
  itkSetObjectMacro(Transform, TransformType);

  /** Get a pointer to the Transform.  */
  itkGetModifiableObjectMacro(Transform, TransformType);

  /** Set the parameters defining the Transform. */
  void SetTransformParameters(const ParametersType & parameters) const;

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE
  { return m_Transform->GetNumberOfParameters(); }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void)
  throw ( ExceptionObject );

protected:
  MeshToMeshMetric();
  virtual ~MeshToMeshMetric() {}
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  FixedMeshConstPointer m_FixedMesh;

  MovingMeshConstPointer m_MovingMesh;

  mutable TransformPointer m_Transform;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MeshToMeshMetric);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshToMeshMetric.hxx"
#endif

#endif
