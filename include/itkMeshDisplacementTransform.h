/*=========================================================================
 *
 *  Copyright NumFOCUS
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

#ifndef itkMeshDisplacementTransform_h
#define itkMeshDisplacementTransform_h

#include "itkTransform.h"
#include "itkMesh.h"
#include "itkMacro.h"
#include "itkMatrix.h"

namespace itk
{

/** \class MeshDisplacementTransform
 *  \brief The class "MeshDisplacementTransformation" defines a finite
 *  dimensional vector space on mesh vertices. Its private member
 *  m_VectorField is a 1D parameter array in the form of
 *  [x_1,y_1,z_1,x_2,y_2,z_2,...], where the subscript denote the index
 *  of the vertex.
 *  A mesh has to be initially associated with a transformation object to
 *  serve as a template. The template essentially designates the number of
 *  vertices, so that m_VectorField can be initialized and allocated with
 *  a correct size (# of vertices * 3)
 *
 * \ingroup ThinShellDemons
 */
template<typename TParametersValueType=double,
           unsigned int NDimensions = 3>
class ITK_TEMPLATE_EXPORT MeshDisplacementTransform :
  public Transform<TParametersValueType, NDimensions, NDimensions>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(MeshDisplacementTransform);

  /** Standard class typedefs. */
  typedef MeshDisplacementTransform                                      Self;
  typedef Transform<TParametersValueType, NDimensions, NDimensions> Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** New macro for creation of through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshDisplacementTransform, Transform);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::FixedParametersType FixedParametersType;
  typedef typename Superclass::ParametersType      ParametersType;

  /** Standard Jacobian container. */
  typedef typename Superclass::JacobianType JacobianType;

  /** The number of parameters defininig this transform. */
  typedef typename Superclass::NumberOfParametersType NumberOfParametersType;

  /** Transform category type. */
  typedef typename Superclass::TransformCategoryType TransformCategoryType;

  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;

  typedef typename Superclass::InverseTransformBasePointer InverseTransformBasePointer;

  typedef itk::Mesh< TParametersValueType, NDimensions >       MeshType;
  typedef typename MeshType::ConstPointer                      MeshConstPointer;
  typedef typename MeshType::PointsContainer::ConstIterator    MeshPointIterator;
  typedef typename MeshType::PointDataContainer::ConstIterator MeshPointDataIterator;

  /** Set/Get the Mesh. */
  itkSetConstObjectMacro(MeshTemplate, MeshType);
  itkGetConstObjectMacro(MeshTemplate, MeshType);

  /** This method sets the parameters for the transform
   * value specified by the user. */
  virtual void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;

  /** Get the Transformation Parameters. */
  virtual const ParametersType & GetParameters() const ITK_OVERRIDE;


  /** Transform by an affine transformation.
   * This method applies the affine transform given by self to a
   * given point or vector, returning the transformed point or
   * vector. */
   OutputPointType TransformPoint(const InputPointType  & point) const ITK_OVERRIDE;

   using Superclass::TransformVector;
   OutputVectorType TransformVector(const InputVectorType & vector) const ITK_OVERRIDE;

   OutputVnlVectorType TransformVector(const InputVnlVectorType & vector) const ITK_OVERRIDE;

   using Superclass::TransformCovariantVector;
   OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType & vector) const ITK_OVERRIDE;

   OutputPointType     TransformNthPoint(const InputPointType  & point, int identifier) const;
  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const ITK_OVERRIDE;

  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual void ComputeJacobianWithRespectToParameters(const InputPointType & point, JacobianType & j) const ITK_OVERRIDE;

  /** Get the jacobian with respect to position, which simply is an identity
   *  jacobian because the transform is position-invariant.
   *  jac will be resized as needed, but it will be more efficient if
   *  it is already properly sized. */
  virtual void ComputeJacobianWithRespectToPosition(const InputPointType & x, JacobianType & jac) const ITK_OVERRIDE;

  /** Set the parameters to the IdentityTransform */
  void SetIdentity();

  /** Create a displace field for the mesh template */
  void Initialize();
  /** Return the number of parameters that completely define the Transfom  */
   virtual NumberOfParametersType GetNumberOfParameters() const ITK_OVERRIDE
   {
     return this->ParametersDimension;
   }

  /** Indicates that this transform is linear. That is, given two
   * points P and Q, and scalar coefficients a and b, then
   *
   * \f[ T( a*P + b*Q ) = a * T(P) + b * T(Q) \f]
   */
  virtual bool IsLinear() const ITK_OVERRIDE
  {
    return true;
  }

  /** Indicates the category transform.
   *  e.g. an affine transform, or a local one, e.g. a deformation field.
   */
  virtual TransformCategoryType GetTransformCategory() const ITK_OVERRIDE
  {
    return Self::Linear;
  }

  /** Set the fixed parameters and update internal transformation.
   * The MeshDisplacement Transform does not require fixed parameters,
   * therefore the implementation of this method is a null operation. */
  virtual void SetFixedParameters(const FixedParametersType &) ITK_OVERRIDE
  {
  }

  /** Get the Fixed Parameters. The MeshDisplacementTransform does not
   * require Fixed parameters, therefore this method returns an
   * parameters array of size zero. */
  virtual const FixedParametersType & GetFixedParameters() const ITK_OVERRIDE
  {
    this->m_FixedParameters.SetSize(0);
    return this->m_FixedParameters;
  }

protected:
  MeshDisplacementTransform();
  ~MeshDisplacementTransform();
  /** Print contents of an MeshDisplacementTransform. */
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:

  unsigned int SpaceDimension;
  unsigned int ParametersDimension;
  //MeshDeformationPointer m_MeshDeformation;
  MeshConstPointer m_MeshTemplate;
  JacobianType     m_IdentityJacobian;
  ParametersType m_VectorField;
};

}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshDisplacementTransform.hxx"
#endif

#endif /* itkMeshDisplacementTransform_h */
