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
#ifndef itkMeshDisplacementTransform_hxx
#define itkMeshDisplacementTransform_hxx

#include "itkMeshDisplacementTransform.h"
#include "itkMath.h"
#include "itkMath.h"

namespace itk
{

template<typename TParametersValueType, unsigned int NDimensions>
MeshDisplacementTransform<TParametersValueType, NDimensions>
	::MeshDisplacementTransform() : Superclass(0)
{
	m_MeshTemplate = ITK_NULLPTR;
	this->SpaceDimension = NDimensions;
	this->ParametersDimension = 0;
}


template<typename TParametersValueType, unsigned int NDimensions>
MeshDisplacementTransform<TParametersValueType, NDimensions>
::~MeshDisplacementTransform()
{
}


template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>
::SetParameters(const ParametersType & parameters)
{
	if( parameters.Size() != this->ParametersDimension )
	{
		itkExceptionMacro( << "Mismatch between parameters size "
			<< parameters.Size() << " and expected number of parameters "
			<< this->ParametersDimension );
	}

	if( &parameters != &( this->m_VectorField ) )
	{
		// Clean up this->m_InternalParametersBuffer because we will
		// use an externally supplied set of parameters as the buffer
		this->m_VectorField = parameters;
	}

	// Modified is always called since we just have a pointer to the
	// parameters and cannot know if the parameters have changed.
	this->Modified();
}

template<typename TParametersValueType, unsigned int NDimensions>
const typename MeshDisplacementTransform<TParametersValueType, NDimensions>::ParametersType &
MeshDisplacementTransform<TParametersValueType, NDimensions>
::GetParameters() const
{

	return this->m_VectorField;
}

template<typename TParametersValueType, unsigned int NDimensions>
void
	MeshDisplacementTransform<TParametersValueType, NDimensions>
	::SetIdentity()
{
	if (!m_MeshTemplate)
	{
		itkExceptionMacro(<< "Mesh template is not present");
	}
	if (ParametersDimension == 0)
	{
		itkExceptionMacro(<< "Mesh template has zero vertex");
	}

	m_VectorField.Fill(0);
}

template<typename TParametersValueType, unsigned int NDimensions>
void
	MeshDisplacementTransform<TParametersValueType, NDimensions>
	::Initialize()
{
	if (!m_MeshTemplate)
	{
		itkExceptionMacro(<< "FixedMesh is not present");
	}

	// the size of the parameters can only be determined after knowing the number of vertices
    // the template mesh should be available before this initialization step
	m_VectorField.SetSize(m_MeshTemplate->GetNumberOfPoints() * SpaceDimension);
	m_VectorField.Fill(0);
	this->ParametersDimension = m_VectorField.GetSize();

}

template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputPointType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformPoint(const InputPointType & point) const
{
	return point;
}

template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputPointType
	MeshDisplacementTransform<TParametersValueType, NDimensions>
	::TransformNthPoint(const InputPointType & point, int identifier) const
{
	InputVectorType vec;
	vec[0] = m_VectorField[identifier];

	return point + vec;
}

template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputVectorType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformVector(const InputVectorType & vect) const
{
  return vect;
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputVnlVectorType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformVector(const InputVnlVectorType & vect) const
{
  return vect;
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::OutputCovariantVectorType
MeshDisplacementTransform<TParametersValueType, NDimensions>
::TransformCovariantVector(const InputCovariantVectorType & vect) const
{
  return vect;
}


template<typename TParametersValueType, unsigned int NDimensions>
typename MeshDisplacementTransform<TParametersValueType, NDimensions>::InverseTransformBasePointer
MeshDisplacementTransform<TParametersValueType, NDimensions>
::GetInverseTransform() const
{
  Pointer inv = New();

  //return GetInverse(inv) ? inv.GetPointer() : ITK_NULLPTR;
  return inv;
}


template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>::ComputeJacobianWithRespectToParameters(
  const InputPointType &,
  JacobianType & jacobian) const
{
  // the Jacobian is constant for this transform, and it has already been
  // initialized in the constructor, so we just need to return it here.
  jacobian = this->m_IdentityJacobian;
}


template<typename TParametersValueType, unsigned int NDimensions>
void
MeshDisplacementTransform<TParametersValueType, NDimensions>
::ComputeJacobianWithRespectToPosition(const InputPointType &,
                                       JacobianType & jac) const
{
  jac.SetSize( NDimensions, NDimensions );
  jac.Fill(0.0);
  for( unsigned int dim = 0; dim < NDimensions; dim++ )
    {
    jac[dim][dim] = 1.0;
    }
}


} // end namespace itk

#endif
