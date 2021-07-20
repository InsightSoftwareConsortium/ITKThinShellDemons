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
#ifndef itkMeshToMeshMetric_hxx
#define itkMeshToMeshMetric_hxx

#include "itkMeshToMeshMetric.h"

namespace itk
{
/** Constructor */
template< typename TFixedMesh, typename TMovingMesh >
MeshToMeshMetric< TFixedMesh, TMovingMesh >
::MeshToMeshMetric()
{
  m_FixedMesh = ITK_NULLPTR;    // has to be provided by the user.
  m_MovingMesh = ITK_NULLPTR; // has to be provided by the user.
  m_Transform = ITK_NULLPTR;    // has to be provided by the user.
}

/** Set the parameters that define a unique transform */
template< typename TFixedMesh, typename TMovingMesh >
void
MeshToMeshMetric< TFixedMesh, TMovingMesh >
::SetTransformParameters(const ParametersType & parameters) const
{
  if ( !m_Transform )
    {
    itkExceptionMacro(<< "Transform has not been assigned");
    }
  m_Transform->SetParameters(parameters);
}

/** Initialize the metric */
template< typename TFixedMesh, typename TMovingMesh >
void
MeshToMeshMetric< TFixedMesh, TMovingMesh >
::Initialize(void)
{
  if ( !m_Transform )
    {
    itkExceptionMacro(<< "Transform is not present");
    }

  if ( !m_MovingMesh )
    {
    itkExceptionMacro(<< "MovingMesh is not present");
    }

  if ( !m_FixedMesh )
    {
    itkExceptionMacro(<< "FixedMesh is not present");
    }

  // If the Mesh is provided by a source, update the source.
  if ( m_MovingMesh->GetSource() )
    {
    m_MovingMesh->GetSource()->Update();
    }

  // If the point set is provided by a source, update the source.
  if ( m_FixedMesh->GetSource() )
    {
    m_FixedMesh->GetSource()->Update();
    }
}

/** PrintSelf */
template< typename TFixedMesh, typename TMovingMesh >
void
MeshToMeshMetric< TFixedMesh, TMovingMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving Mesh: " << m_MovingMesh.GetPointer()  << std::endl;
  os << indent << "Fixed  Mesh: " << m_FixedMesh.GetPointer()   << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
}
} // end namespace itk

#endif
