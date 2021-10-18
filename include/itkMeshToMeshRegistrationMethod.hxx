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
#ifndef itkMeshToMeshRegistrationMethod_hxx
#define itkMeshToMeshRegistrationMethod_hxx

#include "itkMeshToMeshRegistrationMethod.h"
namespace itk
{

template <typename TFixedMesh, typename TMovingMesh>
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::MeshToMeshRegistrationMethod()
{
  this->SetNumberOfRequiredOutputs(1);

  m_InitialTransformParameters = ParametersType(FixedMeshType::PointDimension);
  m_LastTransformParameters = ParametersType(FixedMeshType::PointDimension);

  m_InitialTransformParameters.Fill(0);
  m_LastTransformParameters.Fill(0);

  TransformOutputPointer transformDecorator =
    itkDynamicCastInDebugMode<TransformOutputType *>(this->MakeOutput(0).GetPointer());

  this->ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());
}

template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::SetInitialTransformParameters(const ParametersType & param)
{
  m_InitialTransformParameters = param;
  this->Modified();
}

template <typename TFixedMesh, typename TMovingMesh>
void
  MeshToMeshRegistrationMethod< TFixedMesh, TMovingMesh >
  ::Initialize()
{
  if (!m_FixedMesh)
  {
    itkExceptionMacro(<< "FixedMesh is not present");
  }

  if (!m_MovingMesh)
  {
    itkExceptionMacro(<< "MovingMesh is not present");
  }

  if (!m_Metric)
  {
    itkExceptionMacro(<< "Metric is not present");
  }

  if (!m_Optimizer)
  {
    itkExceptionMacro(<< "Optimizer is not present");
  }

  if (!m_Transform)
  {
    itkExceptionMacro(<< "Transform is not present");
  }

  // Set up the metric
  m_Metric->SetMovingMesh(m_MovingMesh);
  m_Metric->SetFixedMesh(m_FixedMesh);
  m_Metric->SetTransform(m_Transform);
  m_Metric->Initialize();

  // Set up the optimizer
  m_Optimizer->SetCostFunction(m_Metric);

  // Validate initial transform parameters
  if (m_InitialTransformParameters.Size() != m_Transform->GetNumberOfParameters())
  {
    itkExceptionMacro(<< "Size mismatch between initial parameter and transform");
  }

  m_Optimizer->SetInitialPosition(m_InitialTransformParameters);

  // Connect the transform to the Decorator
  TransformOutputType * transformOutput = static_cast<TransformOutputType *>(this->ProcessObject::GetOutput(0));

  transformOutput->Set(m_Transform.GetPointer());
}

template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::GenerateData()
{
  // Initialize the interconnects between components
  try
  {
    this->Initialize();
  }
  catch (ExceptionObject & err)
  {
    m_LastTransformParameters = ParametersType(1);
    m_LastTransformParameters.Fill(0.0f);

    // Pass the  exception to the caller
    throw err;
  }

  // Do the optimization
  try
  {
    m_Optimizer->StartOptimization();
  }
  catch (ExceptionObject & err)
  {
    // An error has occurred in the optimization.
    // Update the parameters
    m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

    // Pass the exception to the caller
    throw err;
  }

  // Get the results
  m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

  m_Transform->SetParameters(m_LastTransformParameters);
}

template <typename TFixedMesh, typename TMovingMesh>
const typename MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::TransformOutputType *
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::GetOutput() const
{
  return static_cast<const TransformOutputType *>(this->ProcessObject::GetOutput(0));
}

template <typename TFixedMesh, typename TMovingMesh>
DataObject::Pointer
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::MakeOutput(DataObjectPointerArraySizeType output)
{
  switch (output)
  {
    case 0:
      return TransformOutputType::New().GetPointer();
      break;
    default:
      itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
      return ITK_NULLPTR;
  }
}

template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::UpdateMovingMesh()
{
  std::cout << "Pranjal Inside UpdateMovingMesh" << std::endl;

  // update the moving mesh with the current transformation
  typedef typename MovingMeshType::PointsContainer OutputPointsContainer;
  typedef typename MovingMeshType::PointsContainer InputPointsContainer;

  const InputPointsContainer *                    inPoints = m_MovingMesh->GetPoints();
  typename MovingMeshType::PointsContainerPointer outPoints = m_MovingMesh->GetPoints();

  typename InputPointsContainer::ConstIterator inputPoint = inPoints->Begin();
  typename InputPointsContainer::ConstIterator inputEnd = inPoints->End();
  typename OutputPointsContainer::Iterator     outputPoint = outPoints->Begin();

  ParametersType m_VectorField = m_Transform->GetParameters();
  int            idx = 0;

  while (inputPoint != inputEnd)
  {
    const typename TMovingMesh::PointType & originalPoint = inputPoint.Value();
    typename TMovingMesh::PointType         displacedPoint;

    for (unsigned int i = 0; i < 3; i++)
    {
      displacedPoint[i] = originalPoint[i] + m_VectorField[idx * 3 + i];
    }
    outputPoint.Value() = displacedPoint;
    ++inputPoint;
    ++outputPoint;
    idx++;
  }
}
} // namespace itk
#endif
