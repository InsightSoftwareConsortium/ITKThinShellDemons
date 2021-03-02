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
#ifndef itkThinShellDemonsMetric_hxx
#define itkThinShellDemonsMetric_hxx

#include "itkThinShellDemonsMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <math.h>

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh >
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ThinShellDemonsMetric()
{
  m_BendWeight = 1;
  m_StretchWeight = 1;
}

/** Initialize the metric */
template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::Initialize(void)
throw ( ExceptionObject )
{
  if ( !this->m_Transform )
    {
    itkExceptionMacro(<< "Transform is not present");
    }

  if ( !this->m_MovingMesh )
    {
    itkExceptionMacro(<< "MovingMesh is not present");
    }

  if ( !this->m_FixedMesh )
    {
    itkExceptionMacro(<< "FixedMesh is not present");
    }

  // If the Mesh is provided by a source, update the source.
  if ( this->m_MovingMesh->GetSource() )
    {
    this->m_MovingMesh->GetSource()->Update();
    }

  // If the point set is provided by a source, update the source.
  if ( this->m_FixedMesh->GetSource() )
    {
    this->m_FixedMesh->GetSource()->Update();
    }

  //generate a VTK copy of the same mesh
  itkMeshTovtkPolyData* dataTransfer = new itkMeshTovtkPolyData();
  dataTransfer->SetInput(this->m_MovingMesh);
  this->movingVTKMesh = dataTransfer->GetOutput();

  // Preprocessing: compute the target position of each vertex in the fixed mesh
  // using Euclidean + Curvature distance
  this->targetMap.Initialize();
  this->ComputeTargetPosition();

  this->ComputeNeighbors();

}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeNeighbors() const
{
  this->neighborMap.Initialize();
  for(int id=0; id<movingVTKMesh->GetNumberOfPoints(); id++)
    {
    //Collect all neighbors
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    movingVTKMesh->GetPointCells(id, cellIdList);
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    for(int i = 0; i < cellIdList->GetNumberOfIds(); i++)
      {
      vtkSmartPointer<vtkIdList> pointIdListTmp = vtkSmartPointer<vtkIdList>::New();
      movingVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdListTmp);
      for(int j=0; j < pointIdListTmp->GetNumberOfIds(); j++)
        {
        if(pointIdListTmp->GetId(j) != id)
          {
          pointIdList->InsertUniqueId (pointIdListTmp->GetId(j) );
          }
        }
      }
    neighborMap[i] = pointIdList;
    }
}


template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeTargetPosition() const
{
  std::cout << "Compute Target Position" << std::endl;
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

  if ( !fixedMesh )
  {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
  }

  MovingMeshConstPointer movingMesh = this->GetMovingMesh();

  if ( !movingMesh )
  {
    itkExceptionMacro(<< "Moving point set has not been assigned");
  }

  MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingMesh->GetPoints()->End();

    // In principal, this part should implement Euclidean + geometric feature similarity
    // Currently, this is simply a closest point search
  int identifier = 0;
  while ( pointItr != pointEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedPoint =
      this->m_Transform->TransformPoint(inputPoint);
    InputPointType targetPoint;

    double minimumDistance = NumericTraits< double >::max();
    FixedPointIterator pointItr2 = fixedMesh->GetPoints()->Begin();
    FixedPointIterator pointEnd2 = fixedMesh->GetPoints()->End();

    while ( pointItr2 != pointEnd2 )
      {
      double dist = pointItr2.Value().SquaredEuclideanDistanceTo(transformedPoint);


      if ( dist < minimumDistance )
        {
        targetPoint.CastFrom( pointItr2.Value() );
        minimumDistance = dist;
        }
      pointItr2++;
      }
    const_cast<TargetMapType*>(&targetMap)->SetElement(identifier, targetPoint);

    ++pointItr;
    identifier++;
    }
}


template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeStretchAndBend( int identifier,
                         const TransformParametersType &parameters,
                         double &stretchEnergy,
                         double &bendEnergy,
                         InputVectorType &stretch,
                         InputVectorType &bend) const
{
  //enumerate all the neighboring vertices (edges) of a given vertex
  //stretching energy : measure the squared derivative along different edge directions
  //bending energy : measure the local laplacian around the local patch using
  //the given vertex and all neighboring vertices
  stretchEnergy = 0;
  bendEnergy = 0;
  stretch.Fill(0);
  bend.Fill(0);

  //Collect all neighbors
  vtkSmartPointer<vtkIdList> pointIdList = neighborMap[identifier];
  int degree = pointIdList->GetNumberOfIds();
  InputVectorType v;
  v[0] = parameters[identifier*3];
  v[1] = parameters[identifier*3+1];
  v[2] = parameters[identifier*3+2];
  InputVectorType bEnergy;
  bEnergy.Fill(0);
  for(int i=0; i < pointIdList->GetNumberOfIds(); i++)
    {
    int neighborIdx = pointIdList->GetId(i);
    InputVectorType vn;
    vn[0] = parameters[neighborIdx*3];
    vn[1] = parameters[neighborIdx*3+1];
    vn[2] = parameters[neighborIdx*3+2];
    InputVectorType dx = v - vn;
    stretchEnergy += dx.GetSquaredNorm();
    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // TODO: should be corrected for different number of neighbors
    int nDegree =  neighborMap[neighborIdx]->GetNumberOfIds();
    stretch += 4 * dx/ (degree+nDegree);
    bEnergy += dx
    bend += 4 * dx / (degree+nDegree);
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;

  // times 4 because edge appears two times in the energy function
  // and the derivative has another factor of 2 from the squared norm
  // TODO: should be corrected for different number of neighbors

  stretchEnergy /= pointIdList->GetNumberOfIds();
  bendEnergy /= pointIdList->GetNumberOfIds();
}

template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetric< TFixedMesh, TMovingMesh >::MeasureType
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::GetValue(const TransformParametersType & parameters) const
{
  std::cout << "Get Value" << std::endl;
  static DerivativeType derivative = DerivativeType(0);
  MeasureType value = 0;
  GetValueAndDerivative(parameters, value, derivative);
  return value;
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::GetDerivative( const TransformParametersType &parameters,
                 DerivativeType &derivative ) const
{
  std::cout << "Get Derivative" << std::endl;
  MeasureType dummy = 0;
  this->GetValueAndDerivative(parameters, dummy, derivative);
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::GetValueAndDerivative(const TransformParametersType &parameters,
                        MeasureType &value, DerivativeType  &derivative) const
{
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();
  if ( !fixedMesh )
    {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
    }

  MovingMeshConstPointer movingMesh = this->GetMovingMesh();
  if ( !movingMesh )
    {
    itkExceptionMacro(<< "Moving point set has not been assigned");
    }

  this->SetTransformParameters(parameters);
  this->ComputeTargetPosition();

  // derivative of data fidelity energy (squared distance to target position)
  MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingMesh->GetPoints()->End();

  if( derivative.GetSize() != movingMesh->GetNumberOfPoints() * 3 )
    {
    derivative = DerivativeType(movingMesh->GetNumberOfPoints() * 3);
    }
  derivative.Fill(0);

  int identifier = 0;
  double functionValue = 0;
  double stretchEnergy = 0;
  double bendEnergy = 0;
  while ( pointItr != pointEnd )
  {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    InputVectorType vec;
    vec[0] = parameters[identifier*3];
    vec[1] = parameters[identifier*3+1];
    vec[2] = parameters[identifier*3+2];
    typename Superclass::OutputPointType transformedPoint = inputPoint + vec;

    InputPointType targetPoint = targetMap.ElementAt(identifier);
    InputVectorType distVec = targetPoint - transformedPoint;

    derivative[identifier*3] = -2 * distVec[0];
    derivative[identifier*3 + 1] = -2 * distVec[1];
    derivative[identifier*3 + 2] = -2 * distVec[2];

    functionValue += distVec.GetSquaredNorm();

    double sE = 0;
    double bE = 0;
    InputVectorType sD;
    InputVectorType bD;
    this->ComputeStretchAndBend(identifier, parameters, sE, bE, sD, bD);
    stretchEnergy += sE ;
    bendEnergy += bE;
    derivative[identifier*3] += m_StretchWeight * sD[0] + m_BendWeight * bD[0];
    derivative[identifier*3+1] += m_StretchWeight * sD[1] + m_BendWeight *bD[1];
    derivative[identifier*3+2] += m_StretchWeight * sD[2] + m_BendWeight * bD[2];

    ++pointItr;
    identifier++;
  }

  std::cout << "Dist:    " << functionValue << std::endl;
  std::cout << "Stretch: " << stretchEnergy << std::endl;
  std::cout << "Bend:    " << bendEnergy << std::endl;

  value = functionValue + m_StretchWeight * stretchEnergy + m_BendWeight * bendEnergy;
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
