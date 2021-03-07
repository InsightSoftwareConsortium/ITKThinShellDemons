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

#include <vtkCurvatures.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

#include <math.h>

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh >
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ThinShellDemonsMetric()
{
  m_BendWeight = 1;
  m_StretchWeight = 1;
  m_GeometricFeatureWeight = 10;
  m_ConfidenceSigma = 3;
  m_UseMaximalDistanceConfidenceSigma = true;
  m_UseConfidenceWeighting = true;
  m_UpdateFeatureMatchingAtEachIteration = false;
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
  itkMeshTovtkPolyData<double>::Pointer dataTransfer =
    itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_MovingMesh);
  this->movingVTKMesh = dataTransfer->GetOutput();

  dataTransfer = itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_FixedMesh);
  this->fixedVTKMesh = dataTransfer->GetOutput();

  this->ComputeNeighbors();
  // Preprocessing: compute the target position of each vertex in the fixed mesh
  // using Euclidean + Curvature distance
  vtkSmartPointer<vtkCurvatures> curvaturesFilter =
  vtkSmartPointer<vtkCurvatures>::New();
  curvaturesFilter->SetInputData(fixedVTKMesh);
  curvaturesFilter->SetCurvatureTypeToGaussian();
  curvaturesFilter->Update();
  fixedCurvature = curvaturesFilter->GetOutput();

  this->targetMap.Initialize();
  this->ComputeTargetPosition();

}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeNeighbors()
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
    neighborMap[id] = pointIdList;
    }
}


template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeTargetPosition() const
{
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();
  MovingMeshConstPointer movingMesh = this->GetMovingMesh();
  MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingMesh->GetPoints()->End();

  double maxDistance = 0;
  vtkPoints *pts = movingVTKMesh->GetPoints();
  int identifier = 0;
  while ( pointItr != pointEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    OutputPointType txPoint = this->m_Transform->TransformPoint(inputPoint);
    pts->SetPoint(identifier, txPoint[0], txPoint[1], txPoint[2]);
    identifier++;
    ++pointItr;
    }

  vtkSmartPointer<vtkCurvatures> curvaturesFilter =
    vtkSmartPointer<vtkCurvatures>::New();
  curvaturesFilter->SetInputData(movingVTKMesh);
  curvaturesFilter->SetCurvatureTypeToGaussian();
  curvaturesFilter->Update();
  vtkSmartPointer<vtkPolyData> movingCurvature = curvaturesFilter->GetOutput();

  vtkSmartPointer<vtkDataArray> movingC = movingCurvature->GetPointData()->GetScalars();
  vtkSmartPointer<vtkDataArray> fixedC = fixedCurvature->GetPointData()->GetScalars();
  pointItr = movingMesh->GetPoints()->Begin();
  pointEnd = movingMesh->GetPoints()->End();
  identifier = 0;
  while ( pointItr != pointEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    OutputPointType transformedPoint = this->m_Transform->TransformPoint(inputPoint);
    InputPointType targetPoint;

    double minimumDistance = NumericTraits< double >::max();
    FixedPointIterator pointItr2 = fixedMesh->GetPoints()->Begin();
    FixedPointIterator pointEnd2 = fixedMesh->GetPoints()->End();
    int fixedIdentifier = 0;
    while ( pointItr2 != pointEnd2 )
      {
      double eDist = pointItr2.Value().SquaredEuclideanDistanceTo(transformedPoint);
      double fDist = m_GeometricFeatureWeight *
        fixedC->GetTuple1(fixedIdentifier) - movingC->GetTuple1(identifier);
      fDist *= fDist;
      double dist = eDist + fDist;
      if ( dist < minimumDistance )
        {
        targetPoint.CastFrom( pointItr2.Value() );
        minimumDistance = dist;
        }
      pointItr2++;
      fixedIdentifier++;
      }
    if( maxDistance < minimumDistance )
      {
      maxDistance = minimumDistance;
      }

    const_cast<TargetMapType*>(&targetMap)->SetElement(identifier, targetPoint);

    ++pointItr;
    identifier++;
    }
  if( m_UseMaximalDistanceConfidenceSigma )
    {
    m_ConfidenceSigma = sqrt( 2*maxDistance  );
    }
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeConfidenceValueAndDerivative(const InputVectorType &v,
    double &confidence,
    InputVectorType &derivative) const
{
  double variance = m_ConfidenceSigma * m_ConfidenceSigma;
  double dist = v.GetSquaredNorm();
  confidence = exp( -dist / variance);
  if( m_UpdateFeatureMatchingAtEachIteration )
    {
    derivative = -confidence * 2 / variance * v;
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
  //bending energy : measure the local Laplacian around the local patch using
  //the given vertex and all neighboring vertices
  stretchEnergy = 0;
  bendEnergy = 0;
  stretch.Fill(0);
  bend.Fill(0);

  //Collect all neighbors
  const vtkSmartPointer<vtkIdList> &pointIdList = neighborMap.ElementAt(identifier);
  int degree = pointIdList->GetNumberOfIds();
  InputVectorType v;
  static unsigned int d = InputVectorType::Dimension;
  for(int i=0; i<d; i++)
    {
    v[i] = parameters[identifier*d+i];
    }
  InputVectorType bEnergy;
  bEnergy.Fill(0);
  for(int i=0; i < pointIdList->GetNumberOfIds(); i++)
    {
    int neighborIdx = pointIdList->GetId(i);
    InputVectorType vn;
    for(int i=0; i<d; i++)
      {
      vn[0+i] = parameters[neighborIdx*d+i];
      }
    InputVectorType dx = v - vn;
    stretchEnergy += dx.GetSquaredNorm();
    bEnergy += dx;

    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // divided by the vertex degrees of current and neighbor vertex
    int nDegree =  neighborMap.ElementAt(neighborIdx)->GetNumberOfIds();
    stretch += dx * 4 / (degree+nDegree);
    bend += dx * 4 / (degree+nDegree);
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= degree;
  bendEnergy /= degree;
}

template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetric< TFixedMesh, TMovingMesh >::MeasureType
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::GetValue(const TransformParametersType & parameters) const
{
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
  if( m_UpdateFeatureMatchingAtEachIteration )
    {
    this->ComputeTargetPosition();
    }

  //derivative of data fidelity energy (squared distance to target position)
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
    double confidence = 1;
    InputVectorType confidenceDerivative;
    if( m_UseConfidenceWeighting )
      {
      this->ComputeConfidenceValueAndDerivative(distVec, confidence, confidenceDerivative);
      }
    double dist = distVec.GetSquaredNorm();
    double sE = 0;
    double bE = 0;
    InputVectorType sD;
    InputVectorType bD;
    this->ComputeStretchAndBend(identifier, parameters, sE, bE, sD, bD);
    stretchEnergy += sE;
    bendEnergy += bE;

    double cost = dist + m_StretchWeight *sE + m_BendWeight * bE;
    InputVectorType dParam = -2.0*distVec + m_StretchWeight * sD + m_BendWeight * bD;
    if(m_UseConfidenceWeighting)
      {
      dParam *= confidence;
      if(m_UpdateFeatureMatchingAtEachIteration)
        {
        dParam -= cost * confidenceDerivative;
        }
      cost *= confidence;
      }

    value += cost;

    derivative[identifier*3]   = dParam[0];
    derivative[identifier*3+1] = dParam[1];
    derivative[identifier*3+2] = dParam[2];

    ++pointItr;
    identifier++;
  }

  value /= movingMesh->GetNumberOfPoints();
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
