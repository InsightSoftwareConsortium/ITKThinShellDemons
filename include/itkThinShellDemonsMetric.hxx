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
  fixedVTKMesh = nullptr;
  movingVTKMesh = nullptr;
  fixedCurvature = nullptr;
  updateConfidenceSigma = true;
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
  this->movingVTKMesh =
    itkMeshTovtkPolyData<MovingMeshType>::Convert(this->m_MovingMesh);
  this->fixedVTKMesh =
    itkMeshTovtkPolyData<FixedMeshType>::Convert(this->m_FixedMesh);

  this->ComputeNeighbors();
  // Preprocessing: compute the target position of each vertex in the fixed mesh
  // using Euclidean + Curvature distance
  vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
  curvaturesFilter->SetInputData(fixedVTKMesh);
  curvaturesFilter->SetCurvatureTypeToGaussian();
  curvaturesFilter->Update();
  this->fixedCurvature = curvaturesFilter->GetOutput();

  this->targetMap.resize(movingVTKMesh->GetNumberOfPoints());
  this->ComputeTargetPosition();

}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeNeighbors()
{
  this->neighborMap.resize(movingVTKMesh->GetNumberOfPoints());
  this->edgeLengthMap.resize(movingVTKMesh->GetNumberOfPoints());
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
          pointIdList->InsertUniqueId( pointIdListTmp->GetId(j) );
          }
        }
      }
    //Store edge lengths;
    edgeLengthMap[id].resize(pointIdList->GetNumberOfIds());
    const InputPointType &p = this->m_MovingMesh->GetPoint(id);
    for(int j=0; j < pointIdList->GetNumberOfIds(); j++)
      {
      const vtkIdType &nid = pointIdList->GetId(j);
      const InputPointType &pn = this->m_MovingMesh->GetPoint(nid);
      edgeLengthMap[id][j] = p.EuclideanDistanceTo(pn);
      //Avoid division by zero
      if( edgeLengthMap[id][j] < itk::NumericTraits<float>::epsilon())
        {
        edgeLengthMap[id][j] = itk::NumericTraits<float>::epsilon();
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
        (fixedC->GetTuple1(fixedIdentifier) - movingC->GetTuple1(identifier));
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

    targetMap[identifier] = targetPoint;

    ++pointItr;
    identifier++;
    }
  if( m_UseMaximalDistanceConfidenceSigma && updateConfidenceSigma)
    {
    this->m_ConfidenceSigma = sqrt( maxDistance  ) / 3;
    updateConfidenceSigma = false;
    }

}


template< typename TFixedMesh, typename TMovingMesh >
double
ThinShellDemonsMetric< TFixedMesh, TMovingMesh >
::ComputeConfidenceValueAndDerivative(const InputVectorType &v,
    InputVectorType &derivative) const
{
  double variance = m_ConfidenceSigma * m_ConfidenceSigma;
  double dist = v.GetSquaredNorm();
  double confidence = exp( -dist / (2*variance));
  if( this->m_UpdateFeatureMatchingAtEachIteration )
    {
    derivative = (-confidence / variance) * v;
    }
  return confidence;
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
  stretchEnergy = 0;
  bendEnergy = 0;
  stretch.Fill(0);
  bend.Fill(0);

  //Collect all neighbors
  const vtkSmartPointer<vtkIdList> &pointIdList = neighborMap[identifier];
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
    int nDegree = neighborMap[neighborIdx]->GetNumberOfIds();
    InputVectorType dx = (v - vn);

    stretchEnergy += dx.GetSquaredNorm();
    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // divided by the vertex degrees of current and neighbor vertex
    stretch += dx * (4 / (degree+nDegree));

    //Normalize bending by edge length
    dx /= edgeLengthMap[identifier][i];
    bEnergy += dx;
    bend += dx * (degree * 4 / (degree+nDegree));
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= degree;
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

  static unsigned int d = InputVectorType::Dimension;
  //derivative of data fidelity energy (squared distance to target position)
  MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingMesh->GetPoints()->End();
  if( derivative.GetSize() != movingMesh->GetNumberOfPoints() * d )
    {
    derivative = DerivativeType(movingMesh->GetNumberOfPoints() * d);
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
    for(int i=0; i<d; i++)
      {
      vec[i] = parameters[identifier*d+i];
      }
    typename Superclass::OutputPointType transformedPoint = inputPoint + vec;

    InputPointType &targetPoint = targetMap[identifier];
    InputVectorType distVec = targetPoint - transformedPoint;
    double confidence = 1;
    InputVectorType confidenceDerivative;
    if( this->m_UseConfidenceWeighting )
      {
      confidence = this->ComputeConfidenceValueAndDerivative(distVec, confidenceDerivative);
      }
    double dist = distVec.GetSquaredNorm();
    double sE = 0;
    double bE = 0;
    InputVectorType sD;
    InputVectorType bD;
    this->ComputeStretchAndBend(identifier, parameters, sE, bE, sD, bD);
    stretchEnergy += sE;
    bendEnergy += bE;

    double cost = confidence*dist + m_StretchWeight *sE + m_BendWeight * bE;
    InputVectorType dParam = -2.0*confidence*distVec + m_StretchWeight * sD + m_BendWeight * bD;
    if(this->m_UseConfidenceWeighting &&
       this->m_UpdateFeatureMatchingAtEachIteration)
      {
      dParam -= dist * confidenceDerivative;
      }
    value += cost;

    for(int i=0; i<d; i++)
      {
      derivative[identifier*d+i]   = dParam[i];
      }
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
