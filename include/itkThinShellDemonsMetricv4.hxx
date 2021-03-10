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
#ifndef itkThinShellDemonsMetricv4_hxx
#define itkThinShellDemonsMetricv4_hxx

#include "itkThinShellDemonsMetricv4.h"

#include <vtkCurvatures.h>
#include <vtkPointData.h>

#include <math.h>

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ThinShellDemonsMetricv4()
{
  m_BendWeight = 1;
  m_StretchWeight = 1;
  m_GeometricFeatureWeight = 10;
  m_ConfidenceSigma = 3;
  m_UseMaximalDistanceConfidenceSigma = true;
  m_UseConfidenceWeighting = true;
  m_UpdateFeatureMatchingAtEachIteration = false;
  m_MovingTransformedFeaturePointsLocator = nullptr;
  fixedVTKMesh = nullptr;
  movingVTKMesh = nullptr;
  fixedCurvature = nullptr;
}

/** Initialize the metric */
template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::Initialize()
{

  if (!this->m_FixedPointSet)
  {
    itkExceptionMacro("Fixed point set is not present");
  }

  if (!this->m_MovingPointSet)
  {
    itkExceptionMacro("Moving point set is not present");
  }

  // If the PointSet is provided by a source, update the source.
  if (this->m_MovingPointSet->GetSource())
  {
    this->m_MovingPointSet->GetSource()->Update();
  }

  // If the point set is provided by a source, update the source.
  if (this->m_FixedPointSet->GetSource())
  {
    this->m_FixedPointSet->GetSource()->Update();
  }

  //generate a VTK copy of the same mesh
  itkMeshTovtkPolyData<double>::Pointer dataTransfer = itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_MovingPointSet);
  this->movingVTKMesh = dataTransfer->GetOutput();

  dataTransfer = itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_FixedPointSet);
  this->fixedVTKMesh = dataTransfer->GetOutput();

  this->ComputeNeighbors();

  Superclass::Initialize();

  //Compute confidence sigma
  if( this->m_UseMaximalDistanceConfidenceSigma )
    {
    this->ComputeMaximalDistanceSigma();
    }
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ComputeNeighbors()
{
  this->neighborMap.resize(fixedVTKMesh->GetNumberOfPoints());
  for(PointIdentifier id=0; id<fixedVTKMesh->GetNumberOfPoints(); id++)
    {
    //Collect all neighbors
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    fixedVTKMesh->GetPointCells(id, cellIdList);
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    for(PointIdentifier i = 0; i < cellIdList->GetNumberOfIds(); i++)
      {
      vtkSmartPointer<vtkIdList> pointIdListTmp = vtkSmartPointer<vtkIdList>::New();
      fixedVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdListTmp);
      for(PointIdentifier j=0; j < pointIdListTmp->GetNumberOfIds(); j++)
        {
        if(pointIdListTmp->GetId(j) != id)
          {
          pointIdList->InsertUniqueId(pointIdListTmp->GetId(j) );
          }
        }
      }
    neighborMap[id] = pointIdList;
    }
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
double
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ComputeConfidenceValueAndDerivative(const VectorType &v, VectorType &derivative) const
{
  double variance = m_ConfidenceSigma * m_ConfidenceSigma;
  double dist = v.GetSquaredNorm();
  double confidence = exp( -dist / (2*variance) );
  if( m_UpdateFeatureMatchingAtEachIteration )
    {
    derivative = (-confidence/variance) * v;
    }
  return confidence;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::VectorType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetMovingDirection(const PointIdentifier &identifier) const
{
  PointType p1 = this->m_FixedPointSet->GetPoint(identifier);
  PointType p2 = this->m_FixedTransformedPointSet->GetPoint(identifier);
  return p2 - p1;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ComputeStretchAndBend( const PointIdentifier &identifier,
                         double &stretchEnergy,
                         double &bendEnergy,
                         VectorType &stretch,
                         VectorType &bend) const
{
  stretchEnergy = 0;
  bendEnergy = 0;
  stretch.Fill(0);
  bend.Fill(0);

  //Collect all neighbors
  const vtkSmartPointer<vtkIdList> &pointIdList = this->neighborMap[identifier];
  int degree = pointIdList->GetNumberOfIds();
  VectorType v = this->GetMovingDirection(identifier);
  VectorType bEnergy;
  bEnergy.Fill(0);
  for(PointIdentifier i=0; i < pointIdList->GetNumberOfIds(); i++)
    {
    PointIdentifier neighborIdx = pointIdList->GetId(i);
    VectorType vn = this->GetMovingDirection(neighborIdx);
    VectorType dx = v - vn;
    stretchEnergy += dx.GetSquaredNorm();
    bEnergy += dx;

    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // divided by the vertex degrees of current and nieghbor vertex
    int nDegree =  this->neighborMap[neighborIdx]->GetNumberOfIds();
    stretch += dx * (4 / (degree+nDegree));
    bend += dx * (degree * 4 / (degree+nDegree));
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= degree;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::MeasureType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetLocalNeighborhoodValue(const PointIdentifier identifier,
                            const PointType &point,
                            const PixelType & pixel) const
{
  MeasureType value = 0;
  LocalDerivativeType derivative;
  this->GetLocalNeighborhoodValueAndDerivative(identifier, point, value, derivative, pixel);
  return value;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetLocalNeighborhoodValueAndDerivative(const PointIdentifier identifier,
                                         const PointType &point,
                                         MeasureType &value,
                                         LocalDerivativeType &derivative,
                                         const PixelType & pixel) const
{
  FeaturePointType fpoint = this->GetFeaturePoint(point, fixedCurvature->GetTuple1(identifier) );
  PointIdentifier mPointId = this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
  PointType closestPoint = this->m_MovingTransformedPointSet->GetPoint(mPointId);

  VectorType direction = closestPoint - point;
  double dist = direction.GetSquaredNorm();
  double confidence = 1;
  VectorType confidenceDerivative;
  if(this->m_UseConfidenceWeighting)
    {
    confidence = this->ComputeConfidenceValueAndDerivative(direction, confidenceDerivative);
    }
  double sE = 0;
  double bE = 0;
  VectorType sD;
  VectorType bD;
  this->ComputeStretchAndBend(identifier, sE, bE, sD, bD);
  VectorType dx = direction * confidence* 2 - m_StretchWeight*sD - bD * m_BendWeight;
  value = confidence * dist + m_StretchWeight * sE + m_BendWeight * bE;
  if(this->m_UseConfidenceWeighting && this->m_UpdateFeatureMatchingAtEachIteration)
    {
    dx += dist * confidenceDerivative;
    }
  derivative = dx;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::InitializePointSets() const
{
  Superclass::InitializePointSets();
  this->InitializeFeaturePointsLocators();
}


template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >::FeaturePointSetPointer
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GenerateFeaturePointSets(bool fixed) const
{

  vtkSmartPointer<vtkPolyData> vMesh;
  //Update meshes according to current transforms
  if(fixed)
    {
    vtkSmartPointer<vtkPoints> pts = fixedVTKMesh->GetPoints();
    for(PointIdentifier i=0; i<this->m_FixedTransformedPointSet->GetNumberOfPoints(); i++ )
      {
      pts->SetPoint(i, this->m_FixedTransformedPointSet->GetPoint(i).data());
      }
    vMesh = fixedVTKMesh;
    }
  else
    {
    vtkSmartPointer<vtkPoints> pts = movingVTKMesh->GetPoints();
    for(PointIdentifier i=0; i<this->m_MovingTransformedPointSet->GetNumberOfPoints(); i++ )
      {
      pts->SetPoint(i, this->m_MovingTransformedPointSet->GetPoint(i).data());
      }
    vMesh = movingVTKMesh;
    }

  vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
  curvaturesFilter->SetInputData(vMesh);
  curvaturesFilter->SetCurvatureTypeToGaussian();
  curvaturesFilter->Update();
  vtkSmartPointer<vtkPolyData> curvaturesOutput = curvaturesFilter->GetOutput();
  vtkSmartPointer<vtkDataArray> curvature = curvaturesOutput->GetPointData()->GetScalars();
  FeaturePointSetPointer features = FeaturePointSetType::New();
  if( fixed )
    {
    fixedCurvature = curvature;
    }
  else
    {
    auto fPoints = features->GetPoints();
    for(PointIdentifier i=0; i<vMesh->GetNumberOfPoints(); i++)
      {
      FeaturePointType point =
        this->GetFeaturePoint(vMesh->GetPoint(i), curvature->GetTuple1(i) );
      fPoints->InsertElement(i, point);
      }
    }
  return features;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ComputeMaximalDistanceSigma()
  const
{
  FeaturePointsContainerPointer mpoints =
    this->m_MovingTransformedFeaturePointsLocator->GetPoints();
  double maximalDistance = 0;
  for(PointIdentifier i = 0; i < fixedVTKMesh->GetNumberOfPoints(); i++)
    {
    FeaturePointType fpoint =
      this->GetFeaturePoint(fixedVTKMesh->GetPoint(i), fixedCurvature->GetTuple1(i));
    PointIdentifier id = this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
    FeaturePointType cpoint = mpoints->GetElement(id);
    double dist = cpoint.SquaredEuclideanDistanceTo(fpoint);
    if( dist > maximalDistance )
    {
      maximalDistance = dist;
      }
    }
  this->m_ConfidenceSigma = sqrt(maximalDistance)/3;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::InitializeFeaturePointsLocators()
  const
{
  //Update fixed curvature
  if(!fixedCurvature || this->m_UpdateFeatureMatchingAtEachIteration){
    this->GenerateFeaturePointSets(true);
  }

  //Update moving curvature feature locator
  if( !this->m_MovingTransformedFeaturePointsLocator
      || this->m_UpdateFeatureMatchingAtEachIteration )
  {
    if (!this->m_MovingTransformedPointSet)
    {
      itkExceptionMacro("The moving transformed point set does not exist.");
    }
    if (!this->m_MovingTransformedFeaturePointsLocator)
    {
      this->m_MovingTransformedFeaturePointsLocator = FeaturePointsLocatorType::New();
    }
    FeaturePointSetPointer features = this->GenerateFeaturePointSets(false);
    this->m_MovingTransformedFeaturePointsLocator->SetPoints(
        features->GetPoints());
    this->m_MovingTransformedFeaturePointsLocator->Initialize();
  }

  //Compute confidence sigma
  /*
  if( this->m_UpdateFeatureMatchingAtEachIteration &&
      this->m_UseMaximalDistanceConfidenceSigma )
    {
    this->ComputeMaximalDistanceSigma();
    }
    */
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::FeaturePointType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetFeaturePoint(const double *v, const double &c) const
{
  FeaturePointType fpoint;
  for(unsigned int i=0; i<PointType::Dimension; i++)
    {
    fpoint[i] = v[i];
    }
  fpoint[PointType::Dimension] = c * m_GeometricFeatureWeight;
  return fpoint;
}

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::FeaturePointType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetFeaturePoint(const PointType &v, const double &c) const
{
  FeaturePointType fpoint;
  for(unsigned int i=0; i<PointType::Dimension; i++)
    {
    fpoint[i] = v[i];
    }
  fpoint[PointType::Dimension] = c * m_GeometricFeatureWeight;
  return fpoint;
}


template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
