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

template< typename TFixedMesh, typename TMovingMesh >
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
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
template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::Initialize()
{

  //generate a VTK copy of the same mesh
  itkMeshTovtkPolyData<double>::Pointer dataTransfer =itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_MovingPointSet);
  this->movingVTKMesh = dataTransfer->GetOutput();

  dataTransfer = itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_FixedPointSet);
  this->fixedVTKMesh = dataTransfer->GetOutput();

  Superclass::Initialize();

  this->ComputeNeighbors();
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::ComputeNeighbors()
{
  this->neighborMap.Initialize();
  for(int id=0; id<fixedVTKMesh->GetNumberOfPoints(); id++)
    {
    //Collect all neighbors
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    fixedVTKMesh->GetPointCells(id, cellIdList);
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    for(int i = 0; i < cellIdList->GetNumberOfIds(); i++)
      {
      vtkSmartPointer<vtkIdList> pointIdListTmp = vtkSmartPointer<vtkIdList>::New();
      fixedVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdListTmp);
      for(int j=0; j < pointIdListTmp->GetNumberOfIds(); j++)
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

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::ComputeConfidenceValueAndDerivative(const VectorType &v,
    double &confidence,
    VectorType &derivative) const
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
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >::VectorType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GetMovingDirection(int identifier) const
{
  PointType p1 = this->m_FixedPointSet->GetPoint(identifier);
  PointType p2 = this->m_FixedTransformedPointSet->GetPoint(identifier);
  return p2 - p1;
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::ComputeStretchAndBend( const PointType &point,
                         double &stretchEnergy,
                         double &bendEnergy,
                         VectorType &stretch,
                         VectorType &bend) const
{
  PointIdentifier identifier = this->m_FixedTransformedPointsLocator->FindClosestPoint(point);

  //enumerate all the neighboring vertices (edges) of a given vertex
  //stretching energy : measure the squared derivative along different edge directions
  //bending energy : measure the local laplacian around the local patch using
  //the given vertex and all neighboring vertices
  stretchEnergy = 0;
  bendEnergy = 0;
  stretch.Fill(0);
  bend.Fill(0);

  //Collect all neighbors
  const vtkSmartPointer<vtkIdList> pointIdList = neighborMap.ElementAt(identifier);
  int degree = pointIdList->GetNumberOfIds();
  VectorType v = GetMovingDirection(identifier);
  VectorType bEnergy;
  bEnergy.Fill(0);
  for(int i=0; i < pointIdList->GetNumberOfIds(); i++)
    {
    int neighborIdx = pointIdList->GetId(i);
    VectorType vn = this->GetMovingDirection(neighborIdx);
    VectorType dx = v - vn;
    stretchEnergy += dx.GetSquaredNorm();
    bEnergy += dx;

    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // divided by the vertex degrees of current and nieghbor vertex
    int nDegree =  neighborMap.ElementAt(neighborIdx)->GetNumberOfIds();
    stretch += dx * 4 / (degree+nDegree);
    //stretch[1] += 4 * dx[1] / (degree+nDegree);
    //stretch[2] += 4 * dx[2] / (degree+nDegree);
    bend += dx * 4 / (degree+nDegree);
    //bend[1] += 4 * dx[1] / (degree+nDegree);
    //bend[2] += 4 * dx[2] / (degree+nDegree);
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= degree;
  bendEnergy /= degree;
}


/**
 * Calculates the local metric value for a single point.
 */
template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >::MeasureType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GetLocalNeighborhoodValue(const PointType &point, const PixelType &pixel) const
{
  LocalDerivativeType derivative;
  MeasureType value;
  this->GetLocalNeighborhoodValueAndDerivative(point, value, derivative);
  return value;
}

/**
 * Calculates the local value and derivative for a single point.
 */
template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GetLocalNeighborhoodValueAndDerivative(const PointType &point,
                                         MeasureType &value,
                                         LocalDerivativeType &derivative,
                                         const PixelType & pixel) const
{
  PointType closestPoint;
  closestPoint.Fill(0.0);

  PointIdentifier identifier =
    this->m_FixedTransformedPointsLocator->FindClosestPoint(point);
  FeaturePointType fpoint =
    this->GetFeaturePoint(point, fixedCurvature->GetTuple1(identifier) );
  PointIdentifier mPointId =
    this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
  closestPoint = this->m_MovingTransformedPointSet->GetPoint(mPointId);

  value = point.SquaredEuclideanDistanceTo(closestPoint);
  VectorType direction = closestPoint - point;
  double confidence = 1;
  VectorType confidenceDerivative;
  if(this->m_UseConfidenceWeighting)
    {
    this->ComputeConfidenceValueAndDerivative(direction, confidence, confidenceDerivative);
    }

  double sE = 0;
  double bE = 0;
  VectorType sD;
  VectorType bD;
  this->ComputeStretchAndBend(point, sE, bE, sD, bD);
  VectorType dvalue = direction * 2 - m_StretchWeight*sD - bD * m_BendWeight;
  value += m_StretchWeight * sE + m_BendWeight * bE;
  if(this->m_UseConfidenceWeighting)
    {
    dvalue *= confidence;
    if(this->m_UpdateFeatureMatchingAtEachIteration)
      {
      dvalue += value * confidenceDerivative;
      }
    value *= confidence;
    }
  derivative =dvalue;
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::InitializePointSets() const
{
  Superclass::InitializePointSets();
  this->InitializeFeaturePointsLocators();
}


template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >::FeaturePointSetPointer
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GenerateFeaturePointSets(bool fixed) const
{

  vtkSmartPointer<vtkPolyData> vMesh;
  //Update meshes according to current transforms
  if(fixed)
    {
    vtkSmartPointer<vtkPoints> pts = fixedVTKMesh->GetPoints();
    for(int i=0; i<this->m_FixedTransformedPointSet->GetNumberOfPoints(); i++ )
      {
      pts->SetPoint(i, this->m_FixedTransformedPointSet->GetPoint(i).data());
      }
    vMesh = fixedVTKMesh;
    }
  else
    {
    vtkSmartPointer<vtkPoints> pts = movingVTKMesh->GetPoints();
    for(int i=0; i<this->m_MovingTransformedPointSet->GetNumberOfPoints(); i++ )
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
    for(int i=0; i<vMesh->GetNumberOfPoints(); i++)
      {
      FeaturePointType point = this->GetFeaturePoint(vMesh->GetPoint(i), curvature->GetTuple1(i) );
      fPoints->InsertElement(i, point);
      }
    }
  return features;
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::ComputeMaximalDistanceSigma()
  const
{
  FeaturePointsContainerPointer mpoints =
    this->m_MovingTransformedFeaturePointsLocator->GetPoints();
  double maximalDistance = 0;
  for(int i = 0; i < fixedVTKMesh->GetNumberOfPoints(); i++)
    {
    FeaturePointType fpoint =
      this->GetFeaturePoint(fixedVTKMesh->GetPoint(i), fixedCurvature->GetTuple1(i));
    int id = this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
    FeaturePointType cpoint = mpoints->GetElement(id);
    double dist = cpoint.SquaredEuclideanDistanceTo(fpoint);
    if( dist > maximalDistance )
    {
      maximalDistance = dist;
      }
    }
  m_ConfidenceSigma = sqrt(maximalDistance*2);
}

template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::InitializeFeaturePointsLocators()
  const
{


  if(!fixedCurvature || this->m_UpdateFeatureMatchingAtEachIteration){
    //Update fixed curvature
    this->GenerateFeaturePointSets(true);
  }

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
  if( this->m_UpdateFeatureMatchingAtEachIteration &&
      this->m_UseMaximalDistanceConfidenceSigma )
    {
    this->ComputeMaximalDistanceSigma();
    }
}

template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >::FeaturePointType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GetFeaturePoint(const double *v, const double &c) const
{
  FeaturePointType fpoint;
  for(int i=0; i< PointType::Dimension; i++)
    {
    fpoint[i] = v[i];
    }
  fpoint[PointType::Dimension+1] = c * m_GeometricFeatureWeight;
  return fpoint;
}

template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >::FeaturePointType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GetFeaturePoint(const PointType &v, const double &c) const
{
  FeaturePointType fpoint;
  for(int i=0; i< PointType::Dimension; i++)
    {
    fpoint[i] = v[i];
    }
  fpoint[PointType::Dimension+1] = c * m_GeometricFeatureWeight;
  return fpoint;
}


template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
