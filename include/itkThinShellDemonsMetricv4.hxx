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
#include <vtkPolyData.h>
#include <vtkPointData.h>

#include <math.h>

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh >
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::ThinShellDemonsMetricv4()
{
  m_ConfidenceSigma = 1;
  m_BendWeight = 1;
  m_StretchWeight = 1;
  m_GeometricFeatureWeight = 100;
}

/** Initialize the metric */
template< typename TFixedMesh, typename TMovingMesh >
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::Initialize()
{
  Superclass::Initialize();
  if ( !this->m_MovingPointSet )
    {
    itkExceptionMacro(<< "MovingMesh is not present");
    }

  if ( !this->m_FixedPointSet )
    {
    itkExceptionMacro(<< "FixedMesh is not present");
    }

  //generate a VTK copy of the same mesh
  itkMeshTovtkPolyData<double>::Pointer dataTransfer =
    itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_MovingPointSet);
  this->movingVTKMesh = dataTransfer->GetOutput();

  //Needs to be in pointLocator
  dataTransfer = itkMeshTovtkPolyData<double>::New();
  dataTransfer->SetInput(this->m_FixedPointSet);
  this->fixedVTKMesh = dataTransfer->GetOutput();
  /*
  vtkSmartPointer<vtkCurvatures> curvaturesFilter =
    vtkSmartPointer<vtkCurvatures>::New();
  curvaturesFilter->SetInputData(fixedVTKMesh);
  curvaturesFilter->SetCurvatureTypeToGaussian();
  curvaturesFilter->Update();
  fixedCurvature = curvaturesFilter->GetOutput();
  */

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
          pointIdList->InsertUniqueId (pointIdListTmp->GetId(j) );
          }
        }
      }
    neighborMap[id] = pointIdList;
    }
}


template< typename TFixedMesh, typename TMovingMesh >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >::VectorType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh >
::GetMovingDirection(int identifier) const
{
  PointType p1 = this->m_FixedPointSet->GetPoint(identifier);
  PointType p2 = this->m_FixedTransformedPointSet->GetPoint(identifier);
  return p1 - p2;
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
  const vtkSmartPointer<vtkIdList> &pointIdList = neighborMap.ElementAt(identifier);
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
    stretch[0] += 4 * dx[0] / (degree+nDegree);
    stretch[1] += 4 * dx[1] / (degree+nDegree);
    stretch[2] += 4 * dx[2] / (degree+nDegree);
    bend[0] += 4 * dx[0] / (degree+nDegree);
    bend[1] += 4 * dx[1] / (degree+nDegree);
    bend[2] += 4 * dx[2] / (degree+nDegree);
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= pointIdList->GetNumberOfIds();
  bendEnergy /= pointIdList->GetNumberOfIds();

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

  //TODO: Need to update / write pointlocator to use geomtric features
  PointIdentifier mPointId = this->m_MovingTransformedPointsLocator->FindClosestPoint(point);
  closestPoint = this->m_MovingPointSet->GetPoint(mPointId);

  value = point.SquaredEuclideanDistanceTo(closestPoint);
  VectorType direction = closestPoint - point;

  double sE = 0;
  double bE = 0;
  VectorType sD;
  VectorType bD;
  this->ComputeStretchAndBend(point, sE, bE, sD, bD);
  derivative = direction + m_StretchWeight*sD + bD * m_BendWeight;
  value += m_StretchWeight * sE + m_BendWeight * bE;
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
