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
#include "itkPointSet.h"
#include <math.h>

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ThinShellDemonsMetricv4()
{
  m_BendWeight = 1;
  m_StretchWeight = 1;
  m_GeometricFeatureWeight = 0;

  m_ConfidenceSigma = 3;
  m_UseMaximalDistanceConfidenceSigma = true;
  m_UseConfidenceWeighting = true;
  m_UpdateFeatureMatchingAtEachIteration = false;
  m_MovingTransformedFeaturePointsLocator = nullptr;
  
  fixedITKMesh = nullptr;
  movingITKMesh = nullptr;
  fixedCurvature = nullptr;
}

/* Set the points and cells for the mesh */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::FillPointAndCell(PointSetPointer &pointset, MeshTypePointer &currentITKMesh)
{
  /* Insert points and cells in the currentITKMesh */
  for (unsigned int n = 0; n < pointset->GetNumberOfPoints(); n++)
  {
    PointType point = pointset->GetPoint(n);
    currentITKMesh->SetPoint(n, point);
  }
  
  for (unsigned int n = 0; n < pointset->GetNumberOfCells(); n++)
  {
    MeshCellAutoPointer tri_cell;
    pointset->GetCell(n, tri_cell);

    // Creating a Cell from the Triangle Cell and inserting it into the Mesh 
    auto * triangleCell = new MeshTriangleCellType;
    
    itk::Array<float> point_ids = tri_cell->GetPointIdsContainer();
    for (unsigned int k = 0; k < 3; ++k)
    {
      triangleCell->SetPointId(k, point_ids[k]);
    }

    MeshCellAutoPointer t_cell;
    t_cell.TakeOwnership(triangleCell);
    currentITKMesh->SetCell(n, t_cell);
  }

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

  this->fixedITKMesh = MeshType::New();
  this->movingITKMesh = MeshType::New();

  /* fill points and cells in fixed mesh */
  FillPointAndCell(this->m_FixedPointSet, this->fixedITKMesh);
  /* fill points and cells in moving mesh */
  FillPointAndCell(this->m_MovingPointSet, this->movingITKMesh);
  
  /* Build the Cell Links for the ITK Mesh for calculating the neighbours*/
  this->fixedITKMesh->BuildCellLinks();
  this->movingITKMesh->BuildCellLinks();
  
  this->curvature_filter = CurvatureFilterType::New();

  /* Compute Neighbors which will be used to calculate the stretch and bend energy*/
  this->ComputeNeighbors();

  Superclass::Initialize();

  //Compute confidence sigma
  if( this->m_UseMaximalDistanceConfidenceSigma )
    {
    this->ComputeMaximalDistanceSigma();
    }
}

/* Iterate over all the cells in which a point belongs and get the points present in those cells*/
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::ComputeNeighbors()
{
  this->neighborMap.resize(fixedITKMesh->GetNumberOfPoints());
  this->edgeLengthMap.resize(fixedITKMesh->GetNumberOfPoints());
  
  for (PointIdentifier id = 0; id < fixedITKMesh->GetNumberOfPoints(); id++)
  {
    /* For iterating over the cells for a given point */
    const std::set<PointIdentifier> link_set = this->fixedITKMesh->GetCellLinks()->ElementAt(id);
    std::set<PointIdentifier> pointIdSet;

    /* Iterate over the cells  and get the neighbouring points */
    for (auto elem : link_set){
        MeshCellAutoPointer tri_cell;
        this->fixedITKMesh->GetCell(elem, tri_cell);
        MeshCellPointIdConstIterator point_ids = tri_cell->GetPointIds();
        for (int ik = 0; ik < 3; ++ik){
          if (point_ids[ik] != id){
            pointIdSet.insert(point_ids[ik]);
          }
        }
    } 

    // Convert Set to Vector for  later use
    std::vector<PointIdentifier> pointIdList( pointIdSet.begin(), pointIdSet.end() );
    
    //Store edge lengths
    edgeLengthMap[id].resize(pointIdList.size());
    
    const PointType & p = this->m_FixedPointSet->GetPoint(id);
    for (unsigned long int j=0; j < pointIdList.size(); ++j)
    {
      PointIdentifier nid = pointIdList[j];
      const PointType &pn = this->m_FixedPointSet->GetPoint(nid);
      edgeLengthMap[id][j] =  p.EuclideanDistanceTo(pn);
      //Avoid division by zero
      if( edgeLengthMap[id][j] < itk::NumericTraits<float>::epsilon())
        {
        edgeLengthMap[id][j] = itk::NumericTraits<float>::epsilon();
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

  // Collect all neighbors
  std::vector<PointIdentifier> pointIdList = this->neighborMap[identifier];
  int degree = pointIdList.size();
  VectorType v = this->GetMovingDirection(identifier);
  VectorType bEnergy;
  bEnergy.Fill(0);

  for (long unsigned int i=0; i < pointIdList.size(); ++i)
  {
    PointIdentifier neighborIdx = pointIdList[i];
    int nDegree = this->neighborMap[neighborIdx].size();

    VectorType vn = this->GetMovingDirection(neighborIdx);
    VectorType dx = (v - vn);

    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // divided by the vertex degrees of current and nieghbor vertex
    stretch += dx * (4 / (degree+nDegree));
    stretchEnergy += dx.GetSquaredNorm();

    //Normalize bending by edge length
    dx /= edgeLengthMap[identifier][i];
    bEnergy += dx;
    bend += dx * (degree * 4 / (degree+nDegree));
    }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= degree;
}

/* Function definition of the original method definition in itkPointSetToPointSetMetric*/
/* Performs the computation in a multi-threaded manner */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::MeasureType
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetLocalNeighborhoodValueWithIndex(const PointIdentifier &identifier,
                            const PointType &point,
                            const PixelType & pixel) const
{
  MeasureType value = 0;
  LocalDerivativeType derivative;
  this->GetLocalNeighborhoodValueAndDerivativeWithIndex(identifier, point, value, derivative, pixel);
  return value;
}

/* Function definition of the original method definition in itkPointSetToPointSetMetric*/
/* Performs the computation in a multi-threaded manner */
/* This method is called inside the CalculateValueAndDerivative in itkPointSetToPointSetMetricWithIndexv4.hxx */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GetLocalNeighborhoodValueAndDerivativeWithIndex(const PointIdentifier &identifier,
                                         const PointType &point,
                                         MeasureType &value,
                                         LocalDerivativeType &derivative,
                                         const PixelType & pixel) const
{
  
  FeaturePointType fpoint = this->GetFeaturePoint(point, fixedCurvature->GetPointData()->ElementAt(identifier));
  
  PointIdentifier mPointId = this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
  PointType closestPoint = this->m_MovingTransformedPointSet->GetPoint(mPointId);

  VectorType direction = closestPoint - point;
  double dist = direction.GetSquaredNorm();
  double confidence = 1;
  VectorType confidenceDerivative{};
  if(this->m_UseConfidenceWeighting)
    {
    confidence = this->ComputeConfidenceValueAndDerivative(direction, confidenceDerivative);
    }
  double sE = 0;
  double bE = 0;
  VectorType sD;
  VectorType bD;
  this->ComputeStretchAndBend(identifier, sE, bE, sD, bD);
  VectorType dx = direction * confidence * 2 - m_StretchWeight*sD - bD * m_BendWeight;

  /* Refer to Equation 2 in the MIUA2015 paper */
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
  /* The call to Superclass initializes the m_MovingTransformedPointSet */
  Superclass::InitializePointSets();
  this->InitializeFeaturePointsLocators();
}


template< typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType >
typename ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >::FeaturePointSetPointer
ThinShellDemonsMetricv4< TFixedMesh, TMovingMesh, TInternalComputationValueType >
::GenerateFeaturePointSets(bool fixed) const
{
  MeshTypePointer currentMesh;

  //Update meshes according to current transforms
  if(fixed)
    {
    for(PointIdentifier i=0; i<this->m_FixedTransformedPointSet->GetNumberOfPoints(); i++ )
      {
      PointType data1 = this->m_FixedTransformedPointSet->GetPoint(i);
      fixedITKMesh->SetPoint(i, data1);
    }

    currentMesh = fixedITKMesh;
  }
  else
  {
    for (PointIdentifier i = 0; i < this->m_MovingTransformedPointSet->GetNumberOfPoints(); i++)
    {
      PointType data1 = this->m_MovingTransformedPointSet->GetPoint(i);
      movingITKMesh->SetPoint(i, data1);
    }

    currentMesh = movingITKMesh;
  }

  curvature_filter->SetTriangleMesh(currentMesh);
  curvature_filter->SetCurvatureTypeToGaussian();
  curvature_filter->Compute();
  auto curvature_output = curvature_filter->GetGaussCurvatureData();
  
  FeaturePointSetPointer        features = FeaturePointSetType::New();

  if( fixed )
    {
      /* Instantiate first time and re-use it for later iterations */
      if(!this->fixedCurvature){
        this->fixedCurvature = MeshType::New();
        PointDataContainerPointer pointData = PointDataContainer::New();
        pointData->Reserve(currentMesh->GetNumberOfPoints());
        this->fixedCurvature->SetPointData(pointData);
      }

      for (PointIdentifier i = 0; i < currentMesh->GetNumberOfPoints(); i++)
      {
        this->fixedCurvature->SetPointData(i, curvature_output->GetElement(i));
      }
    }
  else
    {
    auto fPoints = features->GetPoints();
    for (PointIdentifier i = 0; i < currentMesh->GetNumberOfPoints(); i++)
    {
      FeaturePointType point = this->GetFeaturePoint(currentMesh->GetPoint(i), curvature_output->GetElement(i));
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
  for (PointIdentifier i = 0; i < fixedITKMesh->GetNumberOfPoints(); i++)
  {
    FeaturePointType fpoint = this->GetFeaturePoint(fixedITKMesh->GetPoint(i), fixedCurvature->GetPointData()->ElementAt(i));
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

    // Only for the moving mesh, pass false to the GenerateFeaturePointSets
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

/* returns point with values [x, y, z, feature] */
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
