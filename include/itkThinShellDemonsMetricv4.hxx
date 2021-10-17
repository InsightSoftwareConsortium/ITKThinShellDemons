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
#include "itkMeshTovtkPolyData.h"

// Remove these
#include <vtkCurvatures.h>
#include <vtkPointData.h>
#include <typeinfo>
#include <math.h>

namespace itk
{

template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::ThinShellDemonsMetricv4()
{
  m_BendWeight = 1;
  m_StretchWeight = 1;
  m_GeometricFeatureWeight = 0;

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
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::Initialize()
{

  /* m_FixedPointSet and m_MovingPointSet are defined in the itkPointSetToPointSetMetricWithIndexv4.h */
  /* m_FixedTransformedPointSet and m_MovingTransformedPointSet are defined in the itkPointSetToPointSetMetricWithIndexv4.h */
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


  this->fixedITKMesh1 = MeshType::New();
  this->movingITKMesh1 = MeshType::New();

  /* TODO: Avoid repetition and make a function for this */  
  /* Insert points and cells in the fixedITKMesh1 */
  for (unsigned int n = 0; n < this->m_FixedPointSet->GetNumberOfPoints(); n++)
  {
    PointIdentifier point_id = n;
    PointType point = this->m_FixedPointSet->GetPoint(point_id);
    MeshPointIdentifier id1 = n;
    this->fixedITKMesh1->SetPoint(id1, point);
  }
  
  for (unsigned int n = 0; n < this->m_FixedPointSet->GetNumberOfCells(); n++)
  {
    MeshCellIdentifier cell_id = n;
    MeshCellAutoPointer tri_cell;
    this->m_FixedPointSet->GetCell(cell_id, tri_cell);

    // Creating a Cell from the Triangle Cell and inserting it into the Mesh 
    auto * triangleCell = new MeshTriangleCellType;
    
    itk::Array<float> point_ids = tri_cell->GetPointIdsContainer();
    for (unsigned int k = 0; k < 3; ++k)
    {
      triangleCell->SetPointId(k, point_ids[k]);
    }

    MeshCellAutoPointer t_cell;
    MeshCellIdentifier mesh_cell_id = n;
    t_cell.TakeOwnership(triangleCell);
    this->fixedITKMesh1->SetCell(mesh_cell_id, t_cell);
  }

  /* Insert points and cells in the movingITKMesh1 */
  for (unsigned int n = 0; n < this->m_MovingPointSet->GetNumberOfPoints(); n++)
  {
    PointIdentifier point_id = n;
    PointType point = this->m_MovingPointSet->GetPoint(point_id);
    MeshPointIdentifier id1 = n;
    this->movingITKMesh1->SetPoint(id1, point);
  }
  
  for (unsigned int n = 0; n < this->m_MovingPointSet->GetNumberOfCells(); n++)
  {
    MeshCellIdentifier cell_id = n;
    MeshCellAutoPointer tri_cell;
    this->m_MovingPointSet->GetCell(cell_id, tri_cell);

    // Creating a Cell from the Triangle Cell and inserting it into the Mesh 
    auto * triangleCell = new MeshTriangleCellType;
    
    itk::Array<float> point_ids = tri_cell->GetPointIdsContainer();
    for (unsigned int k = 0; k < 3; ++k)
    {
      triangleCell->SetPointId(k, point_ids[k]);
    }

    MeshCellAutoPointer t_cell;
    MeshCellIdentifier mesh_cell_id = n;
    t_cell.TakeOwnership(triangleCell);
    this->movingITKMesh1->SetCell(mesh_cell_id, t_cell);
  }

  std::cout << "fixedITKMesh1  " <<  this->fixedITKMesh1->GetNumberOfPoints() << " " << this->fixedITKMesh1->GetNumberOfCells() <<std::endl;
  std::cout << "movingITKMesh1  " <<  this->movingITKMesh1->GetNumberOfPoints() << " " <<  this->movingITKMesh1->GetNumberOfCells() <<std::endl;
  
  // generate a VTK copy of the same mesh
  this->movingVTKMesh = itkMeshTovtkPolyData<MovingPointSetType>::Convert(this->m_MovingPointSet);
  this->fixedVTKMesh = itkMeshTovtkPolyData<FixedPointSetType>::Convert(this->m_FixedPointSet);

  this->qeMeshCurvature = QEMeshType::New();
  this->gaussian_curvature_filter = CurvatureFilterType::New();

  /*
  QEWriterTypePointer  PolyDataWriter = QEWriterType::New();
  PolyDataWriter->SetFileName("./qe_curvature_mesh_3.vtk");
  PolyDataWriter->SetInput(output);
  PolyDataWriter->Update();
  */

  //std::cout << "Pranjal Obtained vtk mesh from the itkmesh using the filter " << std::endl;

  /* Compute Neighbors which will be used to calculate the stretch and bend energy*/
  this->ComputeNeighbors();

  Superclass::Initialize();

  std::cout << "Pranjal After the Initialize Method Call " << std::endl;

  // Compute confidence sigma
  if (this->m_UseMaximalDistanceConfidenceSigma)
  {
    this->ComputeMaximalDistanceSigma();
  }
}

/* Iterate over all the cells in which a point belongs and get the points present in those cells*/
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::ComputeNeighbors()
{
  this->neighborMap.resize(fixedVTKMesh->GetNumberOfPoints());
  this->edgeLengthMap.resize(fixedVTKMesh->GetNumberOfPoints());
  
  for (PointIdentifier id = 0; id < fixedVTKMesh->GetNumberOfPoints(); id++)
  {
    // Collect all neighbors
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    
    // Get all the neighboring cells
    fixedVTKMesh->GetPointCells(id, cellIdList);
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    for (PointIdentifier i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {
      vtkSmartPointer<vtkIdList> pointIdListTmp = vtkSmartPointer<vtkIdList>::New();

      // get all the points in the ith cell
      fixedVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdListTmp);
      for (PointIdentifier j = 0; j < pointIdListTmp->GetNumberOfIds(); j++)
      {
        // insert only if it is not the same point
        if (pointIdListTmp->GetId(j) != id)
        {
          pointIdList->InsertUniqueId(pointIdListTmp->GetId(j));
        }
      }
    }

    // Store edge lengths
    edgeLengthMap[id].resize(pointIdList->GetNumberOfIds());
    const PointType & p = this->m_FixedPointSet->GetPoint(id);
    for (PointIdentifier j = 0; j < pointIdList->GetNumberOfIds(); j++)
    {
      const vtkIdType & nid = pointIdList->GetId(j);
      const PointType & pn = this->m_FixedPointSet->GetPoint(nid);
      edgeLengthMap[id][j] = p.EuclideanDistanceTo(pn);
      // Avoid division by zero
      if (edgeLengthMap[id][j] < itk::NumericTraits<float>::epsilon())
      {
        edgeLengthMap[id][j] = itk::NumericTraits<float>::epsilon();
      }
    }

    neighborMap[id] = pointIdList;
  }

  std::cout << "Pranjal ComputeNeighbors done " << std::endl;
}


template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
double
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::ComputeConfidenceValueAndDerivative(
  const VectorType & v,
  VectorType &       derivative) const
{
  double variance = m_ConfidenceSigma * m_ConfidenceSigma;
  double dist = v.GetSquaredNorm();
  double confidence = exp(-dist / (2 * variance));
  if (m_UpdateFeatureMatchingAtEachIteration)
  {
    derivative = (-confidence / variance) * v;
  }
  return confidence;
}

/* Used for calculating the Stretch and Bend cost */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
typename ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::VectorType
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::GetMovingDirection(
  const PointIdentifier & identifier) const
{
  PointType p1 = this->m_FixedPointSet->GetPoint(identifier);
  PointType p2 = this->m_FixedTransformedPointSet->GetPoint(identifier);
  return p2 - p1;
}

template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::ComputeStretchAndBend(
  const PointIdentifier & identifier,
  double &                stretchEnergy,
  double &                bendEnergy,
  VectorType &            stretch,
  VectorType &            bend) const
{
  stretchEnergy = 0;
  bendEnergy = 0;
  stretch.Fill(0);
  bend.Fill(0);

  // Collect all neighbors
  const vtkSmartPointer<vtkIdList> & pointIdList = this->neighborMap[identifier];
  int                                degree = pointIdList->GetNumberOfIds();
  VectorType                         v = this->GetMovingDirection(identifier);
  VectorType                         bEnergy;
  bEnergy.Fill(0);

  for (PointIdentifier i = 0; i < pointIdList->GetNumberOfIds(); i++)
  {
    PointIdentifier neighborIdx = pointIdList->GetId(i);
    int             nDegree = this->neighborMap[neighborIdx]->GetNumberOfIds();

    VectorType vn = this->GetMovingDirection(neighborIdx);
    VectorType dx = (v - vn);

    // times 4 because edge appears two times in the energy function
    // and the derivative has another factor of 2 from the squared norm
    // divided by the vertex degrees of current and nieghbor vertex
    stretch += dx * (4 / (degree + nDegree));
    stretchEnergy += dx.GetSquaredNorm();

    // Normalize bending by edge length
    dx /= edgeLengthMap[identifier][i];
    bEnergy += dx;
    bend += dx * (degree * 4 / (degree + nDegree));
  }

  bendEnergy = bEnergy.GetSquaredNorm() / degree;
  stretchEnergy /= degree;
}

/* Function definition of the original method definition in itkPointSetToPointSetMetric*/
/* Performs the computation in a multi-threaded manner */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
typename ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::MeasureType
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::GetLocalNeighborhoodValueWithIndex(
  const PointIdentifier & identifier,
  const PointType &       point,
  const PixelType &       pixel) const
{
  MeasureType         value = 0;
  LocalDerivativeType derivative;
  this->GetLocalNeighborhoodValueAndDerivativeWithIndex(identifier, point, value, derivative, pixel);
  return value;
}

/* Function definition of the original method definition in itkPointSetToPointSetMetric*/
/* Performs the computation in a multi-threaded manner */
/* This method is called inside the CalculateValueAndDerivative in itkPointSetToPointSetMetricWithIndexv4.hxx */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::
  GetLocalNeighborhoodValueAndDerivativeWithIndex(const PointIdentifier & identifier,
                                                  const PointType &       point,
                                                  MeasureType &           value,
                                                  LocalDerivativeType &   derivative,
                                                  const PixelType &       pixel) const
{
  
  FeaturePointType fpoint = this->GetFeaturePoint(point, fixedCurvature->GetTuple1(identifier));

  // fpoint = [641.261, -15.2199, 471.961, -0.000875806]
  // point = [641.261, -15.2199, 471.961]
  // identifier = 300
  // fixedCurvature->GetTuple1(identifier) = 0.000665184
  
  PointIdentifier  mPointId = this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
  PointType        closestPoint = this->m_MovingTransformedPointSet->GetPoint(mPointId);

  VectorType direction = closestPoint - point;
  double     dist = direction.GetSquaredNorm();
  double     confidence = 1;
  VectorType confidenceDerivative;

  if (this->m_UseConfidenceWeighting)
  {
    confidence = this->ComputeConfidenceValueAndDerivative(direction, confidenceDerivative);
  }
  
  double     sE = 0;
  double     bE = 0;
  VectorType sD;
  VectorType bD;
  this->ComputeStretchAndBend(identifier, sE, bE, sD, bD);
  
  VectorType dx = direction * confidence * 2 - m_StretchWeight * sD - bD * m_BendWeight;

  /* Refer to Equation 2 in the MIUA2015 paper */
  value = confidence * dist + m_StretchWeight * sE + m_BendWeight * bE;
  
  if (this->m_UseConfidenceWeighting && this->m_UpdateFeatureMatchingAtEachIteration)
  {
    dx += dist * confidenceDerivative;
  }
  derivative = dx;
}

template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::InitializePointSets() const
{
  /* The call to Superclass initializes the m_MovingTransformedPointSet */
  Superclass::InitializePointSets();

  //std::cout << "Number of points are m_MovingTransformedPointSet  " <<  this->m_MovingTransformedPointSet->GetNumberOfPoints() << std::endl;

  //std::cout << "Pranjal InitializePointSets ThinShell " << std::endl;
  this->InitializeFeaturePointsLocators();
}

/* It is called in the beginning once to obtain the feature points using the curvature information */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
typename ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::FeaturePointSetPointer
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::GenerateFeaturePointSets(
  bool fixed) const
{
  std::cout << "Pranjal GenerateFeaturePointSets " << std::endl;
  //std::cout << "Pranjal " << TFixedMesh::PointType << "  " << TFixedMesh::PointDimension << std::endl;

  //using PolyDataPointsContainer = typename PolyDataType::PointsContainer;
  //using PolyDataPointIdentifier = typename PolyDataType::PointIdentifier;

  vtkSmartPointer<vtkPolyData> vMesh;
  MeshTypePointer VMesh1;

  // Update meshes according to current transforms
  if (fixed)
  {
    vtkSmartPointer<vtkPoints> pts = fixedVTKMesh->GetPoints();
    for (PointIdentifier i = 0; i < this->m_FixedTransformedPointSet->GetNumberOfPoints(); i++)
    {
      PointType data1 = pts->GetPoint(i);
      pts->SetPoint(i, this->m_FixedTransformedPointSet->GetPoint(i).data());
    }
   
    // Updating the itk Mesh (for calculating the neighbours later)
    for (PointIdentifier i = 0; i < this->m_FixedTransformedPointSet->GetNumberOfPoints(); i++)
    {
      PointType data1 = this->m_FixedTransformedPointSet->GetPoint(i);
      fixedITKMesh1->SetPoint(i, data1);
    }

    std::cout << "Updating 1 done " << std::endl;
    vMesh = fixedVTKMesh;
    VMesh1 = fixedITKMesh1;
  }
  else
  {
    vtkSmartPointer<vtkPoints> pts = movingVTKMesh->GetPoints();
    for (PointIdentifier i = 0; i < this->m_MovingTransformedPointSet->GetNumberOfPoints(); i++)
    {
      pts->SetPoint(i, this->m_MovingTransformedPointSet->GetPoint(i).data());
    }

    // Updating the itk Mesh (for calculating the neighbours later)
    for (PointIdentifier i = 0; i < this->m_MovingTransformedPointSet->GetNumberOfPoints(); i++)
    {
      PointType data1 = this->m_MovingTransformedPointSet->GetPoint(i);
      movingITKMesh1->SetPoint(i, data1);
    }

    std::cout << "Updating 2 done " << std::endl;
    vMesh = movingVTKMesh;
    VMesh1 = movingITKMesh1;
  }

  /* Create a QE Mesh to get the curvature */
  QEMeshTypePointer qeMesh = QEMeshType::New();
  for (unsigned int n = 0; n < VMesh1->GetNumberOfPoints(); n++)
  {
    MeshPointIdentifier point_id = n;
    QEMeshPointType point = VMesh1->GetPoint(point_id);
    QEMeshPointIdentifier id1 = n;
    qeMesh->SetPoint(id1, point);
  }

  /* Clear Re-make the cells in the QE Mesh */
  for (unsigned int n = 0; n < VMesh1->GetNumberOfCells(); n++)
  {
    MeshCellIdentifier cell_id = n;
    MeshCellAutoPointer tri_cell;
    VMesh1->GetCell(cell_id, tri_cell);

    /* Creating a QE Cell from the Triangle Cell and inserting it into the QEMesh */
    auto * triangleCell = new QETriangleCellType;
    QECellAutoPointer qe_cell;

    itk::Array<float> point_ids = tri_cell->GetPointIdsContainer();
    for (unsigned int k = 0; k < 3; ++k)
    {
      triangleCell->SetPointId(k, point_ids[k]);
    }

    QECellIdentifier qe_cell_id = n;
    qe_cell.TakeOwnership(triangleCell);
    qeMesh->SetCell(qe_cell_id, qe_cell);
  }
  
  std::cout << "Pranjal Number of points in the QE Mesh After " << qeMesh->GetNumberOfPoints() << ::endl;
  QEMeshTypePointer curvature_output;
  gaussian_curvature_filter->SetInput(qeMesh);
  gaussian_curvature_filter->Update();
  curvature_output = gaussian_curvature_filter->GetOutput();
  std::cout << "Pranjal QuadEdge Mesh Filter Curvature is done " << curvature_output->GetNumberOfPoints() << " " << curvature_output->GetPointData()->Size() << std::endl;
  
  /* Obtain curvature information at each point */
  vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
  curvaturesFilter->SetInputData(vMesh);
  curvaturesFilter->SetCurvatureTypeToGaussian();
  curvaturesFilter->Update();

  //std::cout << "Number of points in the mesh" << std::endl;
  //std::cout << vMesh->GetNumberOfPoints() << std::endl;

  vtkSmartPointer<vtkPolyData>  curvaturesOutput = curvaturesFilter->GetOutput();
  vtkSmartPointer<vtkDataArray> curvature = curvaturesOutput->GetPointData()->GetScalars();
  FeaturePointSetPointer        features = FeaturePointSetType::New();

  // For printing the results
  for (PointIdentifier i = 0; i < vMesh->GetNumberOfPoints(); i++)
  {
    std::cout << "curvature for " << i << " " << curvature->GetTuple1(i) << " "  << curvature_output->GetPointData()->ElementAt(i) << std::endl;
    //FeaturePointType point = this->GetFeaturePoint(vMesh->GetPoint(i), curvature->GetTuple1(i));
    //fPoints->InsertElement(i, point);
  }

  if (fixed){
    fixedCurvature = curvature;
  }
  else{
    auto fPoints = features->GetPoints();
    for (PointIdentifier i = 0; i < vMesh->GetNumberOfPoints(); i++)
    {
      // std::cout << curvature->GetTuple1(i) << " " << curvature_output->GetPointData()->ElementAt(i) << std::endl;
      FeaturePointType point = this->GetFeaturePoint(vMesh->GetPoint(i), curvature->GetTuple1(i));
      fPoints->InsertElement(i, point);
    }
  }
  return features;
}

template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::ComputeMaximalDistanceSigma() const
{
  FeaturePointsContainerPointer mpoints = this->m_MovingTransformedFeaturePointsLocator->GetPoints();
  double                        maximalDistance = 0;
  for (PointIdentifier i = 0; i < fixedVTKMesh->GetNumberOfPoints(); i++)
  {
    FeaturePointType fpoint = this->GetFeaturePoint(fixedVTKMesh->GetPoint(i), fixedCurvature->GetTuple1(i));
    PointIdentifier  id = this->m_MovingTransformedFeaturePointsLocator->FindClosestPoint(fpoint);
    FeaturePointType cpoint = mpoints->GetElement(id);
    double           dist = cpoint.SquaredEuclideanDistanceTo(fpoint);
    if (dist > maximalDistance)
    {
      maximalDistance = dist;
    }
  }
  this->m_ConfidenceSigma = sqrt(maximalDistance) / 3;
}

template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::InitializeFeaturePointsLocators() const
{

  // Update fixed curvature
  if (!fixedCurvature || this->m_UpdateFeatureMatchingAtEachIteration)
  {
    this->GenerateFeaturePointSets(true);
    std::cout << "Pranjal InitializeFeaturePointsLocators 1 " << std::endl;
  }

  // Update moving curvature feature locator
  if (!this->m_MovingTransformedFeaturePointsLocator || this->m_UpdateFeatureMatchingAtEachIteration)
  {
    if (!this->m_MovingTransformedPointSet)
    {
      itkExceptionMacro("The moving transformed point set does not exist.");
    }
    if (!this->m_MovingTransformedFeaturePointsLocator)
    {
      std::cout << "Pranjal InitializeFeaturePointsLocators 3 " << std::endl;
      this->m_MovingTransformedFeaturePointsLocator = FeaturePointsLocatorType::New();
    }

    std::cout << "Pranjal InitializeFeaturePointsLocators 2 " << std::endl;

    // Only for the moving mesh, pass false to the GenerateFeaturePointSets
    FeaturePointSetPointer features = this->GenerateFeaturePointSets(false);
    this->m_MovingTransformedFeaturePointsLocator->SetPoints(features->GetPoints());
    this->m_MovingTransformedFeaturePointsLocator->Initialize();
  }

  // Compute confidence sigma
  /*
  if( this->m_UpdateFeatureMatchingAtEachIteration &&
      this->m_UseMaximalDistanceConfidenceSigma )
    {
    this->ComputeMaximalDistanceSigma();
    }
    */
}

/* returns point with values [x, y, z, feature] */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
typename ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::FeaturePointType
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::GetFeaturePoint(const double * v,
                                                                                                 const double & c) const
{
  FeaturePointType fpoint;
  for (unsigned int i = 0; i < PointType::Dimension; i++)
  {
    fpoint[i] = v[i];
  }
  fpoint[PointType::Dimension] = c * m_GeometricFeatureWeight;
  return fpoint;
}

/* returns point with values [x, y, z, feature] */
template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
typename ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::FeaturePointType
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::GetFeaturePoint(const PointType & v,
                                                                                                 const double & c) const
{
  FeaturePointType fpoint;
  for (unsigned int i = 0; i < PointType::Dimension; i++)
  {
    fpoint[i] = v[i];
  }
  fpoint[PointType::Dimension] = c * m_GeometricFeatureWeight;
  return fpoint;
}


template <typename TFixedMesh, typename TMovingMesh, typename TInternalComputationValueType>
void
ThinShellDemonsMetricv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>::PrintSelf(std::ostream & os,
                                                                                           Indent         indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
