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
#include <cstdlib>

//#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"

#include "itkCommand.h"
#include "itkThinShellDemonsMetricv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkLBFGS2Optimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkAffineTransform.h"

#include "itkDiscreteGaussianCurvatureQuadEdgeMeshFilter.h"

// Pranjal added these
#include "itkMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"

#include "itkTriangleCell.h"
#include <itkSimplexMesh.h>
#include <itkTriangleMeshToSimplexMeshFilter.h>
#include <itkSimplexMeshAdaptTopologyFilter.h>
#include <itkRegularSphereMeshSource.h>

template <typename TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate(){};

public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    if (typeid(event) != typeid(itk::IterationEvent))
    {
      return;
    }
    const auto * optimizer = dynamic_cast<const TFilter *>(object);

    if (!optimizer)
    {
      itkGenericExceptionMacro("Error dynamic_cast failed");
    }
    std::cout << "It: " << optimizer->GetCurrentIteration();
    std::cout << " metric value: " << optimizer->GetCurrentMetricValue();
    std::cout << std::endl;
  }
};

int
itkThinShellDemonsTestv4_Affine(int args, char ** argv)
{
  const unsigned int Dimension = 3;
  using PointType = float;
  using CoordType = float;

  using TriangleMeshTraits = itk::DefaultDynamicMeshTraits<float, 3, 3, float, float>;
  using TriangleMeshTraitsStatic = itk::DefaultStaticMeshTraits<float, 3, 3, float, float>;

  using SimplexMeshTraits = itk::DefaultDynamicMeshTraits<float, 3, 3, float, float>;
  using SimplexMeshTraitsStatic = itk::DefaultStaticMeshTraits<float, 3, 3, float, float>;

  // Declare the type of the input and output mesh
  using MeshType = itk::Mesh<float, 3>;
  using PointsContainerPointer = MeshType::PointsContainerPointer;
  
  using ReaderType = itk::MeshFileReader<MeshType>;
  using WriterType = itk::MeshFileWriter<MeshType>;

  unsigned int numberOfIterations = 100;

  /*
  Initialize fixed mesh polydata reader
  */
  ReaderType::Pointer           fixedPolyDataReader = ReaderType::New();
  //typedef ReaderType::PointType PointType;
  fixedPolyDataReader->SetFileName(argv[1]);
  //fixedPolyDataReader->SetFileName("/home/pranjal.sahu/ITKThinShellDemons/test/Baseline/fixedMesh.vtk");
  try
  {
    fixedPolyDataReader->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Error during Fixed Mesh Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  MeshType::Pointer fixedMesh = fixedPolyDataReader->GetOutput();

  /*
  Initialize moving mesh polydata reader
  */
  ReaderType::Pointer           movingPolyDataReader = ReaderType::New();
  //typedef ReaderType::PointType PointType;
  movingPolyDataReader->SetFileName(argv[2]);
  // movingPolyDataReader->SetFileName("/home/pranjal.sahu/ITKThinShellDemons/test/Baseline/movingMesh.vtk");
  //movingPolyDataReader->SetFileName("/home/pranjal.sahu/decimate_0.ICP_result.vtk");

  try{
    movingPolyDataReader->Update();
  }
  catch (itk::ExceptionObject & excp){
    std::cerr << "Error during Moving Mesh Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  MeshType::Pointer movingMesh = movingPolyDataReader->GetOutput();


  std::cout << "Fixed mesh points count " << fixedMesh->GetNumberOfPoints() << std::endl;
  std::cout << "Moving mesh points count " << movingMesh->GetNumberOfPoints() << std::endl;

  // using TConvert = itk::TriangleMeshToSimplexMeshFilter<TriangleMeshType, SimplexMeshType>;

  // // Ensure that all cells of the mesh are triangles.
  // for (TriangleMeshType::CellsContainerIterator it = movingMesh->GetCells()->Begin();
  //      it != movingMesh->GetCells()->End();
  //      ++it)
  // {
  //   TriangleMeshType::CellAutoPointer cell;
  //   movingMesh->GetCell(it->Index(), cell);
  //   if (3 != cell->GetNumberOfPoints())
  //   {
  //     std::cerr << "ERROR: All cells must be trianglar." << std::endl;
  //     return EXIT_FAILURE;
  //   }
  // }
















  
  
  // declare the triangle to simplex mesh filter
  // using SimplexFilterType = itk::TriangleMeshToSimplexMeshFilter<TriangleMeshType, SimplexMeshType>;

  // SimplexFilterType::Pointer simplexFilter = SimplexFilterType::New();
  // simplexFilter->SetInput(movingMesh);
  // simplexFilter->Update();
  // std::cout << "before movingMesh Created : " << movingMesh->GetNumberOfPoints() << std::endl;
  // int i = 0;
  // for (unsigned int n = 0; n < simplexMesh->GetNumberOfPoints(); n++)
  // {
  //   SimplexMeshType::PointIdentifier id1 = n;
  //   SimplexMeshType::PixelType point_data = n;
  //   SimplexMeshType::PixelType point_data1;
    
  //   simplexMesh->SetPointData(id1, point_data);
  //   simplexMesh->GetPointData(id1, &point_data1);
    
  //   //std::cout << n << " " << simplexMesh->GetPoint(id1) << " : " << point_data1 << " " << simplexMesh->GetMeanCurvature(id1) << std::endl;
  //   //geometryMap->GetElement(id1)->ComputeGeometry();
  //   //std::cout << n << " " << geometryMap->GetElement(id1)->pos << std::endl;
  // }

  // SimplexMeshType::CellsContainerPointer cells = simplexMesh->GetCells();

  // for (unsigned int n = 0; n < simplexMesh->GetNumberOfCells(); n++)
  // {
  //   SimplexMeshType::CellIdentifier id1 = n;
  //   itk::Array<float> point_ids = cells->GetElement(id1)->GetPointIdsContainer();

  //   //std::cout << "Points are :" << std::endl;
  //   for (unsigned int k = 0; k < point_ids.GetSize(); k++){
  //     SimplexMeshType::PointIdentifier id2 = k;
  //     //std::cout << point_ids[k] << " " << simplexMesh->GetMeanCurvature(id2) << ", ";
  //   }
  //   //std::cout << "--------------------" << std::endl;
  //   //simplexMesh->GetCell(id1, &cell_data);
  //   //std::cout << n << " " << simplexMesh->GetPoint(id1) << " : " << point_data1 << " " << simplexMesh->GetMeanCurvature(id1) << std::endl;
  //   //geometryMap->GetElement(id1)->ComputeGeometry();
  //   //std::cout << n << " " << geometryMap->GetElement(id1)->pos << std::endl;
  // }



  // using FilterType = itk::SimplexMeshAdaptTopologyFilter<SimplexMeshType, SimplexMeshType>;
  // FilterType::Pointer filter = FilterType::New();
  // filter->SetInput(simplexMesh);
  // filter->Update();
  // //filter->Print(std::cout);

  std::cout << "[TEST DONE]" << std::endl;
  





  // /* Creating QuadEdgeMesh from Triangle Mesh by inserting points and cells one by one*/

  // using Traits = itk::QuadEdgeMeshExtendedTraits<CoordType, Dimension, 2, CoordType, CoordType, CoordType, bool, bool>;
  // using QEMeshType = itk::QuadEdgeMesh<CoordType, Dimension, Traits>;
  // using QEPointsContainerPointer = QEMeshType::PointsContainerPointer;

  // QEMeshType::Pointer qe_mesh = QEMeshType::New();

  // std::cout << "Pranjal Number of points in the QE Mesh Before " << qe_mesh->GetNumberOfPoints() << ::endl;
  // for (unsigned int n = 0; n < movingMesh->GetNumberOfPoints(); n++)
  // {
  //   MeshType::PointIdentifier point_id = n;
  //   QEMeshType::PointType point = movingMesh->GetPoint(point_id);

  //   QEMeshType::PointIdentifier id1 = n;
  //   qe_mesh->SetPoint(id1, point);
  // }
  // std::cout << "Pranjal Number of points in the QE Mesh After " << qe_mesh->GetNumberOfPoints() << ::endl;

  // std::cout << "Pranjal Number of cells in the QE Mesh Before " << qe_mesh->GetNumberOfCells() << ::endl;
  // for (unsigned int n = 0; n < movingMesh->GetNumberOfCells(); n++)
  // {
  //   MeshType::CellIdentifier cell_id = n;
  //   MeshType::CellAutoPointer tri_cell;
  //   movingMesh->GetCell(cell_id, tri_cell);

  //   using CellType = typename QEMeshType::CellType;
  //   using TriangleCellType = itk::TriangleCell<CellType>;
  //   using TriangleCellAutoPointer = typename TriangleCellType::SelfAutoPointer;


  //   /* Creating a QE Cell from the Triangle Cell and inserting it into the QEMesh*/
  //   auto * triangleCell = new TriangleCellType;
  //   QEMeshType::CellAutoPointer qe_cell;

  //   itk::Array<float> point_ids = tri_cell->GetPointIdsContainer();
  //   for (unsigned int k = 0; k < 3; ++k)
  //   {
  //     triangleCell->SetPointId(k, point_ids[k]);
  //   }

  //   QEMeshType::CellIdentifier qe_cell_id = n;
  //   qe_cell.TakeOwnership(triangleCell);
  //   qe_mesh->SetCell(qe_cell_id, qe_cell);

  // }
  // std::cout << "Pranjal Number of cells in the QE Mesh After " << qe_mesh->GetNumberOfCells() << ::endl;








  // Convert the triangle mesh to a simplex mesh.
  // TConvert::Pointer convert = TConvert::New();
  // convert->SetInput(movingMesh);
  // convert->Update();

  // SimplexMeshType::Pointer simplex_output = convert->GetOutput();
  // simplex_output->DisconnectPipeline();

  //std::cout << "Pranjal simplex output is " << simplex_output << std::endl;
  
  //using FilterType = itk::SimplexMeshAdaptTopologyFilter<TSimplex, TSimplex>;
  //FilterType::Pointer filter = FilterType::New();
  //filter->SetInput(simplex_output);
  //filter->Update();
  //filter->Print(std::cout);

  //TSimplex::Pointer simplex_mesh_adapt = filter->GetOutput();
  //std::cout << "Pranjal simplex output SimplexMeshAdaptTopologyFilter done " << simplex_mesh_adapt << std::endl;
  //std::cout << "Pranjal simplex output SimplexMeshAdaptTopologyFilter done " << simplex_output << std::endl;

  /* For calculating itkDiscreteGaussianCurvatureQuadEdgeMeshFilterTest */
  
  // qe_reader->SetFileName(argv[2]);
  // qe_reader->Update();
  // QEMeshType::Pointer qe_moving_mesh = qe_reader->GetOutput();

  // qe_reader->SetFileName(argv[1]);
  // qe_reader->Update();
  // QEMeshType::Pointer qe_fixed_mesh = qe_reader->GetOutput();

  // CurvatureFilterType::Pointer gaussian_curvature = CurvatureFilterType::New();
  // //gaussian_curvature->SetInput(qe_fixed_mesh);
  // gaussian_curvature->SetInput(qe_mesh);
  // gaussian_curvature->Update();
  // QEMeshType::Pointer output = gaussian_curvature->GetOutput();

  // std::cout << "Pranjal Gaussian curvature output number of points " << output->GetNumberOfPoints() << std::endl;
  // QEPointsContainerPointer qe_points = qe_moving_mesh->GetPoints();


  // QEWriterType::Pointer  PolyDataWriter = QEWriterType::New();
  // PolyDataWriter->SetFileName("./qe_curvature_mesh1.vtk");
  // PolyDataWriter->SetInput(output);
  // PolyDataWriter->Update();
  
  // PolyDataWriter->SetFileName("./qe_fixed_mesh1.vtk");
  // PolyDataWriter->SetInput(qe_fixed_mesh);
  // PolyDataWriter->Update();
  


  using PixelType = double;
  using FixedImageType = itk::Image<PixelType, Dimension>;
  using MovingImageType = itk::Image<PixelType, Dimension>;

  FixedImageType::SizeType      fixedImageSize;
  FixedImageType::PointType     fixedImageOrigin;
  FixedImageType::DirectionType fixedImageDirection;
  FixedImageType::SpacingType   fixedImageSpacing;

  using PointIdentifier = MeshType::PointIdentifier;
  using BoundingBoxType = itk::BoundingBox<PointIdentifier, Dimension>;
  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
  PointsContainerPointer points = movingMesh->GetPoints();
  boundingBox->SetPoints(points);
  boundingBox->ComputeBoundingBox();

  typename BoundingBoxType::PointType minBounds = boundingBox->GetMinimum();
  typename BoundingBoxType::PointType maxBounds = boundingBox->GetMaximum();

  int    imageDiagonal = 5;
  double spacing = sqrt(boundingBox->GetDiagonalLength2()) / imageDiagonal;
  auto   diff = maxBounds - minBounds;
  fixedImageSize[0] = ceil(1.2 * diff[0] / spacing);
  fixedImageSize[1] = ceil(1.2 * diff[1] / spacing);
  fixedImageSize[2] = ceil(1.2 * diff[2] / spacing);
  fixedImageOrigin[0] = minBounds[0] - 0.1 * diff[0];
  fixedImageOrigin[1] = minBounds[1] - 0.1 * diff[1];
  fixedImageOrigin[2] = minBounds[2] - 0.1 * diff[2];
  fixedImageDirection.SetIdentity();
  fixedImageSpacing.Fill(spacing);

  FixedImageType::Pointer fixedImage = FixedImageType::New();
  fixedImage->SetRegions(fixedImageSize);
  fixedImage->SetOrigin(fixedImageOrigin);
  fixedImage->SetDirection(fixedImageDirection);
  fixedImage->SetSpacing(fixedImageSpacing);
  fixedImage->Allocate();

  using TransformType = itk::AffineTransform<double, Dimension>;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  transform->SetCenter(minBounds + (maxBounds - minBounds) / 2);

  std::cout << "Before creating the ThinShellDemonsMetricv4 " << std::endl;




  
  using PointSetMetricType = itk::ThinShellDemonsMetricv4<MeshType, MeshType>;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetStretchWeight(1);
  metric->SetBendWeight(5);
  metric->SetGeometricFeatureWeight(10);
  metric->UseConfidenceWeightingOn();
  metric->UseMaximalDistanceConfidenceSigmaOn();
  metric->UpdateFeatureMatchingAtEachIterationOff();
  metric->SetMovingTransform(transform);
  // Reversed due to using points instead of an image
  // to keep semantics the same as in itkThinShellDemonsTest.cxx
  // For the ThinShellDemonsMetricv4 the fixed mesh is
  // regularized
  metric->SetFixedPointSet(movingMesh);
  metric->SetMovingPointSet(fixedMesh);
  metric->SetVirtualDomainFromImage(fixedImage);
  metric->Initialize();

  std::cout << "After creating the ThinShellDemonsMetricv4 " << std::endl;

  // Scales estimator
  using ScalesType = itk::RegistrationParameterScalesFromPhysicalShift<PointSetMetricType>;
  ScalesType::Pointer shiftScaleEstimator = ScalesType::New();
  shiftScaleEstimator->SetMetric(metric);
  // Needed with pointset metrics
  shiftScaleEstimator->SetVirtualDomainPointSet(metric->GetVirtualTransformedPointSet());
  

  // optimizer

  // Does currently not support scaling
  // but change requested in:
  // https://github.com/InsightSoftwareConsortium/ITK/pull/2372
  /*
  typedef itk::LBFGS2Optimizerv4 OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetScalesEstimator( shiftScaleEstimator );
*/
  typedef itk::ConjugateGradientLineSearchOptimizerv4 OptimizerType;
  OptimizerType::Pointer                              optimizer = OptimizerType::New();
  optimizer->SetNumberOfIterations(50);
  optimizer->SetScalesEstimator(shiftScaleEstimator);
  optimizer->SetMaximumStepSizeInPhysicalUnits(0.5);
  optimizer->SetMinimumConvergenceValue(0.0);
  optimizer->SetConvergenceWindowSize(10);

  using CommandType = CommandIterationUpdate<OptimizerType>;
  CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver(itk::IterationEvent(), observer);

  using AffineRegistrationType =
    itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, TransformType, FixedImageType, MeshType>;
  AffineRegistrationType::Pointer registration = AffineRegistrationType::New();
  registration->SetNumberOfLevels(1);
  registration->SetObjectName("registration");
  registration->SetFixedPointSet(movingMesh);
  registration->SetMovingPointSet(fixedMesh);
  registration->SetInitialTransform(transform);
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);

  std::cout << "Start Value= " << metric->GetValue() << std::endl;
  try{
    registration->Update();
  }
  catch (itk::ExceptionObject & e){
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
  }

  TransformType::Pointer tx = registration->GetModifiableTransform();
  metric->SetTransform(tx);
  std::cout << "Solution Value= " << metric->GetValue() << std::endl;


//   for (unsigned int n = 0; n < movingMesh->GetNumberOfPoints(); n++)
//   {
//     TriangleMeshType::PointType txMovingPoint = tx->TransformPoint(movingMesh->GetPoint(n));
//     movingMesh->SetPoint(n, txMovingPoint);
//   }

//   WriterType::Pointer writer = WriterType::New();
//   writer->SetInput(movingMesh);
//   writer->SetFileName("affineMovingMesh.vtk");
//   writer->Update();

  return EXIT_SUCCESS;
}
