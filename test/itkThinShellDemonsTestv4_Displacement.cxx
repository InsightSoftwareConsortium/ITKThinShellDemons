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

#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"

#include "itkCommand.h"
#include "itkThinShellDemonsMetricv4.h"
#include "itkImageRegistrationMethodv4.h"
#include <itkDisplacementFieldTransform.h>
//#include "itkMeshDisplacementTransform.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkLBFGS2Optimizerv4.h"

template<typename TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *) caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }
    const auto * optimizer = dynamic_cast< const TFilter * >( object );

    if( !optimizer )
      {
      itkGenericExceptionMacro( "Error dynamic_cast failed" );
      }
    std::cout << "It: " << optimizer->GetCurrentIteration();
    std::cout << " metric value: " << optimizer->GetCurrentMetricValue();
    std::cout << std::endl;
    }
};

int itkThinShellDemonsTestv4_Displacement( int args, char **argv)
{
  const unsigned int Dimension = 3;
  typedef itk::Mesh<double, Dimension>         MeshType;
  using PointsContainerPointer = MeshType::PointsContainerPointer;
  typedef itk::VTKPolyDataReader< MeshType >   ReaderType;

  /*
  Initialize fixed mesh polydata reader
  */
  ReaderType::Pointer fixedPolyDataReader = ReaderType::New();
  typedef ReaderType::PointType PointType;
  fixedPolyDataReader->SetFileName(argv[1]);
  try
  {
    fixedPolyDataReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Error during Fixed Mesh Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  MeshType::Pointer fixedMesh = fixedPolyDataReader->GetOutput();

  /*
  Initialize moving mesh polydata reader
  */
  ReaderType::Pointer  movingPolyDataReader = ReaderType::New();
  typedef ReaderType::PointType PointType;
  movingPolyDataReader->SetFileName(argv[2]);
  try
  {
    movingPolyDataReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Error during Moving Mesh Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  MeshType::Pointer movingMesh = movingPolyDataReader->GetOutput();


  using PixelType = double;
  using FixedImageType = itk::Image<PixelType, Dimension>;
  using MovingImageType = itk::Image<PixelType, Dimension>;


  FixedImageType::SizeType fixedImageSize;
  FixedImageType::PointType fixedImageOrigin;
  FixedImageType::DirectionType fixedImageDirection;
  FixedImageType::SpacingType fixedImageSpacing;

  using PointIdentifier = MeshType::PointIdentifier;
  using BoundingBoxType = itk::BoundingBox<PointIdentifier, Dimension>;
  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
  PointsContainerPointer points = movingMesh->GetPoints();
  boundingBox->SetPoints(points);
  boundingBox->ComputeBoundingBox();
  typename BoundingBoxType::PointType minBounds = boundingBox->GetMinimum();
  typename BoundingBoxType::PointType maxBounds = boundingBox->GetMaximum();

  int imageDiagonal = 200;
  double spacing = sqrt(boundingBox->GetDiagonalLength2()) / imageDiagonal;
  auto diff = maxBounds - minBounds;
  fixedImageSize[0] = ceil( 1.2 * diff[0] / spacing );
  fixedImageSize[1] = ceil( 1.2 * diff[1] / spacing );
  fixedImageSize[2] = ceil( 1.2 * diff[2] / spacing );
  fixedImageOrigin[0] = minBounds[0] - 0.1*diff[0];
  fixedImageOrigin[1] = minBounds[1] - 0.1*diff[1];
  fixedImageOrigin[2] = minBounds[2] - 0.1*diff[2];
  fixedImageDirection.SetIdentity();
  fixedImageSpacing.Fill( spacing );

  FixedImageType::Pointer fixedImage = FixedImageType::New();
  fixedImage->SetRegions( fixedImageSize );
  fixedImage->SetOrigin( fixedImageOrigin );
  fixedImage->SetDirection( fixedImageDirection );
  fixedImage->SetSpacing( fixedImageSpacing );
  fixedImage->Allocate();


  using TransformType = itk::DisplacementFieldTransform<double, Dimension>;
  auto transform = TransformType::New();
  using  DisplacementFieldType = TransformType::DisplacementFieldType;
  DisplacementFieldType::Pointer field = DisplacementFieldType::New();
  field->SetRegions( fixedImageSize );
  field->SetOrigin( fixedImageOrigin );
  field->SetDirection( fixedImageDirection );
  field->SetSpacing( fixedImageSpacing );
  field->Allocate();
  transform->SetDisplacementField(field);

  using PointSetMetricType = itk::ThinShellDemonsMetricv4<MeshType> ;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetStretchWeight(1);
  metric->SetBendWeight(1);
  metric->SetGeometricFeatureWeight(10);
  metric->UseConfidenceWeightingOn();
  metric->UseMaximalDistanceConfidenceSigmaOn();
  metric->UpdateFeatureMatchingAtEachIterationOn();
  metric->SetMovingTransform( transform );
  //Reversed due to using points instead of an image
  //to keep semantics the same as in itkThinShellDemonsTest.cxx
  //For the ThinShellDemonsMetricv4 the fixed mesh is
  //regularized
  metric->SetFixedPointSet( movingMesh );
  metric->SetMovingPointSet( fixedMesh );
  metric->SetVirtualDomainFromImage( fixedImage );
  metric->Initialize();

  // Scales estimator
  using ScalesType = itk::RegistrationParameterScalesFromPhysicalShift< PointSetMetricType >;
  ScalesType::Pointer shiftScaleEstimator = ScalesType::New();
  shiftScaleEstimator->SetMetric( metric );
  // Needed with pointset metrics
  shiftScaleEstimator->SetVirtualDomainPointSet( metric->GetVirtualTransformedPointSet() );


  // optimizer

  // Does currently not support local transform
  // but change requeste in:
  // https://github.com/InsightSoftwareConsortium/ITK/pull/2372
  /*
  typedef itk::LBFGS2Optimizerv4 OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetScalesEstimator( shiftScaleEstimator );
  */

  typedef itk::ConjugateGradientLineSearchOptimizerv4 OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetNumberOfIterations( 50 );
  optimizer->SetScalesEstimator( shiftScaleEstimator );
  optimizer->SetMaximumStepSizeInPhysicalUnits( 0.01 );
  optimizer->SetMinimumConvergenceValue( 0.0 );
  optimizer->SetConvergenceWindowSize( 10 );


  using CommandType = CommandIterationUpdate<OptimizerType>;
  CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  using RegistrationType = itk::ImageRegistrationMethodv4<FixedImageType,
        MovingImageType, TransformType, FixedImageType, MeshType>;
  RegistrationType::Pointer registration = RegistrationType::New();
  registration->SetNumberOfLevels(1);
  registration->SetObjectName("registration");
  registration->SetFixedPointSet(movingMesh);
  registration->SetMovingPointSet(fixedMesh);
  registration->SetInitialTransform(transform);
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);

  std::cout << "Start Value= " << metric->GetValue() << std::endl;
  try
    {
    registration->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }

  TransformType::Pointer tx = registration->GetModifiableTransform();
  metric->SetTransform(tx);
  std::cout << "Solution Value= " << metric->GetValue() << std::endl;
  for (unsigned int n = 0; n < movingMesh->GetNumberOfPoints(); n++)
  {
    PointType txMovingPoint = tx->TransformPoint(movingMesh->GetPoint(n));
    movingMesh->SetPoint(n, txMovingPoint);
  }

  typedef itk::VTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(movingMesh);
  writer->SetFileName( "displacedMovingMesh.vtk" );
  writer->Write();
  return EXIT_SUCCESS;
}
