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

#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkCommand.h"
#include "itkThinShellDemonsMetricv4.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkAffineTransform.h"
#include "itkMeshDisplacementTransform.h"
#include <itkDisplacementFieldTransform.h>
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
    std::cout << "It: " << optimizer->GetCurrentIteration() << " metric value: " << optimizer->GetCurrentMetricValue();
    std::cout << std::endl;
    }
};

int itkThinShellDemonsTestv4_Displacement( int args, char **argv)
{
  std::cout << argv[1] << std::endl;
  std::cout << argv[2] << std::endl;
  std::cout << "Runnng ThinShellDemonsTest" << std::endl;
  const unsigned int Dimension = 3;
  typedef itk::Mesh<double, Dimension>         MeshType;
  using PointsContainerPointer = MeshType::PointsContainerPointer;
  typedef itk::VTKPolyDataReader< MeshType >   ReaderType;

  unsigned int numberOfIterations = 100;

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
  typename BoundingBoxType::PointType maxBounds = boundingBox->GetMinimum();

  fixedImageSize[0] = maxBounds[0]-minBounds[0]+200;
  fixedImageSize[1] = maxBounds[1]-minBounds[1]+200;
  fixedImageSize[2] = maxBounds[2]-minBounds[2]+200;
  fixedImageOrigin[0] = minBounds[0]-100;
  fixedImageOrigin[1] = minBounds[1]-100;
  fixedImageOrigin[2] = minBounds[2]-100;
  fixedImageDirection.SetIdentity();
  fixedImageSpacing.Fill( 1 );

  FixedImageType::Pointer fixedImage = FixedImageType::New();
  fixedImage->SetRegions( fixedImageSize );
  fixedImage->SetOrigin( fixedImageOrigin );
  fixedImage->SetDirection( fixedImageDirection );
  fixedImage->SetSpacing( fixedImageSpacing );
  fixedImage->Allocate();

  //typedef itk::MeshDisplacementTransform<double, Dimension> TransformType;
  //TransformType::Pointer transform = TransformType::New();
  //transform->SetMeshTemplate(movingMesh); // this transformation type needs a mesh as a template
  //transform->Initialize();

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

  using PointSetMetricType = itk::ThinShellDemonsMetricv4<MeshType, MeshType> ;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetStretchWeight(5);
  metric->SetBendWeight(5);
  metric->SetGeometricFeatureWeight(100);
  metric->SetMovingTransform( transform );
  //Reversed due to using points instead of an image
  //to keep semantics the same
  metric->SetFixedMesh( movingMesh );
  metric->SetMovingMesh( fixedMesh );
  metric->SetVirtualDomainFromImage( fixedImage );
  metric->Initialize();

  // optimizer
  /*
  using OptimizerType = itk::GradientDescentOptimizerv4;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetScalesEstimator( shiftScaleEstimator );
  optimizer->SetMaximumStepSizeInPhysicalUnits( 3 );
  optimizer->SetMinimumConvergenceValue( 0.0 );
  optimizer->SetConvergenceWindowSize( 10 );
  */
  typedef itk::ConjugateGradientLineSearchOptimizerv4 OptimizerType;
  //typedef itk::LBFGSOptimizer OptimizerType;
  //typedef itk::GradientDescentOptimizer     OptimizerType;
  //typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  //typedef itk::LBFGS2Optimizerv4 OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetNumberOfIterations( 10 );

  using CommandType = CommandIterationUpdate<OptimizerType>;
  CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  using RegistrationType = itk::ImageRegistrationMethodv4<FixedImageType,
        MovingImageType, TransformType, FixedImageType, MeshType>;
  RegistrationType::Pointer simple = RegistrationType::New();
  simple->SetNumberOfLevels(1);
  simple->SetObjectName("simple");
  simple->SetFixedPointSet(movingMesh);
  simple->SetMovingPointSet(fixedMesh);
  simple->SetInitialTransform(transform);
  simple->SetMetric(metric);
  simple->SetOptimizer(optimizer);

   try
    {
    simple->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }


  //AffineTransformType::InverseTransformBasePointer affineInverseTransform =
  TransformType::Pointer tx = simple->GetModifiableTransform();
  for (unsigned int n = 0; n < movingMesh->GetNumberOfPoints(); n++)
  {
    // compare the points in virtual domain
    //std::cout<< tx->GetInterpolator()->IsInsideBuffer(movingMesh->GetPoint(n)) << std::endl;
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
