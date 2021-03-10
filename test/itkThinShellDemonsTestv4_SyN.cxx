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

#include "itkThinShellDemonsMetricv4.h"
#include "itkAffineTransform.h"
#include "itkDisplacementFieldTransformParametersAdaptor.h"
#include "itkSyNImageRegistrationMethod.h"
#include "itkEuclideanDistancePointSetToPointSetMetricv4.h"

/**
 * The implementation of the thin shell metricin the v4
 * regsitration framwork permits the option
 * to comnine the thin shell regularization with for
 * example the SyN diffeomoprhic transformations.
 */
int itkThinShellDemonsTestv4_SyN( int args, char **argv)
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

  int imageDiagonal = 100;
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

  using VectorType = itk::Vector<double, Dimension>;
  VectorType zeroVector(0.0);

  using DisplacementFieldType = itk::Image<VectorType, Dimension>;
  DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
  displacementField->CopyInformation(fixedImage);
  displacementField->SetRegions(fixedImage->GetBufferedRegion());
  displacementField->Allocate();
  displacementField->FillBuffer(zeroVector);

  DisplacementFieldType::Pointer inverseDisplacementField = DisplacementFieldType::New();
  inverseDisplacementField->CopyInformation(fixedImage);
  inverseDisplacementField->SetRegions(fixedImage->GetBufferedRegion());
  inverseDisplacementField->Allocate();
  inverseDisplacementField->FillBuffer(zeroVector);

  using TransformType = itk::DisplacementFieldTransform<double, Dimension>;

  using DisplacementFieldRegistrationType =
    itk::SyNImageRegistrationMethod<FixedImageType, MovingImageType,
                                    TransformType, FixedImageType, MeshType>;
  DisplacementFieldRegistrationType::Pointer registration =
    DisplacementFieldRegistrationType::New();

  using OutputTransformType = DisplacementFieldRegistrationType::OutputTransformType;
  OutputTransformType::Pointer outputTransform = OutputTransformType::New();
  outputTransform->SetDisplacementField(displacementField);
  outputTransform->SetInverseDisplacementField(inverseDisplacementField);
  registration->SetInitialTransform(outputTransform);
  registration->InPlaceOn();

  using AffineTransformType = itk::AffineTransform<double, MeshType::PointDimension>;
  AffineTransformType::Pointer transform = AffineTransformType::New();
  transform->SetIdentity();
/*
  using PointSetMetricType = itk::EuclideanDistancePointSetToPointSetMetricv4<MeshType>;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
*/
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

  double varianceForUpdateField = spacing*spacing*25;
  double varianceForTotalField = 0.0;
  registration->SetGaussianSmoothingVarianceForTheUpdateField(
      varianceForUpdateField);
  registration->SetGaussianSmoothingVarianceForTheTotalField(
      varianceForTotalField);
  registration->SetFixedPointSet(movingMesh);
  registration->SetMovingPointSet(fixedMesh);
  registration->SetMovingInitialTransform(transform);
  registration->SetMetric(metric);

  unsigned int numberOfLevels = 3;
  registration->SetNumberOfLevels(numberOfLevels);
  DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
  numberOfIterationsPerLevel.SetSize(numberOfLevels);
  numberOfIterationsPerLevel[0] = 5;
  numberOfIterationsPerLevel[1] = 10;
  numberOfIterationsPerLevel[2] = 50;
  registration->SetNumberOfIterationsPerLevel(numberOfIterationsPerLevel);

  DisplacementFieldRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize(numberOfLevels);
  shrinkFactorsPerLevel.Fill(1);

  DisplacementFieldRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize(numberOfLevels);
  smoothingSigmasPerLevel.Fill(0);

  using DisplacementFieldTransformAdaptorType = itk::DisplacementFieldTransformParametersAdaptor<OutputTransformType>;
  DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
  for (unsigned int level = 0; level < numberOfLevels; level++)
  {
    // We use the shrink image filter to calculate the fixed parameters of the virtual
    // domain at each level.  To speed up calculation and avoid unnecessary memory
    // usage, we could calculate these fixed parameters directly.

    using ShrinkFilterType = itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType>;
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetShrinkFactors(shrinkFactorsPerLevel[level]);
    shrinkFilter->SetInput(displacementField);
    shrinkFilter->Update();

    DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor = DisplacementFieldTransformAdaptorType::New();
    fieldTransformAdaptor->SetRequiredSpacing(shrinkFilter->GetOutput()->GetSpacing());
    fieldTransformAdaptor->SetRequiredSize(shrinkFilter->GetOutput()->GetBufferedRegion().GetSize());
    fieldTransformAdaptor->SetRequiredDirection(shrinkFilter->GetOutput()->GetDirection());
    fieldTransformAdaptor->SetRequiredOrigin(shrinkFilter->GetOutput()->GetOrigin());
    fieldTransformAdaptor->SetTransform(outputTransform);

    adaptors.push_back(fieldTransformAdaptor);
  }
  registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
  registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);


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

  OutputTransformType::Pointer tx = registration->GetModifiableTransform();
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
  writer->SetFileName( "synMovingMesh.vtk" );
  writer->Write();
  return EXIT_SUCCESS;
}
