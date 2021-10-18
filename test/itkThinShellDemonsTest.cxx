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
#include "itkThinShellDemonsMetric.h"
#include "itkLBFGSOptimizer.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkMeshDisplacementTransform.h"
#include "itkConjugateGradientOptimizer.h"

int itkThinShellDemonsTest( int args, char **argv)
{
  const unsigned int Dimension = 3;
  typedef itk::Mesh<double, Dimension>       MeshType;
  typedef itk::VTKPolyDataReader< MeshType > ReaderType;

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

  /*
    Initialize Thin Shell Demons transformation
    this transformation type needs a mesh as a template
  */
  typedef itk::MeshDisplacementTransform<double, Dimension> TransformTestType;
  TransformTestType::Pointer transform = TransformTestType::New();
  transform->SetMeshTemplate(movingMesh);
  transform->Initialize();
  transform->SetIdentity();

  /*
    Initialize Thin Shell Demons metric
  */
  typedef itk::ThinShellDemonsMetric<MeshType> MetricType;
  MetricType::Pointer metric = MetricType::New();
  metric->SetStretchWeight(1);
  metric->SetBendWeight(5);
  metric->SetGeometricFeatureWeight(10);
  metric->UseConfidenceWeightingOn();
  metric->UseMaximalDistanceConfidenceSigmaOn();
  metric->UpdateFeatureMatchingAtEachIterationOn();
  metric->SetFixedMesh(fixedMesh);
  metric->SetMovingMesh(movingMesh);
  metric->SetTransform(transform);
  metric->Initialize();


  /*
    Initialize Thin Shell Demons optimizer
  */
  typedef itk::LBFGSOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
/*
  typedef itk::ConjugateGradientOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
*/

  typedef itk::MeshToMeshRegistrationMethod<MeshType, MeshType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform(transform);
  registration->SetInitialTransformParameters(transform->GetParameters());
  registration->SetFixedMesh(fixedMesh);
  registration->SetMovingMesh(movingMesh);


  std::cout << "Start Value= " << metric->GetValue(transform->GetParameters()) << std::endl;
  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & e )
  {
    std::cout << e << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Solution Value= " << metric->GetValue(transform->GetParameters()) << std::endl;


  /*
    output mesh
  */
  typedef itk::VTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  registration->UpdateMovingMesh();
  MeshType::ConstPointer registeredMesh = registration->GetMovingMesh();
  writer->SetInput(registeredMesh);
  writer->SetFileName( "registeredMesh.vtk" );
  writer->Write();
  return EXIT_SUCCESS;
}
