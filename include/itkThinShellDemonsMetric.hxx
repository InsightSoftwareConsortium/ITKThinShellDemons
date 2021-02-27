/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef itkThinShellDemonsMetric_hxx
#define itkThinShellDemonsMetric_hxx

#include "itkThinShellDemonsMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::ThinShellDemonsMetric() :
  m_TargetPositionComputed( false )
{
	m_BendWeight = 1;
	m_StretchWeight = 1;
}
  /** Initialize the metric */
  template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
  void
	  ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
	  ::Initialize(void)
	  throw ( ExceptionObject )
  {
	  if ( !this->m_Transform )
	  {
		  itkExceptionMacro(<< "Transform is not present");
	  }

	  if ( !this->m_MovingMesh )
	  {
		  itkExceptionMacro(<< "MovingMesh is not present");
	  }

	  if ( !this->m_FixedMesh )
	  {
		  itkExceptionMacro(<< "FixedMesh is not present");
	  }

	  // If the Mesh is provided by a source, update the source.
	  if ( this->m_MovingMesh->GetSource() )
	  {
		  this->m_MovingMesh->GetSource()->Update();
	  }

	  // If the point set is provided by a source, update the source.
	  if ( this->m_FixedMesh->GetSource() )
	  {
		  this->m_FixedMesh->GetSource()->Update();
	  }

	  // Preprocessing: compute the target position of each vertex in the fixed mesh
      // using Euclidean + Curvature distance
	  ComputeTargetPosition();
  }

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
	ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
	::ComputeTargetPosition()
{
	FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

	if ( !fixedMesh )
	{
		itkExceptionMacro(<< "Fixed point set has not been assigned");
	}

	MovingMeshConstPointer movingMesh = this->GetMovingMesh();

	if ( !movingMesh )
	{
		itkExceptionMacro(<< "Moving point set has not been assigned");
	}

	this->targetMap.Initialize();

	MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
	MovingPointIterator pointEnd = movingMesh->GetPoints()->End();

    // In principal, this part should implement Euclidean + geometric feature similarity
    // Currently, this is simply a closest point search
	int identifier = 0;
	while ( pointItr != pointEnd )
	{
		InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );
		typename Superclass::OutputPointType transformedPoint =
			this->m_Transform->TransformPoint(inputPoint);
		InputPointType targetPoint;

		double minimumDistance = NumericTraits< double >::max();
		bool   closestPoint = false;

		// If the closestPoint has not been found, go through the list of fixed
		// points and find the closest distance
		if ( !closestPoint )
		{
			FixedPointIterator pointItr2 = fixedMesh->GetPoints()->Begin();
			FixedPointIterator pointEnd2 = fixedMesh->GetPoints()->End();

			while ( pointItr2 != pointEnd2 )
			{
				double dist = pointItr2.Value().SquaredEuclideanDistanceTo(transformedPoint);


				if ( dist < minimumDistance )
				{
					targetPoint.CastFrom( pointItr2.Value() );
					minimumDistance = dist;
				}
				pointItr2++;
			}
		}

		targetMap[identifier] = targetPoint;

		++pointItr;
		identifier++;
	}

	m_TargetPositionComputed = true;

	//generate a VTK copy of the same mesh
	itkMeshTovtkPolyData* dataTransfer = new itkMeshTovtkPolyData();
	dataTransfer->SetInput(movingMesh);
	movingVTKMesh = dataTransfer->GetOutput();
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
typename ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >::MeasureType
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetValue(const TransformParametersType & parameters) const
{
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

  if ( !fixedMesh )
    {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
    }

  MovingMeshConstPointer movingMesh = this->GetMovingMesh();

  if ( !movingMesh )
    {
    itkExceptionMacro(<< "Moving point set has not been assigned");
    }

  MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingMesh->GetPoints()->End();


  this->SetTransformParameters(parameters);

  // data fidelity energy (squared distance to target position)
  int identifier = 0;
  double functionValue = 0;
  while ( pointItr != pointEnd )
    {
	// get the current position of the vertex
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
	InputVectorType vec;
	vec[0] = parameters[identifier*3];
	vec[1] = parameters[identifier*3+1];
	vec[2] = parameters[identifier*3+2];
    // get the transformed vertex
	typename Superclass::OutputPointType transformedPoint = inputPoint + vec;

    // compute squared Euclidean distance to its target position
	InputPointType targetPoint = targetMap.ElementAt (identifier);
	double dist = targetPoint.SquaredEuclideanDistanceTo(transformedPoint);

	functionValue += dist;

    ++pointItr;
    identifier++;
    }

  // stretching energy
  identifier = 0;
  pointItr = movingMesh->GetPoints()->Begin();
  pointEnd = movingMesh->GetPoints()->End();
  while ( pointItr != pointEnd )
  {
	  vtkSmartPointer<vtkIdList> cellIdList =
		  vtkSmartPointer<vtkIdList>::New();
	  movingVTKMesh->GetPointCells(identifier, cellIdList);

	  //enumerate all the neighboring vertices (edges) of a given vertex
	  //stretching energy : measure the squared derivative along different edge directions
	  //bending energy : measure the local laplacian around the local patch using the given vertex and all neighboring vertices
	  int neighborIdx;
	  double lx = 0; //laplacian
	  double ly = 0;
	  double lz = 0;
	  for(int i = 0; i < cellIdList->GetNumberOfIds(); i++)
	  {
		  vtkSmartPointer<vtkIdList> pointIdList =
			  vtkSmartPointer<vtkIdList>::New();
		  movingVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		  if(pointIdList->GetId(0) != identifier)
			  neighborIdx = pointIdList->GetId(0);
		  else
			  neighborIdx = pointIdList->GetId(1);

		  //derivative
		  double dx = parameters[identifier*3] - parameters[neighborIdx*3];
		  double dy = parameters[identifier*3+1] - parameters[neighborIdx*3+1];
		  double dz = parameters[identifier*3+2] - parameters[neighborIdx*3+2];
          // stretching energy associated with an edge
		  functionValue += m_StretchWeight * (dx*dx + dy*dy + dz*dz);

		  lx += dx; ly += dy; lz += dz;
	  }

      //bending energy associated with a vertex-ring stencil
	  functionValue += m_BendWeight * (lx*lx + ly*ly + lz*lz);
	  ++pointItr;
	  identifier++;
  }

  return functionValue;
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{
	FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

	if ( !fixedMesh )
	{
		itkExceptionMacro(<< "Fixed point set has not been assigned");
	}

	MovingMeshConstPointer movingMesh = this->GetMovingMesh();

	if ( !movingMesh )
	{
		itkExceptionMacro(<< "Moving point set has not been assigned");
	}


	if( derivative.GetSize() != movingMesh->GetNumberOfPoints() * 3 )
	{
		derivative = DerivativeType(movingMesh->GetNumberOfPoints() * 3);
	}
	memset( derivative.data_block(),
		0,
		movingMesh->GetNumberOfPoints() * 3 * sizeof( double ) );

	// derivative of data fidelity energy (squared distance to target position)
	MovingPointIterator pointItr = movingMesh->GetPoints()->Begin();
	MovingPointIterator pointEnd = movingMesh->GetPoints()->End();

	int identifier = 0;
	double functionValue = 0;
	while ( pointItr != pointEnd )
	{
		InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );
		InputVectorType vec;
		vec[0] = parameters[identifier*3];
		vec[1] = parameters[identifier*3+1];
		vec[2] = parameters[identifier*3+2];
		typename Superclass::OutputPointType transformedPoint = inputPoint + vec;

		InputPointType targetPoint = targetMap.ElementAt (identifier);

		typename InputPointType::VectorType distVec = targetPoint - inputPoint;
		derivative[identifier*3]     = -2 * distVec[0];
		derivative[identifier*3 + 1] = -2 * distVec[1];
		derivative[identifier*3 + 2] = -2 * distVec[2];

		++pointItr;
		identifier++;
	}

	// derivative of stretching & bending energy
	identifier = 0;
	pointItr = movingMesh->GetPoints()->Begin();
	pointEnd = movingMesh->GetPoints()->End();
	while ( pointItr != pointEnd )
	{
		vtkSmartPointer<vtkIdList> cellIdList =
			vtkSmartPointer<vtkIdList>::New();
		movingVTKMesh->GetPointCells(identifier, cellIdList);

		int neighborIdx;
		double lx = 0;
		double ly = 0;
		double lz = 0;
		for(int i = 0; i < cellIdList->GetNumberOfIds(); i++)
		{
			vtkSmartPointer<vtkIdList> pointIdList =
				vtkSmartPointer<vtkIdList>::New();
			movingVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

			if(pointIdList->GetId(0) != identifier)
				neighborIdx = pointIdList->GetId(0);
			else
				neighborIdx = pointIdList->GetId(1);

			double dx = parameters[identifier*3] - parameters[neighborIdx*3];
			double dy = parameters[identifier*3+1] - parameters[neighborIdx*3+1];
			double dz = parameters[identifier*3+2] - parameters[neighborIdx*3+2];

            // derivative of stretching energy
			derivative[identifier*3]   += 2 * dx * m_StretchWeight;
			derivative[identifier*3+1] += 2 * dy * m_StretchWeight;
			derivative[identifier*3+2] += 2 * dz * m_StretchWeight;
			derivative[neighborIdx*3]   -= 2 * dx * m_StretchWeight;
			derivative[neighborIdx*3+1] -= 2 * dy * m_StretchWeight;
			derivative[neighborIdx*3+2] -= 2 * dz * m_StretchWeight;

			lx += dx; ly += dy; lz += dz;
		}

		for(int i = 0; i < cellIdList->GetNumberOfIds(); i++)
		{
			vtkSmartPointer<vtkIdList> pointIdList =
				vtkSmartPointer<vtkIdList>::New();
			movingVTKMesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

			if(pointIdList->GetId(0) != identifier)
				neighborIdx = pointIdList->GetId(0);
			else
				neighborIdx = pointIdList->GetId(1);

            // derivative of bending energy
			derivative[identifier*3]   += 2 * lx * m_BendWeight;
			derivative[identifier*3+1] += 2 * ly * m_BendWeight;
			derivative[identifier*3+2] += 2 * lz * m_BendWeight;
			derivative[neighborIdx*3]   -= 2 * lx * m_BendWeight;
			derivative[neighborIdx*3+1] -= 2 * ly * m_BendWeight;
			derivative[neighborIdx*3+2] -= 2 * lz * m_BendWeight;
		}

		++pointItr;
		identifier++;
	}
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
  value = this->GetValue(parameters);
  this->GetDerivative(parameters, derivative);
}

template< typename TFixedMesh, typename TMovingMesh, typename TDistanceMap >
void
ThinShellDemonsMetric< TFixedMesh, TMovingMesh, TDistanceMap >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
