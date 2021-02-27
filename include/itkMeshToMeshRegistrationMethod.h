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
#ifndef itkMeshToMeshRegistrationMethod_h
#define itkMeshToMeshRegistrationMethod_h

#include "itkProcessObject.h"
#include "itkMeshToMeshMetric.h"
#include "itkDataObjectDecorator.h"
#include "itkSingleValuedNonLinearOptimizer.h"

namespace itk {
    /** \class MeshToMeshRegistrationMethod
     *  \brief This class is templated over pointset-to-pointset registration method
     *
     */

template< typename TFixedMesh, typename TMovingMesh>
class ITK_TEMPLATE_EXPORT MeshToMeshRegistrationMethod : public ProcessObject
{
public:
	/** Standard class typedefs */
	typedef MeshToMeshRegistrationMethod Self;
	typedef ProcessObject                        Superclass;
	typedef SmartPointer< Self >                 Pointer;
	typedef SmartPointer< const Self >           ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(MeshToMeshRegistrationMethod, ProcessObject);

	/**  Type of the Fixed Mesh. */
	typedef          TFixedMesh                  FixedMeshType;
	typedef typename FixedMeshType::ConstPointer FixedMeshConstPointer;

	/**  Type of the Moving Mesh. */
	typedef          TMovingMesh                  MovingMeshType;
	typedef typename MovingMeshType::Pointer MovingMeshConstPointer;

	/**  Type of the Metric. */
	typedef MeshToMeshMetric< FixedMeshType, MovingMeshType > MetricType;
	typedef typename MetricType::Pointer                                      MetricPointer;

	  /**  Type of the Transform . */
	typedef  typename MetricType::TransformType TransformType;
	typedef  typename TransformType::Pointer    TransformPointer;

	/** Type for the output: Using Decorator pattern for enabling
	*  the Transform to be passed in the data pipeline. */
	typedef  DataObjectDecorator< TransformType >      TransformOutputType;
	typedef typename TransformOutputType::Pointer      TransformOutputPointer;
	typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;

	/**  Type of the Optimizer. */
	typedef	SingleValuedNonLinearOptimizer OptimizerType;
	//typedef	typename	OptimizerType::InternalOptimizerType  vnlOptimizerType;

	/** Type of the Transformation parameters. This is the same type used to
	*  represent the search space of the optimization algorithm. */
	typedef  typename MetricType::TransformParametersType ParametersType;

	/** Set/Get the Fixed Mesh. */
	itkSetConstObjectMacro(FixedMesh, FixedMeshType);
	itkGetConstObjectMacro(FixedMesh, FixedMeshType);

	/** Set/Get the Moving Mesh. */
	itkSetObjectMacro(MovingMesh, MovingMeshType);
	itkGetObjectMacro(MovingMesh, MovingMeshType);

	/** Set/Get the Optimizer. */
	itkSetObjectMacro(Optimizer,  OptimizerType);
	itkGetModifiableObjectMacro(Optimizer, OptimizerType);

	/** Set/Get the Metric. */
	itkSetObjectMacro(Metric, MetricType);
	itkGetModifiableObjectMacro(Metric, MetricType);

	/** Set/Get the Transfrom. */
	itkSetObjectMacro(Transform, TransformType);
	itkGetModifiableObjectMacro(Transform, TransformType);

	/** Set/Get the initial transformation parameters. */
	virtual void SetInitialTransformParameters(const ParametersType & param);

	/** Get the last transformation parameters visited by
	* the optimizer. */
	itkGetConstReferenceMacro(LastTransformParameters, ParametersType);

	/** Initialize by setting the interconnects between the components. */
	void Initialize()
	throw ( ExceptionObject );

	/** Returns the transform resulting from the registration process  */
	const TransformOutputType * GetOutput() const;

	/** Deforms the pointset of the moving mesh using the resulting transformation */
	void UpdateMovingMesh();

	/** Make a DataObject of the correct type to be used as the specified
	* output. */
	typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
	using Superclass::MakeOutput;
	virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;
	int test1;
#ifdef ITKV3_COMPATIBILITY
  /** Method that initiates the registration. This will Initialize and ensure
   * that all inputs the registration needs are in place, via a call to
   * Initialize() will then start the optimization process via a call to
   * StartOptimization()
   * StartRegistration is an old API from before
   * ImageRegistrationMethod was a subclass of ProcessObject.
   * Historically, one could call StartRegistration() instead of
   * calling Update().  However, when called directly by the user, the
   * inputs to ImageRegistrationMethod may not be up to date.  This
   * may cause an unexpected behavior.
   *
   * Since we cannot eliminate StartRegistration for backward
   * compatibility reasons, we check whether StartRegistration was
   * called directly or whether Update() (which in turn called
   * StartRegistration()). */
  void StartRegistration(void) { this->Update(); }
#endif

protected:
	MeshToMeshRegistrationMethod();
	virtual ~MeshToMeshRegistrationMethod() {}

	virtual void GenerateData() ITK_OVERRIDE;

private:

	MetricPointer          m_Metric;
	OptimizerType::Pointer m_Optimizer;

	TransformPointer m_Transform;
	MovingMeshConstPointer m_MovingMesh;
	FixedMeshConstPointer  m_FixedMesh;

	ParametersType m_InitialTransformParameters;
	ParametersType m_LastTransformParameters;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshToMeshRegistrationMethod.hxx"
#endif

#endif
