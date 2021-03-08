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
#ifndef itkThinShellDemonsMetric_h
#define itkThinShellDemonsMetric_h

#include "itkMeshToMeshMetric.h"
#include "itkCovariantVector.h"
#include "itkMesh.h"
#include "itkImage.h"
#include "itkMeshTovtkPolyData.h"
#include "vtkSmartPointer.h"

namespace itk
{
/** \class ThinShellDemonsMetric
 * \brief This Class inherits the base MeshToMeshMetric
 *
 * \brief This Class inherits the basic MeshToMeshMetric.
 * It expects a mesh-to-mesh transformaton to be plugged in.
 * This class computes a metric value, which is a combination
 * of geometric feature matching quality and the Thin Shell
 * deformation Energy. This metric computation part
 * (objective function) is the core of the Thin Shell Demons
 * algorithm. When initializing a metric object of this class
 * with two meshes, the metric object first pre-computes geometric
 * feature matching between the two meshes. The matching results
 * stay the same during the optimization process.
 *
 *  Reference: "Thin Shell Demons"
 *  Zhao Q, Price T, Pizer S, Niethammer M, Alterovitz R, Rosenman J
 *  MIUA 2015
 *
 * \ingroup ThinShellDemons
 */
template< typename TFixedMesh, typename TMovingMesh = TFixedMesh >
class ITK_TEMPLATE_EXPORT ThinShellDemonsMetric:
  public MeshToMeshMetric< TFixedMesh, TMovingMesh >
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(ThinShellDemonsMetric);

  /** Standard class typedefs. */
  typedef ThinShellDemonsMetric                       Self;
  typedef MeshToMeshMetric< TFixedMesh, TMovingMesh > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ThinShellDemonsMetric, MeshToMeshMetric);

  /** Types transferred from the base class. */
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;

  typedef typename Superclass::MeasureType            MeasureType;
  typedef typename Superclass::DerivativeType         DerivativeType;
  typedef typename Superclass::FixedMeshType          FixedMeshType;
  typedef typename Superclass::MovingMeshType         MovingMeshType;
  typedef typename Superclass::FixedMeshConstPointer  FixedMeshConstPointer;
  typedef typename Superclass::MovingMeshConstPointer MovingMeshConstPointer;

  typedef typename Superclass::FixedPointIterator     FixedPointIterator;
  typedef typename Superclass::FixedPointDataIterator FixedPointDataIterator;

  typedef typename Superclass::MovingPointIterator     MovingPointIterator;
  typedef typename Superclass::MovingPointDataIterator MovingPointDataIterator;

  typedef typename Superclass::InputPointType InputPointType;
  typedef typename InputPointType::VectorType InputVectorType;

  typedef typename Superclass::OutputPointType OutputPointType;


  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get the match measure, i.e. the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType & Value,
                             DerivativeType & Derivative) const ITK_OVERRIDE;

  /** Initialize the Metric by computing a target position
   * for each vertex in the fixed mesh using
   * Euclidean + Curvature distance */
  virtual void Initialize(void) throw ( ExceptionObject ) ITK_OVERRIDE;

  /**
   * Stretching penalty weight
   */
  itkSetMacro(StretchWeight, double);
  itkGetConstReferenceMacro(StretchWeight, double);

  /**
   * Bending penalty weight
   */
  itkSetMacro(BendWeight, double);
  itkGetConstReferenceMacro(BendWeight, double);

  /**
   * Weight for curvature match term in geomntric feature matching.
   * Feature distance is:
   * euclidean distance + GeometricFeatureWeight * curvature distance
   */
  itkSetMacro(GeometricFeatureWeight, double);
  itkGetConstReferenceMacro(GeometricFeatureWeight, double);

   /**
   * Update feature match at each iteration.
   *
   * When used in conjunction with UseConfidenceWeighting and
   * UpdateFeatureMatchingAtEachIteration this can lead to unexpected results.
   * If the confidence sigma is smaller than sqrt(2) * the maixmal feature
   * distance the derivate leads to points being pushed away.
   */
  itkSetMacro(ConfidenceSigma, double);
  itkGetConstReferenceMacro(ConfidenceSigma, double);

  /**
   * Update feature match at each iteration.
   *
   * When used in conjunction with UseConfidenceWeighting
   * this can lead to unexpected results. If the confidence sigma is smaller
   * than sqrt(2) * the maixmal feature distance the derivate leads to
   * points being pushed away.
   */
  itkSetMacro(UpdateFeatureMatchingAtEachIteration, bool);
  itkGetConstReferenceMacro(UpdateFeatureMatchingAtEachIteration, bool);
  itkBooleanMacro(UpdateFeatureMatchingAtEachIteration);

  /**
   * Weight cost function by feature distance.
   *
   * When used in conjunction with  UpdateFeatureMatchingAtEachIteration
   * this can lead to unexpected results. If the confidence sigma is smaller
   * than sqrt(2) * the maixmal feature distance the derivate leads to points being
   * pushed away.
   */
  itkSetMacro(UseConfidenceWeighting, bool);
  itkGetConstReferenceMacro(UseConfidenceWeighting, bool);
  itkBooleanMacro(UseConfidenceWeighting);

  /**
   * Automatically compute confidence sigam at each iteration to ensure that
   * no points are being pushed away from each other. The confidence sigma is
   * set to sqrt(2) * the maxmial feature distance
   */
  itkSetMacro(UseMaximalDistanceConfidenceSigma, bool);
  itkGetConstReferenceMacro(UseMaximalDistanceConfidenceSigma, bool);
  itkBooleanMacro(UseMaximalDistanceConfidenceSigma);


protected:
  ThinShellDemonsMetric();
  virtual ~ThinShellDemonsMetric() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:

  typedef std::vector<InputPointType> TargetMap;
  mutable TargetMap targetMap;

  typedef std::vector< vtkSmartPointer<vtkIdList> > NeighborhodMap;
  NeighborhodMap neighborMap;

  vtkSmartPointer<vtkPolyData> movingVTKMesh;
  vtkSmartPointer<vtkPolyData> fixedVTKMesh;
  vtkSmartPointer<vtkPolyData> fixedCurvature;

  mutable bool updateConfidenceSigma;
  double m_StretchWeight;
  double m_BendWeight;
  double m_GeometricFeatureWeight;
  mutable double m_ConfidenceSigma;
  bool m_UseConfidenceWeighting;
  bool m_UpdateFeatureMatchingAtEachIteration;
  bool m_UseMaximalDistanceConfidenceSigma;

  double ComputeConfidenceValueAndDerivative(const InputVectorType &v,
                                             InputVectorType &derivative) const;

  void ComputeStretchAndBend(int identifier,
                             const TransformParametersType &parmaters,
                             double &stretchEnergy,
                             double &bendEnergy,
                             InputVectorType &stretch,
                             InputVectorType &bend) const;

  void ComputeTargetPosition() const;

  void ComputeNeighbors();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThinShellDemonsMetric.hxx"
#endif

#endif
