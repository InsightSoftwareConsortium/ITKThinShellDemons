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
#ifndef itkThinShellDemonsMetricv4_h
#define itkThinShellDemonsMetricv4_h

#include "itkMeshToMeshMetricv4.h"
#include "itkMesh.h"
#include "itkMeshTovtkPolyData.h"
#include "vtkSmartPointer.h"

namespace itk
{
/** \class ThinShellDemonsMetricv4
 *
 * \brief
 * Thin shelll demons metric for use with v4 registration framework.
 * This metric implements a mesh regularization term that penalizes mesh distortions
 * while trying to iteratively match closest points. Currently geometric features matching
 * is not implemented in the closest point computation.
 *
 * Note that in ITK the registration by convection computes a transfrom from
 * the fixed object to the moving object. For image this is natural
 * since to transform the moving image into the fixed domain the the image is
 * pulled back thorugh the transform.
 * For point set registration this can be confusing since certain transforms
 * such a mesh displacements do not readily allow to pull back the moving point
 * set and it is easier to transform the "fixed" point set to the moving point set
 * by tranpsoerting points forward along the transform.
 *
 * For the thin shell demons metric here that means that the regularization is
 * applied to the fixed point set and the transfrom moves the points of the
 * fixed points set to the moving point set
 *
 *  Reference:
 *  "Thin Shell Demons"
 *  Zhao Q, Price T, Pizer S, Niethammer M, Alterovitz R, Rosenman J
 *  MIUA 2015
 *
 * \ingroup ThinShellDemons
 */
template< typename TFixedMesh, typename TMovingMesh >
class ITK_TEMPLATE_EXPORT ThinShellDemonsMetricv4:
  public MeshToMeshMetricv4< TFixedMesh, TMovingMesh >
{
public:

  /** Standard class typedefs. */
  typedef ThinShellDemonsMetricv4                       Self;
  typedef MeshToMeshMetricv4< TFixedMesh, TMovingMesh > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ThinShellDemonsMetricv4, MeshToMeshMetricv4);

  /** Types transferred from the base class. */
  typedef typename Superclass::FixedMeshType          FixedMeshType;
  typedef typename Superclass::MovingMeshType         MovingMeshType;
  typedef typename Superclass::FixedMeshConstPointer  FixedMeshConstPointer;
  typedef typename Superclass::MovingMeshConstPointer MovingMeshConstPointer;

  /** Types transferred from the base class */
  using MeasureType = typename Superclass::MeasureType;
  using LocalDerivativeType = typename Superclass::LocalDerivativeType;
  using PointType = typename Superclass::PointType;
  using PixelType = typename Superclass::PixelType;
  using PointIdentifier = typename Superclass::PointIdentifier;

  using VectorType = typename itk::Vector<double, PointType::Dimension>;

  /**
   * Calculates the local metric value for a single point.
   */
  MeasureType
  GetLocalNeighborhoodValue(const PointType &, const PixelType & pixel = 0) const override;

  /**
   * Calculates the local value and derivative for a single point.
   */
  void
  GetLocalNeighborhoodValueAndDerivative(const PointType &,
                                         MeasureType &,
                                         LocalDerivativeType &,
                                         const PixelType & pixel = 0) const override;

  void Initialize(void) override;

  /** Set/Get algorithm parameters **/
  itkSetMacro(StretchWeight, double);
  itkGetMacro(StretchWeight, double);

  itkSetMacro(BendWeight, double);
  itkGetMacro(BendWeight, double);

  itkSetMacro(ConfidenceSigma, double);
  itkGetMacro(ConfidenceSigma, double);

  itkSetMacro(GeometricFeatureWeight, double);
  itkGetMacro(GeometricFeatureWeight, double);


protected:
  ThinShellDemonsMetricv4();
  virtual ~ThinShellDemonsMetricv4() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ThinShellDemonsMetricv4);

  typedef itk::MapContainer<int, PointType> TargetMapType;
  TargetMapType targetMap;

  typedef itk::MapContainer<int, vtkSmartPointer<vtkIdList>> NeighborhodMapType;
  NeighborhodMapType neighborMap;

  vtkSmartPointer<vtkPolyData> movingVTKMesh;
  vtkSmartPointer<vtkPolyData> fixedVTKMesh;
  //vtkSmartPointer<vtkPolyData> fixedCurvature;

  double m_ConfidenceSigma;
  double m_StretchWeight;
  double m_BendWeight;
  double m_GeometricFeatureWeight;

  void ComputeStretchAndBend(const PointType &point,
                             double &stretchEnergy,
                             double &bendEnergy,
                             VectorType &stretch,
                             VectorType &bend) const;
  void ComputeTargetPosition() const;

  void ComputeNeighbors();

  VectorType GetMovingDirection(int identifier) const;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThinShellDemonsMetricv4.hxx"
#endif

#endif
