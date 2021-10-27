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

#include "itkPointSetToPointSetMetricWithIndexv4.h"

#include <itkMesh.h>
#include <itkQuadEdge.h>
#include <itkQuadEdgeMesh.h>
#include <itkQuadEdgeMeshExtendedTraits.h>
#include <itkDiscreteGaussianCurvatureQuadEdgeMeshFilter.h>

#include "itkMeshFileWriter.h"
#include "itkMeshFileReader.h"

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
 * fixed points set domain to the moving point set domain.
 *
 *  Reference:
 *  "Thin Shell Demons"
 *  Zhao Q, Price T, Pizer S, Niethammer M, Alterovitz R, Rosenman J
 *  MIUA 2015
 *
 * \ingroup ThinShellDemons
 */
template< typename TFixedMesh, typename TMovingMesh = TFixedMesh,
          class TInternalComputationValueType = double >
class ITK_TEMPLATE_EXPORT ThinShellDemonsMetricv4:
  public PointSetToPointSetMetricWithIndexv4< TFixedMesh, TMovingMesh, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(ThinShellDemonsMetricv4);

  /** Standard class typedefs. */
  typedef ThinShellDemonsMetricv4                                        Self;
  typedef PointSetToPointSetMetricWithIndexv4< TFixedMesh, TMovingMesh > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ThinShellDemonsMetricv4, PointSetToPointSetMetricWithIndexv4);

  /** Types transferred from the base class. */
  typedef typename Superclass::FixedPointSetType      FixedPointSetType;
  typedef typename Superclass::MovingPointSetType     MovingPointSetType;

  /** Types transferred from the base class */
  using MeasureType = typename Superclass::MeasureType;
  using LocalDerivativeType = typename Superclass::LocalDerivativeType;
  using PointType = typename Superclass::PointType;
  using PixelType = typename Superclass::PixelType;
  using PointIdentifier = typename Superclass::PointIdentifier;

  using DimensionType = typename Superclass::DimensionType;
  static constexpr DimensionType FixedPointDimension = Superclass::FixedPointDimension;
  static constexpr DimensionType MovingPointDimension = Superclass::MovingPointDimension;

  using VectorType = typename itk::Vector<double, PointType::Dimension>;

  /* For the QE Mesh and to obtain the curvature using it */
  using CoordType = float;
  using QETraits = typename itk::QuadEdgeMeshExtendedTraits<CoordType, PointType::Dimension, 2, CoordType, CoordType, CoordType, bool, bool>;
  using QEMeshType = typename itk::QuadEdgeMesh<CoordType, PointType::Dimension, QETraits>;
  using QEMeshTypePointer = typename QEMeshType::Pointer;
  using QEPointsContainerPointer = typename QEMeshType::PointsContainerPointer;

  using MeshType = TFixedMesh;
  using MeshTypePointer = typename MeshType::Pointer;
  using MeshPointIdentifier = typename MeshType::PointIdentifier;
  using MeshCellType = typename MeshType::CellType;
  using MeshCellIdentifier = typename MeshType::CellIdentifier;
  using MeshCellPointIdConstIterator = typename MeshCellType::PointIdConstIterator;
  using MeshCellAutoPointer = typename MeshCellType::CellAutoPointer;
  using MeshTriangleCellType = itk::TriangleCell<MeshCellType>;
  using MeshCellLinksContainerIterator = typename MeshType::CellLinksContainerIterator;

  using QEMeshPointType = typename QEMeshType::PointType;
  using QEMeshPointIdentifier = typename QEMeshType::PointIdentifier;
  using QECellType = typename QEMeshType::CellType;
  using QECellAutoPointer = typename QECellType::SelfAutoPointer;
  using QECellIdentifier = typename QEMeshType::CellIdentifier;
  using QETriangleCellType = itk::TriangleCell<QECellType>;

  using TriangleCellType = itk::TriangleCell<QECellType>;
  using TriangleCellAutoPointer = typename TriangleCellType::SelfAutoPointer;

  using CurvatureFilterType = typename itk::DiscreteGaussianCurvatureQuadEdgeMeshFilter<QEMeshType, QEMeshType>;
  using CurvatureFilterTypePointer = typename CurvatureFilterType::Pointer;

  using PointSetPointer = typename Superclass::FixedPointSetType::ConstPointer;

  void Initialize(void) override;

  MeasureType
  GetLocalNeighborhoodValueWithIndex(const PointIdentifier &, const PointType &,
                            const PixelType & pixel = 0) const override;

  void
  GetLocalNeighborhoodValueAndDerivativeWithIndex(const PointIdentifier &, const PointType &,
                                         MeasureType &, LocalDerivativeType &,
                                         const PixelType & pixel = 0) const override;

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
  ThinShellDemonsMetricv4();
  virtual ~ThinShellDemonsMetricv4() override = default;

  //Create a points locator for feature matching
  using FeaturePointSetType = PointSet< double, FixedPointDimension+1>;
  using FeaturePointSetPointer = typename FeaturePointSetType::Pointer;
  using FeaturePointType = typename FeaturePointSetType::PointType;
  using FeaturePointsContainer = typename FeaturePointSetType::PointsContainer;
  using FeaturePointsContainerPointer = typename FeaturePointsContainer::Pointer;
  using FeaturePointsLocatorType = PointsLocator<FeaturePointsContainer>;
  using FeaturePointsLocatorPointer = typename FeaturePointsLocatorType::Pointer;

  mutable FeaturePointsLocatorPointer m_MovingTransformedFeaturePointsLocator;

  /**
   * Prepare point sets for use.
   *
   * Override to use geometric features
   */
  virtual void InitializePointSets() const override;
  void InitializeFeaturePointsLocators() const;

  /**
   * This class uses it's own Points locators to
   * accomodate feature matching
   */
  bool RequiresMovingPointsLocator() const override
  {
    return false;
  };

  bool RequiresFixedPointsLocator() const override
  {
    return false;
  };

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
  typedef std::vector<std::vector<PointIdentifier>> NeighborhoodMap;
  NeighborhoodMap neighborMap;

  typedef std::vector< std::vector<double> > EdgeLengthMap;
  EdgeLengthMap edgeLengthMap;

  mutable MeshTypePointer fixedITKMesh;
  mutable MeshTypePointer movingITKMesh;
  mutable QEMeshTypePointer fixedQEMesh;
  mutable QEMeshTypePointer movingQEMesh;
  mutable QEMeshTypePointer fixedCurvature;

  CurvatureFilterTypePointer gaussian_curvature_filter;

  double m_StretchWeight;
  double m_BendWeight;
  double m_GeometricFeatureWeight;
  mutable double m_ConfidenceSigma;
  bool m_UseConfidenceWeighting;
  bool m_UpdateFeatureMatchingAtEachIteration;
  bool m_UseMaximalDistanceConfidenceSigma;

  void FillPointAndCell(PointSetPointer &pointset, MeshTypePointer &currentITKMesh);
  double ComputeConfidenceValueAndDerivative(const VectorType &v,
                                             VectorType &derivative) const;
  void ComputeStretchAndBend(const PointIdentifier &index,
                             double &stretchEnergy,
                             double &bendEnergy,
                             VectorType &stretch,
                             VectorType &bend) const;
  void ComputeNeighbors();
  void ComputeMaximalDistanceSigma() const;
  FeaturePointType GetFeaturePoint(const double *v, const double &c) const;
  FeaturePointType GetFeaturePoint(const PointType &v, const double &c) const;
  VectorType GetMovingDirection(const PointIdentifier &identifier) const;
  FeaturePointSetPointer GenerateFeaturePointSets(bool fixed) const;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThinShellDemonsMetricv4.hxx"
#endif

#endif
