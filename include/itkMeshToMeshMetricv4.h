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
#ifndef itkMeshToMeshMetricv4_h
#define itkMeshToMeshMetricv4_h

#include "itkPointSetToPointSetMetricWithIndexv4.h"

namespace itk
{
/** \class MeshToMeshMetricv4
 *
 * \brief Extension of PointSetToPointSetMetricv4 for use with meshes
 * provides convenience methods to get and set moving and fixed mesh.
 *
 * \ingroup ThinShellDemons
 */
template< typename TFixedMesh,  typename TMovingMesh,
          class TInternalComputationValueType = double >
class ITK_TEMPLATE_EXPORT MeshToMeshMetricv4:
  public PointSetToPointSetMetricWithIndexv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(MeshToMeshMetricv4);

  /** Standard class typedefs. */
  typedef MeshToMeshMetricv4         Self;
  typedef PointSetToPointSetMetricWithIndexv4<TFixedMesh, TMovingMesh, TInternalComputationValueType>
                                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;


  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshToMeshMetricv4, PointSetToPointSetMetricWithIndexv4);

  /**  Type of the moving Mesh. */
  typedef TMovingMesh                           MovingMeshType;
  typedef typename MovingMeshType::ConstPointer MovingMeshConstPointer;

  /**  Type of the fixed Mesh. */
  typedef TFixedMesh                           FixedMeshType;
  typedef typename FixedMeshType::ConstPointer FixedMeshConstPointer;

  /** Constants for the Mesh dimensions */
  itkStaticConstMacro(MovingMeshDimension, unsigned int, TMovingMesh::PointDimension);
  itkStaticConstMacro(FixedMeshDimension, unsigned int, TFixedMesh::PointDimension);


  /** Get/Set the Fixed Mesh.  */
  void SetFixedMesh(FixedMeshConstPointer fixedMesh){
    this->SetFixedPointSet(fixedMesh);
  };

  FixedMeshConstPointer GetFixedMesh() const{
    return this->GetFixedPointSet();
  };

  /** Get/Set the Moving Mesh.  */
  void SetMovingMesh(MovingMeshConstPointer fixedMesh){
    this->SetMovingPointSet(fixedMesh);
  };

  MovingMeshConstPointer GetMovingMesh() const{
    return this->GetMovingPointSet();
  };


protected:
  MeshToMeshMetricv4() = default;
  virtual ~MeshToMeshMetricv4() override = default;
  virtual void PrintSelf(std::ostream &os, Indent indent) const override;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMeshToMeshMetricv4.hxx"
#endif

#endif
