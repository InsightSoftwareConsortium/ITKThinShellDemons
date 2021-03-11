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
#ifndef itkMeshTovtkPolyData_h
#define itkMeshTovtkPolyData_h

#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"

#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include "itkPoint.h"

namespace itk
{
/**
 * \class itkMeshTovtkPolyData
 * \brief
 * \warning
 * \sa
 * \ingroup ThinShellDemons
 */
template<typename TMeshType>
class ITK_TEMPLATE_EXPORT itkMeshTovtkPolyData : public Object
{

 public:
  ITK_DISALLOW_COPY_AND_MOVE(itkMeshTovtkPolyData);

  /** Standard class typedefs. */
  typedef itkMeshTovtkPolyData       Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(itkMeshTovtkPolyData, Object);

  typedef TMeshType MeshType;

  typedef typename MeshType::PointType                PointType;
  typedef typename MeshType::PointsContainer          PointsContainer;
  typedef typename PointsContainer::ConstPointer      PointsContainerPointer;
  typedef typename PointsContainer::ConstIterator     PointsContainerIterator;
  typedef typename MeshType::CellType                 CellType;

  typedef typename MeshType::CellsContainerConstPointer  CellsContainerPointer;
  typedef typename MeshType::CellsContainerConstIterator CellsContainerIterator;

  static vtkSmartPointer<vtkPolyData> Convert(typename MeshType::ConstPointer mesh);

protected:
  itkMeshTovtkPolyData() = default;
  virtual ~itkMeshTovtkPolyData() override = default;

};
}
#ifndef ITK_MANUAL_INSTANTIATION
  #include "itkMeshTovtkPolyData.hxx"
#endif

#endif
