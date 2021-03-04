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
template<typename T>
class ITK_TEMPLATE_EXPORT itkMeshTovtkPolyData : public Object
{

 public:
  /** Standard class typedefs. */
  typedef itkMeshTovtkPolyData       Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(itkMeshTovtkPolyData, Object);

  typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
  typedef itk::Mesh<double,3>                                       TriangleMeshType;

  typedef TriangleMeshType::PointType         PointType;
  typedef TriangleMeshType::PointsContainer   InputPointsContainer;
  typedef InputPointsContainer::ConstPointer  InputPointsContainerPointer;
  typedef InputPointsContainer::ConstIterator InputPointsContainerIterator;
  typedef TriangleMeshType::CellType          CellType;

  typedef TriangleMeshType::CellsContainerConstPointer  CellsContainerPointer;
  typedef TriangleMeshType::CellsContainerConstIterator CellsContainerIterator;
  /**
  The SetInput method provides pointer to the vtkPolyData
  */
  void SetInput(TriangleMeshType::ConstPointer mesh);
  vtkSmartPointer<vtkPolyData> GetOutput();

protected:
  itkMeshTovtkPolyData();
  virtual ~itkMeshTovtkPolyData() override = default;
  TriangleMeshType::ConstPointer m_itkTriangleMesh;
  vtkSmartPointer<vtkPolyData> ConvertitkTovtk();

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(itkMeshTovtkPolyData);
};
}
#ifndef ITK_MANUAL_INSTANTIATION
  #include "itkMeshTovtkPolyData.hxx"
#endif

#endif
