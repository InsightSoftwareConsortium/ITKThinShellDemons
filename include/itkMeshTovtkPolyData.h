#ifndef __itkMeshTovtkPolyData_h__
#define __itkMeshTovtkPolyData_h__

#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"

#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include "itkPoint.h"


/**
  \class itkMeshTovtkPolyData
  \brief
    \warning
  \sa
  */

class itkMeshTovtkPolyData
{

 public:

  itkMeshTovtkPolyData( void );
  virtual ~itkMeshTovtkPolyData( void );

  typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
  typedef itk::Mesh<double,3> TriangleMeshType;
  typedef TriangleMeshType::PointType                       PointType;
  typedef TriangleMeshType::PointsContainer                 InputPointsContainer;
  typedef InputPointsContainer::ConstPointer            InputPointsContainerPointer;
  typedef InputPointsContainer::ConstIterator           InputPointsContainerIterator;
  typedef TriangleMeshType::CellType                        CellType;

  typedef TriangleMeshType::CellsContainerConstPointer           CellsContainerPointer;
  typedef TriangleMeshType::CellsContainerConstIterator          CellsContainerIterator;
  /**
  The SetInput method provides pointer to the vtkPolyData
  */
  void SetInput(TriangleMeshType::ConstPointer mesh);
  vtkPolyData * GetOutput();
  void ConvertitkTovtk();

  TriangleMeshType::ConstPointer m_itkTriangleMesh;

  vtkPoints  * m_Points;
  vtkPolyData * m_PolyData;
  vtkCellArray * m_Polys;

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshTovtkPolyData.hxx"
#endif

#endif
