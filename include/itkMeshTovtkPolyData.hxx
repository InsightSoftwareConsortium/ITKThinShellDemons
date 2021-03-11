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
#ifndef itkMeshTovtkPolyData_hxx
#define itkMeshTovtkPolyData_hxx

#include <iostream>
#include "itkMeshTovtkPolyData.h"

namespace itk{

template<typename TMeshType>
vtkSmartPointer<vtkPolyData>
itkMeshTovtkPolyData<TMeshType>
::Convert(typename MeshType::ConstPointer mesh )
{
  int numPoints = mesh->GetNumberOfPoints();

  PointsContainerPointer  myPoints = mesh->GetPoints();
  PointsContainerIterator points = myPoints->Begin();
  PointType point;

  vtkSmartPointer<vtkPoints>    m_Points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> m_Polys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData>  m_PolyData = vtkSmartPointer<vtkPolyData>::New();

  if( PointType::Dimension > 3)
  {
    return m_PolyData;
  }
  if (numPoints == 0)
  {
    return m_PolyData;
  }
  m_Points->SetNumberOfPoints(numPoints);

  int idx=0;
  double vpoint[3] = {0};
  unsigned int three = 3;
  while( points != myPoints->End() )
  {
    point = points.Value();
    for(unsigned int i =0; i<std::min(PointType::Dimension, three); i++)
    {
      vpoint[i]= point[i];
    }
    m_Points->SetPoint(idx++,vpoint);
    points++;
  }

  m_PolyData->SetPoints(m_Points);

  CellsContainerPointer cells = mesh->GetCells();
  CellsContainerIterator cellIt = cells->Begin();
  while ( cellIt != cells->End() )
  {
    CellType *nextCell = cellIt->Value();

    vtkIdType *pts = new vtkIdType[nextCell->GetNumberOfPoints()];
    typename CellType::PointIdIterator pointIt = nextCell->PointIdsBegin();
    int i=0;
    while (pointIt != nextCell->PointIdsEnd() )
    {
      pts[i++] = *pointIt++;
    }
    m_Polys->InsertNextCell(nextCell->GetNumberOfPoints(), pts);
    delete[] pts;
    cellIt++;
  }

  m_PolyData->SetPolys(m_Polys);
  return m_PolyData;
}
}
#endif
