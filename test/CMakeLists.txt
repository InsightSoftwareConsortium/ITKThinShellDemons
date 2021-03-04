itk_module_test()

set(ThinShellDemonsTests
  itkThinShellDemonsTest.cxx
  itkThinShellDemonsTestv4_Affine.cxx
  itkThinShellDemonsTestv4_Displacement.cxx
  )

CreateTestDriver(ThinShellDemons "${ThinShellDemons-Test_LIBRARIES}" "${ThinShellDemonsTests}")

itk_add_test(NAME itkThinShellDemonsTest
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTest
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )

itk_add_test(NAME itkThinShellDemonsTestv4_Affine
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTestv4_Affine
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )


itk_add_test(NAME itkThinShellDemonsTestv4_Displacement
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTestv4_Displacement
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )
