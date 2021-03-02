itk_module_test()

set(ThinShellDemonsTests
  itkThinShellDemonsTest.cxx
  itkThinShellDemonsTestv4_0.cxx
  )

CreateTestDriver(ThinShellDemons "${ThinShellDemons-Test_LIBRARIES}" "${ThinShellDemonsTests}")

itk_add_test(NAME itkThinShellDemonsTest
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTest
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )

itk_add_test(NAME itkThinShellDemonsTestv4_0
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTestv4_0
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )
