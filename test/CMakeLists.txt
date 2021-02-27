itk_module_test()

set(ThinShellDemonsTests
  itkThinShellDemonsTest.cxx
  )

CreateTestDriver(ThinShellDemons "${ThinShellDemons-Test_LIBRARIES}" "${ThinShellDemonsTests}")

itk_add_test(NAME itkThinShellDemonsTest
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTest
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )
