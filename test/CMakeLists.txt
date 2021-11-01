itk_module_test()

set(ThinShellDemonsTests
  itkThinShellDemonsTestv4_Affine.cxx
  itkThinShellDemonsTestv4_Displacement.cxx
  itkThinShellDemonsTestv4_SyN.cxx
  )

CreateTestDriver(ThinShellDemons "${ThinShellDemons-Test_LIBRARIES}" "${ThinShellDemonsTests}")


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


itk_add_test(NAME itkThinShellDemonsTestv4_SyN
  COMMAND ThinShellDemonsTestDriver
  itkThinShellDemonsTestv4_SyN
  DATA{Baseline/fixedMesh.vtk}
  DATA{Baseline/movingMesh.vtk}
  )
