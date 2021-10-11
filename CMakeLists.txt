cmake_minimum_required(VERSION 3.10.2)
project(ThinShellDemons)

#set(ThinShellDemons_LIBRARIES ThinShellDemons)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
  message(STATUS "ITK_LIBRARIES: ${ITK_LIBRARIES}")
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()
