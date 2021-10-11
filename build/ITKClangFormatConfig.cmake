set(ITK_SOURCE_DIR "")
set(CMAKE_SOURCE_DIR "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons")
set(CLANG_FORMAT_EXECUTABLE "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/ClangFormat/clang-format")

set(WORKING_DIR "")
if(ITK_SOURCE_DIR)
  set(WORKING_DIR "${ITK_SOURCE_DIR}")
else()
  set(WORKING_DIR "${CMAKE_SOURCE_DIR}")
endif()

find_package(Git)
if(GIT_FOUND AND EXISTS "${WORKING_DIR}/.git/config")
  execute_process(COMMAND ${GIT_EXECUTABLE} config clangFormat.binary
    "${CLANG_FORMAT_EXECUTABLE}"
    WORKING_DIRECTORY ${WORKING_DIR})
endif()
