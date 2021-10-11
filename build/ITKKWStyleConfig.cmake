set(ITK_SOURCE_DIR "")
set(CMAKE_SOURCE_DIR "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons")
set(KWSTYLE_EXECUTABLE "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build/KWStyle")

set(WORKING_DIR "")
if(ITK_SOURCE_DIR)
  set(WORKING_DIR "${ITK_SOURCE_DIR}")
else()
  set(WORKING_DIR "${CMAKE_SOURCE_DIR}")
endif()

find_package( Git )
if(GIT_FOUND AND EXISTS "${ITK_SOURCE_DIR}/.git/config")
  execute_process( COMMAND ${GIT_EXECUTABLE} config hooks.KWStyle.path
    "${KWSTYLE_EXECUTABLE}"
    WORKING_DIRECTORY ${WORKING_DIR} )
endif()
