
if(NOT "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-prefix/src/KWStyle-stamp/KWStyle-gitinfo.txt" IS_NEWER_THAN "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-prefix/src/KWStyle-stamp/KWStyle-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-prefix/src/KWStyle-stamp/KWStyle-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone --no-checkout "https://github.com/Kitware/KWStyle.git" "KWStyle"
    WORKING_DIRECTORY "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/Kitware/KWStyle.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout 840d0f63174d6aa21089a0514ee1e733b98c0709 --
  WORKING_DIRECTORY "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '840d0f63174d6aa21089a0514ee1e733b98c0709'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-prefix/src/KWStyle-stamp/KWStyle-gitinfo.txt"
    "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-prefix/src/KWStyle-stamp/KWStyle-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-prefix/src/KWStyle-stamp/KWStyle-gitclone-lastrun.txt'")
endif()

