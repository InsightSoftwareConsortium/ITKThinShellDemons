# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

function(check_file_hash has_hash hash_is_good)
  if("${has_hash}" STREQUAL "")
    message(FATAL_ERROR "has_hash Can't be empty")
  endif()

  if("${hash_is_good}" STREQUAL "")
    message(FATAL_ERROR "hash_is_good Can't be empty")
  endif()

  if("SHA512" STREQUAL "")
    # No check
    set("${has_hash}" FALSE PARENT_SCOPE)
    set("${hash_is_good}" FALSE PARENT_SCOPE)
    return()
  endif()

  set("${has_hash}" TRUE PARENT_SCOPE)

  message(STATUS "verifying file...
       file='/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar'")

  file("SHA512" "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar" actual_value)

  if(NOT "${actual_value}" STREQUAL "b14de32036c48f6c62998e2ebab509e71a0ae71464acb4616484e3a6eb941e1d9fac38559f5d27ea0cbbb512d590279ffb3015fae17779229e1090c2763ebcf3")
    set("${hash_is_good}" FALSE PARENT_SCOPE)
    message(STATUS "SHA512 hash of
    /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar
  does not match expected value
    expected: 'b14de32036c48f6c62998e2ebab509e71a0ae71464acb4616484e3a6eb941e1d9fac38559f5d27ea0cbbb512d590279ffb3015fae17779229e1090c2763ebcf3'
      actual: '${actual_value}'")
  else()
    set("${hash_is_good}" TRUE PARENT_SCOPE)
  endif()
endfunction()

function(sleep_before_download attempt)
  if(attempt EQUAL 0)
    return()
  endif()

  if(attempt EQUAL 1)
    message(STATUS "Retrying...")
    return()
  endif()

  set(sleep_seconds 0)

  if(attempt EQUAL 2)
    set(sleep_seconds 5)
  elseif(attempt EQUAL 3)
    set(sleep_seconds 5)
  elseif(attempt EQUAL 4)
    set(sleep_seconds 15)
  elseif(attempt EQUAL 5)
    set(sleep_seconds 60)
  elseif(attempt EQUAL 6)
    set(sleep_seconds 90)
  elseif(attempt EQUAL 7)
    set(sleep_seconds 300)
  else()
    set(sleep_seconds 1200)
  endif()

  message(STATUS "Retry after ${sleep_seconds} seconds (attempt #${attempt}) ...")

  execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep "${sleep_seconds}")
endfunction()

if("/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if("https://data.kitware.com/api/v1/file/hashsum/sha512/b14de32036c48f6c62998e2ebab509e71a0ae71464acb4616484e3a6eb941e1d9fac38559f5d27ea0cbbb512d590279ffb3015fae17779229e1090c2763ebcf3/download" STREQUAL "")
  message(FATAL_ERROR "REMOTE can't be empty")
endif()

if(EXISTS "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar")
  check_file_hash(has_hash hash_is_good)
  if(has_hash)
    if(hash_is_good)
      message(STATUS "File already exists and hash match (skip download):
  file='/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar'
  SHA512='b14de32036c48f6c62998e2ebab509e71a0ae71464acb4616484e3a6eb941e1d9fac38559f5d27ea0cbbb512d590279ffb3015fae17779229e1090c2763ebcf3'"
      )
      return()
    else()
      message(STATUS "File already exists but hash mismatch. Removing...")
      file(REMOVE "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar")
    endif()
  else()
    message(STATUS "File already exists but no hash specified (use URL_HASH):
  file='/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar'
Old file will be removed and new file downloaded from URL."
    )
    file(REMOVE "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar")
  endif()
endif()

set(retry_number 5)

message(STATUS "Downloading...
   dst='/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar'
   timeout='none'"
)

foreach(i RANGE ${retry_number})
  sleep_before_download(${i})

  foreach(url https://data.kitware.com/api/v1/file/hashsum/sha512/b14de32036c48f6c62998e2ebab509e71a0ae71464acb4616484e3a6eb941e1d9fac38559f5d27ea0cbbb512d590279ffb3015fae17779229e1090c2763ebcf3/download)
    message(STATUS "Using src='${url}'")

    
    
    
    

    file(
        DOWNLOAD
        "${url}" "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar"
        SHOW_PROGRESS
        # no TIMEOUT
        STATUS status
        LOG log
        
        
    )

    list(GET status 0 status_code)
    list(GET status 1 status_string)

    if(status_code EQUAL 0)
      check_file_hash(has_hash hash_is_good)
      if(has_hash AND NOT hash_is_good)
        message(STATUS "Hash mismatch, removing...")
        file(REMOVE "/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/ClangFormat-prefix/src/archive.tar")
      else()
        message(STATUS "Downloading... done")
        return()
      endif()
    else()
      string(APPEND logFailedURLs "error: downloading '${url}' failed
       status_code: ${status_code}
       status_string: ${status_string}
       log:
       --- LOG BEGIN ---
       ${log}
       --- LOG END ---
       "
      )
    endif()
  endforeach()
endforeach()

message(FATAL_ERROR "Each download failed!
  ${logFailedURLs}
  "
)
