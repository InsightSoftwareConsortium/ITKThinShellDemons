# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build

# Include any dependencies generated for this target.
include CMakeFiles/KWStyle.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/KWStyle.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/KWStyle.dir/flags.make

CMakeFiles/KWStyle.dir/kwsStyle.cxx.o: CMakeFiles/KWStyle.dir/flags.make
CMakeFiles/KWStyle.dir/kwsStyle.cxx.o: /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle/kwsStyle.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/KWStyle.dir/kwsStyle.cxx.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/KWStyle.dir/kwsStyle.cxx.o -c /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle/kwsStyle.cxx

CMakeFiles/KWStyle.dir/kwsStyle.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KWStyle.dir/kwsStyle.cxx.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle/kwsStyle.cxx > CMakeFiles/KWStyle.dir/kwsStyle.cxx.i

CMakeFiles/KWStyle.dir/kwsStyle.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KWStyle.dir/kwsStyle.cxx.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle/kwsStyle.cxx -o CMakeFiles/KWStyle.dir/kwsStyle.cxx.s

# Object files for target KWStyle
KWStyle_OBJECTS = \
"CMakeFiles/KWStyle.dir/kwsStyle.cxx.o"

# External object files for target KWStyle
KWStyle_EXTERNAL_OBJECTS =

bin/KWStyle: CMakeFiles/KWStyle.dir/kwsStyle.cxx.o
bin/KWStyle: CMakeFiles/KWStyle.dir/build.make
bin/KWStyle: bin/libKWStyleLib.a
bin/KWStyle: bin/libkwssys.a
bin/KWStyle: CMakeFiles/KWStyle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/KWStyle"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KWStyle.dir/link.txt --verbose=$(VERBOSE)
	/usr/bin/cmake -E copy /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build/bin/KWStyle /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build/KWStyle

# Rule to build all files generated by this target.
CMakeFiles/KWStyle.dir/build: bin/KWStyle

.PHONY : CMakeFiles/KWStyle.dir/build

CMakeFiles/KWStyle.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/KWStyle.dir/cmake_clean.cmake
.PHONY : CMakeFiles/KWStyle.dir/clean

CMakeFiles/KWStyle.dir/depend:
	cd /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build /home/pranjal.sahu/ITK/Modules/External/ITKThinShellDemons/build/KWStyle-build/CMakeFiles/KWStyle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/KWStyle.dir/depend

