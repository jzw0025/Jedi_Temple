# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.9.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.9.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a

# Include any dependencies generated for this target.
include CMakeFiles/main2a.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main2a.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main2a.dir/flags.make

CMakeFiles/main2a.dir/main2a.cc.o: CMakeFiles/main2a.dir/flags.make
CMakeFiles/main2a.dir/main2a.cc.o: main2a.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main2a.dir/main2a.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main2a.dir/main2a.cc.o -c /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a/main2a.cc

CMakeFiles/main2a.dir/main2a.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main2a.dir/main2a.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a/main2a.cc > CMakeFiles/main2a.dir/main2a.cc.i

CMakeFiles/main2a.dir/main2a.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main2a.dir/main2a.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a/main2a.cc -o CMakeFiles/main2a.dir/main2a.cc.s

CMakeFiles/main2a.dir/main2a.cc.o.requires:

.PHONY : CMakeFiles/main2a.dir/main2a.cc.o.requires

CMakeFiles/main2a.dir/main2a.cc.o.provides: CMakeFiles/main2a.dir/main2a.cc.o.requires
	$(MAKE) -f CMakeFiles/main2a.dir/build.make CMakeFiles/main2a.dir/main2a.cc.o.provides.build
.PHONY : CMakeFiles/main2a.dir/main2a.cc.o.provides

CMakeFiles/main2a.dir/main2a.cc.o.provides.build: CMakeFiles/main2a.dir/main2a.cc.o


# Object files for target main2a
main2a_OBJECTS = \
"CMakeFiles/main2a.dir/main2a.cc.o"

# External object files for target main2a
main2a_EXTERNAL_OBJECTS =

main2a: CMakeFiles/main2a.dir/main2a.cc.o
main2a: CMakeFiles/main2a.dir/build.make
main2a: /Users/junchaowei/Desktop/Jedi_Temple/dealii/lib/libdeal_II.g.8.5.0-pre.dylib
main2a: /usr/lib/libbz2.dylib
main2a: /usr/local/lib/libtbb.dylib
main2a: /usr/lib/libz.dylib
main2a: /usr/local/lib/libarpack.dylib
main2a: /usr/local/lib/libhdf5_hl.dylib
main2a: /usr/local/lib/libhdf5.dylib
main2a: CMakeFiles/main2a.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main2a"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main2a.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main2a.dir/build: main2a

.PHONY : CMakeFiles/main2a.dir/build

CMakeFiles/main2a.dir/requires: CMakeFiles/main2a.dir/main2a.cc.o.requires

.PHONY : CMakeFiles/main2a.dir/requires

CMakeFiles/main2a.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main2a.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main2a.dir/clean

CMakeFiles/main2a.dir/depend:
	cd /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a /Users/junchaowei/Desktop/testFEM2/CA2_template/CA2a/CMakeFiles/main2a.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main2a.dir/depend
