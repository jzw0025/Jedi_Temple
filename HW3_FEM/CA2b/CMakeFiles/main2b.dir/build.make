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
CMAKE_SOURCE_DIR = "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b"

# Include any dependencies generated for this target.
include CMakeFiles/main2b.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main2b.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main2b.dir/flags.make

CMakeFiles/main2b.dir/main2b.cc.o: CMakeFiles/main2b.dir/flags.make
CMakeFiles/main2b.dir/main2b.cc.o: main2b.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main2b.dir/main2b.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main2b.dir/main2b.cc.o -c "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/main2b.cc"

CMakeFiles/main2b.dir/main2b.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main2b.dir/main2b.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/main2b.cc" > CMakeFiles/main2b.dir/main2b.cc.i

CMakeFiles/main2b.dir/main2b.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main2b.dir/main2b.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/main2b.cc" -o CMakeFiles/main2b.dir/main2b.cc.s

CMakeFiles/main2b.dir/main2b.cc.o.requires:

.PHONY : CMakeFiles/main2b.dir/main2b.cc.o.requires

CMakeFiles/main2b.dir/main2b.cc.o.provides: CMakeFiles/main2b.dir/main2b.cc.o.requires
	$(MAKE) -f CMakeFiles/main2b.dir/build.make CMakeFiles/main2b.dir/main2b.cc.o.provides.build
.PHONY : CMakeFiles/main2b.dir/main2b.cc.o.provides

CMakeFiles/main2b.dir/main2b.cc.o.provides.build: CMakeFiles/main2b.dir/main2b.cc.o


# Object files for target main2b
main2b_OBJECTS = \
"CMakeFiles/main2b.dir/main2b.cc.o"

# External object files for target main2b
main2b_EXTERNAL_OBJECTS =

main2b: CMakeFiles/main2b.dir/main2b.cc.o
main2b: CMakeFiles/main2b.dir/build.make
main2b: /Users/junchaowei/Desktop/Jedi_Temple/dealii/lib/libdeal_II.g.8.5.0-pre.dylib
main2b: /usr/lib/libbz2.dylib
main2b: /usr/local/lib/libtbb.dylib
main2b: /usr/lib/libz.dylib
main2b: /usr/local/lib/libarpack.dylib
main2b: /usr/local/lib/libhdf5_hl.dylib
main2b: /usr/local/lib/libhdf5.dylib
main2b: CMakeFiles/main2b.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main2b"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main2b.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main2b.dir/build: main2b

.PHONY : CMakeFiles/main2b.dir/build

CMakeFiles/main2b.dir/requires: CMakeFiles/main2b.dir/main2b.cc.o.requires

.PHONY : CMakeFiles/main2b.dir/requires

CMakeFiles/main2b.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main2b.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main2b.dir/clean

CMakeFiles/main2b.dir/depend:
	cd "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/CMakeFiles/main2b.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/main2b.dir/depend

