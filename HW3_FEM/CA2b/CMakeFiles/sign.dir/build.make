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

# Utility rule file for sign.

# Include the progress variables for this target.
include CMakeFiles/sign.dir/progress.make

CMakeFiles/sign:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir="/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Digitally signing main2b"
	/usr/local/Cellar/cmake/3.9.0/bin/cmake -E echo '' && /usr/local/Cellar/cmake/3.9.0/bin/cmake -E echo '***************************************************************************' && /usr/local/Cellar/cmake/3.9.0/bin/cmake -E echo '** Error: No Mac OSX developer certificate specified **' && /usr/local/Cellar/cmake/3.9.0/bin/cmake -E echo '** Please reconfigure with -DOSX_CERTIFICATE_NAME="<...>" **' && /usr/local/Cellar/cmake/3.9.0/bin/cmake -E echo '***************************************************************************' && /usr/local/Cellar/cmake/3.9.0/bin/cmake -E echo ''

sign: CMakeFiles/sign
sign: CMakeFiles/sign.dir/build.make

.PHONY : sign

# Rule to build all files generated by this target.
CMakeFiles/sign.dir/build: sign

.PHONY : CMakeFiles/sign.dir/build

CMakeFiles/sign.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sign.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sign.dir/clean

CMakeFiles/sign.dir/depend:
	cd "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b" "/Users/junchaowei/Desktop/ad_CA2_template (1)/CA2_template/CA2b/CMakeFiles/sign.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/sign.dir/depend
