# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = P:\Linguaggi\CMake-2.8\bin\cmake.exe

# The command to remove a file.
RM = P:\Linguaggi\CMake-2.8\bin\cmake.exe -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = P:\Linguaggi\CMake-2.8\bin\cmake-gui.exe

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "D:\Progetti R\netsim_hg\sim"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "D:\Progetti R\netsim_hg\sim"

# Include any dependencies generated for this target.
include CMakeFiles/mandel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mandel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mandel.dir/flags.make

CMakeFiles/mandel.dir/distrib.c.obj: CMakeFiles/mandel.dir/flags.make
CMakeFiles/mandel.dir/distrib.c.obj: distrib.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mandel.dir/distrib.c.obj"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\mandel.dir\distrib.c.obj   -c "D:\Progetti R\netsim_hg\sim\distrib.c"

CMakeFiles/mandel.dir/distrib.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mandel.dir/distrib.c.i"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\distrib.c" > CMakeFiles\mandel.dir\distrib.c.i

CMakeFiles/mandel.dir/distrib.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mandel.dir/distrib.c.s"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\distrib.c" -o CMakeFiles\mandel.dir\distrib.c.s

CMakeFiles/mandel.dir/distrib.c.obj.requires:
.PHONY : CMakeFiles/mandel.dir/distrib.c.obj.requires

CMakeFiles/mandel.dir/distrib.c.obj.provides: CMakeFiles/mandel.dir/distrib.c.obj.requires
	$(MAKE) -f CMakeFiles\mandel.dir\build.make CMakeFiles/mandel.dir/distrib.c.obj.provides.build
.PHONY : CMakeFiles/mandel.dir/distrib.c.obj.provides

CMakeFiles/mandel.dir/distrib.c.obj.provides.build: CMakeFiles/mandel.dir/distrib.c.obj
.PHONY : CMakeFiles/mandel.dir/distrib.c.obj.provides.build

CMakeFiles/mandel.dir/mandel.c.obj: CMakeFiles/mandel.dir/flags.make
CMakeFiles/mandel.dir/mandel.c.obj: mandel.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mandel.dir/mandel.c.obj"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\mandel.dir\mandel.c.obj   -c "D:\Progetti R\netsim_hg\sim\mandel.c"

CMakeFiles/mandel.dir/mandel.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mandel.dir/mandel.c.i"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\mandel.c" > CMakeFiles\mandel.dir\mandel.c.i

CMakeFiles/mandel.dir/mandel.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mandel.dir/mandel.c.s"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\mandel.c" -o CMakeFiles\mandel.dir\mandel.c.s

CMakeFiles/mandel.dir/mandel.c.obj.requires:
.PHONY : CMakeFiles/mandel.dir/mandel.c.obj.requires

CMakeFiles/mandel.dir/mandel.c.obj.provides: CMakeFiles/mandel.dir/mandel.c.obj.requires
	$(MAKE) -f CMakeFiles\mandel.dir\build.make CMakeFiles/mandel.dir/mandel.c.obj.provides.build
.PHONY : CMakeFiles/mandel.dir/mandel.c.obj.provides

CMakeFiles/mandel.dir/mandel.c.obj.provides.build: CMakeFiles/mandel.dir/mandel.c.obj
.PHONY : CMakeFiles/mandel.dir/mandel.c.obj.provides.build

# Object files for target mandel
mandel_OBJECTS = \
"CMakeFiles/mandel.dir/distrib.c.obj" \
"CMakeFiles/mandel.dir/mandel.c.obj"

# External object files for target mandel
mandel_EXTERNAL_OBJECTS =

mandel.dll: CMakeFiles/mandel.dir/distrib.c.obj
mandel.dll: CMakeFiles/mandel.dir/mandel.c.obj
mandel.dll: P:/Linguaggi/MinGW64/lib32/libglib-2.0.dll.a
mandel.dll: CMakeFiles/mandel.dir/build.make
mandel.dll: CMakeFiles/mandel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C shared library mandel.dll"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\mandel.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mandel.dir/build: mandel.dll
.PHONY : CMakeFiles/mandel.dir/build

CMakeFiles/mandel.dir/requires: CMakeFiles/mandel.dir/distrib.c.obj.requires
CMakeFiles/mandel.dir/requires: CMakeFiles/mandel.dir/mandel.c.obj.requires
.PHONY : CMakeFiles/mandel.dir/requires

CMakeFiles/mandel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\mandel.dir\cmake_clean.cmake
.PHONY : CMakeFiles/mandel.dir/clean

CMakeFiles/mandel.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim\CMakeFiles\mandel.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/mandel.dir/depend

