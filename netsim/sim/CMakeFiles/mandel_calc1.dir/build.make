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
include CMakeFiles/mandel_calc1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mandel_calc1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mandel_calc1.dir/flags.make

CMakeFiles/mandel_calc1.dir/distrib.c.obj: CMakeFiles/mandel_calc1.dir/flags.make
CMakeFiles/mandel_calc1.dir/distrib.c.obj: distrib.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mandel_calc1.dir/distrib.c.obj"
	P:\Linguaggi\C\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\mandel_calc1.dir\distrib.c.obj   -c "D:\Progetti R\netsim_hg\sim\distrib.c"

CMakeFiles/mandel_calc1.dir/distrib.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mandel_calc1.dir/distrib.c.i"
	P:\Linguaggi\C\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\distrib.c" > CMakeFiles\mandel_calc1.dir\distrib.c.i

CMakeFiles/mandel_calc1.dir/distrib.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mandel_calc1.dir/distrib.c.s"
	P:\Linguaggi\C\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\distrib.c" -o CMakeFiles\mandel_calc1.dir\distrib.c.s

CMakeFiles/mandel_calc1.dir/distrib.c.obj.requires:
.PHONY : CMakeFiles/mandel_calc1.dir/distrib.c.obj.requires

CMakeFiles/mandel_calc1.dir/distrib.c.obj.provides: CMakeFiles/mandel_calc1.dir/distrib.c.obj.requires
	$(MAKE) -f CMakeFiles\mandel_calc1.dir\build.make CMakeFiles/mandel_calc1.dir/distrib.c.obj.provides.build
.PHONY : CMakeFiles/mandel_calc1.dir/distrib.c.obj.provides

CMakeFiles/mandel_calc1.dir/distrib.c.obj.provides.build: CMakeFiles/mandel_calc1.dir/distrib.c.obj
.PHONY : CMakeFiles/mandel_calc1.dir/distrib.c.obj.provides.build

CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj: CMakeFiles/mandel_calc1.dir/flags.make
CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj: mandel_calc1.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj"
	P:\Linguaggi\C\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\mandel_calc1.dir\mandel_calc1.c.obj   -c "D:\Progetti R\netsim_hg\sim\mandel_calc1.c"

CMakeFiles/mandel_calc1.dir/mandel_calc1.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mandel_calc1.dir/mandel_calc1.c.i"
	P:\Linguaggi\C\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\mandel_calc1.c" > CMakeFiles\mandel_calc1.dir\mandel_calc1.c.i

CMakeFiles/mandel_calc1.dir/mandel_calc1.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mandel_calc1.dir/mandel_calc1.c.s"
	P:\Linguaggi\C\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\mandel_calc1.c" -o CMakeFiles\mandel_calc1.dir\mandel_calc1.c.s

CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.requires:
.PHONY : CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.requires

CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.provides: CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.requires
	$(MAKE) -f CMakeFiles\mandel_calc1.dir\build.make CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.provides.build
.PHONY : CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.provides

CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.provides.build: CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj
.PHONY : CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.provides.build

# Object files for target mandel_calc1
mandel_calc1_OBJECTS = \
"CMakeFiles/mandel_calc1.dir/distrib.c.obj" \
"CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj"

# External object files for target mandel_calc1
mandel_calc1_EXTERNAL_OBJECTS =

mandel_calc1.dll: CMakeFiles/mandel_calc1.dir/distrib.c.obj
mandel_calc1.dll: CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj
mandel_calc1.dll: P:/Linguaggi/C/glib2.18/lib/libglib-2.0.dll.a
mandel_calc1.dll: P:/Linguaggi/R-2.13.1/bin/i386/R.dll
mandel_calc1.dll: CMakeFiles/mandel_calc1.dir/build.make
mandel_calc1.dll: CMakeFiles/mandel_calc1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C shared library mandel_calc1.dll"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\mandel_calc1.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mandel_calc1.dir/build: mandel_calc1.dll
.PHONY : CMakeFiles/mandel_calc1.dir/build

CMakeFiles/mandel_calc1.dir/requires: CMakeFiles/mandel_calc1.dir/distrib.c.obj.requires
CMakeFiles/mandel_calc1.dir/requires: CMakeFiles/mandel_calc1.dir/mandel_calc1.c.obj.requires
.PHONY : CMakeFiles/mandel_calc1.dir/requires

CMakeFiles/mandel_calc1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\mandel_calc1.dir\cmake_clean.cmake
.PHONY : CMakeFiles/mandel_calc1.dir/clean

CMakeFiles/mandel_calc1.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim\CMakeFiles\mandel_calc1.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/mandel_calc1.dir/depend

