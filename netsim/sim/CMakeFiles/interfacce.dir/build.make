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
include CMakeFiles/interfacce.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/interfacce.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/interfacce.dir/flags.make

CMakeFiles/interfacce.dir/distrib.c.obj: CMakeFiles/interfacce.dir/flags.make
CMakeFiles/interfacce.dir/distrib.c.obj: distrib.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/interfacce.dir/distrib.c.obj"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\interfacce.dir\distrib.c.obj   -c "D:\Progetti R\netsim_hg\sim\distrib.c"

CMakeFiles/interfacce.dir/distrib.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/interfacce.dir/distrib.c.i"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\distrib.c" > CMakeFiles\interfacce.dir\distrib.c.i

CMakeFiles/interfacce.dir/distrib.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/interfacce.dir/distrib.c.s"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\distrib.c" -o CMakeFiles\interfacce.dir\distrib.c.s

CMakeFiles/interfacce.dir/distrib.c.obj.requires:
.PHONY : CMakeFiles/interfacce.dir/distrib.c.obj.requires

CMakeFiles/interfacce.dir/distrib.c.obj.provides: CMakeFiles/interfacce.dir/distrib.c.obj.requires
	$(MAKE) -f CMakeFiles\interfacce.dir\build.make CMakeFiles/interfacce.dir/distrib.c.obj.provides.build
.PHONY : CMakeFiles/interfacce.dir/distrib.c.obj.provides

CMakeFiles/interfacce.dir/distrib.c.obj.provides.build: CMakeFiles/interfacce.dir/distrib.c.obj
.PHONY : CMakeFiles/interfacce.dir/distrib.c.obj.provides.build

CMakeFiles/interfacce.dir/interfacce.c.obj: CMakeFiles/interfacce.dir/flags.make
CMakeFiles/interfacce.dir/interfacce.c.obj: interfacce.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/interfacce.dir/interfacce.c.obj"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\interfacce.dir\interfacce.c.obj   -c "D:\Progetti R\netsim_hg\sim\interfacce.c"

CMakeFiles/interfacce.dir/interfacce.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/interfacce.dir/interfacce.c.i"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\interfacce.c" > CMakeFiles\interfacce.dir\interfacce.c.i

CMakeFiles/interfacce.dir/interfacce.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/interfacce.dir/interfacce.c.s"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\interfacce.c" -o CMakeFiles\interfacce.dir\interfacce.c.s

CMakeFiles/interfacce.dir/interfacce.c.obj.requires:
.PHONY : CMakeFiles/interfacce.dir/interfacce.c.obj.requires

CMakeFiles/interfacce.dir/interfacce.c.obj.provides: CMakeFiles/interfacce.dir/interfacce.c.obj.requires
	$(MAKE) -f CMakeFiles\interfacce.dir\build.make CMakeFiles/interfacce.dir/interfacce.c.obj.provides.build
.PHONY : CMakeFiles/interfacce.dir/interfacce.c.obj.provides

CMakeFiles/interfacce.dir/interfacce.c.obj.provides.build: CMakeFiles/interfacce.dir/interfacce.c.obj
.PHONY : CMakeFiles/interfacce.dir/interfacce.c.obj.provides.build

CMakeFiles/interfacce.dir/globali.c.obj: CMakeFiles/interfacce.dir/flags.make
CMakeFiles/interfacce.dir/globali.c.obj: globali.c
	$(CMAKE_COMMAND) -E cmake_progress_report "D:\Progetti R\netsim_hg\sim\CMakeFiles" $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/interfacce.dir/globali.c.obj"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles\interfacce.dir\globali.c.obj   -c "D:\Progetti R\netsim_hg\sim\globali.c"

CMakeFiles/interfacce.dir/globali.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/interfacce.dir/globali.c.i"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -E "D:\Progetti R\netsim_hg\sim\globali.c" > CMakeFiles\interfacce.dir\globali.c.i

CMakeFiles/interfacce.dir/globali.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/interfacce.dir/globali.c.s"
	P:\Linguaggi\MinGW64\bin\gcc.exe  $(C_DEFINES) $(C_FLAGS) -S "D:\Progetti R\netsim_hg\sim\globali.c" -o CMakeFiles\interfacce.dir\globali.c.s

CMakeFiles/interfacce.dir/globali.c.obj.requires:
.PHONY : CMakeFiles/interfacce.dir/globali.c.obj.requires

CMakeFiles/interfacce.dir/globali.c.obj.provides: CMakeFiles/interfacce.dir/globali.c.obj.requires
	$(MAKE) -f CMakeFiles\interfacce.dir\build.make CMakeFiles/interfacce.dir/globali.c.obj.provides.build
.PHONY : CMakeFiles/interfacce.dir/globali.c.obj.provides

CMakeFiles/interfacce.dir/globali.c.obj.provides.build: CMakeFiles/interfacce.dir/globali.c.obj
.PHONY : CMakeFiles/interfacce.dir/globali.c.obj.provides.build

# Object files for target interfacce
interfacce_OBJECTS = \
"CMakeFiles/interfacce.dir/distrib.c.obj" \
"CMakeFiles/interfacce.dir/interfacce.c.obj" \
"CMakeFiles/interfacce.dir/globali.c.obj"

# External object files for target interfacce
interfacce_EXTERNAL_OBJECTS =

interfacce.dll: CMakeFiles/interfacce.dir/distrib.c.obj
interfacce.dll: CMakeFiles/interfacce.dir/interfacce.c.obj
interfacce.dll: CMakeFiles/interfacce.dir/globali.c.obj
interfacce.dll: P:/Linguaggi/MinGW64/lib32/libglib-2.0.dll.a
interfacce.dll: P:/Linguaggi/MinGW64/lib32/libgsl.dll.a
interfacce.dll: P:/Linguaggi/MinGW64/lib32/libgslcblas.dll.a
interfacce.dll: CMakeFiles/interfacce.dir/build.make
interfacce.dll: CMakeFiles/interfacce.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C shared library interfacce.dll"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\interfacce.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/interfacce.dir/build: interfacce.dll
.PHONY : CMakeFiles/interfacce.dir/build

CMakeFiles/interfacce.dir/requires: CMakeFiles/interfacce.dir/distrib.c.obj.requires
CMakeFiles/interfacce.dir/requires: CMakeFiles/interfacce.dir/interfacce.c.obj.requires
CMakeFiles/interfacce.dir/requires: CMakeFiles/interfacce.dir/globali.c.obj.requires
.PHONY : CMakeFiles/interfacce.dir/requires

CMakeFiles/interfacce.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\interfacce.dir\cmake_clean.cmake
.PHONY : CMakeFiles/interfacce.dir/clean

CMakeFiles/interfacce.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim" "D:\Progetti R\netsim_hg\sim\CMakeFiles\interfacce.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/interfacce.dir/depend

