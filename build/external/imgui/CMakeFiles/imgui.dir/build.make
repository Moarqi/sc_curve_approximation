# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/markus/dev/sc_curve_approximation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/markus/dev/sc_curve_approximation/build

# Include any dependencies generated for this target.
include external/imgui/CMakeFiles/imgui.dir/depend.make

# Include the progress variables for this target.
include external/imgui/CMakeFiles/imgui.dir/progress.make

# Include the compile flags for this target's objects.
include external/imgui/CMakeFiles/imgui.dir/flags.make

external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o: external/imgui/CMakeFiles/imgui.dir/flags.make
external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o: ../external/imgui/imgui.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/markus/dev/sc_curve_approximation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui.cpp.o -c /home/markus/dev/sc_curve_approximation/external/imgui/imgui.cpp

external/imgui/CMakeFiles/imgui.dir/imgui.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui.cpp.i"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/markus/dev/sc_curve_approximation/external/imgui/imgui.cpp > CMakeFiles/imgui.dir/imgui.cpp.i

external/imgui/CMakeFiles/imgui.dir/imgui.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui.cpp.s"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/markus/dev/sc_curve_approximation/external/imgui/imgui.cpp -o CMakeFiles/imgui.dir/imgui.cpp.s

external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.requires:

.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.requires

external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.provides: external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.requires
	$(MAKE) -f external/imgui/CMakeFiles/imgui.dir/build.make external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.provides.build
.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.provides

external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.provides.build: external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o


external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o: external/imgui/CMakeFiles/imgui.dir/flags.make
external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o: ../external/imgui/imgui_demo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/markus/dev/sc_curve_approximation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui_demo.cpp.o -c /home/markus/dev/sc_curve_approximation/external/imgui/imgui_demo.cpp

external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui_demo.cpp.i"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/markus/dev/sc_curve_approximation/external/imgui/imgui_demo.cpp > CMakeFiles/imgui.dir/imgui_demo.cpp.i

external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui_demo.cpp.s"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/markus/dev/sc_curve_approximation/external/imgui/imgui_demo.cpp -o CMakeFiles/imgui.dir/imgui_demo.cpp.s

external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.requires:

.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.requires

external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.provides: external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.requires
	$(MAKE) -f external/imgui/CMakeFiles/imgui.dir/build.make external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.provides.build
.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.provides

external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.provides.build: external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o


external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o: external/imgui/CMakeFiles/imgui.dir/flags.make
external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o: ../external/imgui/imgui_draw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/markus/dev/sc_curve_approximation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui_draw.cpp.o -c /home/markus/dev/sc_curve_approximation/external/imgui/imgui_draw.cpp

external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui_draw.cpp.i"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/markus/dev/sc_curve_approximation/external/imgui/imgui_draw.cpp > CMakeFiles/imgui.dir/imgui_draw.cpp.i

external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui_draw.cpp.s"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/markus/dev/sc_curve_approximation/external/imgui/imgui_draw.cpp -o CMakeFiles/imgui.dir/imgui_draw.cpp.s

external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.requires:

.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.requires

external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.provides: external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.requires
	$(MAKE) -f external/imgui/CMakeFiles/imgui.dir/build.make external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.provides.build
.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.provides

external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.provides.build: external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o


external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o: external/imgui/CMakeFiles/imgui.dir/flags.make
external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o: ../external/imgui/imgui_glfw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/markus/dev/sc_curve_approximation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui_glfw.cpp.o -c /home/markus/dev/sc_curve_approximation/external/imgui/imgui_glfw.cpp

external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui_glfw.cpp.i"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/markus/dev/sc_curve_approximation/external/imgui/imgui_glfw.cpp > CMakeFiles/imgui.dir/imgui_glfw.cpp.i

external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui_glfw.cpp.s"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/markus/dev/sc_curve_approximation/external/imgui/imgui_glfw.cpp -o CMakeFiles/imgui.dir/imgui_glfw.cpp.s

external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.requires:

.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.requires

external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.provides: external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.requires
	$(MAKE) -f external/imgui/CMakeFiles/imgui.dir/build.make external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.provides.build
.PHONY : external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.provides

external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.provides.build: external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o


# Object files for target imgui
imgui_OBJECTS = \
"CMakeFiles/imgui.dir/imgui.cpp.o" \
"CMakeFiles/imgui.dir/imgui_demo.cpp.o" \
"CMakeFiles/imgui.dir/imgui_draw.cpp.o" \
"CMakeFiles/imgui.dir/imgui_glfw.cpp.o"

# External object files for target imgui
imgui_EXTERNAL_OBJECTS =

libimgui.a: external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o
libimgui.a: external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o
libimgui.a: external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o
libimgui.a: external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o
libimgui.a: external/imgui/CMakeFiles/imgui.dir/build.make
libimgui.a: external/imgui/CMakeFiles/imgui.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/markus/dev/sc_curve_approximation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library ../../libimgui.a"
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean_target.cmake
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imgui.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/imgui/CMakeFiles/imgui.dir/build: libimgui.a

.PHONY : external/imgui/CMakeFiles/imgui.dir/build

external/imgui/CMakeFiles/imgui.dir/requires: external/imgui/CMakeFiles/imgui.dir/imgui.cpp.o.requires
external/imgui/CMakeFiles/imgui.dir/requires: external/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o.requires
external/imgui/CMakeFiles/imgui.dir/requires: external/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o.requires
external/imgui/CMakeFiles/imgui.dir/requires: external/imgui/CMakeFiles/imgui.dir/imgui_glfw.cpp.o.requires

.PHONY : external/imgui/CMakeFiles/imgui.dir/requires

external/imgui/CMakeFiles/imgui.dir/clean:
	cd /home/markus/dev/sc_curve_approximation/build/external/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean.cmake
.PHONY : external/imgui/CMakeFiles/imgui.dir/clean

external/imgui/CMakeFiles/imgui.dir/depend:
	cd /home/markus/dev/sc_curve_approximation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/markus/dev/sc_curve_approximation /home/markus/dev/sc_curve_approximation/external/imgui /home/markus/dev/sc_curve_approximation/build /home/markus/dev/sc_curve_approximation/build/external/imgui /home/markus/dev/sc_curve_approximation/build/external/imgui/CMakeFiles/imgui.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/imgui/CMakeFiles/imgui.dir/depend

