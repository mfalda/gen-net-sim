cmake_minimum_required (VERSION 2.6)
project (pars)

#set(CMAKE_C_COMPILER "R CMD SHLIB")

# parametri aggiuntivi
add_definitions("-m32 -std=c99 -fmessage-length=0 -D_FORTIFY_SOURCE=2 -Wconversion -Wfloat-equal")

# librerie aggiuntive
find_library(M_LIB m)
find_package(GLIB2)

if(APPLE)
	set(CARBON "-framework Carbon -liconv -lintl")
else()
	set(CARBON "")
endif()

# piattaforma
if(UNIX)
  if(APPLE)
    set(NOME "pars_OSX")
  else(APPLE)
    set(NOME "pars_linux")
  endif(APPLE)
else(UNIX)
  if(WIN32)
    set(NOME "pars_win32")
    add_definitions("-m32")
    set(CMAKE_MODULE_LINKER_FLAGS -m32)
  else(WIN32)
    error("OS unknown")
  endif(WIN32)
endif(UNIX)

message("GLIB2 libraries: " ${GLIB2_LIBRARIES})

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")
include_directories(${GLIB2_INCLUDE_DIR})

# add the binary tree to the search path for library files
link_directories(".")

# add the libraries
add_library(${NOME} pars.c)
target_link_libraries(${NOME} ${GLIB2_LIBRARIES})
set_target_properties(${NOME} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")

# add the target
add_executable(prova prova.c)
target_link_libraries(prova ${NOME} ${GLIB2_LIBRARIES} ${CARBON})
set_target_properties(prova PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
