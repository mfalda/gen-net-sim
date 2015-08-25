cmake_minimum_required (VERSION 2.6)
project (hmm)

#set(CMAKE_C_COMPILER "R CMD SHLIB")

set(DEBG OFF)
set(FINALE ON)
set(DET OFF)

if(NOT DEBG AND NOT FINALE AND NOT DET)
	message(FATAL_ERROR "Definire almeno un tipo di versione")
endif()

# debugging (Release o Debug)
if(NOT FINALE)
	set(CMAKE_BUILD_TYPE Debug)
   add_definitions("-std=c99 -fPIC -Wall -Wextra -Wformat=2 -Wswitch-default -Wcast-align -Wpointer-arith -Wbad-function-cast -Wstrict-prototypes -Winline -Wundef -Wnested-externs -Wcast-qual -Wshadow -Wwrite-strings -Wconversion -Wunreachable-code -Wstrict-aliasing=2 -pedantic -O0 -ggdb3 -g -fPIC")
else()
	set(CMAKE_BUILD_TYPE Release)
   add_definitions("-pedantic -O3 -DNDEBUG -std=gnu99 -fPIC -g")
endif()
message("build type system is ${CMAKE_BUILD_TYPE}")

# piattaforma
if(UNIX)
  if(APPLE)
    set(OS "OSX")
  else(APPLE)
    set(OS "linux")
  endif(APPLE)
else(UNIX)
  if(WIN32)
    set(OS "win32")
  else(WIN32)
    error("OS unknown")
  endif(WIN32)
endif(UNIX)

if(NOT FINALE)
  message("Versione di debug")
  if (DET)
	  	message("det ${DET}")
		add_definitions("-DDET")
		set(R2C "r2c_${OS}_moduli")
	else()
  		set(R2C "r2c_${OS}_debug")
  endif()
else()
  message("Versione di release")
  set(R2C "r2c_${OS}")
endif()

# librerie aggiuntive
if (WIN32)
  set(M_LIB "")
else()
  find_library(M_LIB m)
endif()
find_package(GLIB2)
find_package(R)

if(APPLE)
	set(CARBON "-framework Carbon -liconv -lintl")
	set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
else()
	set(CARBON "")
endif()

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")
include_directories(${GLIB2_INCLUDE_DIR})
include_directories(${R_INCLUDE_DIR})

message("GLIB2 libraries are ${GLIB2_LIBRARIES}")

# add the binary tree to the search path for library files
link_directories(".")

# specifica i file componenti

SET(OBJS_COMM distrib.c)

SET(OBJS_create_mod globali.c  create_mod.c cluster_coeff.c hubs.c read_table.c pathlength.c)

SET(OBJS_HMM HMM.c globali.c hist.c scoremodular.c sampleB.c sample.c write_table.c probmod.c cluster_coeff.c connetti_scalefree.c assign_nodes.c check_conn.c)

SET(OBJS_HMM_und HMM_und.c globali.c hist.c scoremodular.c sampleB.c sample.c write_table.c probmod_und.c cluster_coeff.c connetti_scalefree.c assign_nodes_und.c check_conn.c)

# add the libraries
add_library(create_mod SHARED ${OBJS_COMM} ${OBJS_create_mod})
add_library(hmm SHARED ${OBJS_COMM} ${OBJS_HMM})
add_library(hmm_und SHARED ${OBJS_COMM} ${OBJS_HMM_und})

set_target_properties(create_mod PROPERTIES PREFIX "")
set_target_properties(hmm PROPERTIES PREFIX "")
set_target_properties(hmm_und PROPERTIES PREFIX "")

target_link_libraries(create_mod ${R2C} ${M_LIB} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
target_link_libraries(hmm ${R2C} ${M_LIB} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
target_link_libraries(hmm_und ${R2C} ${M_LIB} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
