cmake_minimum_required (VERSION 2.6)
project (netsim)

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
   add_definitions("-std=c99 -fPIC -Wall -Wextra -Wformat=2 -Wswitch-default -Wcast-align -Wpointer-arith -Wbad-function-cast -Wstrict-prototypes -Winline -Wundef -Wnested-externs -Wcast-qual -Wshadow -Wwrite-strings -Wconversion -Wunreachable-code -Wstrict-aliasing=2 -ffloat-store -fno-common -fstrict-aliasing -pedantic -O0 -ggdb3 -pg")
else()
	set(CMAKE_BUILD_TYPE Release)
   add_definitions("-pedantic -O3 -DNDEBUG -funroll-loops")
endif()
message("build type system is ${CMAKE_BUILD_TYPE}")

# piattaforma
if(UNIX)
  if(APPLE)
    set(OS "OSX")    
  else(APPLE)
    set(OS "linux")
  endif(APPLE)
  set(muPARS "muparser")
else(UNIX)
  if(WIN32)
    set(OS "win32")
    add_definitions("-m32")
    set(CMAKE_SHARED_LINKER_FLAGS "-m32 -shared -LP:/Linguaggi/R-2.13.1/bin/i386 -lR")
    set(muPARS "muparser32")
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

set(PARS "pars_${OS}")

# librerie aggiuntive
#if(WIN32)
  #set(M_LIB "G:/MinGW/lib/libm.a")
#else()
 #find_library(M_LIB m)
#endif()
 find_package(GLIB2)
 find_package(GSL)
 find_package(R)

message("R libraries are ${R_LIBRARY_BASE}")

if(APPLE)
	set(CARBON "-framework Carbon -liconv -lintl")
	set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
else()
	set(CARBON "")
endif()

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")
include_directories(${GLIB2_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${R_INCLUDE_DIR})

message("GLIB2 libraries are ${GLIB2_LIBRARIES}")

# add the binary tree to the search path for library files
link_directories(".")

# specifica i file componenti

SET(OBJS_COMM distrib.c)

SET(OBJS_LSODA lsoda_oscillatore.c)

SET(OBJS_INTRF interfacce.c globali.c)

#SET(OBJS_TRPOL2 trpol2_R.c)

SET(OBJS_HIST hist.c)

SET(OBJS_CALC calc.c ext_parser.c)

SET(OBJS_SN simulatenet.c globali.c createRules.c create_logicRule.c sample.c connectivity_modular.c module1.c module2.c module3.c mod1.c mod2.c mod3.c triangola.c probmod.c scoremodular.c cluster_coeff.c assign_nodes.c hist.c connectivity_geometric.c connectivity_random.c sampleB.c connectivity_scalefree.c score_sf.c createNEG.c dinamica.c target.c boole_result.c next_op.c write_table.c lsoda.c read_table.c ext_parser.c)

SET(OBJS_SP simulateprofiles.c globali.c createRules.c create_logicRule.c sample.c dinamica.c target.c boole_result.c next_op.c write_table.c lsoda.c read_table.c ext_parser.c)

SET(OBJS_MANDEL mandel.c)

SET(OBJS_MCALC mandel_calc.c)

SET(OBJS_MCALC1 mandel_calc1.c)

# add the libraries
add_library(interfacce SHARED ${OBJS_COMM} ${OBJS_INTRF})
add_library(lsoda_oscillatore SHARED ${OBJS_COMM} ${OBJS_LSODA})
#add_library(trpol2 SHARED ${OBJS_COMM} ${OBJS_TRPOL2})
add_library(simulatenet SHARED ${OBJS_COMM} ${OBJS_SN})
add_library(simulateprofiles SHARED ${OBJS_COMM} ${OBJS_SP})
add_library(hist SHARED ${OBJS_COMM} ${OBJS_HIST})
add_library(calc SHARED ${OBJS_COMM} ${OBJS_CALC})
add_library(mandel SHARED ${OBJS_COMM} ${OBJS_MANDEL})
add_library(mandel_calc SHARED ${OBJS_COMM} ${OBJS_MCALC})
#add_library(mandel_calc1 SHARED ${OBJS_COMM} ${OBJS_MCALC1})

set_target_properties(interfacce PROPERTIES PREFIX "")
set_target_properties(lsoda_oscillatore PROPERTIES PREFIX "")
#set_target_properties(trpol2 PROPERTIES PREFIX "")
set_target_properties(simulatenet PROPERTIES PREFIX "")
set_target_properties(simulateprofiles PROPERTIES PREFIX "")
set_target_properties(hist PROPERTIES PREFIX "")
set_target_properties(calc PROPERTIES PREFIX "")
set_target_properties(mandel PROPERTIES PREFIX "")
set_target_properties(mandel_calc PROPERTIES PREFIX "")
#set_target_properties(mandel_calc1 PROPERTIES PREFIX "")

  target_link_libraries(calc ${muPARS} ${R_LIBRARIES} ${R2C} ${GLIB2_LIBRARIES} ${CARBON})
  target_link_libraries(simulatenet ${muPARS} ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${GSL_LIBRARIES} ${CARBON})
  target_link_libraries(simulateprofiles ${muPARS} ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${GSL_LIBRARIES} ${CARBON})
  target_link_libraries(hist ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
  #target_link_libraries(trpol2 ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
  target_link_libraries(lsoda_oscillatore ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${GSL_LIBRARIES} ${CARBON})
  target_link_libraries(interfacce ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${GSL_LIBRARIES} ${CARBON})
  target_link_libraries(mandel ${R2C} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
  target_link_libraries(mandel_calc ${R2C} ${muPARS} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
  #target_link_libraries(mandel_calc1 ${R2C} ${PARS} ${GLIB2_LIBRARIES} ${R_LIBRARIES} ${CARBON})
if(NOT FINALE)
  set_target_properties(interfacce PROPERTIES COMPILE_FLAGS "-pg" LINK_FLAGS "-pg")
endif()
