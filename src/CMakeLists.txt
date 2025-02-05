SET(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow -std=c++11")
SET(CMAKE_C_FLAGS_DEBUG "-g -Wall -Weffc++ -Wshadow")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow -DNDEBUG -std=c++11")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Weffc++ -Wshadow -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -w -DNDEBUG -std=c++11")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -w -DNDEBUG")
add_compile_definitions(BOOST_TIMER_ENABLE_DEPRECATED)

SET(exODT_SRCS
  ALE.cpp
  ALE_util.cpp
  exODT.cpp
  model_scaled.cpp
  traceback_scaled.cpp
  ALE_util.cpp
  sample_scaled.cpp
  undated.cpp
  fractionMissing.cpp
)

add_library(exODT STATIC ${exODT_SRCS})
target_compile_features(exODT PUBLIC cxx_std_11)
target_link_libraries(exODT PUBLIC ${BPP_LIBRARIES} 
  Boost::algorithm Boost::mpi Boost::serialization Boost::timer)

IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF(NOT NO_VIRTUAL_COV)

IF(NO_VIRTUAL_COV)
  MESSAGE(STATUS "Covariant return with virtual inheritance disabled.")
  target_compile_definitions(exODT PUBLIC NO_VIRTUAL_COV=1)
ENDIF(NO_VIRTUAL_COV)

add_executable(ALEml ALEml.cpp)
target_link_libraries(ALEml PUBLIC exODT)

add_executable(ALEml_undated ALEml_undated.cpp)
target_link_libraries(ALEml_undated PUBLIC exODT)

add_executable(ALEmcmc_undated ALEmcmc_undated.cpp)
target_link_libraries(ALEmcmc_undated PUBLIC exODT)

add_executable(ALEevaluate_undated ALEevaluate_undated.cpp)
target_link_libraries(ALEevaluate_undated PUBLIC exODT)


add_executable(ALEsample ALEsample.cpp)
target_link_libraries(ALEsample PUBLIC exODT)

add_executable(ALEobserve ALEobserve.cpp)
target_link_libraries(ALEobserve PUBLIC exODT)

add_executable(ALEadd ALEadd.cpp)
target_link_libraries(ALEadd PUBLIC exODT)

add_executable(ALEcount ALEcount.cpp)
target_link_libraries(ALEcount PUBLIC exODT)

add_executable(CCPscore CCPscore.cpp)
target_link_libraries(CCPscore PUBLIC exODT)

add_executable(ls_leaves ls_leaves.cpp)
target_link_libraries(ls_leaves PUBLIC exODT)

#add_executable(parse_maxtic parse_maxtic.cpp)
#target_link_libraries(parse_maxtic exODT
# ${BPP_LIBRARIES}
#)

#add_executable(bls bls.cpp)
#target_link_libraries(bls exODT
# ${BPP_LIBRARIES}
#)

add_executable(ALEprune ALEprune.cpp)
target_link_libraries(ALEprune PUBLIC ${BPP_LIBRARIES} 
  Boost::algorithm Boost::timer Boost::math)

add_executable(times_undated times_undated.cpp)
target_link_libraries(times_undated PUBLIC exODT)

add_executable(computeALEcomplexity computeALEcomplexity.cpp)
target_link_libraries(computeALEcomplexity PUBLIC exODT)

# needs OMP
# add_executable(mlsampler exODT_sim.h mlsampler.cpp)
# target_link_libraries(mlsampler exODT
#   ${Boost_MPI_LIBRARY}
#   ${MPI_LIBRARIES}
#   ${BPP_LIBRARIES}
# )

add_executable(mlresampler_undated mlresampler_undated.cpp)
target_link_libraries(mlresampler_undated PUBLIC exODT)


# do not link with recent BPP
# add_executable(ALE_tutorial ALE_tutorial.cpp)
# target_link_libraries(ALE_tutorial exODT
#   ${Boost_MPI_LIBRARY}
#   ${MPI_LIBRARIES}
#   ${BPP_LIBRARIES}
# )


add_executable(simulation simulation.cpp exODT_sim.cpp)
target_link_libraries(simulation PUBLIC exODT)

#
#add_executable(wALE_ml_sample wALE_ml_sample.cpp)
#target_link_libraries(wALE_ml_sample exODT
# ${Boost_MPI_LIBRARY}
#  ${Boost_SERIALIZATION_LIBRARY}
#  ${BPP_LIBRARIES}
#)

#add_executable(wALE_ml_sample_undated wALE_ml_sample_undated.cpp)
#target_link_libraries(wALE_ml_sample_undated exODT
#  ${BPP_LIBRARIES}
#)

#add_executable(mpi_ml mpi_tree.cpp mpi_ml.cpp)
#target_link_libraries(mpi_ml exODT
#  ${Boost_MPI_LIBRARY}
#  ${Boost_SERIALIZATION_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)
   
#add_executable(mpi_ml_undated mpi_tree.cpp mpi_ml_undated.cpp)
#target_link_libraries(mpi_ml_undated exODT
#  ${Boost_MPI_LIBRARY}
#  ${Boost_SERIALIZATION_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)   

#add_executable(mpi_ml-bw_undated mpi_tree.cpp mpi_ml-bw_undated.cpp)
#target_link_libraries(mpi_ml-bw_undated exODT
#  ${Boost_MPI_LIBRARY}
#  ${Boost_SERIALIZATION_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)

#add_executable(mpi_S_ml mpi_tree.h mpi_tree.cpp mpi_S_ml.cpp)
#target_link_libraries(mpi_S_ml exODT
#  ${Boost_MPI_LIBRARY}
#  ${Boost_SERIALIZATION_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)

#add_executable(wol_host wol_host.cpp)
#target_link_libraries(wol_host exODT
#  ${Boost_MPI_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)

#add_executable(wol_paras wol_paras.cpp)
#target_link_libraries(wol_paras exODT
#  ${Boost_MPI_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)

#add_executable(wol_tree wol_tree.cpp)
#target_link_libraries(wol_tree exODT
#  ${Boost_MPI_LIBRARY}
#  ${MPI_LIBRARIES}
#  ${BPP_LIBRARIES}
#)

#add_executable(dates dates.cpp)
#target_link_libraries(dates
#  ${BPP_LIBRARIES}
#)
