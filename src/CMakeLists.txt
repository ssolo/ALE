 SET(exODT_SRCS
  ALE.h ALE.cpp
  ALE_util.h ALE_util.cpp
  exODT.h exODT.cpp
  model_scaled.cpp
  traceback_scaled.cpp
  ALE_util.h ALE_util.cpp
  sample_scaled.cpp
  undated.cpp
  pairHasher.h
  fractionMissing.h fractionMissing.cpp
  )

add_library(exODT ${exODT_SRCS})
target_link_libraries(exODT 
  ${BPP_LIBRARIES}
)

add_executable(ALEml ALEml.cpp)
target_link_libraries(ALEml exODT
  ${BPP_LIBRARIES}
)

add_executable(ALEml_undated ALEml_undated.cpp)
target_link_libraries(ALEml_undated exODT
  ${BPP_LIBRARIES}
)

add_executable(ALEmcmc_undated ALEmcmc_undated.cpp)
target_link_libraries(ALEmcmc_undated exODT
  ${BPP_LIBRARIES}
)

add_executable(ALEevaluate_undated ALEevaluate_undated.cpp)
target_link_libraries(ALEevaluate_undated exODT
  ${BPP_LIBRARIES}
)


add_executable(ALEsample ALEsample.cpp)
target_link_libraries(ALEsample exODT
  ${BPP_LIBRARIES}
)

add_executable(ALEobserve ALEobserve.cpp)
target_link_libraries(ALEobserve exODT
  ${BPP_LIBRARIES}
)

add_executable(ALEadd ALEadd.cpp)
target_link_libraries(ALEadd exODT
  ${BPP_LIBRARIES}
)

add_executable(ALEcount ALEcount.cpp)
target_link_libraries(ALEcount exODT
  ${BPP_LIBRARIES}
)

add_executable(CCPscore CCPscore.cpp)
target_link_libraries(CCPscore exODT
  ${BPP_LIBRARIES}
)

add_executable(ls_leaves ls_leaves.cpp)
target_link_libraries(ls_leaves exODT
  ${BPP_LIBRARIES}
)

#add_executable(parse_maxtic parse_maxtic.cpp)
#target_link_libraries(parse_maxtic exODT
# ${BPP_LIBRARIES}
#)

#add_executable(bls bls.cpp)
#target_link_libraries(bls exODT
# ${BPP_LIBRARIES}
#)


add_executable(times_undated times_undated.cpp)
target_link_libraries(times_undated exODT
  ${BPP_LIBRARIES}
)

add_executable(computeALEcomplexity computeALEcomplexity.cpp)
target_link_libraries(computeALEcomplexity exODT
  ${BPP_LIBRARIES}
)

# needs OMP
# add_executable(mlsampler exODT_sim.h mlsampler.cpp)
# target_link_libraries(mlsampler exODT
#   ${Boost_MPI_LIBRARY}
#   ${MPI_LIBRARIES}
#   ${BPP_LIBRARIES}
# )

add_executable(mlresampler_undated mlresampler_undated.cpp)
target_link_libraries(mlresampler_undated exODT
  ${BPP_LIBRARIES}
)


# do not link with recent BPP
# add_executable(ALE_tutorial ALE_tutorial.cpp)
# target_link_libraries(ALE_tutorial exODT
#   ${Boost_MPI_LIBRARY}
#   ${MPI_LIBRARIES}
#   ${BPP_LIBRARIES}
# )


add_executable(simulation simulation.cpp exODT_sim.cpp)
target_link_libraries(simulation exODT
  ${BPP_LIBRARIES}
)

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
   
add_executable(mpi_ml_undated mpi_tree.cpp mpi_ml_undated.cpp)
target_link_libraries(mpi_ml_undated exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)   

add_executable(mpi_ml-bw_undated mpi_tree.cpp mpi_ml-bw_undated.cpp)
target_link_libraries(mpi_ml-bw_undated exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

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
