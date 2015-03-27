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
  )

add_library(exODT ${exODT_SRCS})


target_link_libraries(exODT 
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(ALEml ALEml.cpp)
target_link_libraries(ALEml exODT
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(ALEml_undated ALEml_undated.cpp)
target_link_libraries(ALEml_undated exODT
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(ALEsample ALEsample.cpp)
target_link_libraries(ALEsample exODT
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(ALEobserve ALEobserve.cpp)
target_link_libraries(ALEobserve exODT
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(ALEcount ALEcount.cpp)
target_link_libraries(ALEcount exODT
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)


add_executable(times times.cpp)
target_link_libraries(times exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(computeALEcomplexity computeALEcomplexity.cpp)
target_link_libraries(computeALEcomplexity exODT
  ${Boost_MPI_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

# needs OMP
# add_executable(mlsampler exODT_sim.h mlsampler.cpp)
# target_link_libraries(mlsampler exODT
#   ${Boost_MPI_LIBRARY}
#   ${MPI_LIBRARIES}
#   ${BPP_LIBRARIES}
# )

add_executable(wALE_ml_sample wALE_ml_sample.cpp)
target_link_libraries(wALE_ml_sample exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

# do not link with recent BPP
# add_executable(ALE_tutorial ALE_tutorial.cpp)
# target_link_libraries(ALE_tutorial exODT
#   ${Boost_MPI_LIBRARY}
#   ${MPI_LIBRARIES}
#   ${BPP_LIBRARIES}
# )

add_executable(simulateSpAndGeneTrees simulateSpAndGeneTrees.cpp exODT_sim.cpp)
target_link_libraries(simulateSpAndGeneTrees exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)

add_executable(mpi_ml mpi_tree.cpp mpi_ml.cpp)
target_link_libraries(mpi_ml exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)
   
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

add_executable(mpi_S_ml mpi_tree.h mpi_tree.cpp mpi_S_ml.cpp)
target_link_libraries(mpi_S_ml exODT
  ${Boost_MPI_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${MPI_LIBRARIES}
  ${BPP_LIBRARIES}
)
