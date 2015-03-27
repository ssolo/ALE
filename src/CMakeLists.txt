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