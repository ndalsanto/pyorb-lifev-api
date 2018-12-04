add_executable(loadlifev loadlifev.cpp)
target_link_libraries(loadlifev pyorb-lifev-api ${LINKED_LIBRARIES} )

add_executable(use_lifev_wrappers use_lifev_wrappers.cpp)
target_link_libraries(use_lifev_wrappers pyorb-lifev-api ${LINKED_LIBRARIES} )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/4cube1.mesh
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/SolverParamList.xml
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
