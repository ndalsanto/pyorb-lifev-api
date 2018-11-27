add_executable(loadlifev loadlifev.cpp)
# target_include_directories(loadlifev PUBLIC ${INCLUDED_HEADERS})
target_link_libraries(loadlifev ${PROJECT_BINARY_DIR}/libtessst.so ${LINKED_LIBRARIES} )
add_test(loadlifev_test loadlifev)
