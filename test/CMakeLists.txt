add_executable(loadlifev loadlifev.cpp)
target_link_libraries(loadlifev ${PROJECT_BINARY_DIR}/libtessst.so)
# target_link_libraries(loadlifev ${TPL_LifeV_LIBRARY_DIRS}/lib*.so)
add_test(loadlifev_test loadlifev)
