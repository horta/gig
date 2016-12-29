include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_gig test_gig.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_gig LINK_PUBLIC gig)
else()
  target_link_libraries (test_gig LINK_PUBLIC gig m)
endif()
