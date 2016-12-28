cmake_minimum_required (VERSION 3.2.2)

project(Gig C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

file(STRINGS "VERSION" GIG_VERSION)

include_directories(.)

add_subdirectory(gig)

add_library(gig SHARED $<TARGET_OBJECTS:gig_objs>)

add_library(gig_static STATIC $<TARGET_OBJECTS:gig_objs>)

set_target_properties(gig
  PROPERTIES
    VERSION ${GIG_VERSION}
  PUBLIC_HEADER "gig/gig.h")

install(TARGETS gig
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
      PUBLIC_HEADER DESTINATION include/gig)

install(TARGETS gig_static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
      PUBLIC_HEADER DESTINATION include/gig)

# add_subdirectory(test)
#
# enable_testing()
# add_test(NAME TestLinking COMMAND test_linking)
