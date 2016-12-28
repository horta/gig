add_library(gig_objs OBJECT gig.cpp)
set_target_properties(gig_objs
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
