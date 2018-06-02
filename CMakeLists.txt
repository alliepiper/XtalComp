cmake_minimum_required (VERSION 2.6)
project (XtalComp)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library (XtalComp STATIC xtalcomp.cpp xctransform.cpp xcvector.h xcmatrix.h
  stablecomparison.h)

add_library (XtalCompShared SHARED xtalcomp.cpp xctransform.cpp xcvector.h xcmatrix.h
stablecomparison.h)
set_target_properties(XtalCompShared PROPERTIES OUTPUT_NAME XtalComp)

target_link_libraries(XtalComp spglib)
target_link_libraries(XtalCompShared spglib)

add_subdirectory(spglib)

option( BUILD_CGI
       "Whether to compile the CGI handler as well as the XtalComp code."
       OFF )
if(BUILD_CGI)
  add_executable (xtalcomp.cgi cgi/xtalcomp-cgi.cpp)
  target_link_libraries (xtalcomp.cgi XtalComp)
endif(BUILD_CGI)

add_executable (test test.cpp)
target_link_libraries (test XtalComp)

install(TARGETS XtalComp XtalCompShared LIBRARY DESTINATION "lib" ARCHIVE DESTINATION "lib")
install(FILES "xcmatrix.h" "xctransform.h" "xcvector.h" "xtalcomp.h" DESTINATION "include/xtalcomp")
