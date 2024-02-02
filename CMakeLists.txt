cmake_minimum_required (VERSION 3.14) # Recommend a newer version for FetchContent
project (XtalComp)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library (XtalComp STATIC xtalcomp.cpp xctransform.cpp xcvector.h xcmatrix.h
  stablecomparison.h)

add_library (XtalCompShared SHARED xtalcomp.cpp xctransform.cpp xcvector.h xcmatrix.h
stablecomparison.h)
set_target_properties(XtalCompShared PROPERTIES OUTPUT_NAME XtalComp)

include(FetchContent)
FetchContent_Declare(Spglib
    GIT_REPOSITORY https://github.com/spglib/spglib
    GIT_TAG v2.3.0
)
FetchContent_MakeAvailable(Spglib)

target_include_directories(XtalComp PRIVATE ${spglib_SOURCE_DIR}/src)
target_include_directories(XtalCompShared PRIVATE ${spglib_SOURCE_DIR}/src)

target_link_libraries(XtalComp PRIVATE Spglib::symspg)
target_link_libraries(XtalCompShared PRIVATE Spglib::symspg)

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
