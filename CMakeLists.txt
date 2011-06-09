cmake_minimum_required (VERSION 2.6)
project (XtalComp)

add_library (XtalComp xtalcomp.cpp xctransform.cpp)

add_executable (test test.cpp)
target_link_libraries (test XtalComp)

option( BUILD_CGI
 	"Whether to compile the CGI handler as well as the XtalComp code."
 	OFF )
if(BUILD_CGI)
  add_executable (xtalcomp.cgi cgi/xtalcomp-cgi.cpp)
  target_link_libraries (xtalcomp.cgi XtalComp)
endif(BUILD_CGI)
