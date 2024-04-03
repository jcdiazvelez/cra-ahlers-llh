cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -DBOOST_LOG_DYN_LINK")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

project(iter-lhreco-proj CXX C)

include(pretty)
include(config)
include(tools)
include(tooldef)
include(CheckIncludeFile)

FIND_PACKAGE(HEALPIX REQUIRED)
FIND_PACKAGE(CFITSIO REQUIRED)
FIND_PACKAGE(Boost 1.36.0 COMPONENTS system filesystem program_options iostreams REQUIRED)

include_directories(
                ${Boost_INCLUDE_DIR} 
                ${CFITSIO_INCLUDE_DIRS} 
                ${HEALPIX-CXX_INCLUDE_DIRS} include)


set(LIBS ${LIBS} 
    ${HEALPIX-CXX_LIBRARIES} 
    ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES} ) 



add_executable(multi-llh
            private/iter-lhreco/multi-llh.cc 
            private/iter-lhreco/pickle.cc
            private/iter-lhreco/illh-utils.cc
            )

add_executable(rotation-test
            private/test/rotation.cc
            private/iter-lhreco/illh-utils.cc
            )




target_link_libraries(multi-llh 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})

target_link_libraries(rotation-test PRIVATE
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})


enable_testing()

# define tests
add_test(
  NAME rotation_test_128
  COMMAND $<TARGET_FILE:rotation-test> 128
  )
add_test(
  NAME rotation_test_64
  COMMAND $<TARGET_FILE:rotation-test> 64
  )
add_test(
  NAME rotation_test_32
  COMMAND $<TARGET_FILE:rotation-test> 32
  )
add_test(
  NAME rotation_test_16
  COMMAND $<TARGET_FILE:rotation-test> 16
  )




if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
       install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
         DESTINATION ${CMAKE_BINARY_DIR}/scripts
         PATTERN ".svn" EXCLUDE
         PATTERN ".git" EXCLUDE
         PATTERN "*.py"
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples)
       install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples
         DESTINATION ${CMAKE_BINARY_DIR}/examples
         PATTERN ".svn" EXCLUDE
         PATTERN ".git" EXCLUDE
         PATTERN "*.py"
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples)


