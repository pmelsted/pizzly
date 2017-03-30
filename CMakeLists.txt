# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (pizzly CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
find_package( Threads REQUIRED )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories ( . )
include_directories( ${ZLIB_INCLUDE_DIRS} )

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -std=c++14  -DSEQAN_HAS_ZLIB=1" )
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )

# Add executable and link against SeqAn dependencies.
add_executable (pizzly main.cpp)
target_link_libraries (pizzly pthread ${ZLIB_LIBRARIES} )
