
add_subdirectory(openfast)

if(BUILD_OPENFAST_CPP_API)
   add_subdirectory(openfast-cpp)
endif()

if(BUILD_FASTFARM)
   add_subdirectory(fast-farm)
endif()

if(BUILD_OPENFAST_PYTHONBINDINGS)
   add_subdirectory(openfast-pythonbindings)
endif()