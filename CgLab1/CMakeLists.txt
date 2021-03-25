# Top-Level CMakeList.txt

cmake_minimum_required (VERSION 2.6)
project ( FirstLab )

add_executable( FirstLab skeleton.cpp)
add_executable( StarField starfield.cpp)

set_target_properties(StarField
    PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)

find_package (SDL)

if ( NOT SDL_FOUND )
   message ( FATAL_ERROR "SDL not found!" )
endif ( NOT SDL_FOUND )

if(SDL_FOUND)
	include_directories(
		${SDL_INCLUDE_DIR}
		${PROJECT_SOURCE_DIR}/glm
	)
	#link_libraries(${SDL_LIBRARY})
	target_link_libraries(FirstLab ${SDL_LIBRARY})
	target_link_libraries(StarField ${SDL_LIBRARY})
endif(SDL_FOUND)

