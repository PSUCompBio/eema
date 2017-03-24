if(${CMAKE_ARGC} GREATER 1)
	if ("${CMAKE_ARGV3}" STREQUAL "ENONLY")
		message("#define ENGLISH_ONLY")
	else()
		math(EXPR UPTO ${CMAKE_ARGC}-1)
		foreach(i RANGE 3 ${UPTO})
			message("#define LANG_${CMAKE_ARGV${i}}")
		endforeach()
	endif()
endif()

