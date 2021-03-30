function(link_open_ephys_lib target libname)
	add_library(${libname} SHARED IMPORTED)

	if(MSVC)
		set(LIBLOC IMPORTED_IMPLIB)
	else()
		set(LIBLOC IMPORTED_LOCATION)
	endif()

	# try to find release version first
	set(LIB_CONFIGS ${CMAKE_CONFIGURATION_TYPES})
	list(INSERT LIB_CONFIGS 0 Release)
	list(REMOVE_DUPLICATES LIB_CONFIGS)

	foreach(config ${LIB_CONFIGS})
		string(TOUPPER ${config} CONFIG)

		set(LIB_PATH_${CONFIG} ${libname}-NOTFOUND) # force re-search

		find_library(LIB_PATH_${CONFIG} NAMES ${libname}
			PATHS ${GUI_COMMONLIB_DIR}/${config}/lib/${CMAKE_LIBRARY_ARCHITECTURE})

		if (LIB_PATH_${CONFIG})
			set_target_properties(${libname} PROPERTIES ${LIBLOC}_${CONFIG} ${LIB_PATH_${CONFIG}})

			# use Release as default
			if (${config} STREQUAL "Release")
				set_target_properties(${libname} PROPERTIES ${LIBLOC} ${LIB_PATH_${CONFIG}})
			endif()

		elseif(LIB_PATH_RELEASE)
			message(WARNING "${libname} ${config} library not found - will link to Release version instead.")
		else()
			message(WARNING "${libname} ${config} library not found - this configuration will not be available.")
		endif()
	endforeach()

	target_link_libraries(${target} ${libname})
endfunction(link_open_ephys_lib)