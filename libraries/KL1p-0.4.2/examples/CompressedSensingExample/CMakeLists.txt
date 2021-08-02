CMAKE_MINIMUM_REQUIRED(VERSION 2.8.1)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Avoid CMake regeneration within generated project.
SET(CMAKE_SUPPRESS_REGENERATION "1")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Project declaration.
PROJECT(CompressedSensingExample)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Compiler configuration
IF(WIN32)
	# External headers folder.
	INCLUDE_DIRECTORIES(../../libs/Armadillo/include)
	INCLUDE_DIRECTORIES(../../include)

	#External libraries folder.
	LINK_DIRECTORIES(../../bin/win/)
	
	
	# Compiler options definition.
	STRING(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG_INIT ${CMAKE_CXX_FLAGS_DEBUG_INIT})
    STRING(REPLACE "/MDd" "/MTd" CMAKE_C_FLAGS_DEBUG_INIT   ${CMAKE_C_FLAGS_DEBUG_INIT})
    STRING(REPLACE "/MD"  "/MT"  CMAKE_C_FLAGS_RELEASE      ${CMAKE_C_FLAGS_RELEASE})
    STRING(REPLACE "/MDd" "/MTd" CMAKE_C_FLAGS_DEBUG        ${CMAKE_C_FLAGS_DEBUG})
    STRING(REPLACE "/MD"  "/MT"  CMAKE_CXX_FLAGS_RELEASE    ${CMAKE_CXX_FLAGS_RELEASE})
    STRING(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG      ${CMAKE_CXX_FLAGS_DEBUG})
		
	# Output path.	
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../../../../../bin/win/)	
	
	FOREACH(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})	# For multi-config builds.
		STRING(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
		SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ../../../../../../bin/win/)
	ENDFOREACH(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)
ENDIF(WIN32)

IF(UNIX)
	IF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)
		IF(${FORCE_32BITS_BUILD} EQUAL 1)
			# Only build 32-bit binaries
			ADD_DEFINITIONS(-m32)
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32")
			SET(CMAKE_SHARED_LIBRARY_C_FLAGS "${CMAKE_SHARED_LIBRARY_C_FLAGS} -m32")
			SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -m32")
			SET(CMAKE_SIZEOF_VOID_P 4)
		ENDIF(${FORCE_32BITS_BUILD} EQUAL 1)
	ENDIF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)
	
	# External headers folder.
	INCLUDE_DIRECTORIES(../../libs/Armadillo/include)
	INCLUDE_DIRECTORIES(../../include)
	
	# External libraries folder.
	LINK_DIRECTORIES(../../bin/unix/)	
	
	# Output path.
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../../../bin/unix/)
	
	FOREACH(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})	# For multi-config builds.
		STRING(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
		SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ../../../../bin/unix/)
	ENDFOREACH(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)
ENDIF(UNIX)

IF(APPLE)
	IF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)
		IF(${FORCE_32BITS_BUILD} EQUAL 1)
			# Only build 32-bit binaries
			ADD_DEFINITIONS(-m32)
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32")
			SET(CMAKE_SHARED_LIBRARY_C_FLAGS "${CMAKE_SHARED_LIBRARY_C_FLAGS} -m32")
			SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -m32")
			SET(CMAKE_SIZEOF_VOID_P 4)
		ENDIF(${FORCE_32BITS_BUILD} EQUAL 1)
	ENDIF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)
	
	# External headers folder.
	INCLUDE_DIRECTORIES(../../libs/Armadillo/include)
	INCLUDE_DIRECTORIES(../../include)
	
	# External libraries folder.
	LINK_DIRECTORIES(../../bin/macos/)	
	
	# Output path.
	SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../../../bin/macos/)
	
	FOREACH(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})	# For multi-config builds.
		STRING(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
		SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ../../../../bin/macos/)
	ENDFOREACH(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)
ENDIF(APPLE)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Generation of the source files.
IF(WIN32)
	FILE(
		GLOB
		SourceFiles
		*.h
		*.cpp 
		win/*.h
		win/*.cpp      
	)
ENDIF(WIN32)

IF(UNIX)
	FILE(
		GLOB
		SourceFiles
		*.h
		*.cpp 
		unix/*.h
		unix/*.cpp 
	)
ENDIF(UNIX)

IF(APPLE)
	FILE(
		GLOB
		SourceFiles
		*.h
		*.cpp 
		macos/*.h
		macos/*.cpp 
	)
ENDIF(APPLE)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Project filters definition.
SOURCE_GROUP("Headers" REGULAR_EXPRESSION .*[.]h)
SOURCE_GROUP("Sources" REGULAR_EXPRESSION .*[.]cpp)

IF(WIN32)
	SOURCE_GROUP("Sources" FILES win/main.cpp)
ENDIF(WIN32)

IF(UNIX)
	SOURCE_GROUP("Sources" FILES unix/main.cpp)
ENDIF(UNIX)

IF(APPLE)
	SOURCE_GROUP("Sources" FILES macos/main.cpp)
ENDIF(APPLE)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Declaration of the executable
IF(CMAKE_BUILD_TYPE MATCHES "Debug")	
	IF(CMAKE_SIZEOF_VOID_P EQUAL 4)						# 32 bits generation.
		ADD_EXECUTABLE(CompressedSensingExample_d ${SourceFiles})
		IF(WIN32)
			TARGET_LINK_LIBRARIES(CompressedSensingExample_d KLab_d)
		ENDIF(WIN32)
		IF(UNIX)
			TARGET_LINK_LIBRARIES(CompressedSensingExample_d KLab_d pthread dl)
		ENDIF(UNIX)
		IF(APPLE)
			TARGET_LINK_LIBRARIES(CompressedSensingExample_d KLab_d pthread dl "-framework CoreServices")
		ENDIF(APPLE)
	ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)		
	IF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)		# 64bits generation.
		ADD_EXECUTABLE(CompressedSensingExample64_d ${SourceFiles})
		IF(WIN32)
			TARGET_LINK_LIBRARIES(CompressedSensingExample64_d KLab64_d)
		ENDIF(WIN32)
		IF(UNIX)
			TARGET_LINK_LIBRARIES(CompressedSensingExample64_d KLab64_d pthread dl)
		ENDIF(UNIX)
		IF(APPLE)
			TARGET_LINK_LIBRARIES(CompressedSensingExample64_d KLab64_d pthread dl "-framework CoreServices")
		ENDIF(APPLE)
	ENDIF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")	
IF(CMAKE_BUILD_TYPE MATCHES "Release")
	IF(CMAKE_SIZEOF_VOID_P EQUAL 4)						# 32 bits generation.
		ADD_EXECUTABLE(CompressedSensingExample ${SourceFiles})	
		IF(WIN32)
			TARGET_LINK_LIBRARIES(CompressedSensingExample KLab)
		ENDIF(WIN32)
		IF(UNIX)
			TARGET_LINK_LIBRARIES(CompressedSensingExample KLab pthread dl)
		ENDIF(UNIX)
		IF(APPLE)
			TARGET_LINK_LIBRARIES(CompressedSensingExample KLab pthread dl "-framework CoreServices")
		ENDIF(APPLE)
	ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
	IF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)		# 64bits generation.
		ADD_EXECUTABLE(CompressedSensingExample64 ${SourceFiles})
		IF(WIN32)
			TARGET_LINK_LIBRARIES(CompressedSensingExample64 KLab64)
		ENDIF(WIN32)
		IF(UNIX)
			TARGET_LINK_LIBRARIES(CompressedSensingExample64 KLab64 pthread dl)
		ENDIF(UNIX)
		IF(APPLE)
			TARGET_LINK_LIBRARIES(CompressedSensingExample64 KLab64 pthread dl "-framework CoreServices")
		ENDIF(APPLE)
	ENDIF(CMAKE_CL_64 OR CMAKE_SIZEOF_VOID_P EQUAL 8)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Release")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
