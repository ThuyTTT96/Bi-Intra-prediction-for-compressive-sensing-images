ECHO OFF

SET CURRENT_PATH=%CD%\
SET BUILD_PATH=..\..\libs\KLab\main\win\


cd %BUILD_PATH%
call "Visual2005_x86Debug_Make.bat"

cd %CURRENT_PATH%
cd %BUILD_PATH%
cd Visual2005_x86Debug
call "build.bat"

cd %CURRENT_PATH%
cd %BUILD_PATH%
call "Visual2005_x86Release_Make.bat"

cd %CURRENT_PATH%
cd %BUILD_PATH%
cd Visual2005_x86Release
call "build.bat"
