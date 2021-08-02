ECHO OFF

SET CURRENT_PATH=%CD%\
SET BUILD_PATH=..\..\libs\KLab\main\win\


cd %BUILD_PATH%
call "Visual2005_x64Debug_Make.bat"

cd %CURRENT_PATH%
cd %BUILD_PATH%
cd Visual2005_x64Debug
call "build.bat"

cd %CURRENT_PATH%
cd %BUILD_PATH%
call "Visual2005_x64Release_Make.bat"

cd %CURRENT_PATH%
cd %BUILD_PATH%
cd Visual2005_x64Release
call "build.bat"
