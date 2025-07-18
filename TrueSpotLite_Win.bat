::	https://www.mathworks.com/matlabcentral/answers/1666774-creating-shortcut-for-matlab-app-mlapp
::	https://stackoverflow.com/questions/3827567/how-to-get-the-path-of-the-batch-script-in-windows

@echo off

setlocal MYDIR=%~dp0
setlocal MLAPP_PATH=%MY_DIR%\src\gui\TrueSpotMinimal.mlapp

matlab -minimize -nosplash -sd "%MY_DIR%\src" -r "%MLAPP_PATH%;"