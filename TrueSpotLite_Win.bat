::	https://www.mathworks.com/matlabcentral/answers/1666774-creating-shortcut-for-matlab-app-mlapp
::	https://stackoverflow.com/questions/3827567/how-to-get-the-path-of-the-batch-script-in-windows

@echo off

set MYDIR=%~dp0
:: echo %MYDIR%
set MLAPP_PATH=%MYDIR%src\gui\TrueSpotMinimal.mlapp

matlab -minimize -nosplash -sd "%MYDIR%src" -r "addpath('./gui'); TrueSpotMinimal;"