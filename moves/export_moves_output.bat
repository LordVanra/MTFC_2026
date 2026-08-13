@echo off
REM ===========================================================================
REM  export_moves_output.bat
REM  MOVES writes results into MariaDB databases, not files. This finds the
REM  bundled MariaDB client and exports what we need to tab-separated files.
REM
REM  Put this file anywhere and double-click it, or run from a Command Prompt.
REM  Output lands next to this file as *.txt (tab-separated - open in Excel).
REM ===========================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo.
echo === Locating the MariaDB client shipped with MOVES ===

set MYSQL=
where mysql.exe >nul 2>&1 && set MYSQL=mysql.exe

if "!MYSQL!"=="" (
  for %%P in (
    "C:\Program Files\MariaDB 10.6\bin\mysql.exe"
    "C:\Program Files\MariaDB 10.11\bin\mysql.exe"
    "C:\Program Files\MariaDB 11.4\bin\mysql.exe"
    "C:\Program Files (x86)\MariaDB 10.6\bin\mysql.exe"
    "C:\MOVES5.0.1\MariaDB\bin\mysql.exe"
    "%USERPROFILE%\MOVES5.0.1\MariaDB\bin\mysql.exe"
    "C:\MOVES5\MariaDB\bin\mysql.exe"
  ) do (
    if exist %%P set MYSQL=%%P
  )
)

if "!MYSQL!"=="" (
  echo.
  echo Could not find mysql.exe automatically.
  echo Search for it with this in a Command Prompt, then edit the MYSQL line above:
  echo     dir /s /b "C:\Program Files\mysql.exe"
  echo     dir /s /b "C:\MOVES*\mysql.exe"
  pause
  exit /b 1
)

echo Using: !MYSQL!

REM MOVES default credentials
set CRED=-u moves -pmoves
!MYSQL! %CRED% -e "SELECT 1" >nul 2>&1
if errorlevel 1 (
  echo moves/moves rejected, trying root with no password...
  set CRED=-u root
  !MYSQL! !CRED! -e "SELECT 1" >nul 2>&1
  if errorlevel 1 (
    echo.
    echo Could not connect. Check Settings inside MOVES for the database user.
    pause
    exit /b 1
  )
)

echo.
echo === Databases on this server ===
!MYSQL! !CRED! -e "SHOW DATABASES" > 00_databases.txt
type 00_databases.txt

for %%D in (bos_2025_out bos_2050_out bos_2026_out) do (
  !MYSQL! !CRED! -e "USE %%D" >nul 2>&1
  if not errorlevel 1 (
    echo.
    echo === Exporting %%D ===

    REM What tables exist - tells us whether this was Rates or Inventory
    !MYSQL! !CRED! -e "SHOW TABLES" %%D > %%D_00_tables.txt

    REM Rates output: the numbers we actually want
    !MYSQL! !CRED! --batch --raw %%D -e ^
      "SELECT MOVESScenarioID, yearID, monthID, hourID, fuelTypeID, processID, pollutantID, avgSpeedBinID, roadTypeID, sourceTypeID, SUM(ratePerDistance) AS ratePerDistance FROM rateperdistance WHERE pollutantID IN (110,116,117) GROUP BY MOVESScenarioID,yearID,monthID,hourID,fuelTypeID,processID,pollutantID,avgSpeedBinID,roadTypeID,sourceTypeID" ^
      > %%D_01_rateperdistance.txt 2>%%D_01_error.txt

    REM Speed bin reference
    !MYSQL! !CRED! --batch --raw %%D -e "SELECT * FROM avgspeedbin" > %%D_02_speedbins.txt 2>nul

    REM Run log - records any warnings MOVES raised
    !MYSQL! !CRED! --batch --raw %%D -e "SELECT * FROM movesrun" > %%D_03_movesrun.txt 2>nul
    !MYSQL! !CRED! --batch --raw %%D -e "SELECT * FROM moveserror" > %%D_04_errors.txt 2>nul

    echo Wrote %%D_*.txt
  )
)

echo.
echo === Done. Send the *.txt files back. ===
pause
