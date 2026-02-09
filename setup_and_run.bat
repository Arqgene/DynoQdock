@echo off
setlocal enabledelayedexpansion

:: ArqGene - Ultimate Automated Setup & Start Script
:: This script strictly follows the Step-by-Step Guide for Windows installation.

title ArqGene - Automated Molecular Docking Platform Setup

echo ============================================================
echo           ArqGene - Molecular Docking Platform
echo           Complete Automated Installation ^& Setup
echo ============================================================
echo.

:: STEP 1: INSTALL PYTHON (Strictly following Guide Step 1)
echo [+] Step 1: Checking for Python 3.11...
python --version 2>nul | findstr "3.11" >nul
if %errorlevel% neq 0 (
    echo [!] Python 3.11 not found. Downloading Python 3.11.x...
    powershell -Command "[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri 'https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe' -OutFile 'python_installer.exe'"
    if not exist "python_installer.exe" (
        echo [!] Failed to download Python installer. Please check your internet connection.
        pause
        exit /b 1
    )
    echo [+] Installing Python 3.11 and checking 'Add Python to PATH'...
    echo [!] Starting installer... Please follow the on-screen instructions if a window appears.
    
    :: Using /passive instead of /quiet to show progress while still being non-interactive
    start /wait "" "python_installer.exe" /passive InstallAllUsers=1 PrependPath=1 Include_test=0
    
    del python_installer.exe
    echo [+] Python installation complete.
    
    :: Force path refresh for current session using multiple possible install locations
    set "PY_PATH=C:\Program Files\Python311"
    if not exist "!PY_PATH!\python.exe" set "PY_PATH=%LocalAppData%\Programs\Python\Python311"
    
    set "PATH=!PY_PATH!;!PY_PATH!\Scripts;%PATH%"
    
    :: Final check
    python --version >nul 2>&1
    if %errorlevel% neq 0 (
        echo [!] Python was installed but is still not in the PATH. 
        echo [!] Please restart this script or your computer.
        pause
        exit /b 1
    )
) else (
    echo [OK] Python 3.11 is already installed.
)

:: STEP 2: INSTALL OPENBABEL (Strictly following Guide Step 2)
echo [+] Step 2: Checking for OpenBabel 3.1.1...
obabel -V 2>nul >nul
if %errorlevel% neq 0 (
    echo [!] OpenBabel not found. Downloading OpenBabel-3.1.1-x64.exe...
    powershell -Command "[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri 'https://github.com/openbabel/openbabel/releases/download/openbabel-3-1-1/OpenBabel-3.1.1-x64.exe' -OutFile 'babel_installer.exe'"
    if not exist "babel_installer.exe" (
        echo [!] Failed to download OpenBabel installer.
        pause
        exit /b 1
    )
    echo [+] Installing OpenBabel 3.1.1...
    :: Using /passive for better visual feedback
    start /wait "" "babel_installer.exe" /S
    del babel_installer.exe
    echo [+] Adding OpenBabel to PATH (C:\Program Files\OpenBabel-3.1.1)...
    set "PATH=%PATH%;C:\Program Files\OpenBabel-3.1.1"
    setx PATH "%PATH%;C:\Program Files\OpenBabel-3.1.1" /M >nul
) else (
    echo [OK] OpenBabel is already installed.
)

:: STEP 3: INSTALL SMINA (Strictly following Guide Step 3)
echo [+] Step 3: Checking for Smina docking engine...
if not exist "smina.exe" (
    echo [!] Smina not found. Downloading smina.static.exe...
    powershell -Command "[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri 'https://sourceforge.net/projects/smina/files/smina.static.exe/download' -OutFile 'smina.exe'"
    echo [OK] Smina executable prepared in project root.
) else (
    echo [OK] Smina is already present.
)

:: STEP 4: SETUP PROJECT (Strictly following Guide Step 4)
echo [+] Step 4: Setting up Project Environment...
if not exist "venv" (
    echo [+] Creating Virtual Environment (venv)...
    python -m venv venv
)

echo [+] Activating Virtual Environment...
call "venv\Scripts\activate.bat"

echo [+] Installing/Updating Python Dependencies from requirements.txt...
python -m pip install --upgrade pip
pip install -r requirements.txt

:: Initialize required data directories
if not exist "data" mkdir data
if not exist "data\poses" mkdir data\poses

:: STEP 5: RUN APPLICATION (Strictly following Guide Step 5)
echo.
echo ============================================================
echo   Setup Complete! Launching ArqGene...
echo ============================================================
echo.

:: Start the Flask application in its own window
echo [+] Starting ArqGene Server...
start "ArqGene Server" cmd /c "python main.py"

:: Wait for server initialization
echo [+] Waiting for server to initialize...
timeout /t 15 /nobreak >nul

:: Automatically open the web browser
echo [+] Launching Web Browser (http://localhost:5000)...
:: Using explorer to ensure it opens in default browser correctly
explorer "http://localhost:5000"

echo.
echo ============================================================
echo   Server is running. Close the background window to stop.
echo ============================================================
echo.

pause
