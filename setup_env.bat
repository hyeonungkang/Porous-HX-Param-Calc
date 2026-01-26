@echo off
setlocal

set VENV_DIR=.venv

if not exist %VENV_DIR% (
    echo [INFO] Creating virtual environment...
    python -m venv %VENV_DIR%
)

call %VENV_DIR%\Scripts\activate

echo [INFO] Installing requirements...
pip install -r requirements.txt

echo [INFO] Virtual environment is ready.
cmd /k
