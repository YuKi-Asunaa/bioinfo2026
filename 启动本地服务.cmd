@echo off
cd /d "%~dp0"
title bioinfo2026 preview server
echo.
echo   Starting local preview on http://localhost:8000/
echo   Close this window to stop.
echo.
where node >nul 2>nul
if not errorlevel 1 (
  start "" /min cmd /c "timeout /t 2 >nul & start http://localhost:8000/02-Imaging/index.html"
  node serve.cjs
) else (
  echo Node.js not found, fallback to Python http.server ...
  start "" /min cmd /c "timeout /t 2 >nul & start http://localhost:8000/02-Imaging/index.html"
  python -m http.server 8000
)
pause
