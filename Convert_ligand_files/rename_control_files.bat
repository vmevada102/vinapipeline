@echo off
setlocal enabledelayedexpansion

for %%f in (*_mmff94_E=*.pdbqt) do (
    set "fname=%%~nf"
    rem Extract part before "_mmff94"
    for /f "tokens=1 delims=_" %%a in ("!fname!") do (
        if /I not "%%a"=="Control" (
            ren "%%f" "Control_%%a.pdbqt"
        )
    )
)

endlocal
