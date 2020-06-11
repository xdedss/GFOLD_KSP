@echo off
echo process.bat
if "%1"=="" (
echo "Empty input!"
pause
goto :end
)
@echo on
if exist GFOLD_%1 (rd /s/q GFOLD_%1)
python ./GFOLD_codegen.py %1
python ./postprocess.py GFOLD_%1 gfold_solver_%1
cd GFOLD_%1
python ./setup.py install
@echo -
@echo -
@echo OK
:end
