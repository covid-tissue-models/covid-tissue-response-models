@ECHO OFF

REM Set PREFIX_CC3D to the root directory of your CC3D installation!
@SET PREFIX_CC3D=C:\CompuCell3D-py3-64bit


REM Environment setup

@SET PYTHON_INSTALL_PATH=%PREFIX_CC3D%\Python36
@SET PYTHONPATH=%PREFIX_CC3D%\lib\site-packages

REM Run it!

@SET exit_code=0
"%PYTHON_INSTALL_PATH%\python" "%CD%\batch_run.py"
@SET exit_code= %errorlevel%
