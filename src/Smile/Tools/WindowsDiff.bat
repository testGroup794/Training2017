@echo off
:loop
    j-rand > j.in
    j-ans < j.in > j-ans.out
    j-me < j.in > j.out
    fc j.out j-ans.out
if %errorlevel%==0 goto loop
pause