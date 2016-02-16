start /affinity F PPTT.exe
wmic process where name="PPTT.exe" CALL setpriority 256