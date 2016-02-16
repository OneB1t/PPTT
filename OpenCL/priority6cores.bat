start /affinity 3F PPTT.exe
wmic process where name="PPTT.exe" CALL setpriority 256