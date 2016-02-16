start /affinity 1 PPTT.exe
wmic process where name="PPTT.exe" CALL setpriority 256