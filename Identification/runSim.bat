REM echo "Starting" >> ..\debug.txt
..\dymosim dsin.txt 
C:\Users\filip\Anaconda3\python.exe "..\..\post_process.py" > post_process_log.txt
REM echo "Stopping" >> ..\debug.txt

