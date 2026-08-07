#!/bin/zsh
# Drive ch3_leantall_drive to completion, restarting MATLAB when it dies.
# The driver is resumable (Results/ch3_leantall_state.mat), so each restart
# continues at the first unfinished stage.

cd "/Users/mohammadmalek/Desktop/Research Files/RABBIT-Bipedal-Robot" || exit 1
SP="$(pwd)/Chapter3/Optimization"
LOG="Results/ch3_leantall.log"
OUT="Results/ch3_leantall.stdout"

for i in $(seq 1 40); do
    if grep -qE "MARKER_ALLDONE|MARKER_STOPPED" "$LOG" 2>/dev/null; then break; fi
    echo "[wrapper] attempt $i at $(date +%H:%M:%S)" >> "$LOG"
    /Applications/MATLAB_R2021b.app/bin/matlab -batch \
        "addpath('$SP'); ch3_leantall_drive" >> "$OUT" 2>&1
done

if grep -q MARKER_ALLDONE "$LOG" 2>/dev/null; then
    echo "[wrapper] COMPLETE at $(date +%H:%M:%S)" >> "$LOG"
elif grep -q MARKER_STOPPED "$LOG" 2>/dev/null; then
    echo "[wrapper] STOPPED by a failed verify at $(date +%H:%M:%S)" >> "$LOG"
else
    echo "[wrapper] GAVE UP after 40 attempts at $(date +%H:%M:%S)" >> "$LOG"
fi
