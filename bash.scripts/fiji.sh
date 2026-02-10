fiji=/home/user/Downloads/fiji-linux64/Fiji.app/ImageJ-linux64
taskset -c $(shuf -i 0-$(($(nproc --all)-1)) -n $1 | paste -sd ,) $fiji
