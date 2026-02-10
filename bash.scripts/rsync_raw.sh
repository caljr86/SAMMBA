#!/bin/bash

# Usage: ./rsync_raw.sh PATTERN PLATE

PATTERN="$1"
PLATE="$2"

TMP_LIST=$(mktemp)

# Remote host details
REMOTE_USER="user"
REMOTE_HOST="10.XX.X.YY" 
REMOTE_BASE="/home/user/$PLATE"

# Get matching files from remote "croped/" folders
ssh $REMOTE_USER@$REMOTE_HOST "find '$REMOTE_BASE' -type f -path '*/raw/*' -name '*$PATTERN*'" > "$TMP_LIST"

sed -i "s|/home/user/||" "$TMP_LIST"

# Rsync matching files to current directory, preserving relative paths
rsync -av --files-from="$TMP_LIST" --relative "$REMOTE_USER@$REMOTE_HOST:/home/user/" .

# Clean up
rm "$TMP_LIST"
