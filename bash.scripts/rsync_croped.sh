#!/bin/bash

# Usage: ./rsync_croped.sh PATTERN PLATE

PATTERN="$1"
PLATE="$2"

TMP_LIST=$(mktemp)

# Remote host details
REMOTE_USER="user"
REMOTE_HOST="server IP"
REMOTE_BASE="/home/Disk1/user/$PLATE"

# Get matching files from remote "croped/" folders
ssh $REMOTE_USER@$REMOTE_HOST "find '$REMOTE_BASE' -type f -path '*/croped/*' -name '*$PATTERN*'" > "$TMP_LIST"

sed -i "s|/home/Disk1/user/||" "$TMP_LIST"

# Rsync matching files to current directory, preserving relative paths
rsync -av --files-from="$TMP_LIST" --relative "$REMOTE_USER@$REMOTE_HOST:/home/Disk1/user/" .

# Clean up
rm "$TMP_LIST"
