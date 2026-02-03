#!/bin/bash

# Usage: ./rsync_croped.sh PATTERN PLATE

PATTERN="$1"
PLATE="$2"

TMP_LIST=$(mktemp)

# Remote host details
REMOTE_USER="gareth"
REMOTE_HOST="10.36.4.90"
REMOTE_BASE="/media/Disk2/cicero/$PLATE"

# Get matching files from remote "croped/" folders
ssh $REMOTE_USER@$REMOTE_HOST "find '$REMOTE_BASE' -type f -path '*/croped/*' -name '*$PATTERN*'" > "$TMP_LIST"

# Adjust the paths to be relative from $PLATE
# For example: /media/Disk1/cicero/plate58/250706/croped/img_C13.zip
# becomes: plate58/250706/croped/img_C13.zip
sed -i "s|/media/Disk2/cicero/||" "$TMP_LIST"

# Rsync matching files to current directory, preserving relative paths
rsync -av --files-from="$TMP_LIST" --relative "$REMOTE_USER@$REMOTE_HOST:/media/Disk2/cicero/" .

# Clean up
rm "$TMP_LIST"
