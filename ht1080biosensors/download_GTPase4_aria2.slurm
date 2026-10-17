#!/usr/bin/env bash
#SBATCH --job-name=g4-aria2
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --output=/home/4473331/projects/scaling_sandbox/ht1080biosensors/GTPase4.%j.out
#SBATCH --error=/home/4473331/projects/scaling_sandbox/ht1080biosensors/GTPase4.%j.err

set -euo pipefail

destination_dir=/home/4473331/projects/scaling_sandbox/ht1080biosensors
destination_file="$destination_dir/GTPase4.rar"
download_url=https://ndownloader.figshare.com/files/31467194
expected_size=21056813510
expected_md5=c468b3f693e13e33ba44f4902d002738
aria2c=/home/4473331/.conda/envs/scaling_sandbox/bin/aria2c

echo "Starting at $(date --iso-8601=seconds) on $(hostname)"
"$aria2c" \
    --no-want-digest-header=true \
    --use-head=false \
    --user-agent=curl/7.61.1 \
    --continue=true \
    --always-resume=true \
    --max-connection-per-server=16 \
    --split=16 \
    --min-split-size=1M \
    --file-allocation=none \
    --max-tries=20 \
    --retry-wait=5 \
    --connect-timeout=30 \
    --timeout=60 \
    --check-certificate=true \
    --ca-certificate=/etc/pki/tls/certs/ca-bundle.crt \
    --console-log-level=notice \
    --summary-interval=60 \
    --dir="$destination_dir" \
    --out=GTPase4.rar \
    "$download_url"

actual_size=$(stat --printf='%s' "$destination_file")
if (( actual_size != expected_size )); then
    echo "ERROR: Size mismatch: ${actual_size} != ${expected_size}" >&2
    exit 1
fi

printf '%s  %s\n' "$expected_md5" "$destination_file" | md5sum --check
echo "Size and MD5 verified at $(date --iso-8601=seconds)"
