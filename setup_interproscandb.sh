#!/usr/bin/env bash
echo "=== Installing InterProScan DB ==="
mkdir -p databases/interproscan
cd databases/interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz
tar -pxvzf interproscan-5.78-109.0-64-bit.tar.gz
cd interproscan-5.78-109.0
python3 setup.py -f interproscan.properties
echo "Done! In Nextflow, use: --interproscan_db $(pwd)"