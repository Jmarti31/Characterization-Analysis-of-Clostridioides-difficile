#!/usr/bin/env bash
# getNGS.sh


# Retrieve the Clostridioides difficile NGS reads.
fasterq-dump --split-3 $SRR_ID -O data/
