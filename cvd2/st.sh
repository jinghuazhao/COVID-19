#!/usr/bin/bash

format.sh

metal ACE2.metal 2>&1 | \
tee ACE2-1.tbl.log

select.sh ACE2
