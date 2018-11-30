#!/bin/bash

./wrmssgenc res .enc 1 2 1 1e-6 0
./wrmssgdec res .enc dec 1 2 1 0

find . -type f -name '*.enc' -exec du -ch {} + | grep total$
