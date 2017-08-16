#!/bin/bash

ls -1 [SED]RR[0-9]*stringtie* > | stringtie --merge -G $1 -o $2
