#!/bin/bash

nohup mpiexec -np 8 ./roms roms.in  >roms.log 2> roms.err &
