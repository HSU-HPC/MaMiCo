#!/bin/bash

for((i=16;i<17; i++)) do 
	rm *.dat
	rm *.txt
	rm CheckpointSimpleMD_10000_periodic_0_0_0__0_0.checkpoint
	./test

	mv CheckpointSimpleMD_10000_periodic_0_0_0__10000_0.checkpoint Checkpoints/CheckpointSimpleMD_10000_periodic_${i}.checkpoint
done
	
