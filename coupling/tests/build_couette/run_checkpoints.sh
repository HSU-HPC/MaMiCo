#!/bin/bash

for((i=1;i<10; i++)) do 
	rm *.dat
	rm *.txt
	rm CheckpointSimpleMD_10000_periodic_0_0_0__0_0.checkpoint
	./test

	mv CheckpointSimpleMD_10000_periodic_0_0_0__10000_0.checkpoint Checkpoints2/CheckpointSimpleMD_10000_periodic_${i}.checkpoint
done
	
