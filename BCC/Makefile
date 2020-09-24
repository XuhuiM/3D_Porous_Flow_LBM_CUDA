CC=gcc

All: flowmain

common.o: common.c
	$(CC) -o $@ -c -O3 $<

flowmain: *.cu common.o
	nvcc flowmain.cu common.o -o flowmain -gencode arch=compute_35,code=compute_35 -I ~/NVIDIA_CUDA-8.0_Samples/common/inc 
	
clean:
	rm -rf flowmain *.o *.dat

