all: javac run

javac:
	cd src; javac changchun/DetectEdits.java -cp .:../lib/htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar -d ../bin/

run:
	cd src; java changchun.DetectEdits ../edit_positions.bed /corgi/henriksson/changchun/p2/P19764_101_S1_L001_R1_001.fastq.gz.out.bam
