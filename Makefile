all: javac run

javac:
	cd src; javac changchun/DetectEdits.java -cp .:../lib/htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar:../lib/sqlite-jdbc-3.34.0.jar -d ../bin/
	#cd src; javac changchun/FormatAlignment.java -cp .:../lib/htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar:../lib/sqlite-jdbc-3.34.0.jar -d ../bin/


run:
	cd bin; java -cp .:../lib/htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar:../lib/sqlite-jdbc-3.34.0.jar \
		changchun.DetectEdits \
		../edit_positions.bed \
		/corgi/henriksson/changchun/p2/P19764_101_S1_L001_R1_001.fastq.gz.out.bam \
		out.stat

toupp:
	scp -r . rackham.uppmax.uu.se:/home/mahogny/mystore/dataset/changchun/

tows:
	scp alignments.sqlite beagle.henlab.org:/data/henlab/project/changchun_mutant/website/data
