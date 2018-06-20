INDEX_DIR = /var/data/centrifuge-indexes

test:
	./run_centrifuge.py \
		-q $(HOME)/in \
		-I $(INDEX_DIR) -i p_compressed+h+v -f fastq \
		-o $(HOME)/centrifuge-out \
		-X 10000 -t 4 -P 4

test_paired:
	./run_centrifuge.py \
		-1 $(HOME)/in/DNA_1_CGATGT_L002_R1_001.fastq,$(HOME)/in/DNA_1_CGATGT_L002_R1_002.fastq,$(HOME)/in/DNA_1_CGATGT_L002_R1_003.fastq,$(HOME)/in/DNA_2_TGACCA_L002_R1_001.fastq,$(HOME)/in/DNA_2_TGACCA_L002_R1_002.fastq,$(HOME)/in/DNA_2_TGACCA_L002_R1_003.fastq,$(HOME)/in/DNA_3_ACAGTG_L002_R1_001.fastq,$(HOME)/in/DNA_3_ACAGTG_L002_R1_002.fastq,$(HOME)/in/DNA_3_ACAGTG_L002_R1_003.fastq,$(HOME)/in/DNA_4_GCCAAT_L002_R1_001.fastq,$(HOME)/in/DNA_4_GCCAAT_L002_R1_002.fastq,$(HOME)/in/DNA_4_GCCAAT_L002_R1_003.fastq \
		-2 $(HOME)/in/DNA_1_CGATGT_L002_R2_001.fastq,$(HOME)/in/DNA_1_CGATGT_L002_R2_002.fastq,$(HOME)/in/DNA_1_CGATGT_L002_R2_003.fastq,$(HOME)/in/DNA_2_TGACCA_L002_R2_001.fastq,$(HOME)/in/DNA_2_TGACCA_L002_R2_002.fastq,$(HOME)/in/DNA_2_TGACCA_L002_R2_003.fastq,$(HOME)/in/DNA_3_ACAGTG_L002_R2_001.fastq,$(HOME)/in/DNA_3_ACAGTG_L002_R2_002.fastq,$(HOME)/in/DNA_3_ACAGTG_L002_R2_003.fastq,$(HOME)/in/DNA_4_GCCAAT_L002_R2_001.fastq,$(HOME)/in/DNA_4_GCCAAT_L002_R2_002.fastq,$(HOME)/in/DNA_4_GCCAAT_L002_R2_003.fastq \
		-I $(INDEX_DIR) -i p_compressed+h+v -f fastq \
		-o $(HOME)/centrifuge-out \
		-x 9606,32630 -t 4 -P 4


