IMG="../stampede/radcot.img"
INPUT_DIR=$(WORK)/in-smaller
OUT_DIR=$(WORK)/radcot_test
GENOME_DIR=$(WORK)/genomes

test:
	singularity exec $(IMG) run_centrifuge.py \
		-q $(WORK)/in-smaller \
		-i p_compressed+h+v -f fastq \
		-o $(WORK)/radcot-test \
		-X 10000 -t 4 -P 4

test_paired:
	singularity exec $(IMG) run_centrifuge.py \
		-1 $(WORK)/in-smaller/DNA_1_CGATGT_L002_R1_001.fastq,$(WORK)/in-smaller/DNA_1_CGATGT_L002_R1_002.fastq,$(WORK)/in-smaller/DNA_1_CGATGT_L002_R1_003.fastq,$(WORK)/in-smaller/DNA_2_TGACCA_L002_R1_001.fastq,$(WORK)/in-smaller/DNA_2_TGACCA_L002_R1_002.fastq,$(WORK)/in-smaller/DNA_2_TGACCA_L002_R1_003.fastq,$(WORK)/in-smaller/DNA_3_ACAGTG_L002_R1_001.fastq,$(WORK)/in-smaller/DNA_3_ACAGTG_L002_R1_002.fastq,$(WORK)/in-smaller/DNA_3_ACAGTG_L002_R1_003.fastq,$(WORK)/in-smaller/DNA_4_GCCAAT_L002_R1_001.fastq,$(WORK)/in-smaller/DNA_4_GCCAAT_L002_R1_002.fastq,$(WORK)/in-smaller/DNA_4_GCCAAT_L002_R1_003.fastq \
		-2 $(WORK)/in-smaller/DNA_1_CGATGT_L002_R2_001.fastq,$(WORK)/in-smaller/DNA_1_CGATGT_L002_R2_002.fastq,$(WORK)/in-smaller/DNA_1_CGATGT_L002_R2_003.fastq,$(WORK)/in-smaller/DNA_2_TGACCA_L002_R2_001.fastq,$(WORK)/in-smaller/DNA_2_TGACCA_L002_R2_002.fastq,$(WORK)/in-smaller/DNA_2_TGACCA_L002_R2_003.fastq,$(WORK)/in-smaller/DNA_3_ACAGTG_L002_R2_001.fastq,$(WORK)/in-smaller/DNA_3_ACAGTG_L002_R2_002.fastq,$(WORK)/in-smaller/DNA_3_ACAGTG_L002_R2_003.fastq,$(WORK)/in-smaller/DNA_4_GCCAAT_L002_R2_001.fastq,$(WORK)/in-smaller/DNA_4_GCCAAT_L002_R2_002.fastq,$(WORK)/in-smaller/DNA_4_GCCAAT_L002_R2_003.fastq \
		-i p_compressed+h+v -f fastq \
		-o $(WORK)/radcot-test \
		-x 9606,32630 -t 4 -P 4

test_runall:
	./runAll.py \
		-i $(INPUT_DIR) \
		-o $(OUT_DIR) \
		-g $(GENOME_DIR) \
		-m $(WORK)/subset_crc_mouse_metadata.txt \
		-d -t 12 -P 4

test_bowtie:
	./runAll.py \
		-o $(OUT_DIR) \
		-i $(INPUT_DIR) \
		-g $(GENOME_DIR) \
		-m $(WORK)/subset_crc_mouse_metadata.txt \
		--debug --procs 4 \
		--skip-centrifuge

