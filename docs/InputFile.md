Table format for providing read files
-------------------------------------

This pipeline is designed to run parallel on a large number of samples. Since every sample can contain several read files in fastq format the input is provided in the form of a tabular input file. It must have the followin format 

|Sample_Id|Long_Reads|Short_Reads_1|Short_reads_2|
----------|-------------|-------------|----------|
|testSet |testReads/testReads_nanopore.fastq.gz|testReads/testReads_illumina_S1.fastq.gz|testReads/testReads_illumina_S2.fastq.gz|

Content of the columns:
   1) *Sample_Id* Unique identifier to this sample, will be to label output files.
   2) *Short_Reads_1* Path to .fastq.gz file with part 1 of Illumina paired end reads.
   3) *Short_Reads_2* Path to .fastq.gz file with part 2 of Illumina paired end reads.
   4) *Long_Reads*  Path to .fastq.gz file containing long reads

See also the includede Testfile `samples_test.tsv` for reference. 

### Multiple Nanopore ReadFiles:
- If you have several FastQ files for one sample, which is sometimes the case for ONT basecalling output you need to combine them before using them in this pipeline. (`zcat file1.fastq.gz file2.fastq.gz > combined.fastq.gz`)



