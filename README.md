# genome_size_estimation_from_coverage GSEC.py

A script to trim and estimate genome sequencing coverage from a base coverage file. Base coverage file should be a tab delimited file with a sequence name, position number, and coverage. These can be generated by many genome manipulation tools, including samtools (using the depth command), BBMap/pileup.sh (with the "basecov=" option), bedtools (with the "-d" option).

Trimming the ends of coverage pileups is often necessary since read mapping programs tend to struggle when extending contigs beyond reference sequences. This decrease in coverage can cause the apparent overall coverage of the contig to be artificially low. GSEC.py removes these low-coverage positions by removing bases from both ends and/or retaining only a certain percentage of positions from the middle of the contig, and calculates a corrected coverage for each contig.