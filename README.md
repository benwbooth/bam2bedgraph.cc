bam2bedgraph
==============

Convert .bam alignment files to bedgraph or bigwig format. 

Requires:

- a modern C++ compiler with C++14 features
- bamtools (https://github.com/pezmaster31/bamtools)
- boost

To build:

Run ```make```

Help documentation:

```
Usage: ./bam2bedgraph [ options ] bamfile
  --help                  produce help
  --bamfile FILE          Input BAM filename
  --split                 Use CIGAR string to split alignment into separate
                          exons (default)
  --nosplit
  --autostrand ANNOT_FILE Attempt to determine the strandedness of the input
                          data using an annotation file. Must be a BAM file.
  --strand [TYPE]         Split output bedgraph by strand: Possible values: u s
                          r uu us ur su ss sr ru rs rr, first char is read1,
                          second is read2, u=unstranded, s=stranded, r=reverse
  --read                  Split output bedgraph by read number
  --noread                (default)
  --zero                  Pad output bedgraph with zeroes
  --nozero                (default)
  --fixchr                Transform chromosome names to be UCSC-compatible
  --nofixchr              (default)
  --paired                Only output paired read alignments
  --nopaired              (default)
  --proper                Only output proper-paired read alignments
  --noproper              (default)
  --primary               Only output primary alignments
  --noprimary             (default)
  --bigwig                Output bigwig files (requires bedGraphToBigWig in
                          $PATH)
  --nobigwig              (default)
  --uniq                  Keep only unique alignments (NH:i:1)
  --nouniq                (default)
  --out FILE              Output file prefix
  --trackline             Output a UCSC track line (default)
  --notrackline
  --trackname TRACKNAME   Name of track for the track line
  --5prime                Only report the 5' ends
  --no5prime              (default)
  --3prime                Only report the 3' ends
  --no3prime              (default)
  --edge LENGTH (=1)      Length of edges for --5prime and --3prime options
                          (default: 1)
```
