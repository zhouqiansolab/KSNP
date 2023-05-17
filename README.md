# KSNP: k-mer based haplotype assembly

## Download and install

**Dependencies**
- cmake >= 3.10
- gcc >= 6.4
- htslib 1.15.1

KSNP can be complied on Linux or MAC OS.
Before KSNP installation, keep sure the dependency htslib has been correctly installed.
If it is not in default system path, please add its directory to environmental variable `LD_LIBRARY_PATH`.

```
git clone https://github.com/zhouqiansolab/KSNP
cd KSNP; mkdir build; cd build
cmake ..; make
```

If the installation is successful, the build subdirectory will contain the `ksnp` binary.
## Usage
```
ksnp -k <k-mer size> -b <BAM> -r <FASTA> -v <VCF> -o <output file>
  ## The k-mer size supports INT value from 2 to 5. Default value is set to 2.
  ## BAM file contains the aligned reads. It must be sorted and indexed.
  ## FASTA file is the reference sequence used for reads alignment and variants calling. It should be indexed by 'samtools faidx'
  ## VCF file contains the heterozygous variants to phase. It should be properly filtered before phasing.
  ## The output file keeps all varinats in input VCF file but with phased information. Without specifying it, the results will be print to stdout.
  ## Sample usage: ksnp -b aln.bam -r ref.fa -v variants.vcf -o phased.vcf
```
## For testing
```
cd test
ksnp -k 2 -b aln.bam -r ref.fa -v variants.vcf -o test_ksnp.vcf
  ##expect_output.vcf is an expected output file that can be used for comparing with the test results.
```
## Contact
```
zhouqian_solab@163.com
```
