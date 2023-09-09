# KSNP: k-mer based haplotype assembly

## Download and install

### Dependencies
- cmake >= 3.10
- gcc >= 6.4
- htslib >= 1.15
- zlib >= 1.2

KSNP can be installed on Linux or MAC OS.

### Install from conda

We recommend users install KSNP using conda command.

```
# New environment is recommended for avoiding dependency conflicts.
conda create -n ksnp_env
conda install -c conda-forge -c bioconda ksnp
```

### Install from source

Before compile the KSNP code locally, keep sure the dependency htslib has been correctly installed.
Please add the subdirectory containing `libhts.so` to the default environmental variable for
searching dynamic-link libraries (such as `LD_LIBRARY_PATH` environment variable in Linux OS) 
or use the pre-compile option `-DHTSLIB=` when running `cmake` to specify the location.

```
git clone https://github.com/zhouqiansolab/KSNP.git
cd KSNP; mkdir build; cd build
cmake -DHSTLIB=<directory of libhts.so> ..
make
```
If the installation is successful, the build subdirectory will contain the executable file `ksnp`.

## Usage
```
ksnp -k <k-mer size> -b <BAM> -r <FASTA> -v <VCF> -o <output file>  [-c chr]
  ## The k-mer size supports INT value from 2 to 5. Default value is set to 2.
  ## BAM file contains the aligned reads. It must be sorted and indexed.
  ## FASTA file is the reference sequence used for reads alignment and variants calling. It should be indexed by 'samtools faidx'
  ## Users may specify a chromosome using option -c <chr-name> to manually parallize KSNP across chromosomes.
  ## VCF file contains the heterozygous variants to phase. It should be properly filtered before phasing.
  ## The output file keeps all varinats in input VCF file but with phased information. Without specifying it, the results will be print to stdout.
  ## Sample usage: ksnp -b aln.bam -r ref.fa -v variants.vcf -o phased.vcf
  ## -c allows users to specify individual chromosome for phasing, enabling chromosome-level parallel processing without the need for splitting input files.
```
## For testing
```
cd test
ksnp -k 2 -b aln.bam -r ref.fa -v variants.vcf -o test_ksnp.vcf
```
expect_output.vcf is an expected output file that can be used for comparing with the test results.

## Contact
```
zhouqian_solab@163.com
```
