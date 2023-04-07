# KSNP: k-mer based haplotype assembly

## Download and install

**Dependencies**
- cmake >= 3.10
- gcc >= 6.4
- htslib 1.15.1

Before KSNP installation, keep sure the dependency htslib (https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2) 
has been successfully installed.

```
git clone https://github.com/zhouqiansolab/KSNP
cd KSNP; mkdir build; cd build
cmake .. -DHTSLIB=#htslib installtion path 
make
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

## Contact
```
zhouqian_solab@163.com
```
