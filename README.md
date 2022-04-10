# KSNP: k-mer based haplotype assembly

## Download and install

**Dependencies**
- cmake >= 3.10
- gcc >= 6.4
- zlib >= 1.2
- htslib 1.15.1

```
git clone https://github.com/zhouqiansolab/KSNP-k-mer-based-haplotype-assembly.git
cd KSNP-k-mer-based-haplotype-assembly
./install.sh
```
If the installation is successful, the build subdirectory will contain the `ksnp` binary.
## Usage
```
ksnp -k <k-mer size> -b <BAM file> -v <VCF file> -o <output file>
  ## K-mer size supports INT value from 2 to 10 (default value is 2)
  ## The in.bam is sorted and indexed. The chr_name is not required when the in.bam only includes one chromosome.
  ## The VCF file should contain SNP positions on one chromosome. Please split the VCF file into small files according to the chr_name. The VCF file can be in be in compressed gz format or decompressed format.
  ## The output file includes all of the variants in input VCF file with phased information.
  ## Without specifying the output file, the phased VCF is print to stdout.
  
```
## Contact
```
zhouqian_solab@163.com
```
