# KSNP: k-mer based haplotype assembly

## Download and install

**Dependencies**
- cmake >= 3.10
- gcc >= 6.4
- zlib
- make

```
git clone https://github.com/zhouqiansolab/KSNP-k-mer-based-haplotype-assembly.git
mkdir build; cd build
cmake ..; make
```
## Usage
```
samtools view in.bam chr_name | ksnp -k <k-mer size> -v <VCF file> -o <output file>
  ## K-mer size supports INT value from 2 to 10 (default value is 2)
  ## The in.bam is sorted and indexed. The chr_name is not required when the in.bam only includes one chromosome.
  ## The VCF file should contain SNP positions on one chromosome. Please split the VCF file into small files according to the chr_name. The VCF file can be in be in compressed gz format or decompressed format.
  ## The output file includes all of the variants in input VCF file with phased information.
  
```
## Contact
```
zhouqian_solab@163.com
```
