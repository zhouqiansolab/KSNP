cmake_minimum_required(VERSION 3.10)
project(ksnp)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-g -O3)
include_directories(htslib-1.15.1/include)
link_directories(htslib-1.15.1/lib)

add_executable(ksnp main.cpp
        vcf_reader.cpp vcf_reader.h
        snp_dbg.cpp snp_dbg.h
        phased_block.cpp phased_block.h)

target_link_libraries(ksnp -lhts -lz)