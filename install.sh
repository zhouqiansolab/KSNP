wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
tar -xjf htslib-1.15.1.tar.bz2
rm htslib-1.15.1.tar.bz2

cd htslib-1.15.1 || exit
autoreconf -i
./configure --prefix=$(pwd)
make
make install

if [ ! -d lib ] || [ ! -e include ]; then
	echo -e 'htslib installation failed. Please check your network.'
	exit
fi

cd ..
mkdir build
cd build || exit
cmake ..; make