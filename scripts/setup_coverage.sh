# install screen
sudo apt-get update
sudo apt-get install screen
# install parallel
sudo apt-get install parallel
# install samtools
mkdir ~/repos
cd ~/repos
git clone https://github.com/samtools/samtools.git
git clone git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
sudo apt-get update
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
autoheader
autoconf -Wno-syntax
./configure
make
make install
# setup volume vdc for data
# data storage mounted at /pvol
sudo mkfs -t ext4 /dev/vdc
sudo mkdir /pvol
sudo mount /dev/vdc /pvol -t auto
dufo mount --all
sudo chown ubuntu:ubuntu /pvol
# copy data to server
# example destination: ubuntu@123.123.255.255:/pvol
scp -rp <src> <destination>
# setup volume vde for results
# result storage mounted at /evol
sudo mkfs -t ext4 /dev/vde
sudo mkdir /evol
sudo mount /dev/vde /evol -t auto
dufo mount --all
sudo chown ubuntu:ubuntu /evol
# setup directories
mkdir /evol/pvol
mkdir /evol/pvol/Tamsin
mkdir /evol/pvol/Tamsin/BAM
# launch screen and run samtools
ls /pvol/Tamsin/BAM/*.bam | parallel "samtools depth {} > /evol{}.txt" 
# run script to count coverage histogram
cd ~/repos/phylonco/scripts
chmod +x count_coverage.sh
./count_coverage.sh /evol/pvol/Tamsin/BAM/*.txt  