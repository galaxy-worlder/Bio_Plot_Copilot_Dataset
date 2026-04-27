#
# Linux packages dependencies
# sudo add-apt-repository ppa:ubuntugis/ppa
# sudo apt-get update
# sudo apt install -y libfontconfig1-dev build-essential checkinstall gdal-bin libgdal-dev
# sudo apt install -y libssl-dev zlib1g-dev libxml2-dev curl libcurl4-openssl-dev libfribidi-dev libharfbuzz-dev  libjpeg8-dev libtiff-dev libpng-dev libudunits2-dev libbz2-dev
# sudo apt install -y liblapack-dev libopenblas-dev gfortran
# sudo apt install -y python python3-venv python3-numpy python3-pip tensorflow-gpu
# sudo apt install -y ncbi-blast+ bowtie2 samtools bedtools
# wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2.tar.gz
# tar -zxvf cmake-3.23.2.tar.gz
# cd cmake-3.23.2
# ./bootstrap
# make
# sudo make install


install.packages(c( "usethis", "pkgdown", "rcmdcheck", "roxygen2", "rversions", "urlchecker", "devtools"))
# On windows install Rtools package seprately
# https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html
devtools::install_cran(c("readr", "dplyr", "data.table", "fitdistrplus", "rstatix", "bedr", "reshape2", "ggseqlogo", "randomColorR", "foreach"))
devtools::install_cran(c("igraph", "baseline", "smoother", "dbscan", "units"))
devtools::install_cran(c("ggvenn", "randomcoloR", "ggbeeswarm", "ggpmisc", "ggridges", "ggrepel", "units", "gridpattern", "ggpattern", "ggprism", "ggpubr"))
devtools::install_cran(c("reticulate", "tensorflow", "keras"))
devtools::install_bioc(c("Rsamtools", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "rtracklayer", "Biostrings", "ComplexHeatmap"))
devtools::install_deps("breaktools/")

#
# Download genomes
#
# python -m pip install genome_downloader
# python -m genome_downloader mm10 genomes


#
# Create Python environment for Keras and Tensorflow
# To run Tensorboard (for validating model training) use command
#
# tensorboard --logdir logs/
#
reticulate::virtualenv_create("r-tensorflow", python="/usr/bin/python3")
reticulate::py_install("tensorflow-gpu", envname="r-tensorflow")
tensorflow::install_tensorflow(envname="r-tensorflow")
keras::install_keras(envname="r-tensorflow")

# Install appropriate versions of cuDNN and CUDA. On WSL2 you will need to install Docker Desktop (https://www.docker.com/products/docker-desktop/)
# Version table: https://www.tensorflow.org/install/source#gpu
 tensorflow_gpu-2.11.0

# Follow this instruction
# https://medium.com/mlearning-ai/how-to-install-the-nvidia-cuda-toolkit-11-2-a94d86a45d38

wget -qO - https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/nvidia.gpg
wget -qO - https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/3bf863cc.pub | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/nvidia-wsl2.gpg
sudo sh -c 'echo "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" > /etc/apt/sources.list.d/cuda.list'

sudo apt-get update
sudo apt-get --yes install cuda-toolkit-11-2

#
# Install CUDA
#
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-wsl-ubuntu.pin
sudo mv cuda-wsl-ubuntu.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.2.0/local_installers/cuda-repo-wsl-ubuntu-11-2-local_11.2.0-1_amd64.deb
sudo dpkg -i cuda-repo-wsl-ubuntu-11-2-local_11.2.0-1_amd64.deb
#sudo cat /var/cuda-repo-wsl-ubuntu-11-2-local/7fa2af80.pub | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/nvidia-wsl2.gpg
wget -qO - https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/3bf863cc.pub | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/nvidia-wsl2.gpg
sudo add-apt-repository 'deb https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/ /'
sudo apt-get update
sudo apt-get -y install cuda-11-2
# Might need to fix conflicts
sudo apt --fix-broken install -o Dpkg::Options::="--force-overwrite"

#
# Install cuDNN
#
# Download cuDNN version 8.1 from https://developer.nvidia.com/rdp/cudnn-archive
sudo dpkg -i libcudnn8_8.1.1.33-1+cuda11.2_amd64.deb

# Test installation
git clone https://github.com/nvidia/cuda-samples && cd cuda-samples/ && make && cd ../../../
cuda-samples/Samples/1_Utilities/deviceQuery



#sudo cp /var/cuda-repo-wsl-ubuntu-11-2-local/cuda-*-keyring.gpg /usr/share/keyrings/

#wget -qO - https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/nvidia.gpg
wget -qO - https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/3bf863cc.pub | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/nvidia-wsl2.gpg
#sudo cp /var/cuda-repo-wsl-ubuntu-11-2-local/cuda-*-keyring.gpg /usr/share/keyrings/




sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/3bf863cc.pub

wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-wsl-ubuntu.pin
sudo mv cuda-wsl-ubuntu.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/3bf863cc.pub
wget https://developer.download.nvidia.com/compute/cuda/11.2.0/local_installers/cuda-repo-wsl-ubuntu-11-2-local_11.2.0-1_amd64.deb
#sudo add-apt-repository 'deb https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/ /'
sudo dpkg -i cuda-repo-wsl-ubuntu-11-2-local_11.2.0-1_amd64.deb
sudo cp /var/cuda-repo-wsl-ubuntu-11-2-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda

# Check CUDA version
#  dpkg -l | grep -i nvidia
#
# Test GPU in tensorflow
from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())

#
# OS-X
# ===============================
# Install blast, bowtie2 and dustmasker if you don't have them. The easiest way is to use Anaconda installer. You might need to add Anaconda binaries folder to your PATH environment variable.
# If you are using OSX+Rstudion also add PATH to ~/.Renviron file:
#
# Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ":/path/to/anaconda/bin"))
#
# https://anaconda.org/bioconda/blast
# https://anaconda.org/bioconda/bowtie2
# https://anaconda.org/bioconda/blat
#

#
# Windows Native (Doesn't work)
# ===============================
#
# On Windows you will have to install Perl before continuing with bowtie2 instllation.
# conda install -c conda-forge perl
#
# On Windows you will need to download BLAST (contains dustmasker) and bowtie2 manually from
# https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
# https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/
#
# The problem is that Perl (used by bowtie2) can't handle spaces in paths properly and has problems accessing "c:/Program Files"
#

#
# Windows WSL2 (works)
# ===============================
#
# On Windows you will have to install Perl before continuing with bowtie2 instllation.
# conda install -c conda-forge perl
#
# On Windows you will need to download BLAST (contains dustmasker) and bowtie2 manually from
# https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
# https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/
#
# The problem is that Perl (used by bowtie2) can't handle spaces in paths properly and has problems accessing "c:/Program Files"
#

dir.create("reports", recursive=T, showWarnings=F)
dir.create("genomes", recursive=T, showWarnings=F)
dir.create("tmp", recursive=T, showWarnings=F)

#
# Download TLX files from NCBI
# https://drive.google.com/u/0/uc?id=1gUVUePDl89nnYBTb4ZjL03l8NaLdytdB&export=download
#
unzip("data/data.zip", exdir="data/")
dir.create("data/TLX", recursive=T, showWarnings=F)
file.copy(Sys.glob("data/TLX_paper/*"), "data/TLX", overwrite=T, recursive=T)
file.copy(Sys.glob("data/TLX_public/*"), "data/TLX", overwrite=T, recursive=T)
# file.copy(Sys.glob("reports/00-upload_ncbi/TLX/*.tlx"), "data/TLX", overwrite=TRUE)


