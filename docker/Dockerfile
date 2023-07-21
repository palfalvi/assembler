# Set the base image to debian bionic
FROM ubuntu:20.04
USER root
ENV TZ=Asia/Tokyo
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR $HOME/miniconda3
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# File Author / Maintainer
MAINTAINER Gergo Palfalvi <palfalvi.gergo@gmail.com>

RUN apt-get update && apt-get upgrade -y &&\
    apt-get install -y sudo \
    autoconf \
    automake \
    bamtools \
    bedtools \
    bcftools \
    bowtie2 \
    build-essential \
    cdbfasta \
    cmake \
    cpanminus \
    csh \
    default-jre \
    diamond-aligner \
    git \
    gnuplot \
    last-align \
    libatlas-base-dev \
    libbamtools-dev \
    libboost-all-dev \
    libdbd-mysql-perl \
    libfile-which-perl \
    libgsl-dev \
    libhts-dev \
    liblpsolve55-dev \
    libmysql++-dev \
    libparallel-forkmanager-perl \
    libpthread-stubs0-dev \
    libsqlite3-dev \
    libsuitesparse-dev \
    libyaml-perl \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    make \
    mc \
    mummer \
    python3-pip \
    rsync \
    samtools \
    subversion \
    software-properties-common \
    time \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    jellyfish




# Start working in /opt
WORKDIR /opt
ENV PASA_CO eebd6223e1a7aa008c8aa63c8969296a0285bd8b

# BRAKER installation
# HTSlib, braker, GeneMark-ES, Augustus
# GENEMARK-EX need a new link in every week or so
# Perl modules
RUN cpanm \
  File::Spec::Functions \
  Hash::Merge \
  List::Util \
  Logger::Simple \
  Module::Load::Conditional \
  Parallel::ForkManager \
  POSIX \
  Scalar::Util::Numeric \
  YAML \
  Math::Utils \
  MCE::Mutex \
  DBD::mysql \
  Thread::Queue \
  URI::Escape \
  threads &&\
  sudo ln -s /usr/bin/python3 /usr/bin/python &&\
  wget -c https://www.python.org/ftp/python/3.9.1/Python-3.9.1.tgz -O - | tar -xz &&\
  python3 -m pip install matplotlib numpy pandas biopython mikado &&\
  git clone https://github.com/gatech-genemark/ProtHint.git &&\
  #wget -c https://github.com/gatech-genemark/ProtHint/releases/download/v2.4.0/ProtHint-2.4.0.tar.gz -O - | tar -xz &&\
  wget -c https://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz -O - | tar -xz && mv gth-1.7.3-Linux_x86_64-64bit gth &&\
  wget -c http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz -O - | tar -xz &&\
  wget -c https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && tar -xf htslib-1.11.tar.bz2 && cd htslib-1.11 && ./configure && make && sudo make install && cd .. &&\
  #wget -c https://github.com/Gaius-Augustus/BRAKER/archive/v2.1.5.tar.gz -O - | tar -xz &&\
  git clone https://github.com/Gaius-Augustus/BRAKER.git &&\
  wget -c http://topaz.gatech.edu/GeneMark/tmp/GMtool_tJdsf/gmes_linux_64.tar.gz -O - | tar -xz &&\
  git clone https://github.com/Gaius-Augustus/Augustus && cd ./Augustus && sed -i 's/	cd bam2wig; make/	#cd bam2wig; make/g' ./auxprogs/Makefile && make clean && make && make install && cd .. &&\
  git clone https://github.com/splicebox/psiclass.git && cd psiclass && sudo make && cd .. &&\
  wget -c https://github.com/ruolin/strawberry/releases/download/v1.1.2/strawberry && chmod +x strawberry &&\
  git clone https://github.com/gpertea/stringtie && cd stringtie && make release && cd .. &&\
  wget -c https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.11.0/trinityrnaseq-v2.11.0.FULL.tar.gz -O - | tar -xz && cd trinityrnaseq-v2.11.0 && make && make install && cd .. &&\
  wget -c https://github.com/alexdobin/STAR/archive/2.7.7a.tar.gz -O - | tar -xz && cd STAR-2.7.7a/source && make STAR && cd ../.. &&\
  git clone https://github.com/lh3/minimap2 && cd minimap2 && make && cd .. &&\
  git clone https://github.com/iMetOsaka/UNAGI && chmod a+x UNAGI/unagi &&\
  sed -i 's/minimap_path=.\/tools\/minimap2\/minimap2/minimap_path=minimap2/g' /opt/UNAGI/app/conf.ini &&\
  sed -i 's/samtools_path=.\/tools\/samtools\/samtools/samtools_path=samtools/g' /opt/UNAGI/app/conf.ini &&\
  sed -i 's/bedtools_path=.\/tools\/bedtools\/bin\/bedtools/bedtools_path=bedtools/g' /opt/UNAGI/app/conf.ini &&\
  wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && chmod 755 ./blat &&\
  wget -c http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz -O - | tar -xz && cd ./fasta-36.3.8g/src && make -f ../make/Makefile.linux_sse2 all && cp ../bin/fasta36 /usr/local/bin/fasta &&\
  git clone https://github.com/PASApipeline/PASApipeline.git && cd PASApipeline && git checkout $PASA_CO && git submodule init && git submodule update && make &&\
  cpanm File::HomeDir &&\
  chmod -R a+w /opt/Augustus/config/ &&\
  rm -rf /opt/Augustus/config && rm -rf /opt/augustus-3.4.0/config &&\
  cp /opt/ProtHint/bin/* /opt/gth/bin/ && cp -r /opt/ProtHint/dependencies/ /opt/gth/ &&\
  wget -c https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz -O - | tar -xz && mv salmon-latest_linux_x86_64 salmon

  #tama
  #taco

  #busco
RUN wget -c https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz -O - | tar -xz
RUN git clone https://github.com/GenomeRIK/tama

RUN sudo apt-get install build-essentials &&\
    wget -c https://github.com/alekseyzimin/masurca/releases/download/v4.0.1/MaSuRCA-4.0.1.tar.gz -O - | tar -xz && cd MaSuRCA-4.0.1 && ./install.sh


ENV PATH=/opt/MaSuRCA-4.0.1/:/opt/tama:/opt/tama/tama_go/format_converter/:/opt/salmon-latest_linux_x86_64/bin:/opt/trinityrnaseq-v2.11.0/:/opt/gth/bin/:/opt/ProtHint/bin:/opt/gmes_linux_64/:/opt/gmes_linux_64/bin/:/opt/ProtHint/dependencies/:/opt/diamond-linux64/:/opt/PASApipeline/:/opt/PASApipeline/bin/:/opt/UNAGI/:/opt/minimap2/:/opt/psiclass/:/opt/augustus-3.4.0/auxprogs/bam2hints:/opt/augustus-3.4.0/bin:/opt/augustus-3.4.0/scipts:/opt/canu/Linux-amd64/bin:/opt/BRAKER/scripts/:$PATH
ENV AUGUSTUS_BIN_PATH=/opt/Augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts
ENV AUGUSTUS_CONFIG_PATH=/opt/Augustus/config
ENV BSSMDIR=/opt/gth/bin/bssm
ENV GTHDATADIR=/opt/gth/bin/gthdata
ENV GTHNOFLOCK=yes
