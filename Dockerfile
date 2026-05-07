ARG BASE_IMAGE_PLATFORM=linux/amd64
FROM --platform=${BASE_IMAGE_PLATFORM} ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# ubuntu environment setup
RUN apt-get update -y && apt-get upgrade -y &&\
	apt-get install -y \
	wget \
	ca-certificates \
	unzip \
	tree \
	build-essential checkinstall zlib1g-dev \
	libssl-dev \
	git \
	dc \
	ffmpeg \
	libsm6 \
	libxext6 \
	python3 \
	python3-pip \
	python3-venv \
        bc

RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2.tar.gz &&\
	tar -xvzf cmake-3.23.2.tar.gz &&\
	rm cmake-3.23.2.tar.gz &&\
	cd cmake-3.23.2 &&\
	./bootstrap &&\
	make &&\
	make install

# create and switch to working directory
RUN mkdir /aida/
WORKDIR /aida/

# NiftyReg preparation and installation
RUN apt-get update && apt-get install -y --no-install-recommends \
    libpng-dev libjpeg-dev libtiff-dev \
 && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /aida/NiftyReg/niftyreg_source /aida/NiftyReg/niftyreg_build /aida/NiftyReg/niftyreg_install
WORKDIR /aida/NiftyReg

RUN git clone https://github.com/KCL-BMEIS/niftyreg.git niftyreg_source && \
    cd niftyreg_source && \
    git reset --hard 83d8d1182ed4c227ce4764f1fdab3b1797eecd8d

WORKDIR /aida/NiftyReg/niftyreg_build
RUN cmake -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_INSTALL_PREFIX=/aida/NiftyReg/niftyreg_install \
          -D CMAKE_C_COMPILER=/usr/bin/gcc \
          ../niftyreg_source && \
    make -j"$(nproc)" && \
    make install

ENV NIFTYREG_INSTALL=/aida/NiftyReg/niftyreg_install
ENV PATH="${PATH}:${NIFTYREG_INSTALL}/bin"
ENV LD_LIBRARY_PATH="/aida/dsi_studio_ubuntu2204/dsi-studio:/usr/local/lib:/usr/lib:/lib:${NIFTYREG_INSTALL}/lib"
WORKDIR /aida

# download DSI studio
# https://github.com/frankyeh/DSI-Studio/releases/download/2023.12.06/dsi_studio_ubuntu2204.zip
RUN wget https://github.com/frankyeh/DSI-Studio/releases/download/2025.04.16/dsi_studio_ubuntu2204.zip &&\
	unzip dsi_studio_ubuntu2204.zip -d /aida/dsi_studio_ubuntu2204 &&\
	rm dsi_studio_ubuntu2204.zip

# Install ANTs (if no 22.04 binary, keep 18.04 version)
# https://github.com/ANTsX/ANTs/releases/download/v2.6.2/ants-2.6.2-ubuntu-22.04-X64-gcc.zip
RUN wget https://github.com/ANTsX/ANTs/releases/download/v2.6.2/ants-2.6.2-ubuntu-22.04-X64-gcc.zip &&\
	unzip ants-2.6.2-ubuntu-22.04-X64-gcc.zip -d ants-2.6.2 &&\
	rm ants-2.6.2-ubuntu-22.04-X64-gcc.zip
ENV PATH=$PATH:/aida/ants-2.6.2/ants-2.6.2/bin
RUN N4BiasFieldCorrection --version

# Python setup (Default for Ubuntu 22.04 is Python 3.10)
ENV VIRTUAL_ENV=/opt/env
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN python3 -m pip install --upgrade pip setuptools
COPY requirements.txt requirements.txt
RUN pip install --upgrade pip &&\
	pip install -r requirements.txt

# installation of FSL 5.0.11 with modified installer 
COPY fslinstaller_mod.py ./
RUN python3 fslinstaller_mod.py -V 5.0.11

# Configure environment
ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ

# Configure environment
ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH=${FSLDIR}/bin:${PATH}

# copy bin/ and lib/ from AIDAmri into image
COPY bin/ bin/
RUN chmod u+x bin/3.2_DTIConnectivity/dsi_main.py
ENV PATH=/aida/bin:/aida/bin/3.2_DTIConnectivity:$PATH
RUN cp bin/3.2_DTIConnectivity/dsi_main.py dsi_main
COPY lib/ lib/
# make install_immv executable and run it
RUN chmod +x /aida/bin/install_immv.sh
RUN /aida/bin/install_immv.sh
RUN echo "/aida/dsi_studio_ubuntu2204/dsi-studio/dsi_studio" > /aida/bin/3.2_DTIConnectivity/dsi_studioPath.txt
RUN test -x /aida/dsi_studio_ubuntu2204/dsi-studio/dsi_studio

RUN pip install dipy scikit-learn
RUN pip install fslpy
RUN wget -O /aida/bin/bet4animal "https://git.fmrib.ox.ac.uk/fsl/bet2/-/raw/master/bet4animal?ref_type=heads&inline=false" && \
    chmod +x /aida/bin/bet4animal
RUN cd /aida/bin && bash -c "/aida/bin/install_immv.sh"
