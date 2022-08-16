FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# ubuntu environment setup
RUN apt-get update -y && apt-get upgrade -y &&\
	apt-get install -y \
	wget \
	unzip \
	build-essential checkinstall zlib1g-dev \
	libssl-dev \
	git \
	dc \
	ffmpeg \
	libsm6 \
	libxext6

RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2.tar.gz &&\
	tar -xvzf cmake-3.23.2.tar.gz &&\
	rm cmake-3.23.2.tar.gz &&\
	cd cmake-3.23.2 &&\
	./bootstrap &&\
	make &&\
	make install

# create and switch to working directory
RUN mkdir aida/
WORKDIR /aida/

# Python setup
RUN apt install -y python3 python3-pip &&\
	python3 -m pip install --user --upgrade pip &&\
	apt-get install -y python3-venv
ENV VIRTUAL_ENV=/opt/env
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN	python3 -m pip install --upgrade setuptools
COPY requirements.txt requirements.txt
RUN pip install --upgrade pip &&\
	pip install -r requirements.txt

# installation of FSL 5.0.11 with modified installer 
# (disabling interactive allocation query)
COPY fslinstaller_mod.py ./
#RUN python3 fslinstaller_mod.py -V 5.0.11
RUN python3 fslinstaller_mod.py

# Configure environment
ENV FSLDIR=/usr/local/fsl
RUN . ${FSLDIR}/etc/fslconf/fsl.sh
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH=${FSLDIR}/bin:${PATH}
RUN export FSLDIR PATHs

# Niftyreg preparation and installation
RUN mkdir -p NiftyReg/niftyreg_source/
WORKDIR /aida/NiftyReg
RUN git clone git://git.code.sf.net/p/niftyreg/git niftyreg_source &&\
cd niftyreg_source &&\
git reset --hard 83d8d1182ed4c227ce4764f1fdab3b1797eecd8d &&\
 	mkdir niftyreg_install niftyreg_build && cd .. &&\
 	cmake niftyreg_source &&\
 	cmake -D CMAKE_BUILD_TYPE=Release niftyreg_source &&\
 	cmake -D CMAKE_INSTALL_PREFIX=niftyreg_source/niftyreg_build niftyreg_source &&\
 	cmake -D CMAKE_C_COMPILER=/usr/bin/gcc-7 niftyreg_source &&\
 	make && make install
RUN export NIFTYRREG_INSTALL=../niftyreg_install
ENV PATH=${PATH}:${NIFTYREG_INSTALL}/bin
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NIFTYREG_INSTALL}/lib
RUN export PATH && export LD_LIBRARY_PATH

WORKDIR /aida
# download DSI studio
RUN wget https://github.com/frankyeh/DSI-Studio/releases/download/2022.08.03/dsi_studio_ubuntu_1804.zip &&\
	unzip dsi_studio_ubuntu_1804.zip -d dsi_studio_ubuntu_1804 &&\
	rm dsi_studio_ubuntu_1804.zip

# copy bin/ and lib/ from AIDAmri into image
COPY bin/ bin/
COPY lib/ lib/
RUN echo "/aida/bin/dsi_studio_ubuntu_1804/dsi-studio/dsi_studio" > bin/3.2_DTIConnectivity/dsi_studioPath.txt