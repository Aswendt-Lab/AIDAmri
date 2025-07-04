![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/Aswendt-Lab/AIDAmri/docker-image.yml) ![Static Badge](https://img.shields.io/badge/Docker_image-11.97_GB-blue) [![Static Badge](https://img.shields.io/badge/data_structure-BIDS-yellow)](https://bids.neuroimaging.io/news.html) [![Static Badge](https://img.shields.io/badge/Niftyreg-CBSI-orange)](https://github.com/KCL-BMEIS/niftyreg) [![Static Badge](https://img.shields.io/badge/DSI--Studio-2023-orange)](https://dsi-studio.labsolver.org/) [![Static Badge](https://img.shields.io/badge/FSL-5.0.11-orange)]([https://dsi-studio.labsolver.org/](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)) ![Static Badge](https://img.shields.io/badge/Python-3.7-orange)

[1.2]: http://i.imgur.com/wWzX9uB.png
[1]: http://www.twitter.com/AswendtMarkus
<!--social icon from https://github.com/carlsednaoui/gitsocial -->

<img align="left" src="https://github.com/maswendt/AIDAmri/blob/master/AIDA_Logo.png" width="120">
<h1>AIDA<i>mri</i></h1>

Atlas-based Imaging Data Analysis Pipeline (AIDA) for structural and functional MRI of the mouse brain
<br/>
## Key Features of AIDA<em>mri</em>

1. **Automated Preprocessing**  
   It performs tasks like image re-orientation, bias-field correction, and brain extraction with minimal user input required.

2. **Atlas-Based Registration**  
   AIDA<em>mri</em> uses the **Allen Mouse Brain Reference Atlas** for accurate region-based analysis of MRI data, allowing researchers to compare results across different studies efficiently. A modified atlas version with larger labels to better match MRI resolution is provided. Users can define specific **regions of interest (ROIs)** for analysis, such as stroke lesions.

3. **Modular Design**  
   The pipeline is developed in Python, making it cross-platform and open-source, allowing for easy integration and modification.

4. **Validation**  
   The pipeline was validated with different MRI datasets, including those involving stroke models, demonstrating its robustness even in the presence of significant brain deformations.

5. **Functional and Structural Connectivity Analysis**  
   The output of the pipeline includes connectivity matrices that can be used for further analysis of brain network changes in health and disease.

<p align="center">
  <img src="https://github.com/maswendt/AIDAmri/blob/master/AIDAmri_drawing.png" style="max-width: 100%; height: auto;">
</p>

Pipeline overview from [Pallast et al.](https://doi.org/10.3389/fninf.2019.00042)



## Version history

[Information latest Version 2.0](https://github.com/maswendt/AIDAmri/releases/tag/v2.0)

[**Manual**](https://github.com/maswendt/AIDA/blob/master/manual.pdf)

[Information about Version 1.2 (Docker stable release)](https://github.com/maswendt/AIDAmri/releases/tag/v1.2)
<br/>
[Information about Version 1.1.1 (Docker pre-release)](https://github.com/maswendt/AIDAmri/releases/tag/1.1.1)
<br/>
[Information about Version 1.1 (Stable)](https://github.com/maswendt/AIDAmri/releases/tag/v1.1)
<br/>
[Information about Version 1.0](https://github.com/maswendt/AIDAmri/releases/tag/v1.0)

<h3><b>Important note: read this before you install AIDAmri for the first time</h3></b>

We fully moved to the containerized version of AIDAmri via [Docker](https://docs.docker.com/get-docker/). All information can be found in the manual above. Please report issues and bugs directly in the issue section of this repository or at gitter (Link below in the contact section).

<details>
<summary>Note for Linux Users</summary></b>
When building the AIDAmri Docker image on a Linux system, you may encounter warning messages related to undefined environment variables. Specifically, the following warnings may appear:

```
3 warnings found (use docker --debug to expand):
 - UndefinedVar: Usage of undefined variable '$NIFTYREG_INSTALL' (line 44)
 - UndefinedVar: Usage of undefined variable '$NIFTYREG_INSTALL' (line 43)
 - UndefinedVar: Usage of undefined variable '$LD_LIBRARY_PATH' (line 44)
```

These warnings indicate that certain environment variables referenced in the Dockerfile are either not defined or not properly initialized during the build process.

### Recommended Solution:
To ensure compatibility and suppress these warnings, edit the `Dockerfile` in the AIDAmri directory. Replace lines **30–93** with the corrected version below, ensuring that all relevant environment variables are explicitly declared and exported: 

```
# NiftyReg preparation and installation
RUN apt update && apt install -y gcc-7 g++-7

RUN mkdir -p /aida/NiftyReg/niftyreg_source /aida/NiftyReg/niftyreg_build /aida/NiftyReg/niftyreg_install

WORKDIR /aida/NiftyReg

RUN git clone git://git.code.sf.net/p/niftyreg/git niftyreg_source && \
    cd niftyreg_source && \
    git reset --hard 83d8d1182ed4c227ce4764f1fdab3b1797eecd8d

WORKDIR /aida/NiftyReg/niftyreg_build

RUN cmake -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_INSTALL_PREFIX=/aida/NiftyReg/niftyreg_install \
          -D CMAKE_C_COMPILER=/usr/bin/gcc-7 \
          ../niftyreg_source && \
    make -j$(nproc) && \
    make install

ENV NIFTYREG_INSTALL=/aida/NiftyReg/niftyreg_install
ENV PATH="${PATH}:${NIFTYREG_INSTALL}/bin"
ENV LD_LIBRARY_PATH="/usr/local/lib:/usr/lib:/lib:${NIFTYREG_INSTALL}/lib"

WORKDIR /aida
# download DSI studio
RUN wget https://github.com/frankyeh/DSI-Studio/releases/download/2023.07.08/dsi_studio_ubuntu1804.zip &&\
	unzip dsi_studio_ubuntu1804.zip -d dsi_studio_ubuntu1804 &&\
	rm dsi_studio_ubuntu1804.zip

# Python setup
RUN apt install -y python3.7 python3-pip &&\
	python3 -m pip install --user --upgrade pip &&\
	apt-get install -y python3.7-venv &&\
	apt clean &&\
	rm -rf /var/lib/apt/lists/*
ENV VIRTUAL_ENV=/opt/env
RUN python3.7 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN	python3 -m pip install --upgrade setuptools
COPY requirements.txt requirements.txt
RUN pip install --upgrade pip &&\
	pip install -r requirements.txt

# installation of FSL 5.0.11 with modified installer 
# (disabling interactive allocation query)
COPY fslinstaller_mod.py ./
RUN python3 fslinstaller_mod.py -V 5.0.11

# Configure environment
ENV FSLDIR=/usr/local/fsl
RUN . ${FSLDIR}/etc/fslconf/fsl.sh
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH=${FSLDIR}/bin:${PATH}
RUN export FSLDIR PATHs


# copy bin/ and lib/ from AIDAmri into image
COPY bin/ bin/
RUN chmod u+x bin/3.2_DTIConnectivity/dsi_main.py
ENV PATH=/aida/bin/3.2_DTIConnectivity:$PATH
RUN cp bin/3.2_DTIConnectivity/dsi_main.py dsi_main
COPY lib/ lib/
RUN echo "/aida/bin/dsi_studio_ubuntu_1804/dsi-studio/dsi_studio" > bin/3.2_DTIConnectivity/dsi_studioPath.txt
```

</details>

## EXAMPLE FILES

Download [**here**](https://gin.g-node.org/Aswendt_Lab/testdata_AIDA) (you probably have to clone the dataset from the gin repo. The files are annexed files, also use the raw_data folder as the test data).\
Mouse MRI data, acquired with Bruker 9.4T - cryo coil setup: adult C57BL7/6 mouse, 
T2-weighted (anatomical scan),
DTI (structural connectivity scan),
rs-fMRI (functional connectivity scan).

## ARA CREATOR
[Matlab script](https://github.com/maswendt/AIDAmri/ARA) to generate a custom version of the Allen Mouse Brain Atlas.

[<h3><b>CONTACT</h3></b>](https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/)
If you encounter problems, report directly in [![Gitter](https://badges.gitter.im/AIDA_tools/community.svg)](https://gitter.im/AIDA_tools/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

or 

join our Open Office Hour - each Thursday 3:00 pm (UTC+2) [![Zoom](https://img.shields.io/badge/Zoom-2D8CFF?style=for-the-badge&logo=zoom&logoColor)](https://uni-koeln.zoom.us/meeting/register/tJYsceyorDoqGdX4H8Z7c86_qxoaq6yOdFGM)

For all other inquiries: Markus Aswendt (markus.aswendt@uk-koeln.de)

<h3><b>LICENSE/CITATION</h3></b>
GNU General Public License v3.0
<br/>
<br/>
If you use our software or modify parts of it and use it in other ways, please cite: 
<br/>
<br/>

*Pallast N, Diedenhofen M, Blaschke S, Wieters F, Wiedermann D, Hoehn M, Fink GR, Aswendt M. Processing Pipeline for Atlas-Based Imaging Data Analysis of Structural and Functional Mouse Brain MRI (AIDAmri). Front Neuroinform. 2019 Jun 4;13:42.[doi: 10.3389/fninf.2019.00042.](https://doi.org/10.3389/fninf.2019.00042)*
___
<details>
<summary>REFERENCES</summary></b>

+ Brain Connectivity Toolbox
    + [M. Rubinov and O. Sporns (2010). Complex Network Measures of Brain Connectivity: Uses 
and Interpretations. NeuroImage 52 (3), 1059–69.](https://www.sciencedirect.com/science/article/abs/pii/S105381190901074X)
+ Allen Mouse Brain Reference Atlas
    + [Wang et al. (2020). The Allen Mouse Brain Common Coordinate Framework: A 3D Reference Atlas. Cell 181 (4), 936-953.](https://pubmed.ncbi.nlm.nih.gov/32386544/)
+ Niftyreg
    + [Ourselin, et al. (2001). Reconstructing a 3D structure from serial
histological sections. Image and Vision Computing, 19(1-2), 25–31.](https://www.sciencedirect.com/science/article/pii/S0262885600000524)
    + [Modat, et al. (2014). Global image registration using a symmetric block-
matching approach. Journal of Medical Imaging, 1(2), 024003–024003.](https://www.ncbi.nlm.nih.gov/pubmed/26158035)
    + [Rueckert, et al.. (1999). Nonrigid registration using free-form
deformations: Application to breast MR images. IEEE Transactions on Medical
Imaging, 18(8), 712–721.](https://ieeexplore.ieee.org/document/796284)
    + [Modat, et al. (2010). Fast free-form deformation using graphics processing
units. Computer Methods And Programs In Biomedicine,98(3), 278–284.](https://www.ncbi.nlm.nih.gov/pubmed/19818524)
+ FSL
    + [M.W. Woolrich, S. Jbabdi, B. Patenaude, M. Chappell, S. Makni, T. Behrens, C. Beckmann, M. Jenkinson, S.M. Smith. Bayesian analysis of neuroimaging data in FSL. NeuroImage, 45:S173-86, 2009](https://www.ncbi.nlm.nih.gov/pubmed/19059349)
    + [S.M. Smith, M. Jenkinson, M.W. Woolrich, C.F. Beckmann, T.E.J. Behrens, H. Johansen-Berg, P.R. Bannister, M. De Luca, I. Drobnjak, D.E. Flitney, R. Niazy, J. Saunders, J. Vickers, Y. Zhang, N. De Stefano, J.M. Brady, and P.M. Matthews. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 23(S1):208-19, 2004](https://www.sciencedirect.com/science/article/pii/S1053811904003933?via%3Dihub)
    + [M. Jenkinson, C.F. Beckmann, T.E. Behrens, M.W. Woolrich, S.M. Smith. FSL. NeuroImage, 62:782-90, 2012](https://www.sciencedirect.com/science/article/pii/S1053811911010603?via%3Dihub) 
+ DSIstudio
    + [Yeh, Fang-Cheng, et al. Deterministic diffusion fiber tracking improved by quantitative anisotropy. (2013): e80713. PLoS ONE 8(11)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080713)
</details>
