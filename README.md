[1.2]: http://i.imgur.com/wWzX9uB.png
[1]: http://www.twitter.com/AswendtMarkus
<!--social icon from https://github.com/carlsednaoui/gitsocial -->

<img align="left" src="https://github.com/maswendt/AIDAmri/blob/master/AIDA_Logo.png" width="120">
<h1>AIDA<i>mri</i></h1>

Atlas-based Imaging Data Analysis Pipeline (AIDA) for structural and functional MRI of the mouse brain
<br/>
<br/>

[Information latest Version 1.3](https://github.com/maswendt/AIDAmri/releases/tag/v1.3)

<h3><b>Manual v1.3</h3></b>

[**Link**](https://github.com/maswendt/AIDA/blob/master/manual.pdf)

<h3><b>Important note: read this before you install AIDAmri for the first time</h3></b>

We fully moved to the containerized version of AIDAmri via [Docker](https://docs.docker.com/get-docker/). All information can be found in the manual above. Please report issues and bugs directly in the issue section of this repository or at gitter (Link below in the contact section).

<h3><b>EXAMPLE FILES</h3></b>

Download [**here**](https://doid.gin.g-node.org/70e11fe472242e2d4f96c53ac9b0a556/) (you may need to browse the archive or the repository and directly download the test-data file from there. Do not download the zipped archive).\
Mouse MRI data, acquired with Bruker 9.4T - cryo coil setup: adult C57BL7/6 mouse, 
T2-weighted (anatomical scan),
DTI (structural connectivity scan),
rs-fMRI (functional connectivity scan).

<h3><b>VERSION HISTORY</h3></b>

[Information about Version 1.2 (Docker stable release)](ttps://github.com/maswendt/AIDAmri/releases/tag/v1.2)
[Information about Version 1.1.1 (Docker pre-release)](https://github.com/maswendt/AIDAmri/releases/tag/1.1.1)
<br/>
[Information about Version 1.1 (Stable)](https://github.com/maswendt/AIDAmri/releases/tag/v1.1)
<br/>
[Information about Version 1.0](https://github.com/maswendt/AIDAmri/releases/tag/v1.0)

[<h3><b>CONTACT</h3></b>](https://neurologie.uk-koeln.de/forschung/ag-neuroimaging-neuroengineering/)
If you encounter problems, report directly in [![Gitter](https://badges.gitter.im/AIDA_tools/community.svg)](https://gitter.im/AIDA_tools/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

or 

join our Open Office Hour - each Thursday 3:00 pm (UTC+2) [![Zoom](https://img.shields.io/badge/Zoom-2D8CFF?style=for-the-badge&logo=zoom&logoColor)](https://uni-koeln.zoom.us/meeting/register/tJYsceyorDoqGdX4H8Z7c86_qxoaq6yOdFGM)\

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
