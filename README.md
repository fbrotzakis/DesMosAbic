# DesMosAbic

This repository provides with custom code to sample protein aminoacid sequences of a predetermined scaffold, able to exhibit slower unbinding kinetics from targets. This is achieved by a reweighting procedure using proteinMPNN, acting on series of input unbinding trajectories, which can be generated e.g by Enhanced Sampling Molecular Dynamics. In notebook  `bayesian_opt.ipynb` we use as example the optimization of CDR3 of a Nanobody scafold (H11) for slow unbinding rates from the main Receptor Binding Domain epitope of the SARS-CoV-2 spike antigen.

## Installation 
```
git clone https://gitlab.developers.cam.ac.uk/ch/sormanni/abnativ.git
git clone "https://github.com/dina-lab3D/NanoNet"
git clone https://github.com/tommyhuangthu/FASPR.git
conda env create --name DesMosAbic --file environment_test.yml
cd abnativ 
conda activate DesMosAbic
pip install .
cd ..
cd FASPR
g++ -O3  -o FASPR src/*.cpp
cd ../
jupyter-notebook bayesian_opt.ipynb
```
## Software requirements

The code can be used in Linux or macOS. 
