# gpr_demos
Easy-to-run versions of Gaussian Process Regression.

## Purpose of this Repository
Gaussian Process Regression (GPR) often feels esoteric due to high-level mathematics involved. This repository is dedicated to simple GPR problems and visualization thereof. It is meant to not only aid my personal journey as I learn GPR, but also be a tool for others to easily download and run basic versions of GPR. Hence, barriers of entry are minimzed (free and thoroughly commented code, minimal system requirements, compact datsets) and the MATLAB language is used.

## Purpose of the Project 
Multi-beam sonar systems use acoustics to measure distance from system to points along the seafloor. These measurements are converted into depths and subsequent bathymetric (seafloor) maps are generated. Traditional regression methodologies only provide depths at prediction points not directly measured by multi-beam sonar; however, GPR provides both depths and uncertainty at a given prediction point. Multi-beam sonar surveys generate datasets with millions of points and it is well known that the GPR algorithm is intractable with massive datasets (the algorithm scales with O(N^3) due to matrix inversion). Hence, previous work (Krasnosky, 2021) demonstrated computation speed-up from the approach of memory storage, cholesky factor updates, and processing with GPUs/CUDA. What has not be implemented by our team (but has been in literature), is the use of informative downsampling of a massive dataset with N datapoints to generate a sparse representative subset containing M datapoints that is able to be processed efficiently.

The use of a representative subset (M datapoints) with Kranosky's computational resource management should yield drastically accelerated calculation time while maintaining accuracy comparable to full GPR (all N datapoints).

### About the Author
My name is Phil Parisi and I am working toward my MS/PhD in Ocean Engineering. Developing a skillset to work on autonomous marine robotics. Feel free to connect! <br>
LinkedIn: https://www.linkedin.com/in/phil-parisi/ <br>
Portfolio & Travel: https://philparisi.weebly.com/ <br>
GitHub: https://github.com/PhilParisi <br>
YouTube Channel: https://www.youtube.com/PhilsBeginnerCode <br>
Patreon: https://patreon.com/PhilsBeginnerCode <br>

### Acknowledgements
This work would not be possible without my PhD Advisor, Dr. Chris Roman (URI), and preceeding work done by Dr. Kristopher Krasnosky (URI). Shoutout to the brilliant minds behind the development of classic and modern GPR for Machine Learning: Rasmussen, Williams, Seeger, Lawrence, Csato, Qi, and many others.

