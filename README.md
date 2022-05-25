# gpr_demos
Easy-to-run versions of Gaussian Process Regression.

## Purpose of this Repository
Gaussian Process Regression (GPR) often feels esoteric due to the nature of high-level mathematics involved in derivations. This Repository is dedicated to the visualization of GPR. It is meant to aid my personal journey of working with GPR as well as serve as a tool for other to easily download and run basic versions of GPR. Hence, barriers of entry are minimzed (system requirements are minimal, datasets are small and generated natively).

## Purpose of the Project 
Traditional regression methodologies only provide depths at prediction points not directly measured by multi-beam sonar; however, GPR provides both depths and uncertainty at a given prediction point. Multi-beam sonar surveys generate datasets with millions of points and it is well known that GPR is intractable with massive datasets (the algorithm scales with O(N^3) due to matrix inversion). Hence, previous work (Krasnosky, 2021) demonstrates computation speed-up from the approach of memory storage, cholesky factor updates, and processing with GPUs and CUDA. What has not be implemented, is the use to informative downsampling of massive dataset of N terms to generate a sparse representative subset containing M datapoints that is able to be processed efficiently.

The use of a representative subset (M datapoints) with Kranosky's computational resource management should yield drastically accelerated calculation time while maintaining accuracy comparable to full GPR (all N datapoints).

### About the Author
My name is Phil Parisi and I am working toward my MS/PhD in Ocean Engineering. Developing a skillset to work on autonomous marine robotics.
LinkedIn: https://www.linkedin.com/in/phil-parisi/
Portfolio & Travel: https://philparisi.weebly.com/
GitHub: https://github.com/PhilParisi
YouTube Channel: https://www.youtube.com/PhilsBeginnerCode

### Acknowledgements
This work would not be possible without my PhD Advisor, Dr. Chris Roman (URI), and preceeding work done by Dr. Kristopher Krasnosky (URI). 

