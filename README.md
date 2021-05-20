Dynammo: a suite of tools for multiple co-evolving time series forecasting and mining.
=======

Linear dynamical systems (aka Kalman filters) 

It provides matlab packages to learn parameters for linear dynamical systems (kalman filters), and variants under various contraints and missing values. 
It contains the following package:

DynaMMo: Mining, Compression, and Missing Value Imputation
---
1. dynammo
  Expectation-Maximization algorithm to learn the parameters of a linear dynamical system, with or without missing values.  
  There also also additional code for learning with missing observations, recovering the missing values, and compressing data.
  located in the dynammo subdirectory
  check out make demo
  
  ```
    DynaMMo: Mining and Summarization of Coevolving Sequences with Missing Values,
    Lei Li, James McCann, Nancy Pollard, and Christos Faloutsos.
    In Proceeding of the 15th ACM SIGKDD international conference on Knowledge discovery and data mining (KDD) , ACM , New York, NY, USA, 2009.
  ```
  
  The original paper is [here](https://dl.acm.org/doi/10.1145/1557019.1557078). 

PLiF: Feature extraction and clustering
---
2. PLiF
  A package for extracting shift-invariant and harmonics features from time series. These features are useful for clustering and classification. 
  
  ```
    Parsimonious linear fingerprinting for time series,
    Lei Li, B. Aditya Prakash, and Christos Faloutsos.
    In the Proceedings of the Very Large Data Bases Endowment (VLDB) , VLDB Endowment , 3, pages 385-396, 2010.
  ```
  The original paper is [here](https://dl.acm.org/doi/10.14778/1920841.1920893).
  
BoLeRO: Missing value prediction for motion capture data, guided by human body knowledge
---
3. Bolero
  A package for missing value imputation for human motion capture sequence, utilizing the human body skeleton. The decoding is essentially dynammo method plus the additional bone length constraints. 
  
  ```
    BoLeRO: a principled technique for including bone length constraints in motion capture occlusion filling,
    Lei Li, James McCann, Nancy Pollard, and Christos Faloutsos.
    In Proceedings of the 2010 ACM SIGGRAPH/Eurographics Symposium on Computer Animation (SCA) , Eurographics Association , Aire-la-Ville, Switzerland, Switzerland, pages 179-188, 2010.
  ```
  The original paper is [here](https://dl.acm.org/doi/10.5555/1921427.1921454).
  
Complex-valued Linear Dynamical System
--- 
4. CLDS
  A package for learning complex-valued dynamical systems. This is a generalization to traditional real-valued Kalman filter. It includes discrete Fourier transform as a special case. The original paper is [here](https://icml.cc/Conferences/2011/papers/159_icmlpaper.pdf).
  
  ```
    Time Series Clustering: Complex is Simpler!,
    Lei Li and B. Aditya Prakash.
    In Proceedings of the 28th International Conference on Machine Learning (ICML) , 2011.
  ```
