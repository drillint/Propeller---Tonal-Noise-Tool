# Hanson's Model in Frequency Domain - Tonal Noise of Rotors in Uniform Inflow
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?host=gitlab.tudelft.nl&repo=jatindergoyal/hanson-model-helicoidal-theory)

This code is based on Hanson's Paper titled "Helicoidal Surface Theory for Harmonic Noise of Propellers in the Far Field" with DOI [https://doi.org/10.2514/3.50873](https://doi.org/10.2514/3.50873). This code is meant to be used for tonal noise prediction of rotors with uniform inflow with propeller axis aligned with the inflow. This code is based on the frequency domain formulation and is only valid for far-field noise predictions.


## Requirements

This code has been developed and tested in Matlab R2021a and R2023a.  


## Installation

No additional installation is required for the code except an activated installation of Matlab (preferably 2021a).  


## Structure

```
.
├── CITATION.cff
├── data
│   └── Xprop
│       ├── FA_rR.csv
│       ├── Loading
│       │    ├── BEM_Cd_rR_J1.60_idx71.csv
│       │    ├── BEM_Cl_rR_J1.60_idx71.csv
│       │    ├── BEM_dCd_dxc_J1.60_idx71.csv
│       │    ├── BEM_dCl_dxc_J1.60_idx71.csv
│       │    ├── BEM_dQ_NdrR_J1.60_idx71_pitch15B3.csv
│       │    └── BEM_dT_NdrR_J1.60_idx71_pitch15B3.csv
│       ├── MCA_rR.csv
│       ├── bd_rR.csv
│       ├── t_xc.csv
│       └── tb_rR.csv
├── LICENSES  
│   ├── Apache-License-v2.0.txt
│   ├── CC-BY-4.0.txt
│   ├── license_cosspace.txt
│   └── license_simps.txt
├── README.md
└── src
    ├── Mainfile.m   
    ├── PsiLoad.m
    ├── PsiV.m
    ├── SPLprocessing.m
    ├── TQLD.m
    ├── UserInput.m
    ├── cosspace.m
    ├── hanson.m
    ├── inputHanson.m
    ├── interpData.m
    ├── simps.m
    └── timeSignal.m

```          

The `data/` directory contains the propeller geometry files determining the blade properties. The `Loading/` directory contains the example blade loading distributions, which are used to calculate loading noise.

The **UserInput.m** file in `src/` is the file where ambient conditions, propeller parameters, blade properties, and operating conditions are given.   

The **Mainfile.m** in `src/` is the script to run. There, you can set the number of harmonics to be calculated by adjusting the parameter `m`.  

With the default settings in the code, you should get `P`, `t`, and `SPL` as the output structure arrays (known as a struct). 

The `P` struct contains the output from Hanson's model in the frequency domain. The first dimension of the variables in struct `P` corresponds to different harmonics, and the second dimension corresponds to the different observer positions. The meaning of various variables contained in the struct `P` is explained below:

- `P.Vm` refers to the thickness noise (volume displacement).  
- `P.Lm` refers to lift noise.  
- `P.Dm` refers to drag noise.  
- `P.Tm` refers to thrust noise.  
- `P.Qm` refers to torque noise.  

Either lift and drag combined or thrust and torque combined give the total loading noise.

The `t` struct contains the pressure signal in Pascals in the variable `t.p`. The first dimension of the variable `t.p` refers to different observer locations, whereas the second dimension corresponds to the time in seconds contained in the variable `t.range`.

The `SPL` struct contains the processed outputs in decibels. The variables `Vm`, `Lm`, `Dm`, `Tm`, and `Qm` are the same for struct `P`. The remaining variables are explained below:

- `SPL.loadm` refers to **loading noise** in dB calculated for different harmonics, expected size = (no. of harmonics X no. of observer locations).  
- `SPL.totalm` refers to **total tonal noise** (sum of thickness and loading noise) in dB calculated for different harmonics, expected size = (no. of harmonics X no. of observer locations).  
- `SPL.V` refers to OverAll Sound Pressure Level (OASPL) in dB for **thickness noise** calculated by summing the contribution of all the harmonics, expected size = (1 X no. of observer locations).  
- `SPL.L` refers to OASPL in dB for **lift noise** calculated by summing the contribution of all the harmonics, expected size = (1 X no. of observer locations).  
- `SPL.D` refers to OASPL in dB for **drag noise** calculated by summing the contribution of all the harmonics, expected size = (1 X no. of observer locations).  
- `SPL.T` refers to OASPL in dB for **thrust noise** calculated by summing the contribution of all the harmonics, expected size = (1 X no. of observer locations).  
- `SPL.Q` refers to OASPL in dB for **torque noise** calculated by summing the contribution of all the harmonics, expected size = (1 X no. of observer locations).  
- `SPL.loadingLD` refers to OASPL in dB for **loading noise** calculated by summing the contribution of all the harmonics of lift and drag noise, expected size = (1 X no. of observer locations).  
- `SPL.total` refers to OASPL in dB for **total tonal noise** calculated by summing the contribution of all the harmonics of loading and thickness noise, expected size = (1 X no. of observer locations).  
- `SPL.loadingTQ` refers to OASPL in dB for **loading noise** calculated by summing the contribution of all the harmonics of thrust and torque noise, expected size = (1 X no. of observer locations).  

If the code ran OK, then `SPL.loadingLD` and `SPL.loadingTQ` should have identical values.  


## Usage 

To run this code:

- Setup the case in **UserInput.m** file.  
- Set the number of harmonics to be calculated in **Mainfile.m** file by adjusting the parameter `m`.  
- Evaluate your case by running the **Mainfile.m** file in Matlab.  


## Documentation

The theory and mathematical expressions behind the code can be found in the paper "Helicoidal Surface Theory for Harmonic Noise of Propellers in the Far Field", [https://doi.org/10.2514/3.50873](https://doi.org/10.2514/3.50873). The paper presents a mathematical model to predict and analyze the noise generated by spinning propellers when observed from a distance, or the "far field." The model takes into account the propeller's geometry, aerodynamic forces, and the concept of a helicoidal surface to understand how the changes in air pressure and velocity around the propeller generate noise.
 
## Author(s)

This software has been developed by **Jatinder Goyal** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0001-7448-8617](https://orcid.org/0000-0001-7448-8617), Technische Universiteit Delft


## License

The function `simps.m` is licensed under **license_simps** (see [license_simps](LICENSES/license_simps.txt) file) with copyright &copy; 2009, Damien Garcia.   

The function `cosspace.m` is licensed under **license_cosspace** (see [license_cosspace](LICENSES/license_cosspace.txt) file) with copyright &copy; 2016, The MathWorks, Inc. *All rights reserved*.  

The contents in the `InputData/` directory are licensed under a **CC-BY 4.0** license (see [CC-BY-4.0](LICENSES/CC-BY-4.0.txt) file). The source code and any other file in this repository are licensed under an **Apache License v2.0** (see [Apache-License-v2.0](LICENSES/Apache-License-v2.0.txt) file).

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program “Hanson's Model in Frequency Domain - Tonal Noise of Rotors in Uniform Inflow” (meaning the source code files licensed under **Apache License v2.0** as explained above). It is a Matlab code to calculate the tonal noise of rotors operating in a uniform inflow condition written by the Author(s).  
Henri Werij, Dean of Faculty of Aerospace Engineering, Technische Universiteit Delft.

&copy; 2023, J. Goyal


## References

- Hanson, D. B. (1980). Helicoidal surface theory for harmonic noise of propellers in the far field. AIAA journal, 18(10), 1213-1220.  
- Magliozzi, B., Hanson, D. B., & Amiet, R. K. (1991). Propeller and propfan noise. Aeroacoustics of flight vehicles: theory and practice, 1, 1-64.  
- Cole Stephens (2023). cosspace [https://www.mathworks.com/matlabcentral/fileexchange/5491-cosspace](https://www.mathworks.com/matlabcentral/fileexchange/5491-cosspace), MATLAB Central File Exchange. Retrieved November 6, 2023.  
- Damien Garcia (2023). Simpson's rule for numerical integration [https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration](https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration), MATLAB Central File Exchange. Retrieved November 6, 2023.  

## Cite this repository

If you use this software, please cite it as below or check out the **CITATION.cff** file.

**How to cite this repository**: Goyal, J., 2024. Hanson's Model in Frequency Domain - Tonal Noise of Rotors in Uniform Inflow. 4TU.ResearchData. Software. [https://doi.org/10.4121/7da5aa45-e44b-4fa3-9407-8bf61e835d99](https://doi.org/10.4121/7da5aa45-e44b-4fa3-9407-8bf61e835d99) 


## Would you like to contribute?

If you have any comments, feedback, or recommendations, feel free to **reach out** by sending an email to jatindergoel777@gmail.com

If you want to contribute directly, you are welcome to **fork** this repository.

Thank you, and enjoy!
