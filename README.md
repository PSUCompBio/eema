# PSU Compbio FEM
Finite element code with support to embedded elements.

This is the property of The Penn State Computational Biomechanics Group.
This code is being developed by Harsha Teja Garimella and Jesse Gerber under the supervision of Dr. Reuben H Kraft.

## Motivation:
Computational Brain Biomechanics

## Acknowledgements:
Funding from CFD Research Corporation and ARL

## Contact Details:

Harsha T Garimella, <br />
Ph.D. Candidate, Mechanical Engineering, <br /> 
The Pennsylvania State University, <br />
University Park, Pennsylvania, USA. <br />
Email: harshatejagarimella@gmail.com <br />

Jesse Gerber, <br />
M.S. Student, Mechanical Engineering, <br />
The Pennsylvania State University, <br />
University Park, Pennsylvania, USA. <br />
Email: jig6@psu.edu <br />

Reuben H. Kraft, Ph.D. <br />
Shuman Asst. Professor, <br />
Department of Mechanical Engineering, <br />
Department of Biomedical Engineering, <br />
The Pennsylvania State University, <br />
University Park, Pennsylvania, USA. <br />
Email: reuben.kraft@psu.edu <br />


## TO INSTALL USING CMAKE
1. cd eema
2. mkdir build
3. cd build
4. ccmake ..
5. make 

## How to run the examples?
1. you need to compile with examples ON in the 
cmake configuration (step 4 above)
2. navigate to ~/eema/build/examples/example-1
3. ./eema_example1 . input.inp

## How to use doxygen?
1. Please edit the Doxyfile in the 'Documentation' folder to enter correct location of the project.
2. Enter the 'Documentation' folder.
2. Enter the following command: "./../third-party-libs/doxygen/build/bin/doxygen \<DOXYGEN-CONFIGURATION-FILE>".

