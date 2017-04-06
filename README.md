# PSU Compbio FEM
Finite element code with support to embedded elements.

This is the property of The Penn State Computational Biomechanics Group.
This code is developed by Harsha Teja Garimella under the supervision of Dr. Reuben H Kraft.

## Motivation:
Computational Brain Biomechanics

## Acknowledgements:
Funding from CFDRC and Funding from ARL

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
1. cd embeddedfem
2. mkdir build
3. cd build
4. ccmake ..
5. make 

## How to run the code:
1. you need to compile with examples ON in the 
cmake configuration (step 4 above)
2. navigate to ~/embeddedfem/build/examples/example-1
3. ./example1 . input.inp

## How to run the code?

1. Build the code using the command "make" or "make all". <br />
   If the system you are running the code has multiple processors, 
   one can use the command "make -j \<NUMBER-OF-PROCESSORS>". <br />
   The above step reduces the build time.

2. Run the code using the command "./main \<JOB-HOME-FOLDER> \<JOB-INPUT-FILE-NAME> <br />
   for example: ./main examples/example-1 input.inp

## How to use doxygen?

1. Please edit the Doxyfile in the 'Documentation' folder to enter correct location of the project.
2. Enter the 'Documentation' folder.
2. Enter the following command: "./../third-party-libs/doxygen/build/bin/doxygen \<DOXYGEN-CONFIGURATION-FILE>".

