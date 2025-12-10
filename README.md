# **NewtonX-QChem Interface for Non-Adiabatic Dynamics**

This repository contains a refactored and robust Perl script, run-qchem-eom.pl, that serves as a communication layer between the [NewtonX](https://newtonx.org/) dynamics software and the [Q-Chem](https://www.q-chem.com/) quantum chemistry package.

The interface is designed to handle non-adiabatic dynamics simulations using the Equation-of-Motion Coupled-Cluster (EOM-CCSD) methods available in Q-Chem.

## **Compatibility**

This interface has been tested and is confirmed to work with the following software versions:

* **NewtonX:** 2.4  
* **Q-Chem:** 5.4 (Should be compatible with later versions)

**Supported Methods:**

* **EOM-IP-CCSD** (Equation-of-Motion Ionization Potential Coupled-Cluster Singles and Doubles)

## **Key Features**

* **EOM-CC Support:** Capable of running dynamics with EOM-IP-CCSD to calculate energies, gradients, and non-adiabatic couplings (NACs) between target states.  
* **Automatic Archiving:** At each time step, the script automatically archives the full Q-Chem output files to INFO\_RESTART/qchem\_outputs/ for detailed post-analysis.  
* **Robust Parsing:** Specifically designed to read the "State A / State B" gradient and NAC blocks produced by Q-Chem's EOM-CC module, ensuring correct phase handling and state assignment.  
* **Detailed Output:** Prints formatted energies, gradients, and NACs to standard output at each time step for easy monitoring of the trajectory progress.

## **Installation and Prerequisites**

Before using the script, you must configure NewtonX to recognize it as a valid program interface.

1. **Locate colib\_perl.pm:** Open the colib\_perl.pm file, which is located at $NX/../lib/colib\_perl.pm, where $NX is the environment variable pointing to your NewtonX installation's bin directory.  
2. **Add Program Definition:** Add the following Perl code block to the colib\_perl.pm file. This defines "program 22.0" as the Q-Chem interface.
   
   '''perl
   if (($prog \>= 21.95) and ($prog \< 22.05)){  
          %progconf=(progname      \= "qchem",  
                     methodname \=\> "qchem",  
                     ic         \=\> "y",  
                     dyn        \=\> "y",  
                     hyb        \=\> "n",  
                     ip         \=\> "n",  
                     key        \=\> sprintf("%4.1f",22.0),  
                     label      \= "QChem",  
                     method     \= "EOM-IP-CCSD",  
                     parfile    \= "qchem.par",  
                     nad\_exec   \= "y",  
                     lvprt\_d    \=\> 1,  
                     vdoth\_d    \=\> 0,  
                     never\_state\_d  \=\> 0,  
                     cio\_options\_d  \=\> "NULL",  
                     cisc\_options\_d \=\> "NULL",  
                     cprog\_d    \=\> 0,  
                     progic     \=\> "",  
                     progdyn    \=\> "run-qchem-eom.pl");  
   }
   '''

4. **Place the Script:** Copy the run-qchem-eom.pl script into your NewtonX bin directory ($NX/bin/) and ensure it is executable (chmod \+x run-qchem-eom.pl).

## **Running a Trajectory**

To start a non-adiabatic dynamics simulation, you need to set up a trajectory directory with a specific structure and set of input files.

1. **Set prog in control.dyn:** In your main dynamics control file (control.dyn), set the prog parameter to 22.0 to tell NewtonX to use this Q-Chem interface.  
   prog \= 22.0  
2. **Directory Structure:** Your trajectory directory (e.g., TRAJ1) must contain the following files and subdirectories:  
   TRAJ1/    
   ├── control.dyn      \# Main NewtonX dynamics parameters    
   ├── sh.inp           \# Surface hopping parameters    
   ├── geom             \# Initial geometry (in Bohr)    
   ├── veloc            \# Initial velocities (in atomic units)    
   └── JOB\_NAD/         \# Directory for non-adiabatic calculations    
       └── qchem.inp    \# Q-Chem input template

3. **Q-Chem Template (qchem.inp):** The JOB\_NAD/qchem.inp file is your Q-Chem input template. You must specify the EOM-CCSD settings, basis set, and standard rem variables here.  
   **Important Template Rules:**  
   * Include $molecule and $end tags for the geometry section, but **do not** include the actual coordinates. The interface will automatically insert the updated geometry and the charge/multiplicity line found in the template.  
   * Ensure CC\_STATE\_TO\_OPT is present if you want the interface to verify/set the target state optimization list (though the interface may simply pass your template through depending on configuration).

**Example Template Snippet:**$molecule  
0 2  
$end

$rem  
JOBTYPE             FORCE  
METHOD              EOM-CCSD  
BASIS               cc-pVDZ  
IP\_STATES           \[3\]  
CC\_STATE\_TO\_OPT     \[1,3\]  
CALC\_NAC            2  
...  
$end

## **Credits**

**Authors:** Dakshitha Abeygunewardane and Spiridoula Matsika

**Affiliation:** Matsika Lab, Temple University, USA

**Year:** 2025
