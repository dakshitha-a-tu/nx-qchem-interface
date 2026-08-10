# **NewtonX-QChem Interface for Non-Adiabatic Dynamics**

This repository contains a robust Perl script (`run-qchem-eom.pl`) that serves as an automated communication layer between the [NewtonX](https://newtonx.org/) dynamics software and the [Q-Chem](https://www.q-chem.com/) quantum chemistry package.

Specifically designed to handle non-adiabatic molecular dynamics simulations (trajectory surface hopping), this interface orchestrates on-the-fly execution of Q-Chem to extract potential energies, analytical gradients, and non-adiabatic coupling matrix elements (NACMEs) using Equation-of-Motion Coupled-Cluster (EOM-CC) methods.

## **Compatibility**

This interface has been tested and is confirmed to work with the following software versions:

* **NewtonX:** 2.4  
* **Q-Chem:** 5.4 (Should be compatible with later versions)

**Supported Methods:**
* **EOM-IP-CCSD** (Equation-of-Motion Ionization Potential Coupled-Cluster Singles and Doubles)

---

## **Key Features**

* **EOM-CC Support:** Capable of running dynamics with EOM-IP-CCSD to calculate energies, gradients, and NACs between target states.  
* **Automatic Archiving:** At each time step, the script automatically archives the full Q-Chem output files to `INFO_RESTART/qchem_outputs/` for detailed post-analysis and easy restarts.  
* **Robust Parsing:** Specifically designed to read the "State A / State B" gradient and NAC blocks produced by Q-Chem's EOM-CC module, ensuring correct phase handling and state assignment.  
* **Detailed Output:** Prints formatted energies, gradients, and NACs to standard output at each time step for easy monitoring of the trajectory progress.

---

## **How It Works**

During a NewtonX dynamics run, this script is called at every time step to perform the necessary quantum chemical calculations. Its execution follows a strict multi-step workflow:

1. **Initialization and Status Tracking:** The script reads the current dynamics status from NewtonX (`control.d`), identifying the active state, the current time step, and the total number of states. It verifies that the calculation is a non-adiabatic dynamics job (adiabatic runs are deliberately aborted).
2. **Input Preparation:** It reads the current molecular geometry from the NewtonX `geom` file (provided in Bohr) and converts the coordinates to Angstroms. The script then safely injects this updated geometry into your template Q-Chem input file (`JOB_NAD/qchem.inp`) while preserving the `$molecule` block structure. It also dynamically updates the `CC_STATE_TO_OPT` keyword so Q-Chem explicitly targets the active electronic state for gradient calculations.
3. **Execution and Archiving:** Q-Chem is executed via the command line (defaulting to 4 threads via `-nt 4`). The output of each time step is then archived into the `../INFO_RESTART/qchem_outputs/` directory, named according to the current trajectory time.
4. **Data Extraction:** The script parses the resulting `qchem.out` file to extract three critical sets of data:
   * **Energies:** Scans for `Total energy` entries associated with EOM transitions, validates the state count, and writes them to `epot` (in Hartrees).
   * **Gradients:** Locates `G_I` and `G_J` gradient blocks. It maps the forces to the correct state index, outputting the active state gradient to `grad` and all gradients to `grad.all`. (For single-state calculations, it falls back to parsing the `Final gradient`).
   * **Couplings (NACMEs):** Extracts the spatial derivative couplings from `NAC d^x_IJ (CI part)` blocks. The parser strictly manages Q-Chem's output order mapping $I$ and $J$ indices to NewtonX's expected lower-triangular matrix format, adjusting mathematical signs ($d_{JI} = -d_{IJ}$) as necessary.
5. **Phase Continuity Correction:** Because quantum mechanical wavefunctions carry an arbitrary phase sign, the script interfaces with NewtonX's internal `escalar` utility. It calculates the dot product between the coupling vectors at $t$ and $t-\Delta t$, smoothly correcting any unphysical phase inversions before writing the final `nad_vectors`.

---

## **Installation and Prerequisites**

Before using the script, you must configure NewtonX to recognize it as a valid program interface.

1. **Locate `colib_perl.pm`:** Open the `colib_perl.pm` file, located at `$NX/lib/colib_perl.pm` (where `$NX` is the environment variable pointing to your NewtonX installation's bin directory).  
2. **Add Program Definition:** Add the following Perl code block to the `colib_perl.pm` file. This defines "program 22.0" as the Q-Chem interface.
   
   ```perl
   if (($prog >= 21.95) and ($prog < 22.05)){  
          %progconf=(progname      => "qchem",  
                     methodname => "qchem",  
                     ic         => "y",  
                     dyn        => "y",  
                     hyb        => "n",  
                     ip         => "n",  
                     key        => sprintf("%4.1f",22.0),  
                     label      => "QChem",  
                     method     => "EOM-IP-CCSD",  
                     parfile    => "qchem.par",  
                     nad_exec   => "y",  
                     lvprt_d    => 1,  
                     vdoth_d    => 0,  
                     never_state_d  => 0,  
                     cio_options_d  => "NULL",  
                     cisc_options_d => "NULL",  
                     cprog_d    => 0,  
                     progic     => "",  
                     progdyn    => "run-qchem-eom.pl");  
   }

4. **Place the Script:** Copy the `run-qchem-eom.pl` script into your NewtonX bin directory (`$NX/bin/`) and ensure it is executable (`chmod +x run-qchem-eom.pl`).

## **Running a Trajectory**
To start a non-adiabatic dynamics simulation, you need to set up a trajectory directory with a specific structure and set of input files.

1. **Set `prog` in `control.dyn`:** In your main dynamics control file (`control.dyn`), set the `prog` parameter to `22.0` to tell NewtonX to use this Q-Chem interface.

```prog = 22.0```

2. **Directory Structure:** Your trajectory directory (e.g., `TRAJ1`) must contain the following files and subdirectories:

```
TRAJ1/    
├── control.dyn      # Main NewtonX dynamics parameters    
├── sh.inp           # Surface hopping parameters    
├── geom             # Initial geometry (in Bohr)    
├── veloc            # Initial velocities (in atomic units)    
└── JOB_NAD/         # Directory for non-adiabatic calculations    
    └── qchem.inp    # Q-Chem input template
```
3. **Q-Chem Template (`qchem.inp`):** The `JOB_NAD/qchem.inp` file is your Q-Chem input template. Specify the EOM-CCSD settings, basis set, and standard `$rem` variables here.

   **Important Template Rules:**
   * Include `$molecule` and `$end` tags for the geometry section. The actual atomic coordinates inside this block are optional (they will be overwritten by the script), but the **charge and multiplicity line is mandatory**. You must include the charge/multiplicity line immediately following the `$molecule` tag.
   * Ensure `CC_STATE_TO_OPT` is present in your $rem section. This will also be overwritten by the script dynamically.

Example Template Snippet:

```
$molecule 
1 2
O  0.000000  0.000000  0.000000
H  0.758602  0.000000  0.504284
H  0.758602  0.000000 -0.504284
$end

$rem   
METHOD              EOM-CCSD  
BASIS               cc-pVDZ  
IP_STATES           [3]  
CC_STATE_TO_OPT     [1,3]  
CALC_NAC            2  
...  
$end
```
---
## **Credits**

**Authors:** Dakshitha Abeygunewardane and Spiridoula Matsika

**Affiliation:** Matsika Lab, Temple University, USA

Year: 2025
