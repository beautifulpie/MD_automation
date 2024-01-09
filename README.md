# MD_automation
Molecular dynamics automation program for GROMACS
Made by Jung Youngwoo
From Spider_core 

## Usage
1. Install the GROMACS (https://tutorials.gromacs.org/), VMD(https://www.ks.uiuc.edu/Research/vmd/) or Access the port 51000, cuop account
2. Check the name in input.txt file (please do not write the file path, Just writing a name of file)
3. Run the code below
```shell
python main.py
```
4. Waiting (Let's change the code to use the GPU. I replaced the mdrun code with '-nb gpu' in the documentation, but I get an error. So I'm going to fix it.)
5. Enter the number of the operation you want to perform. Multiple inputs separated by spaces are also supported.