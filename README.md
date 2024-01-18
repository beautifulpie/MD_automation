# MD_automation
Molecular dynamics automation program for GROMACS
Made by Jung Youngwoo
From Spider_core 

## Usage
1. Install the GROMACS (https://tutorials.gromacs.org/), VMD(https://www.ks.uiuc.edu/Research/vmd/) or Access the port 51000, cuop account
2. Run the code below
```shell
python GUI.py
```
3. Click the each button, and add the files. (You can check the uploaded files in text box)
4. Check the name in input_file_path.txt file (please do not write the file path, Just writing a name of file)
5. Enter the number of the operation you want to perform. Multiple inputs separated by spaces are also supported.
6. Waiting (I will change the code to use the GPU. Document say that we can use gpu using the mdrun code with '-nb gpu', but I wasn't work. So I'm going to fix it.)
