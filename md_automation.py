import subprocess
import time
import sys
from datetime import datetime

def read_lines_from_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.readlines()

def extract_itp_files_names(itp_files):
    itp_files_name = []
    output = []
    for itp_file in itp_files:
        itp_files_name.append(itp_file.split('/')[-1].strip())
    return output

def run_gmx_command(command, message):
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print(message)
        print ("====================================")
    except subprocess.CalledProcessError as e:
        print(f"오류 발생: {e}")
        print ("====================================")

def sleep(seconds):
    time.sleep(seconds)

def cleanup_files():
    subprocess.run('rm *.top', shell=True, check=False)
    subprocess.run('rm mdout.mdp', shell=True, check=False)
    subprocess.run('rm posre.itp', shell=True, check=False)
    #subprocess.run('rm *.gro', shell=True, check=False)
    subprocess.run('rm *.top', shell=True, check=False)
    #subprocess.run('rm *.trr', shell=True, check=False)
    #   subprocess.run('rm *.tpr', shell=True, check=False)
    subprocess.run('rm *.log', shell=True, check=False)

def print_rotating_bar(duration_seconds=3):
    rotation_symbols = ['\\', '|', '/', '-']
    start_time = time.time()

    while time.time() - start_time < duration_seconds:
        for symbol in rotation_symbols:
            sys.stdout.write('\r' + symbol)
            sys.stdout.flush()
            time.sleep(0.25)

    sys.stdout.write('\r')  # Clear the line
    sys.stdout.flush()

def update_top_file(top_file_path, itp_files_name):
    try:
        # top 파일 읽기
        with open(top_file_path, 'r') as f:
            lines = f.readlines()

        # [ molecules ] 섹션 찾기
        molecules_start = lines.index('[ moleculetype ]\n')
        txt_end = len(lines)

        print(f"txt_end :{txt_end}")

        for j in range(len(itp_files_name)):
            lines.append('#')

        print(f"len(itp_files_name) :{len(itp_files_name)}")
        print(len(lines))

        for k in range( txt_end - molecules_start  ):
            lines[ txt_end + len(itp_files_name) - k -1 ] = lines[ txt_end - k - 1 ]

        # [ molecules ] 섹션 업데이트
        for i in range(len(itp_files_name)):
            lines[molecules_start - i ] ="#include " +  '\"' + itp_files_name[i] + '\"\n'
            print("#include " +  '\"' + itp_files_name[i] + '\"\n')
        # top 파일 쓰기
        with open(top_file_path, 'w') as f:
            f.writelines(lines)

        print("topol.top 파일 업데이트 완료")
    except Exception as e:
        print(f"오류 발생: {e}")

def account()

def command_runner(command, message):
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print(message)
        print ("====================================")
    except subprocess.CalledProcessError as e:
        print(f"오류 발생: {e}")
        print ("====================================")


def additional_job():
        # 사용자 입력 및 추가 계산 실행
    user_input = input("계산할 것을 선택하세요. (1: Energy, 2: RMSD, 3: RMSF, 4: Gyrate) : ")
    input_list = user_input.split()

    for option in input_list:
        if option == '1':
            command = f'gmx energy -f md_0_1.edr -o potential.xvg'
            run_gmx_command(command, "에너지 계산 완료")
            
        elif option == '2':
            command = f'gmx rms -s md_0_1.tpr -f md_0_1.xtc -o rmsd.xvg'
            run_gmx_command(command, "RMSD 계산")
            
        elif option == '3':
            command = f'gmx rmsf -s md_0_1.tpr -f md_0_1.xtc -o rmsf.xvg'
            run_gmx_command(command, "RMSF 계산 완료")
            
        elif option == '4':
            command = f'gmx gyrate -s md_0_1.tpr -f md_0_1.xtc -o gyrate.xvg'
            run_gmx_command(command, "Gyrate 계산 완료")
            

def make_ndx(input_file_name = 'index'):
    try:
        command = f'echo "q" | gmx make_ndx -f {input_file_name}.gro -o {input_file_name}.ndx '
        run_gmx_command(command, "ndx 파일 생성 완료")
    except:
        print("ndx 파일 생성 오류 발생"

def grompp(mdp_file_name = 'minim', gro_file_name = 'md_0_1', top_file_name = 'topol', output_file_name = 'md_0_1.tpr'):
    try:
        make_ndx(input_file_name = gro_file_name + ".gro")
        command = f'gmx grompp -f {mdp_file_name}.mdp -c {gro_file_name}.gro -p {top_file_name}.top -o {output_file_name}.tpr -n {gro_file_name}.ndx -maxwarn 3'
        run_gmx_command(command, "grompp 완료")
    except:
        print("grompp 오류 발생")

def mdrun(input_file_name = 'md_0_1'):
    try:
        command = f'gmx mdrun -v -deffnm {input_file_name} -nb gpu -ntmpi 1'
        run_gmx_command(command, "mdrun 완료")
    except:
        print("mdrun 오류 발생")


def generate_mdp_file(job_name = 'Molecular_dynamics', minim_step = 50000, nvt_step = 50000, npt_step = 50000, md_step = 500000, calculation_molecules = ['Protein', 'Non-Protein']):
    try : 
        with open('minim.mdp', 'w') as f:
            f.write("; minim.mdp - used as input into grompp to generate em.tpr\n")
            f.write("; Made by Youngwoo Jung\n")
            f.write("; Parameters describing what to do, when to stop and what to save\n")
            f.write(f"integrator  = steep         ; Algorithm (steep = steepest descent minimization)\n")
            f.write("emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n")
            f.write("emstep      = 0.01          ; Energy step size\n")
            f.write("nsteps      = " + str(minim_step) + "         ; Maximum number of (minimization) steps to perform\n")
            f.write("\n")
            f.write("; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n")
            f.write("nstlist         = 1         ; Frequency to update the neighbor list and long range forces\n")
            f.write("cutoff-scheme   = Verlet    ; Buffered neighbor searching\n")
            f.write("ns_type         = grid      ; Method to determine neighbor list (simple, grid)\n")
            f.write("coulombtype     = PME       ; Treatment of long range electrostatic interactions\n")
            f.write("rcoulomb        = 1.0       ; Short-range electrostatic cut-off\n")
            f.write("rvdw            = 1.0       ; Short-range Van der Waals cut-off\n")
            f.write("pbc             = xyz       ; Periodic Boundary Conditions (yes/no)\n")

        with open('nvt.mdp', 'w') as f:
            f.write(";Made by Jung\n")
            f.write(f"title = {job_name} NVT ensemble \n")
            f.write("define = -DPOSRES ; position restrain the protein \n")
            f.write("; Run parameters \n")
            f.write("integrator = md ; leap-frog integrator \n")
            f.write(f"nsteps = {str(nvt_step)} ; 2 * {str(nvt_step)} = {str(nvt_step/1000)} ps \n")
            f.write("dt = 0.002 ; 2 fs \n")
            f.write("; Output control \n")
            f.write("nstxout = 500 ; save coordinates every 1.0 ps \n")
            f.write("nstvout = 500 ; save velocities every 1.0 ps \n")
            f.write("nstenergy = 500 ; save energies every 1.0 ps \n")
            f.write("nstlog = 500 ; update log file every 1.0 ps \n")
            f.write("; Bond parameters \n")
            f.write("continuation = no ; first dynamics run \n")
            f.write("constraint_algorithm = lincs ; holonomic constraints \n")
            f.write("constraints = h-bonds ; bonds involving H are constrained \n")
            f.write("lincs_iter = 1 ; accuracy of LINCS \n")
            f.write("lincs_order = 4 ; also related to accuracy \n")
            f.write("; Nonbonded setting \n")
            f.write("cutoff-scheme = Verlet ; Buffered neighbor searching \n")
            f.write("ns_type = grid ; search neighboring grid cells \n")
            f.write("nstlist = 10 ; 20 fs, largely irrelevant with Verlet \n")
            f.write("rcoulomb = 1.0 ; short-range electrostatic cutoff (in nm) \n")
            f.write("rvdw = 1.0 ; short-range van der Waals cutoff (in nm) \n")
            f.write("DispCorr = EnerPres ; account for cut-off vdW scheme \n")
            f.write("; Electrostatics \n")
            f.write("coulombtype = PME ; Particle Mesh Ewald for long-range electrostatics \n")
            f.write("pme_order = 4 ; cubic interpolation \n")
            f.write("fourierspacing = 0.16 ; grid spacing for FFT \n")
            f.write("; Temperature coupling is on \n")
            f.write("tcoupl = V-rescale ; modified Berendsen thermostat \n")
            f.write(f"tc-grps = {str(calculation_molecules[0])} {str(calculation_molecules[1])} ; two coupling groups - more accurate \n")
            f.write("tau_t = 0.1 0.1 ; time constant, in ps \n")
            f.write("ref_t = 300 300 ; reference temperature, one for each group, in K \n")
            f.write("; Pressure coupling is off \n")
            f.write("pcoupl = no ; no pressure coupling in NVT \n")
            f.write("; Periodic boundary conditions \n")
            f.write("pbc = xyz ; 3-D PBC \n")
            f.write("; Velocity generation \n")
            f.write("gen_vel = yes ; assign velocities from Maxwell distribution \n")
            f.write("gen_temp = 300 ; generate initial velocities for 300 K \n")
            f.write("gen_seed = -1 ; generate a random seed \n")

        with open('npt.mdp', 'w') as f:
            f.write("; Made by Youngwoo Jung\n")
            f.write(f'title                   = {job_name} NPT ensemble \n')
            f.write('define                  = -DPOSRES  ; position restrain the protein\n')
            f.write('; Run parameters\n')
            f.wrtie('integrator              = md        ; leap-frog integrator\n')
            f.write('nsteps                  = ' + str(npt_step) + ' ; 2 * ' + str(npt_step) + ' = ' + int(npt_step)/1000 + ' ps\n')
            f.write('dt                      = 0.002     ; 2 fs\n')
            f.write('; Output control\n')
            f.write('nstxout                 = 500      ; save coordinates every 1.0 ps\n')
            f.write('nstvout                 = 500      ; save velocities every 1.0 ps\n')
            f.write('nstenergy               = 500      ; save energies every 1.0 ps\n')
            f.write('nstlog                  = 500      ; update log file every 1.0 ps\n')
            f.write('; Bond parameters\n')
            f.write('continuation            = yes       ; Restarting after NVT \n')
            f.write('constraint_algorithm    = lincs     ; holonomic constraints \n')
            f.write('constraints             = h-bonds   ; bonds involving H are constrained \n')
            f.write('lincs_iter              = 1         ; accuracy of LINCS \n')
            f.write('lincs_order             = 4         ; also related to accuracy \n')
            f.write('; Nonbonded settings \n')
            f.write('cutoff-scheme           = Verlet    ; Buffered neighbor searching \n')
            f.write('ns_type                 = grid      ; search neighboring grid cells \n')
            f.write('nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet \n')
            f.write('rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm) \n')
            f.write('rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm) \n')
            f.write('DispCorr                = EnerPres  ; account for cut-off vdW scheme \n')
            f.write('; Electrostatics \n')
            f.write('coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics \n')
            f.write('pme_order               = 4         ; cubic interpolation \n')
            f.write('fourierspacing          = 0.16      ; grid spacing for FFT \n')
            f.write('; Temperature coupling is on \n')
            f.write('tcoupl                  = V-rescale             ; modified Berendsen thermostat \n')
            f.write(f'tc-grps                 = {str(calculation_molecules[0])}     {str(calculation_molecules[1])}   ; two coupling groups - more accurate \n')
            f.write('tau_t                   = 0.1         0.1       ; time constant, in ps \n')
            f.write('ref_t                   = 300         300       ; reference temperature, one for each group, in K \n')
            f.write('; Pressure coupling is on \n')
            f.write('pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT \n')
            f.write('pcoupltype              = isotropic             ; uniform scaling of box vectors \n')
            f.write('tau_p                   = 2.0                   ; time constant, in ps \n')
            f.write('ref_p                   = 1.0                   ; reference pressure, in bar \n')
            f.write('compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1 \n')
            f.write('refcoord_scaling        = com \n')
            f.write('; Periodic boundary conditions \n')
            f.write('pbc                     = xyz       ; 3-D PBC \n')
            f.write('; Velocity generation \n')
            f.write('gen_vel                 = no        ; Velocity generation is off \n')

        with open('md.mdp', 'w') as f:
            f.write("; Made by Youngwoo Jung\n")
            f.write(f'title                   = {job_name} MD run \n')
            f.write('; Run parameters\n')
            f.write('integrator              = md        ; leap-frog integrator\n')
            f.write('nsteps                  = ' + str(md_step) + ' ; 2 * ' + str(md_step) + ' = ' + int(md_step)/1000 + ' ps\n')
            f.write('dt                      = 0.002     ; 2 fs\n')
            f.write('; Output control\n')
            f.write('nstxout                 = 0      ; suppress bulky .trr file by specifying\n')
            f.write('nstvout                 = 0      ; 0 for output frequency of nstxout\n')
            f.write('nstfout                 = 0      ; nstvout, and nstfout\n')
            f.write('nstenergy               = 5000    ; save energies every 10.0 ps\n')
            f.write('nstlog                  = 5000    ; update log file every 10.0 ps\n')
            f.write('nstxout-compressed      = 5000     ; save compressed coordinates every 10.0 ps\n')
            f.write('compressed-x-grps       = System  ; save the whole system\n')
            f.write('; Bond parameters\n')
            f.write('continuation            = yes       ; Restarting after NPT \n')
            f.write('constraint_algorithm    = lincs     ; holonomic constraints \n')
            f.write('constraints             = h-bonds   ; bonds involving H are constrained \n')
            f.write('lincs_iter              = 1         ; accuracy of LINCS \n')
            f.write('lincs_order             = 4         ; also related to accuracy \n')
            f.write('; Neighborsearching\n')
            f.write('cutoff-scheme           = Verlet    ; Buffered neighbor searching \n')
            f.write('ns_type                 = grid      ; search neighboring grid cells \n')
            f.write('nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet \n')
            f.write('rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm) \n')
            f.write('rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm) \n')
            f.write('; Electrostatics\n')
            f.write('coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics \n')
            f.write('pme_order               = 4         ; cubic interpolation \n')
            f.write('fourierspacing          = 0.16      ; grid spacing for FFT \n')
            f.write('; Temperature coupling is on\n')
            f.write('tcoupl                  = V-rescale             ; modified Berendsen thermostat \n')
            f.write(f'tc-grps                 = {str(calculation_molecules[0])}     {str(calculation_molecules[1])}   ; two coupling groups - more accurate \n')
            f.write('tau_t                   = 0.1         0.1       ; time constant, in ps \n')
            f.write('ref_t                   = 300         300       ; reference temperature, one for each group, in K \n')
            f.write('; Pressure coupling is on\n')
            f.write('pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT \n')
            f.write('pcoupltype              = isotropic             ; uniform scaling of box vectors \n')
            f.write('tau_p                   = 2.0                   ; time constant, in ps \n')
            f.write('ref_p                   = 1.0                   ; reference pressure, in bar \n')
            f.write('compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1 \n')
            f.write('; Periodic boundary conditions\n')
            f.write('pbc                     = xyz       ; 3-D PBC \n')
            f.write('; Dispersion correction\n')
            f.write('DispCorr                = EnerPres  ; account for cut-off vdW scheme \n')
            f.write('; Velocity generation\n')
            f.write('gen_vel                 = no        ; Velocity generation is off \n')
        print("모든 mdp 파일 생성 완료")

    except : 
        print("mdp generation 오류 발생")