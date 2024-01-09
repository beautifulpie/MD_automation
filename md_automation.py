import subprocess
import time

def read_lines_from_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.readlines()

def extract_itp_files_names(itp_files):
    itp_files_name = []
    for itp_file in itp_files:
        itp_files_name.append(itp_file.split('/')[-1].strip())
    return itp_files_name

def print_header():
    print("====================================")
    print('\033[1m'+"MD_automation"+'\0333')
    print("Made by Youngwoo Jung")
    print("Final update : 2024.01.09")
    print("====================================")

def run_gmx_command(command, message):
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print(message)
    except subprocess.CalledProcessError as e:
        print(f"오류 발생: {e}")

def sleep(seconds):
    time.sleep(seconds)

def cleanup_files():
    subprocess.run('rm topol.top*', shell=True, check=True)
    subprocess.run('rm mdout.mdp*', shell=True, check=True)
    subprocess.run('rm posre.itp*', shell=True, check=True)
    subprocess.run('rm conf.gro*', shell=True, check=True)
    subprocess.run('rm out.gro*', shell=True, check=True)
