import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QPushButton, QFileDialog, QPlainTextEdit, QWidget
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtCore import QThread, pyqtSignal
import time

import md_automation as md
import subprocess
import os
from pyfiglet import Figlet
import platform

current_os = platform.system()

f = Figlet(font='slant')

def get_file_names_in_current_directory():
    current_directory = os.getcwd()  # 현재 디렉토리 경로
    files = [f for f in os.listdir(current_directory) if os.path.isfile(os.path.join(current_directory, f))]
    return files

def remove_path_from_file_names(file_names):
    # 파일 이름에서 경로 제거
    file_names_only = [os.path.splitext(os.path.basename(file))[0] for file in file_names]
    return file_names_only

def generate_input_path(pdb_file, minim_mdp, npt_mdp, nvt_mdp, md_mdp, itp_files):
    # pdb 파일 이름
    pdb_file_name = os.path.splitext(os.path.basename(pdb_file))[0] + '.pdb'
    # mdp 파일 이름
    minim_mdp_name = os.path.splitext(os.path.basename(minim_mdp))[0] + '.mdp'
    npt_mdp_name = os.path.splitext(os.path.basename(npt_mdp))[0] + '.mdp'
    nvt_mdp_name = os.path.splitext(os.path.basename(nvt_mdp))[0] + '.mdp'
    md_mdp_name = os.path.splitext(os.path.basename(md_mdp))[0] + '.mdp'
    # itp 파일 이름
    itp_files_name = remove_path_from_file_names(itp_files) 
    # input_path.txt 파일 생성
    with open('임시.txt', 'w') as f:
        f.write("====================================\n")
        f.write("MD_automation\n")
        f.write("\n")
        f.write("Made by Youngwoo Jung\n")
        f.write("====================================\n")
        f.write(f"structure file name : \n {pdb_file_name}\n")
        f.write("\n")
        f.write(f"Simulator parameter file name : \n {md_mdp_name}\n {minim_mdp_name}\n {nvt_mdp_name}\n {npt_mdp_name}\n")
        f.write("\n")
        f.write("Itp file names : \n")
        if itp_files_name == []:
            f.write("posre.itp")
        for name in itp_files_name:
            f.write(f"{name}\n")

    print(f"pdb_file_path : {pdb_file}")
    print(f'minim_mdp_path : {minim_mdp}')
    print(f'npt_mdp_path : {npt_mdp}')
    print(f'nvt_mdp_path : {nvt_mdp}')
    print(f'md_mdp_path : {md_mdp}')
    print(f'pdb_file_name : {pdb_file_name}')
    print(f'minim_mdp_name : {minim_mdp_name}')
    print(f'itp_files_path : {itp_files}')

class MDGui(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle('MD Automation GUI')
        self.setGeometry(100, 100, 880, 600)

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout(self.central_widget)

        self.label = QLabel('MD Automation', self)
        self.label.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(self.label)

        # PDB File
        self.add_file_selection_button('Select PDB File', self.get_pdb_file)

        # MDP Files
        self.add_file_selection_button('Select md.mdp', self.get_md_mdp_file)
        self.add_file_selection_button('Select minim.mdp', self.get_minim_mdp_file)
        self.add_file_selection_button('Select npt.mdp', self.get_npt_mdp_file)
        self.add_file_selection_button('Select nvt.mdp', self.get_nvt_mdp_file)

        # ITP Files
        self.add_file_selection_button('Select ITP Files', self.get_itp_files)

        # Output Text
        self.output_text = QPlainTextEdit(self)
        self.layout.addWidget(self.output_text)
        
        self.output_text.appendPlainText("MD_automation")
        self.output_text.appendPlainText("Made by Jung Youngwoo")
        self.output_text.appendPlainText("Final update : 2024.01.11")
        
        print("==============================================================================================")
        print(f.renderText('Spider Core MD Automation'))
        print(f.renderText('Made by J.YW'))
        print("==============================================================================================")
        # Run Button    
        self.run_button = QPushButton('Run MD Automation', self)
        self.run_button.clicked.connect(self.run_md_automation)
        self.layout.addWidget(self.run_button)

        # Selected File Variables
        self.selected_pdb_file = ""
        self.selected_md_mdp_file = ""
        self.selected_minim_mdp_file = ""
        self.selected_npt_mdp_file = ""
        self.selected_nvt_mdp_file = ""
        self.selected_itp_files = []

    def add_file_selection_button(self, button_text, click_function):
        button = QPushButton(button_text, self)
        button.clicked.connect(click_function)
        self.layout.addWidget(button)

    def get_pdb_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("PDB Files (*.pdb);;All Files (*)")
        file_dialog.setDefaultSuffix("pdb")
        file_path, _ = file_dialog.getOpenFileName(self, "Select PDB File", "", "PDB Files (*.pdb);;All Files (*)")

        if file_path:
            self.selected_pdb_file = file_path
            self.output_text.appendPlainText(f"Selected PDB file: {self.selected_pdb_file}")
        
        print(f"Get pdb file : {self.selected_pdb_file}")

    def get_itp_files(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("ITP Files (*.itp);;All Files (*)")
        file_dialog.setDefaultSuffix("itp")
        file_paths, _ = file_dialog.getOpenFileNames(self, "Select ITP Files", "", "ITP Files (*.itp);;All Files (*)")

        if file_paths:
            self.selected_itp_files = file_paths
            self.output_text.appendPlainText(f"Selected ITP files: {', '.join(self.selected_itp_files)}")
        
        print(f"Get itp files: {self.selected_itp_files}")

    def get_md_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (*.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select md.mdp File", "", "MDP Files (*.mdp);;All Files (*)")

        if file_path:
            self.selected_md_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected md.mdp file: {self.selected_md_mdp_file}")
        
        print(f"Get md.mdp file :{self.selected_md_mdp_file}")

    def get_minim_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (minim.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select minim.mdp File", "", "MDP Files (minim.mdp);;All Files (*)")

        if file_path:
            self.selected_minim_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected minim.mdp file: {self.selected_minim_mdp_file}")
        print(f"Get minim.mdp file {self.selected_minim_mdp_file}")

    def get_npt_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (npt.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select npt.mdp File", "", "MDP Files (npt.mdp);;All Files (*)")

        if file_path:
            self.selected_npt_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected npt.mdp file: {self.selected_npt_mdp_file}")
        print(f"Get npt.mdp file : {self.selected_npt_mdp_file}")

    def get_nvt_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (nvt.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select nvt.mdp File", "", "MDP Files (nvt.mdp);;All Files (*)")

        if file_path:
            self.selected_nvt_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected nvt.mdp file: {self.selected_nvt_mdp_file}")
        print(f"Get nvt.mdp file : {self.selected_nvt_mdp_file}")

    def run_md_automation(self):
        inst = 0
        if not self.selected_pdb_file :
            self.output_text.appendPlainText("Please add pdb file.")
            inst = 1
        if not self.selected_md_mdp_file :
            self.output_text.appendPlainText("Please add md.mdp file.")
            inst = 1
        if not self.selected_minim_mdp_file :
            self.output_text.appendPlainText("Please add minim.mdp file.")
            inst = 1
        if not self.selected_npt_mdp_file :
            self.output_text.appendPlainText("Please add npt.mdp file.")
            inst = 1
        if not self.selected_nvt_mdp_file :
            self.output_text.appendPlainText("Please add nvt.mdp file.")
            inst = 1
        if inst == 1 :
            return

        # Run MD automation with the selected files
        self.label.setText("MD Automation is running... Please wait.")
        self.output_text.appendPlainText("====================================")
        self.output_text.appendPlainText("Generate the input_path_file")
        generate_input_path(self.selected_pdb_file, self.selected_minim_mdp_file, self.selected_npt_mdp_file, self.selected_nvt_mdp_file, self.selected_md_mdp_file, self.selected_itp_files)

        self.output_text.appendPlainText("====================================")
        self.output_text.appendPlainText("Run the main.py")

        if current_os == 'Linux':
            subprocess.run('python main.py', shell=True, check=True)
            time.sleep(0.5)
            print ("Run the main.py")

        # Display success message
        self.label.setText("MD Automation completed successfully!")
        self.output_text.appendPlainText("====================================")
        self.output_text.appendPlainText("Job Success!")    

    def update_progress(self):
        self.progress_bar.setValue(self.progress_value)
        self.timer.stop()

def main():
    app = QApplication(sys.argv)
    mainWindow = MDGui()
    mainWindow.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()