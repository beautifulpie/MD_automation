import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QPushButton, QFileDialog, QPlainTextEdit, QWidget, QProgressBar
from PyQt5.QtCore import Qt, QTimer
import md_automation as md
import subprocess

class MDGui(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle('MD Automation GUI')
        self.setGeometry(100, 100, 800, 600)

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

        # Progress Bar
        self.progress_bar = QProgressBar(self)
        self.layout.addWidget(self.progress_bar)

        # Output Text
        self.output_text = QPlainTextEdit(self)
        self.layout.addWidget(self.output_text)

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

        # Timer for Progress Bar Update
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_progress)
        self.progress_value = 0

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
        
        print("Get pdb file")

    def get_itp_files(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("ITP Files (*.itp);;All Files (*)")
        file_dialog.setDefaultSuffix("itp")
        file_paths, _ = file_dialog.getOpenFileNames(self, "Select ITP Files", "", "ITP Files (*.itp);;All Files (*)")

        if file_paths:
            self.selected_itp_files = file_paths
            self.output_text.appendPlainText(f"Selected ITP files: {', '.join(self.selected_itp_files)}")
        
        print("Get itp files")

    def get_md_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (*.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select md.mdp File", "", "MDP Files (*.mdp);;All Files (*)")

        if file_path:
            self.selected_md_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected md.mdp file: {self.selected_md_mdp_file}")
        
        print("Get md.mdp file")

    def get_minim_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (minim.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select minim.mdp File", "", "MDP Files (minim.mdp);;All Files (*)")

        if file_path:
            self.selected_minim_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected minim.mdp file: {self.selected_minim_mdp_file}")
        print("Get minim.mdp file")

    def get_npt_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (npt.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select npt.mdp File", "", "MDP Files (npt.mdp);;All Files (*)")

        if file_path:
            self.selected_npt_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected npt.mdp file: {self.selected_npt_mdp_file}")
        print("Get npt.mdp file")

    def get_nvt_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (nvt.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select nvt.mdp File", "", "MDP Files (nvt.mdp);;All Files (*)")

        if file_path:
            self.selected_nvt_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected nvt.mdp file: {self.selected_nvt_mdp_file}")
        print("Get nvt.mdp file")


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
        self.output_text.appendPlainText("MD_automation")
        self.output_text.appendPlainText("Maed by Jung Youngwoo")
        self.output_text.appendPlainText("====================================")

        commands = [
            # Use self.selected_files to get the paths of the selected files
            # Example: self.selected_files['.pdb'], self.selected_files['.mdp'], etc.
        ]

        total_commands = len(commands)
        self.progress_bar.setRange(0, total_commands - 1)
        self.progress_bar.setValue(0)

        # Run the MD automation logic
        for i, command in enumerate(commands):
            self.output_text.appendPlainText(f"Running command: {command}")
            try:
                subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.output_text.appendPlainText("작업 완료")
            except subprocess.CalledProcessError as e:
                self.output_text.appendPlainText(f"Error: {e}")

            # Update progress bar value
            self.progress_value = i
            self.timer.start(100)  # Start timer to update progress bar with a slight delay

        # Display success message
        self.label.setText("MD Automation completed successfully!")
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
