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
        self.setGeometry(100, 100, 600, 400)

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout(self.central_widget)

        self.label = QLabel('MD Automation', self)
        self.label.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(self.label)

        self.pdb_file_button = QPushButton('Select PDB File', self)
        self.pdb_file_button.clicked.connect(self.get_pdb_file)
        self.layout.addWidget(self.pdb_file_button)

        self.mdp_file_button = QPushButton('Select MDP File', self)
        self.mdp_file_button.clicked.connect(self.get_mdp_file)
        self.layout.addWidget(self.mdp_file_button)

        self.progress_bar = QProgressBar(self)
        self.layout.addWidget(self.progress_bar)

        self.output_text = QPlainTextEdit(self)
        self.layout.addWidget(self.output_text)

        self.run_button = QPushButton('Run MD Automation', self)
        self.run_button.clicked.connect(self.run_md_automation)
        self.layout.addWidget(self.run_button)

        self.selected_pdb_file = ""
        self.selected_mdp_file = ""

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_progress)
        self.progress_value = 0

    def get_pdb_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("PDB Files (*.pdb);;All Files (*)")
        file_dialog.setDefaultSuffix("pdb")
        file_path, _ = file_dialog.getOpenFileName(self, "Select PDB File", "", "PDB Files (*.pdb);;All Files (*)")

        if file_path:
            self.selected_pdb_file = file_path
            self.output_text.appendPlainText(f"Selected PDB file: {self.selected_pdb_file}")

    def get_mdp_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("MDP Files (*.mdp);;All Files (*)")
        file_dialog.setDefaultSuffix("mdp")
        file_path, _ = file_dialog.getOpenFileName(self, "Select MDP File", "", "MDP Files (*.mdp);;All Files (*)")

        if file_path:
            self.selected_mdp_file = file_path
            self.output_text.appendPlainText(f"Selected MDP file: {self.selected_mdp_file}")

    def run_md_automation(self):
        if not self.selected_pdb_file or not self.selected_mdp_file:
            self.output_text.appendPlainText("Please select both PDB and MDP files.")
            return

        # Run MD automation with the selected files
        self.label.setText("MD Automation is running... Please wait.")
        self.output_text.appendPlainText("====================================")
        self.output_text.appendPlainText('\033[1m' + "MD_automation" + '\0333')

        commands = [
            # ... (replace with your existing commands using self.selected_pdb_file and self.selected_mdp_file)
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
