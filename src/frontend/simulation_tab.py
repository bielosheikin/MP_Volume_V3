from PyQt5.QtWidgets import QWidget, QFormLayout, QDoubleSpinBox, QPushButton

class SimulationParamsTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()

        self.time_step = QDoubleSpinBox()
        self.time_step.setDecimals(3)
        self.time_step.setRange(1e-6, 1.0)
        self.time_step.setValue(0.001)
        layout.addRow("Time Step (s):", self.time_step)

        self.total_time = QDoubleSpinBox()
        self.total_time.setDecimals(1)
        self.total_time.setRange(0.0, 10000.0)
        self.total_time.setValue(1000.0)
        layout.addRow("Total Simulation Time (s):", self.total_time)

        self.run_button = QPushButton("Run")
        layout.addWidget(self.run_button)

        self.setLayout(layout)

    def get_data(self):
        return {
            "time_step": self.time_step.value(),
            "total_time": self.total_time.value(),
        }