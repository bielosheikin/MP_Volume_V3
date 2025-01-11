from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class ResultsTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        self.setLayout(layout)

    def plot_results(self, histories_dict):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        plot_every = 10
        ax.plot(histories_dict['simulation_time'][::plot_every], histories_dict['Vesicle_pH'][::plot_every])
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Vesicle pH')
        ax.set_title('Simulation Results: Vesicle pH Over Time')
        self.canvas.draw()