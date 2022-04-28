from qtpy.QtWidgets import (QWidget, QVBoxLayout,QFileDialog, QPushButton,
QSpinBox, QLabel, QGridLayout, QHBoxLayout, QGroupBox, QComboBox)
from .folder_list_widget import FolderList
from .serial_analysis import run_cellpose

from pathlib import Path
import skimage.io
from skimage.measure import regionprops_table
from cellpose import models
from napari_skimage_regionprops._table import TableWidget

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)



class SerialWidget(QWidget):
    # your QWidget.__init__ can optionally request the napari viewer instance
    # in one of two ways:
    # 1. use a parameter called `napari_viewer`, as done here
    # 2. use a type annotation of 'napari.viewer.Viewer' for any parameter
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer

        self.cellpose_model_path = None
        self.cellpose_model = None
        self.output_folder = None
        self.props_table = None

        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)

        self.main_layout.addWidget(QLabel("List of images"))

        self.file_list = FolderList(napari_viewer)
        self.main_layout.addWidget(self.file_list)

        self.folder_group = VHGroup('Folder selection')
        self.main_layout.addWidget(self.folder_group.gbox)

        self.btn_select_file_folder = QPushButton("Select data folder")
        self.folder_group.glayout.addWidget(self.btn_select_file_folder)

        self.btn_select_output_folder = QPushButton("Select output folder")
        self.folder_group.glayout.addWidget(self.btn_select_output_folder)

        self.qcbox_model_choice = QComboBox(visible=True)
        self.qcbox_model_choice.addItems(['custom', 'cyto', 'cyto2', 'nuclei'])
        self.folder_group.glayout.addWidget(self.qcbox_model_choice)

        self.btn_select_cellpose_model = QPushButton("Select cellpose model")
        self.folder_group.glayout.addWidget(self.btn_select_cellpose_model)

        self.run_group = VHGroup('Run analysis')
        self.main_layout.addWidget(self.run_group.gbox)

        self.btn_run_on_current = QPushButton("Run on current image")
        self.run_group.glayout.addWidget(self.btn_run_on_current)

        self.btn_run_on_folder = QPushButton("Run on folder")
        self.run_group.glayout.addWidget(self.btn_run_on_folder)

        self.options_group = VHGroup('Options', orientation='G')
        self.main_layout.addWidget(self.options_group.gbox)

        self.options_group.glayout.addWidget(QLabel("Batch size"), 0, 0, 1, 1)
        self.spinbox_batch_size = QSpinBox()
        self.spinbox_batch_size.setValue(3)
        self.options_group.glayout.addWidget(self.spinbox_batch_size, 0, 1, 1, 1)

        self.options_group.glayout.addWidget(QLabel("Rescaling factor"), 1, 0, 1, 1)
        self.spinbox_rescaling = QSpinBox()
        self.spinbox_rescaling.setValue(1)
        self.options_group.glayout.addWidget(self.spinbox_rescaling, 1, 1, 1, 1)

        self.diameter_label = QLabel("Diameter", visible=False)
        self.options_group.glayout.addWidget(self.diameter_label, 2, 0, 1, 1)
        self.spinbox_diameter = QSpinBox(visible=False)
        self.spinbox_diameter.setValue(30)
        self.options_group.glayout.addWidget(self.spinbox_diameter, 2, 1, 1, 1)

        self.add_connections()

        self.plot_group = VHGroup('Plots')
        self.main_layout.addWidget(self.plot_group.gbox)

        self.sc = MplCanvas(self, width=5, height=4, dpi=100)
        self.plot_group.glayout.addWidget(self.sc)
        self.toolbar = NavigationToolbar(self.sc, self)

        self.plot_group.glayout.addWidget(self.toolbar)
        self.plot_group.glayout.addWidget(self.sc)
        

    def add_connections(self):

        self.btn_select_file_folder.clicked.connect(self._on_click_select_file_folder)
        self.btn_select_cellpose_model.clicked.connect(self._on_click_select_cellpose_model)
        self.btn_select_output_folder.clicked.connect(self._on_click_select_output_folder)
        self.file_list.currentItemChanged.connect(self._on_select_file)
        self.btn_run_on_current.clicked.connect(self._on_click_run_on_current)
        self.btn_run_on_folder.clicked.connect(self._on_click_run_on_folder)
        self.qcbox_model_choice.currentTextChanged.connect(self._on_change_modeltype)
    

    def open_file(self):
        """Open file selected in list. Returns True if file was opened."""
        
        # clear existing layers.
        #self.clear_layers()
        self.viewer.layers.clear()

        # if file list is empty stop here
        if self.file_list.currentItem() is None:
            return False
        
        # open image
        image_name = self.file_list.currentItem().text()
        image_path = self.file_folder.joinpath(image_name)
        self.viewer.open(image_path)

        if self.output_folder is not None:
            mask_path = Path(self.output_folder).joinpath(image_path.stem+'_mask.tif')
            if mask_path.exists():
                mask = skimage.io.imread(mask_path)
                self.viewer.add_labels(mask, name='mask')
                self.add_table_props()

    def _on_click_select_file_folder(self):
        """Interactively select folder to analyze"""

        self.file_folder = Path(str(QFileDialog.getExistingDirectory(self, "Select Directory")))
        self.file_list.update_from_path(self.file_folder)

    def _on_click_select_output_folder (self):
        """Interactively select folder to analyze"""

        self.output_folder = Path(str(QFileDialog.getExistingDirectory(self, "Select Directory")))

    def _on_click_select_cellpose_model(self):
        """Interactively select cellpose model"""

        self.cellpose_model_path = QFileDialog.getOpenFileName(self, "Select model file")[0]
    
        
    def _on_select_file(self, current_item, previous_item):
        success = self.open_file()
        if not success:
            return

    def _on_click_run_on_current(self):
        """Run cellpose on current image"""

        model_type = self.qcbox_model_choice.currentText()

        if self.output_folder is None:
            self._on_click_select_output_folder()
        if (self.cellpose_model_path is None) and (model_type == 'custom'):
            self._on_click_select_cellpose_model()

        image_path = self.file_list.folder_path.joinpath(self.file_list.currentItem().text())
        
        diameter = None
        if self.qcbox_model_choice.currentText() == 'custom':
            self.cellpose_model = models.CellposeModel(pretrained_model=self.cellpose_model_path)
        else:
            self.cellpose_model = models.Cellpose(model_type=model_type)
            diameter = self.spinbox_diameter.value()
        
        # run cellpose
        segmented = run_cellpose(
            image_path=image_path,
            cellpose_model=self.cellpose_model,
            output_path=self.output_folder,
            scaling_factor=self.spinbox_rescaling.value(),
            diameter=diameter
        )
        self.viewer.add_labels(segmented, name='mask')
        self.add_table_props()

    def _on_click_run_on_folder(self):
        """Run cellpose on all images in folder"""

        if self.output_folder is None:
            self._on_click_select_output_folder()
        if self.cellpose_model_path is None:
            self._on_click_select_cellpose_model()
        
        file_list = [self.file_list.item(x).text() for x in range(self.file_list.count())]
        file_list = [f for f in file_list if f[0] != '.']
        file_list = [self.file_folder.joinpath(x) for x in file_list]

        n = self.spinbox_batch_size.value()
        file_list_partition = [file_list[i:i + n] for i in range(0, len(file_list), n)]

        self.cellpose_model = models.CellposeModel(pretrained_model=self.cellpose_model_path)

        for batch in file_list_partition:
            run_cellpose(
                image_path=batch,
                cellpose_model=self.cellpose_model,
                output_path=self.output_folder,
                scaling_factor=1.0,
            )

    def _on_change_modeltype(self):
        "if selecting non-custom model, show diameter box"

        if self.qcbox_model_choice.currentText() != 'custom':
            self.diameter_label.setVisible(True)
            self.spinbox_diameter.setVisible(True)
        else:
            self.diameter_label.setVisible(False)
            self.spinbox_diameter.setVisible(False)
        
    def add_table_props(self):
        """Add table with properties of segmented cells"""

        props = regionprops_table(self.viewer.layers['mask'].data,
                          properties=('label', 'area', 'eccentricity', 'solidity', 'feret_diameter_max'))
        self.viewer.layers['mask'].properties = props
        if self.props_table is None:
            self.props_table = TableWidget(layer=self.viewer.layers['mask'])
            self.main_layout.addWidget(self.props_table)
        else:
            self.props_table._layer = self.viewer.layers['mask'] 
            self.props_table.update_content()

        self.sc.axes.clear()       
        self.sc.axes.hist(props['feret_diameter_max'])
        self.sc.axes.figure.canvas.draw()



    

class VHGroup():
    """Group box with specific layout
    Parameters
    ----------
    name: str
        Name of the group box
    orientation: str
        'V' for vertical, 'H' for horizontal, 'G' for grid
    """

    def __init__(self, name, orientation='V'):
        self.gbox = QGroupBox(name)
        if orientation=='V':
            self.glayout = QVBoxLayout()
        elif orientation=='H':
            self.glayout = QHBoxLayout()
        elif orientation=='G':
            self.glayout = QGridLayout()
        else:
            raise Exception(f"Unknown orientation {orientation}") 

        self.gbox.setLayout(self.glayout)