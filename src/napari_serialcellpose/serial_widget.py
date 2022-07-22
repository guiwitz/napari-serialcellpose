from qtpy.QtWidgets import (QWidget, QVBoxLayout,QFileDialog, QPushButton,
QSpinBox, QDoubleSpinBox, QLabel, QGridLayout, QHBoxLayout, QGroupBox, QComboBox, QTabWidget,
QCheckBox, QLayout)
from qtpy.QtCore import Qt
import magicgui.widgets
from napari.layers import Image

from .folder_list_widget import FolderList
from .serial_analysis import run_cellpose, load_props, load_allprops

from pathlib import Path
import skimage.io
import numpy as np
from cellpose import models
from napari_skimage_regionprops._table import TableWidget

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


class SerialWidget(QWidget):
    """
    Implementation of a napari widget allowing to select a folder filled with images and 
    segment them using cellpose.
    """
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer

        self.cellpose_model_path = None
        self.cellpose_model = None
        self.output_folder = None
        self.props_table = None
        self.allprops = None

        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)

        self.tabs = QTabWidget()
        self.main_layout.addWidget(self.tabs)

        # segmentation tab
        self.segmentation = QWidget()
        self._segmentation_layout = QVBoxLayout()
        self.segmentation.setLayout(self._segmentation_layout)
        self.tabs.addTab(self.segmentation, 'Segmentation')

        # channels tab
        self.options_tab = QWidget()
        self._options_tab_layout = QVBoxLayout()
        self.options_tab.setLayout(self._options_tab_layout)
        self.tabs.addTab(self.options_tab, 'Options')

        # properties tab
        self.properties = QWidget()
        self._properties_layout = QVBoxLayout()
        self.properties.setLayout(self._properties_layout)
        self.tabs.addTab(self.properties, 'Properties')

        # summary tab
        self.summary = QWidget()
        self._summary_layout = QVBoxLayout()
        self.summary.setLayout(self._summary_layout)
        self.tabs.addTab(self.summary, 'Summary')

        self._segmentation_layout.addWidget(QLabel("List of images"))

        self.file_list = FolderList(napari_viewer)
        self._segmentation_layout.addWidget(self.file_list)

        self.folder_group = VHGroup('Folder selection')
        self._segmentation_layout.addWidget(self.folder_group.gbox)

        self.btn_select_file_folder = QPushButton("Select data folder")
        self.folder_group.glayout.addWidget(self.btn_select_file_folder)

        self.btn_select_output_folder = QPushButton("Select output folder")
        self.folder_group.glayout.addWidget(self.btn_select_output_folder)

        self.qcbox_model_choice = QComboBox(visible=True)
        self.qcbox_model_choice.addItems([
            'custom', 'cyto', 'cyto2', 'nuclei', 'tissuenet'])
        self.folder_group.glayout.addWidget(self.qcbox_model_choice)

        self.btn_select_cellpose_model = QPushButton("Select custom cellpose model file")
        self.folder_group.glayout.addWidget(self.btn_select_cellpose_model)

        self.run_group = VHGroup('Run analysis')
        self._segmentation_layout.addWidget(self.run_group.gbox)

        self.btn_run_on_current = QPushButton("Run on current image")
        self.run_group.glayout.addWidget(self.btn_run_on_current)

        self.btn_run_on_folder = QPushButton("Run on folder")
        self.run_group.glayout.addWidget(self.btn_run_on_folder)

        self.check_usegpu = QCheckBox('Use GPU')
        self.run_group.glayout.addWidget(self.check_usegpu)

        self._options_tab_layout.setAlignment(Qt.AlignTop)
        self.options_group = VHGroup('Options', orientation='G')
        self._options_tab_layout.addWidget(self.options_group.gbox)

        self.options_group.glayout.addWidget(QLabel("Batch size"), 0, 0, 1, 1)
        self.spinbox_batch_size = QSpinBox()
        self.spinbox_batch_size.setValue(3)
        self.options_group.glayout.addWidget(self.spinbox_batch_size, 0, 1, 1, 1)

        # Not yet implemented
        # self.options_group.glayout.addWidget(QLabel("Rescaling factor"), 1, 0, 1, 1)
        # self.spinbox_rescaling = QSpinBox()
        # self.spinbox_rescaling.setValue(1)
        # self.options_group.glayout.addWidget(self.spinbox_rescaling, 1, 1, 1, 1)

        self.flow_threshold_label = QLabel("Flow threshold")
        self.options_group.glayout.addWidget(self.flow_threshold_label, 3, 0, 1, 1)
        self.flow_threshold = QDoubleSpinBox()
        self.flow_threshold.setSingleStep(0.1)
        self.flow_threshold.setValue(0.4)
        self.options_group.glayout.addWidget(self.flow_threshold, 3, 1, 1, 1)

        self.cellprob_threshold_label = QLabel("Cell probability threshold")
        self.options_group.glayout.addWidget(self.cellprob_threshold_label, 4, 0, 1, 1)
        self.cellprob_threshold = QDoubleSpinBox()
        self.cellprob_threshold.setSingleStep(0.1)
        self.cellprob_threshold.setValue(0.0)
        self.options_group.glayout.addWidget(self.cellprob_threshold, 4, 1, 1, 1)

        self.check_clear_border = QCheckBox('Clear labels on border')
        self.check_clear_border.setChecked(True)
        self.options_group.glayout.addWidget(self.check_clear_border) 

        self.mainoptions_group = VHGroup('Main options', orientation='G')
        self._segmentation_layout.addWidget(self.mainoptions_group.gbox)

        self.qcbox_channel_to_segment = QComboBox()
        self.mainoptions_group.glayout.addWidget(QLabel('Channel to segment'), 0, 0, 1, 1)
        self.mainoptions_group.glayout.addWidget(self.qcbox_channel_to_segment, 0,1,1,1)
        self.mainoptions_group.glayout.addWidget(QLabel('Helper channel'), 1, 0, 1, 1)
        self.qcbox_channel_helper = QComboBox()
        self.mainoptions_group.glayout.addWidget(self.qcbox_channel_helper, 1,1,1,1)
        self.diameter_label = QLabel("Diameter", visible=False)
        self.mainoptions_group.glayout.addWidget(self.diameter_label, 2, 0, 1, 1)
        self.spinbox_diameter = QSpinBox(visible=False)
        self.spinbox_diameter.setValue(30)
        self.spinbox_diameter.setMaximum(1000)
        self.mainoptions_group.glayout.addWidget(self.spinbox_diameter, 2, 1, 1, 1)

        self.plot_group = VHGroup('Plots')
        self._properties_layout.addWidget(self.plot_group.gbox)

        self.sc = MplCanvas(self, row=1, col=2, width=6, height=4, dpi=100)
        self.toolbar = NavigationToolbar(self.sc, self)
        self.plot_group.glayout.addWidget(self.toolbar)
        self.plot_group.glayout.addWidget(self.sc)

        self.sc_sum = MplCanvas(self, row=1, col=2, width=6, height=4, dpi=100)
        self.toolbar_sum = NavigationToolbar(self.sc_sum, self)
        self._summary_layout.addWidget(self.toolbar_sum)
        self._summary_layout.addWidget(self.sc_sum)
        self.btn_load_summary = QPushButton("Load summary")
        self._summary_layout.addWidget(self.btn_load_summary)

        self.eccentricity_slider = magicgui.widgets.FloatSlider(min=0, max=1, step=0.01, value=1)
        self._summary_layout.addWidget(self.eccentricity_slider.native)

        self.add_connections()

    def add_connections(self):
        """Add callbacks"""

        self.btn_select_file_folder.clicked.connect(self._on_click_select_file_folder)
        self.btn_select_cellpose_model.clicked.connect(self._on_click_select_cellpose_model)
        self.btn_select_output_folder.clicked.connect(self._on_click_select_output_folder)
        self.file_list.currentItemChanged.connect(self._on_select_file)
        self.btn_run_on_current.clicked.connect(self._on_click_run_on_current)
        self.btn_run_on_folder.clicked.connect(self._on_click_run_on_folder)
        self.qcbox_model_choice.currentTextChanged.connect(self._on_change_modeltype)
        self.btn_load_summary.clicked.connect(self._on_click_load_summary)
        self.eccentricity_slider.changed.connect(self.update_eccentricity)
        self.viewer.layers.events.connect(self._on_change_layers)

    def open_file(self):
        """Open file selected in list. Returns True if file was opened."""

        # clear existing layers.
        self.viewer.layers.clear()
        
        # if file list is empty stop here
        if self.file_list.currentItem() is None:
            return False
        
        # open image
        image_name = self.file_list.currentItem().text()
        image_path = self.file_list.folder_path.joinpath(image_name)

        self.viewer.open(image_path)

        if self.output_folder is not None:
            mask_path = Path(self.output_folder).joinpath(image_path.stem+'_mask.tif')
            if mask_path.exists():
                mask = skimage.io.imread(mask_path)
                self.viewer.add_labels(mask, name='mask')
                props = load_props(self.output_folder, image_name)
                self.add_table_props(props)

        return True

    def _on_click_select_file_folder(self):
        """Interactively select folder to analyze"""

        file_folder = Path(str(QFileDialog.getExistingDirectory(self, "Select Directory")))
        self.file_list.update_from_path(file_folder)

    def _on_click_select_output_folder (self):
        """Interactively select folder where to save results"""

        self.output_folder = Path(str(QFileDialog.getExistingDirectory(self, "Select Directory")))

    def _on_click_select_cellpose_model(self):
        """Interactively select cellpose model"""

        self.cellpose_model_path = QFileDialog.getOpenFileName(self, "Select model file")[0]
    
    def _on_select_file(self, current_item, previous_item):
        
        success = self.open_file()
        if not success:
            return False

    def _on_click_run_on_current(self):
        """Run cellpose on current image"""

        model_type = self.qcbox_model_choice.currentText()
        self.output_and_model_check(choose_output=False)

        image_path = self.file_list.folder_path.joinpath(self.file_list.currentItem().text())
        
        self.cellpose_model, diameter = self.get_cellpose_model(model_type=model_type)
        
        channel_to_segment, channel_helper = self.get_channels_to_use()

        # run cellpose
        segmented = run_cellpose(
            image_path=image_path,
            cellpose_model=self.cellpose_model,
            output_path=self.output_folder,
            diameter=diameter,
            flow_threshold=self.flow_threshold.value(),
            cellprob_threshold=self.cellprob_threshold.value(),
            clear_border=self.check_clear_border.isChecked(),
            channel_to_segment=channel_to_segment,
            channel_helper=channel_helper,
        )
        self.viewer.add_labels(segmented, name='mask')
        if self.output_folder is not None:
            props = load_props(self.output_folder, image_path)
            self.add_table_props(props)

    def _on_click_run_on_folder(self):
        """Run cellpose on all images in folder"""

        model_type = self.qcbox_model_choice.currentText()
        self.output_and_model_check()
        
        file_list = [self.file_list.item(x).text() for x in range(self.file_list.count())]
        file_list = [f for f in file_list if f[0] != '.']
        file_list = [self.file_list.folder_path.joinpath(x) for x in file_list]

        n = self.spinbox_batch_size.value()
        file_list_partition = [file_list[i:i + n] for i in range(0, len(file_list), n)]

        channel_to_segment, channel_helper = self.get_channels_to_use()
        self.cellpose_model, diameter = self.get_cellpose_model(model_type=model_type)

        for batch in file_list_partition:
            run_cellpose(
                image_path=batch,
                cellpose_model=self.cellpose_model,
                output_path=self.output_folder,
                diameter=diameter,
                flow_threshold=self.flow_threshold.value(),
                cellprob_threshold=self.cellprob_threshold.value(),
                clear_border=self.check_clear_border.isChecked(),
                channel_to_segment=channel_to_segment,
                channel_helper=channel_helper
            )

    def get_channels_to_use(self):
        """Translate selected channels in QCombox into indices.
        As the first choice is None, channels are already incremented by one
        as expected by cellpose"""
        
        channel_to_segment = 0
        channel_helper = 0
        if self.qcbox_channel_to_segment.currentText() != 'None':
            print(f'self.qcbox_channel_to_segment.currentIndex(): {self.qcbox_channel_to_segment.currentIndex()}')
            channel_to_segment = self.qcbox_channel_to_segment.currentIndex()
        if self.qcbox_channel_helper.currentText() != 'None':
            channel_helper = self.qcbox_channel_helper.currentIndex()
        
        return channel_to_segment, channel_helper

    def output_and_model_check(self, choose_output=True):
        """Check if output folder and model are set"""

        model_type = self.qcbox_model_choice.currentText()
        if (self.output_folder is None) and (choose_output):
            self._on_click_select_output_folder()
        if (self.cellpose_model_path is None) and (model_type == 'custom'):
            self._on_click_select_cellpose_model()

    def get_cellpose_model(self, model_type):
        """Get cellpose model. For non-custom model provide a model name
        in model_type. For custom models the returned diameter is None, otherwise
        it is the GUI value."""

        diameter = None
        if self.qcbox_model_choice.currentText() == 'custom':
            cellpose_model = models.CellposeModel(
                gpu=self.check_usegpu.isChecked(),
                pretrained_model=self.cellpose_model_path)
        else:
            cellpose_model = models.Cellpose(
                gpu=self.check_usegpu.isChecked(),
                model_type=model_type)
            diameter = self.spinbox_diameter.value()
        
        return cellpose_model, diameter

    def _on_change_layers(self):

        self.qcbox_channel_to_segment.clear()
        self.qcbox_channel_to_segment.addItems(['None']+[x.name for x in self.viewer.layers if isinstance(x, Image)])
        self.qcbox_channel_helper.clear()
        self.qcbox_channel_helper.addItems(['None']+[x.name for x in self.viewer.layers if isinstance(x, Image)])

    def _on_click_load_summary(self):
        """Load summary from folder"""

        self.allprops = load_allprops(self.output_folder)
        props = self.allprops[self.allprops.eccentricity < self.eccentricity_slider.value]
        prop_names = ['eccentricity', 'feret_diameter_max']
        for i in range(len(prop_names)):
            self.sc_sum.ax[0,i].clear()
            self.sc_sum.ax[0,i].hist(props[prop_names[i]], rwidth=0.85)
            self.sc_sum.ax[0,i].figure.canvas.draw()
            self.sc_sum.ax[0,i].tick_params(colors='black',labelsize=12)
            self.sc_sum.ax[0,i].set_title(prop_names[i], fontsize=15, color='black')

    def _on_change_modeltype(self):
        "if selecting non-custom model, show diameter box"

        if self.qcbox_model_choice.currentText() != 'custom':
            self.diameter_label.setVisible(True)
            self.spinbox_diameter.setVisible(True)
        else:
            self.diameter_label.setVisible(False)
            self.spinbox_diameter.setVisible(False)
        
    def add_table_props(self, props):
        """Add table with properties of segmented cells"""

        self.viewer.layers['mask'].properties = props
        if self.props_table is None:
            self.props_table = TableWidget(layer=self.viewer.layers['mask'])
            self._properties_layout.addWidget(self.props_table)
        else:
            self.props_table._layer = self.viewer.layers['mask'] 
            self.props_table.update_content()

        #self.sc.axes.clear()       
        #self.sc.axes.hist(props['feret_diameter_max'], rwidth=0.85)
        #self.sc.axes.figure.canvas.draw()

        prop_names = ['eccentricity', 'feret_diameter_max']
        for i in range(len(prop_names)):
            self.sc.ax[0,i].clear()
            self.sc.ax[0,i].hist(props[prop_names[i]], rwidth=0.85)
            self.sc.ax[0,i].figure.canvas.draw()
            self.sc.ax[0,i].tick_params(colors='black',labelsize=12)
            self.sc.ax[0,i].set_title(prop_names[i], fontsize=15, color='black')

    def update_eccentricity(self, vakue):
        """Update eccentricity plot"""

        if self.allprops is None:
            self._on_click_load_summary()
        else:
            props = self.allprops[self.allprops.eccentricity < self.eccentricity_slider.value]
            prop_names = ['eccentricity', 'feret_diameter_max']
            for i in range(len(prop_names)):
                self.sc_sum.ax[0,i].clear()
                self.sc_sum.ax[0,i].hist(props[prop_names[i]], rwidth=0.85)
                self.sc_sum.ax[0,i].figure.canvas.draw()
                self.sc_sum.ax[0,i].tick_params(colors='black',labelsize=12)
                self.sc_sum.ax[0,i].set_title(prop_names[i], fontsize=15, color='black')


class VHGroup():
    """Group box with specific layout.

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

class MplCanvas(FigureCanvasQTAgg):
    """Matplotlib QT canvas for grid of plots.
    
    Parameters
    ----------
    col: int
        number of columns
    row: int
        number of rows
    width: int
        width of the canvas
    height: int
        height of the canvas
    dpi: int
        dpi of the canvas
    """

    def __init__(self, parent=None, col=1, row=1, width=6, height=4, dpi=100):
        self.ax=np.array([[None]*col]*row)
        fig = Figure(figsize=(width, height), dpi=dpi, facecolor='white')
        count = 1
        for i in range(row):
            for j in range(col):
                self.ax[i,j] = fig.add_subplot(row, col, count)
                count+=1
        super(MplCanvas, self).__init__(fig)