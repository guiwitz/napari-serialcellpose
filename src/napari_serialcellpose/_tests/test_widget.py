from napari_serialcellpose import SerialWidget
import numpy as np

from pathlib import Path
import shutil

def test_load_widget(make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    assert type(widget) == SerialWidget