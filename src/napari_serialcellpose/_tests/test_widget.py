from napari_serialcellpose import SerialWidget
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
import time
import os
import tempfile
import shutil

def test_load_single_image(make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)
    
    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_singlechannel/')
              
    widget.file_list.update_from_path(mypath)
    assert len(viewer.layers) == 0 
    widget.file_list.setCurrentRow(0)
    assert len(viewer.layers) == 1

def test_analyse_single_image_no_save(qtbot, make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)
    
    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_singlechannel/')
              
    widget.file_list.update_from_path(mypath)
    widget.file_list.setCurrentRow(0)

    # Check that selecting cyto2 displays diameter choice
    assert widget.spinbox_diameter.isVisible() is False
    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2')
    )
    # not working for some reason
    #assert widget.spinbox_diameter.isVisible() is True
    
    # set diameter and run segmentation
    widget.spinbox_diameter.setValue(70)
    widget._on_click_run_on_current()
    
    # check that segmentatio has been added, named 'mask' and results in 33 objects
    def check_layers():
        assert len(viewer.layers) == 2
        
    qtbot.waitUntil(check_layers, timeout=30000)
    assert viewer.layers[1].name == 'mask'
    assert viewer.layers[1].data.max() == 33

def test_analyse_single_image_save(qtbot, make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_multichannel')
              
    output_dir = Path(tempfile.mkdtemp())
    
    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)
    # Check that selecting cyto2 displays diameter choice
    assert widget.spinbox_diameter.isVisible() is False
    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2')
    )
    # not working for some reason
    #assert widget.spinbox_diameter.isVisible() is True
    
    # set diameter and run segmentation
    widget.spinbox_diameter.setValue(60)
    widget.qcbox_channel_to_segment.setCurrentIndex(2)
    widget.qcbox_channel_helper.setCurrentIndex(1)
    widget.check_clear_border.setChecked(False)

    widget.check_props['size'].setChecked(True)
    widget.check_props['intensity'].setChecked(True)
    widget.qcbox_channel_analysis.setCurrentRow(1)
    widget._on_click_run_on_current()

    def check_outputs():
        assert len(list(output_dir.glob('*mask.tif'))) == 1

    qtbot.waitUntil(check_outputs, timeout=30000)

    assert len(list(output_dir.joinpath('tables').glob('*_props.csv'))) == 1
    shutil.rmtree(output_dir)

def test_analyse_multi_image(qtbot, make_napari_viewer):
    """Test analysis of multiple images in a folder. No properties are analyzed."""
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)


    mypath = Path('src/napari_serialcellpose/_tests/data/multifile/')

    output_dir = Path(tempfile.mkdtemp())


    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.spinbox_diameter.setValue(70)
    widget._on_click_run_on_folder()

    def check_output():
        assert len(list(output_dir.glob('*mask.tif'))) == 4

    qtbot.waitUntil(check_output, timeout=30000)
    shutil.rmtree(output_dir)

def test_analyse_multi_image_props(qtbot, make_napari_viewer):

    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/multifile/')
    output_dir = Path(tempfile.mkdtemp())


    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.spinbox_diameter.setValue(70)
    widget.qcbox_channel_to_segment.setCurrentIndex(2)
    widget.qcbox_channel_helper.setCurrentIndex(1)
    widget.check_props['size'].setChecked(True)
    widget.check_props['intensity'].setChecked(True)
    widget.qcbox_channel_analysis.setCurrentRow(1)

    widget._on_click_run_on_folder()

    def check_outputs():
        assert len(list(output_dir.glob('*mask.tif'))) == 4
    
    qtbot.waitUntil(check_outputs, timeout=30000)
     # check that the properties are correct
    df = pd.read_csv(output_dir.joinpath(
      'tables',
        Path(widget.file_list.currentItem().text()).stem + '_props.csv')
    )
     # check number of columns in df
    assert df.shape[1] == 8   
     

def test_analyse_multichannels(qtbot, make_napari_viewer):
    """Test that multiple channels can be used for intensity measurements"""
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_multichannel/')
    output_dir = Path(tempfile.mkdtemp())

    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.spinbox_diameter.setValue(70)
    widget.qcbox_channel_to_segment.setCurrentIndex(2)
    widget.qcbox_channel_helper.setCurrentIndex(1)
    widget.check_props['size'].setChecked(True)
    widget.check_props['intensity'].setChecked(True)
    widget.qcbox_channel_analysis.item(0).setSelected(True)
    widget.qcbox_channel_analysis.item(1).setSelected(True)

    widget._on_click_run_on_folder()

    def check_outputs():
        assert len(list(output_dir.glob('*mask.tif'))) == 1

    qtbot.waitUntil(check_outputs, timeout=30000)

    # check that the properties are correct
    df = pd.read_csv(output_dir.joinpath(
        'tables',
        Path(widget.file_list.currentItem().text()).stem + '_props.csv'
    )
    )
    # check number of columns in df
    assert df.shape[1] == 11

def test_mask_loading(qtbot, make_napari_viewer):

    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/multifile/')
    output_dir = Path(tempfile.mkdtemp())


    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.spinbox_diameter.setValue(70)
    widget._on_click_run_on_current()

    # check that segmentation has been added
    def check_layers():
        assert len(viewer.layers) == 3
        
    qtbot.waitUntil(check_layers, timeout=30000)
    
    # check that when selecting the second file, we get only 2 channels and no mask
    widget.file_list.setCurrentRow(1)
    assert len(viewer.layers) == 2

    # check that when selecting the first file, we get  2 channels and a mask
    widget.file_list.setCurrentRow(0)
    assert len(viewer.layers) == 3

def test_analyse_single_image_options_yml(qtbot, make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_multichannel')
    yml_path = Path('src/napari_serialcellpose/_tests/my_options.yml')
              
    widget.file_list.update_from_path(mypath)
    widget.file_list.setCurrentRow(0)
    
    # set diameter and run segmentation
    widget.spinbox_diameter.setValue(70)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.qcbox_channel_to_segment.setCurrentIndex(2)
    widget.qcbox_channel_helper.setCurrentIndex(1)

    # overwrite options from yml file
    widget.options_file_path = yml_path

    widget._on_click_run_on_current()

    # check that segmentation has been added
    def check_layers():
        assert len(viewer.layers) == 3
        
    qtbot.waitUntil(check_layers, timeout=30000)
    # check that because of small diameter from yml file, we get only 7 elements
    assert viewer.layers[2].data.max() == 7