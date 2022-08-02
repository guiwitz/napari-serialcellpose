from napari_serialcellpose import SerialWidget
import numpy as np
import pandas as pd

from pathlib import Path
import shutil

def test_load_single_image(make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_singlechannel/')
              
    widget.file_list.update_from_path(mypath)
    assert len(viewer.layers) == 0 
    widget.file_list.setCurrentRow(0)
    assert len(viewer.layers) == 1

def test_analyse_single_image_no_save(make_napari_viewer):
    
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
    assert len(viewer.layers) == 2
    assert viewer.layers[1].name == 'mask'
    assert viewer.layers[1].data.max() == 33

def test_analyse_single_image_save(make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/single_file_multichannel')
              
    output_dir = Path('src/napari_serialcellpose/_tests/data/analyzed_single')
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

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
    widget.qcbox_channel_analysis.setCurrentIndex(2)

    widget._on_click_run_on_current()

    assert len(list(output_dir.glob('*mask.tif'))) == 1
    assert len(list(output_dir.joinpath('tables').glob('*_props.csv'))) == 1

def test_analyse_multi_image(make_napari_viewer):

    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)


    mypath = Path('src/napari_serialcellpose/_tests/data/multifile/')
    output_dir = Path('src/napari_serialcellpose/_tests/data/analyzed_multiple')
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.spinbox_diameter.setValue(70)
    widget._on_click_run_on_current()

    assert len(list(output_dir.glob('*mask.tif'))) == 1

    widget._on_click_run_on_folder()
    assert len(list(output_dir.glob('*mask.tif'))) == 4

def test_analyse_multi_image_props(make_napari_viewer):

    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/multifile/')
    output_dir = Path('src/napari_serialcellpose/_tests/data/analyzed_multiple3')
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

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
    widget.qcbox_channel_analysis.setCurrentIndex(2)

    widget._on_click_run_on_folder()
    assert len(list(output_dir.glob('*mask.tif'))) == 4

    # check that the properties are correct
    df = pd.read_csv(output_dir.joinpath(
        'tables',
        Path(widget.file_list.currentItem().text()).stem + '_props.csv'
    )
    )
    assert df.shape[1] == 9

def test_mask_loading(make_napari_viewer):

    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/multifile/')
    output_dir = Path('src/napari_serialcellpose/_tests/data/analyzed_multiple2')
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

    widget.file_list.update_from_path(mypath)
    widget.output_folder = output_dir
    widget.file_list.setCurrentRow(0)

    widget.qcbox_model_choice.setCurrentIndex(
        [widget.qcbox_model_choice.itemText(i) for i in range(widget.qcbox_model_choice.count())].index('cyto2'))
    widget.spinbox_diameter.setValue(70)
    widget._on_click_run_on_current()

    # check that when selecting the second file, we get only 2 channels and no mask
    widget.file_list.setCurrentRow(1)
    assert len(viewer.layers) == 2

    # check that when selecting the first file, we get  2 channels and a mask
    widget.file_list.setCurrentRow(0)
    assert len(viewer.layers) == 3

def test_analyse_single_image_options_yml(make_napari_viewer):
    
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

    # check that because of small diameter from yml file, we get only 5 elements
    assert viewer.layers[2].data.max() == 7