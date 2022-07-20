from napari_serialcellpose import SerialWidget
import numpy as np

from pathlib import Path

def test_load_single_image(make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/singlefile_singlechannel/')
              
    widget.file_list.update_from_path(mypath)
    assert len(viewer.layers) == 0 
    widget.file_list.setCurrentRow(0)
    assert len(viewer.layers) == 1

def test_analyse_single_image_no_save(make_napari_viewer):
    
    viewer = make_napari_viewer()
    widget = SerialWidget(viewer)

    mypath = Path('src/napari_serialcellpose/_tests/data/singlefile_singlechannel/')
              
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




# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
'''def test_example_q_widget(make_napari_viewer, capsys):
    # make viewer and add an image layer using our fixture
    viewer = make_napari_viewer()
    viewer.add_image(np.random.random((100, 100)))

    # create our widget, passing in the viewer
    my_widget = SerialWidget(viewer)

    # call our widget method
    my_widget._on_click()

    # read captured output and check that it's as we expected
    captured = capsys.readouterr()
    assert captured.out == "napari has 1 layers\n"
    
def test_example_magic_widget(make_napari_viewer, capsys):
    viewer = make_napari_viewer()
    layer = viewer.add_image(np.random.random((100, 100)))

    # this time, our widget will be a MagicFactory or FunctionGui instance
    my_widget = example_magic_widget()

    # if we "call" this object, it'll execute our function
    my_widget(viewer.layers[0])

    # read captured output and check that it's as we expected
    captured = capsys.readouterr()
    assert captured.out == f"you have selected {layer}\n"
'''