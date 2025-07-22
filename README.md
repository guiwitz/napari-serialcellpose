# napari-serialcellpose

[![License](https://img.shields.io/pypi/l/napari-serialcellpose.svg?color=green)](https://github.com/guiwitz/napari-serialcellpose/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/napari-serialcellpose.svg?color=green)](https://pypi.org/project/napari-serialcellpose)
[![Python Version](https://img.shields.io/pypi/pyversions/napari-serialcellpose.svg?color=green)](https://python.org)
[![tests](https://github.com/guiwitz/napari-serialcellpose/workflows/tests/badge.svg)](https://github.com/guiwitz/napari-serialcellpose/actions)
[![codecov](https://codecov.io/gh/guiwitz/napari-serialcellpose/branch/main/graph/badge.svg)](https://codecov.io/gh/guiwitz/napari-serialcellpose)
[![napari hub](https://img.shields.io/endpoint?url=https://api.napari-hub.org/shields/napari-serialcellpose)](https://napari-hub.org/plugins/napari-serialcellpose)

This napari plugin allows you to segment single images or series of images using built-in or custom Cellpose models as well as to analyze the properties of these segmented regions ("region properties"). Properties can be visualized for a single image or a complete experiment in the form of histograms that can also be filtered (e.g. based on area size, mean intensity etc.) Thanks to the [napari-skimage-regionprops](https://github.com/haesleinhuepf/napari-skimage-regionprops) plugin, properties of segmented objects can be interactively explored at a single object level.

## Main goal

The main goal of this plugin is to simplify the classical image processing pipeline of image segmentation followed by region analysis via Cellpose. It allows to quickly get a quantification of a set of images without the need for any scripting.

## Installation

In order to use this plugin, whe highly recommend to create a specific environment and to install the required software in it. You can create a conda environment using:

    conda create -n serialcellpose python=3.12 napari pyqt -c conda-forge

Then activate it and install the plugin:
    
    conda activate serialcellpose
    pip install napari-serialcellpose

Alternatively, you can install the plugin via the napari plugin manager.

## Note on versions

Versions older than 0.3.0 depend on the aicsimageio importer which is not supported anymore and which will likely lead to errors when trying to run the plugin. Version 0.3.0 moved to the [bioio](https://github.com/bioio-devs/bioio) importer and some improvements were made to the data handling but the core functionality is the same. **Version 0.3.0 is the last version that will support pre-SAM Cellpose models i.e. Cellpose=4.0. If you want to use Cellpose=4.0, please use version napari-serialcellpose=0.4.0 or later.** You can pick versions in the napari plugin manger or directly via pip e.g.:

    pip install napari-serialcellpose==0.3.0

### Potential issue with PyTorch

Cellpose and therefore the plugin and napari can crash without warning in some cases with ```torch==1.12.0```. This can be fixed by reverting to an earlier version using:
    
    pip install torch==1.11.0

### GPU

In order to use a GPU:

1. Uninstall the PyTorch version that gets installed by default with Cellpose:

        pip uninstall torch

2. Make sure your have up-to-date drivers for your NVIDIA card installed.

3. Re-install a GPU version of PyTorch via conda using a command that you can find [here](https://pytorch.org/get-started/locally/) (this takes care of the cuda toolkit, cudnn etc. so **no need to install manually anything more than the driver**). The command will look like this:

        conda install pytorch torchvision cudatoolkit=11.3 -c pytorch

### Plugin Updates

To update the plugin, you only need to activate the existing environment and install the new version:

    conda activate serialcellpose
    pip install napari-serialcellpose -U

## Usage: segmentation

The main interface is shown below. The sequence of events should be the following:

1. Select a folder containing images. The list of files within that folder will appear in the area above. You can also just drag and drop a folder or an image in that area. When selecting an image, it gets displayed in the viewer. Images are opened via [bioio](https://github.com/bioio-devs/bioio) but not all bioio plugins come pre-installed (only bioio-tifffile and bioio-imageio) so check the [plugin list](https://github.com/bioio-devs/bioio?tab=readme-ov-file#plug-in-registry) to open other file formats. You can use grayscale images, RGB images or multi-channel images. In the latter case, **make sure each channel opens as a separate layer when you open them**.
2. If you want to save the segmentation and tables with properties, select a folder that will contain the output.
3. Select the type of cellpose model.
4. If you use a custom model, select its location.
5. Run the analysis on the currently selected image or on all files in the folder.

### Options

6. Select if you want to use a GPU or not.
7. If you are using multi-channel images, you can specify which channel to segment and optionally which to use as "nuclei" channel to help cell segmentation.
8. In case you are using one of the built-in models, you can set the estimated diameter of your objects.
9. In the Options tab you will find a few more options for segmentation, including the two thresholds ```flow_threshold``` and ```cellprob_threshold```. You can also decide to discard objects touching the border. Using the ```Select options yml file``` you can select a ```.yml``` file which contains a list of additional options to pass to the ```eval``` method of the Cellpose model. **Note that options specified in the yml file will override options set in the GUI**. The file [my_options.yml](https://raw.githubusercontent.com/guiwitz/napari-serialcellpose/main/src/napari_serialcellpose/_tests/my_options.yml) is an example of such a file where for example the ```diameter``` (also available in the GUI) and ```resample``` (not available in the GUI) options are set. 

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui1.png" alt="image" width="500">
<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui1b.png" alt="image" width="500">

### Properties

10. After segmentation, properties of the objects can automatically be computed. You can select which properties should be computed in the Options tab. As defined in ```napari-skimage-regionprops``` properties are grouped by types. If you want to measure intensity properties such as mean intensity, you have to specify which channel (```Analysis channel```) you want to perform the measurement on.

### Output

The results of the analysis are saved in the folder chosen in #2. The segmentation mask is saved with the same name as the original image with the suffix ```_mask.tif```. A table with properties is saved in the subfolder ```tables``` also with the same name as the image with the suffix ```props.csv```. If you run the plugin on multiple files in a folder, a ```summary.csv``` file is also generated which compiles all the data.

## Usage: post-processing

After the analysis is done, when you select an image, the corresponding segmentation mask is shown on top of the image as shown below. This also works for saved segmentations: in that case you just select a folder with data and the corresponding output folder.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui2.png" alt="image" width="500">

### Properties

If you head to the **Properties** tab, you will find there two histograms showing the distribution of two properties that you can choose from a list at the top of the window. Below the plot you find the table containing information for each cell (each line is a cell).

As shown below, if you select the box ```show selected```, you can select items in the properties table and it will highlight the corresponding cell in the viewer. If you select the pipet tool, you can also select a cell and see the corresponding line in the table highlighted.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui3.png" alt="image" width="500">

### Summary

Finally if you select the **Summary** tab, and click on ```Load summary```, it will load all data of the current output folder and create histograms of two properties that can be selected. An additional property can be used for filtering the data. Using the sliders, one can set a minimum and maximum threshold on the "filtering property", which will create a sub-selection of the data.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui4.png" alt="image" width="500">

## Data

Sample data were acquired by Fabian Blank at the DBMR, University of Bern.

## License

Distributed under the terms of the [BSD-3] license,
"napari-serialcellpose" is free and open source software

## Issues

If you encounter any problems, please [file an issue] along with a detailed description.

[napari]: https://github.com/napari/napari
[Cookiecutter]: https://github.com/audreyr/cookiecutter
[@napari]: https://github.com/napari
[MIT]: http://opensource.org/licenses/MIT
[BSD-3]: http://opensource.org/licenses/BSD-3-Clause
[GNU GPL v3.0]: http://www.gnu.org/licenses/gpl-3.0.txt
[GNU LGPL v3.0]: http://www.gnu.org/licenses/lgpl-3.0.txt
[Apache Software License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[Mozilla Public License 2.0]: https://www.mozilla.org/media/MPL/2.0/index.txt
[cookiecutter-napari-plugin]: https://github.com/napari/cookiecutter-napari-plugin

[file an issue]: https://github.com/guiwitz/napari-serialcellpose/issues

[napari]: https://github.com/napari/napari
[tox]: https://tox.readthedocs.io/en/latest/
[pip]: https://pypi.org/project/pip/
[PyPI]: https://pypi.org/
