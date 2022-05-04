# napari-serialcellpose

[![License](https://img.shields.io/pypi/l/napari-serialcellpose.svg?color=green)](https://github.com/guiwitz/napari-serialcellpose/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/napari-serialcellpose.svg?color=green)](https://pypi.org/project/napari-serialcellpose)
[![Python Version](https://img.shields.io/pypi/pyversions/napari-serialcellpose.svg?color=green)](https://python.org)
[![tests](https://github.com/guiwitz/napari-serialcellpose/workflows/tests/badge.svg)](https://github.com/guiwitz/napari-serialcellpose/actions)
[![codecov](https://codecov.io/gh/guiwitz/napari-serialcellpose/branch/main/graph/badge.svg)](https://codecov.io/gh/guiwitz/napari-serialcellpose)
[![napari hub](https://img.shields.io/endpoint?url=https://api.napari-hub.org/shields/napari-serialcellpose)](https://napari-hub.org/plugins/napari-serialcellpose)

This napari plugin allows to segment all images within a folder using cellpose. The cellpose model can be either a custom model or an official cellpose model. In addition, a set of geometrical properties can be visualized as histograms for each image or for the entire folder. The set of properties is currently fix, but we plan to make those selectable.

This plugin uses the great [napari-skimage-regionprops](https://github.com/haesleinhuepf/napari-skimage-regionprops) plugin to show properties as interactive tables.

## Installation

In order to let users choose whether they want to use the GPU or not for segmentation, cellpose is not added as a dependency of this package. Therefore, to use this plugin you need to create an environment and install a few packages manually. First create the environment:

    conda create -n serialcellpose python=3.8.5

Then activate it and install napari and the plugin:
    
    conda activate serialcellpose
    pip install "napari[all]"
    pip install git+https://github.com/guiwitz/napari-serialcellpose.git

### CPU

Then for CPU work, just install cellpose in the regular way:
    
    pip install cellpose

### GPU

For GPU work, you just need to make sure that you have an NVIDIA card and the appropriate driver for it. Then you can install PyTorch via conda using a command that you can find [here](https://pytorch.org/get-started/locally/) (this takes care of the cuda toolkit, cudnn etc. so **no need to install manually anything more than the driver**). The command will look like this:

    conda install pytorch torchvision cudatoolkit=11.3 -c pytorch

Finally you can install cellpose:
    
    pip install cellpose

**Note that it is important to install pytorch before cellpose, otherwise cellpose will install non-GPU dependencies.**

### Plugin Updates

To update the plugin, you only need to activate the existing environment and install the new version:

    conda activate serialcellpose
    pip install git+https://github.com/guiwitz/napari-serialcellpose.git -U

## Usage: segmentation

The main interface is shown below. The sequence of events should be the following:
1. Select a folder containing images. Either gray scale images or RGB that are then converted to gray scale. Once selected, the list of files will appear in the window above. When selecting an image, it gets displayed in the viewer.
2. Select a folder that will contain the output (segmentations and tables)
3. Select the type of cellpose model.
4. If you use a custom model, select it.
5. Set parameters. If you use an official cellpose model, you can e.g. set the diameter. You can also set how many images are sent together for analysis as a batch (for GPU, this will depend on the memory size. For CPU the effect is marginal).
6. Select if you want to use a GPU or not.
7. Run the analysis on the currently selected image or on all files in the folder.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui1.png" alt="image" width="500">

The results of the analysis are saved in the folder chosen in #2. The segmentation mask is saved with the same name as the original image with the suffix ```_mask.tif```. A table with properties is saved in the subfolder ```tables``` also with the same name as the image with the suffix ```props.csv```.
## Usage: post-processing

After the analysis is done, when you select an image, the corresponding segmentation mask is shown on top of the image as shown below. This also works for saved segmentations: in that case you just select a folder with data and the corresponding output folder.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui2.png" alt="image" width="500">

If you head to the **Properties** tab, you will find there two histograms showing the distribution of eccentricity and maximum Feret diameter. Below you find the table containing information for each cell (each line is a cell).

As shown below, if you select the box ```show selected```, you can select items in the properties table and it will highlight the corresponding cell in the viewer. If you select the pipet tool, you can also select a cell and see the corresponding line in the table highlighted.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui3.png" alt="image" width="500">

Finally if you select the **Summary** tab, and click on ```Load summary```, it will load all data of the given folder and create histograms of eccentricity and Feret diameter. The slider allows to put a threshold on eccentricity in order to select only round cells.

<img src="https://github.com/guiwitz/napari-serialcellpose/raw/main/illustrations/napari_serialcellpose_gui4.png" alt="image" width="500">

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
