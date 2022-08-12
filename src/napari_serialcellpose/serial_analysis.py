from pathlib import Path
import warnings
import skimage.io
import skimage.segmentation
from skimage.measure import regionprops_table as sk_regionprops_table
from napari_skimage_regionprops._regionprops import regionprops_table
import pandas as pd
import numpy as np
from aicsimageio import AICSImage
import yaml

def run_cellpose(image_path, cellpose_model, output_path, scaling_factor=1,
                 diameter=None, flow_threshold=0.4, cellprob_threshold=0.0,
                 clear_border=True, channel_to_segment=0, channel_helper=0,
                 channel_measure=None, channel_measure_names=None, properties=None,
                 options_file=None, force_no_rgb=False):
    """Run cellpose on image.
    
    Parameters
    ----------
    image_path : str or Path
        path to image
    cellpose_model : cellpose model instance
    output_path : str or Path
        path to output folder
    scaling_factor : int
        scaling factor for image (not implemented)
    diameter : int
        diameter of cells to segment, only useful for native cellpose models
    flow_threshold : float
        cellpose setting: maximum allowed error of the flows for each mask
    cellprob_threshold : float
        cellpose setting: pixels greater than the cellprob_threshold are used to run dynamics and determine ROIs
    clear_border : bool
        remove cells touching border
    channel_to_segment : int, default 0
        index of channel to segment, if image is multi-channel
    channel_helper : int, default 0
        index of helper nucleus channel for models using both cell and nucleus channels
    channel_measure: int or list of int, default None
        index of channel(s) in which to measure intensity
    channel_measure_names: list of str, default None
        names of channel(s) in which to measure intensity
    properties = list of str, default None
        list of types of properties to compute. Any of 'intensity', 'perimeter', 'shape', 'position', 'moments'
    options_file: str or Path, default None
        path to yaml options file for cellpose
    force_no_rgb: bool, default False
        if image is RGB convert it to multi-channel

    Returns
    -------
    cellpose_output : list of arrays
        list of segmented images
    props : pandas dataframe
        properties of segmented cells for the last analyzed image
    """

    if not isinstance(image_path, list):
        image_path = [image_path]

    if properties is None:
        properties = []

    channels = [0, 0]
    image_aics = [AICSImage(x) for x in image_path]
    if (len(image_aics[0].dims.shape) == 6) and (not force_no_rgb):
        image = [x.get_image_data('YXS', C=0, T=0, Z=0) for x in image_aics]
        is_rgb = True
        image_measure = [None]*len(image)
    else:
        #!!! Note that there is a bug in aicsimageio that causes it to place the S dimension in the wrong place.
        #!!! This is why currently we import S planes separately and stack them.
        if force_no_rgb:
            if channel_helper == 0:
                image = [x.get_image_data('YX', S=np.max([0,channel_to_segment-1]) , T=0, Z=0, C=0) for x in image_aics]
            else:
                im1 = [x.get_image_data('YX', S=np.max([0,channel_to_segment-1]) , T=0, Z=0, C=0) for x in image_aics]
                im2 = [x.get_image_data('YX', S=np.max([0,channel_helper-1]) , T=0, Z=0, C=0) for x in image_aics]
                image = [np.stack([im1[i], im2[i]], axis=0) for i in range(len(im1))]
                #image = [x.get_image_data('SYX', S=[np.max([0,channel_to_segment-1]), np.max([0,channel_helper-1])] , T=0, Z=0, C=0) for x in image_aics]
                channels = [1, 2]
        else:
            if channel_helper == 0:
                image = [x.get_image_data('YX', C=np.max([0,channel_to_segment-1]) , T=0, Z=0) for x in image_aics]
            else:
                image = [x.get_image_data('CYX', C=[np.max([0,channel_to_segment-1]), np.max([0,channel_helper-1])] , T=0, Z=0) for x in image_aics]
                channels = [1, 2]

        image_measure=None
        if channel_measure is not None:
            if force_no_rgb:
                image_measure = [np.stack([x.get_image_data('YX', S=s, T=0, Z=0, C=0) for s in channel_measure], axis=2) for x in image_aics]
                #image_measure = [x.get_image_data('YXS', S=channel_measure, T=0, Z=0, C=0) for x in image_aics]
            else:
                image_measure = [x.get_image_data('YXC', C=channel_measure, T=0, Z=0) for x in image_aics]
        else:
            image_measure = [None]*len(image)
        is_rgb = False

    for i in range(len(image)):
        if image[i].ndim == 3:
            if is_rgb:
                image_gray = skimage.color.rgb2gray(image[i])
                image_gray = skimage.util.img_as_ubyte(image_gray)
                image[i] = image_gray
        if scaling_factor != 1:
            image[i] = image[i][::scaling_factor, ::scaling_factor]
    
    # handle yaml options file
    default_options = {'diameter': diameter, 'flow_threshold': flow_threshold, 'cellprob_threshold': cellprob_threshold}
    options_yml = {}
    if options_file is not None:
        with open(options_file) as file:
            options_yml = yaml.load(file, Loader=yaml.FullLoader)
        list_of_cellpose_options = cellpose_model.eval.__code__.co_varnames
        for k in options_yml.keys():
            if k not in list_of_cellpose_options:
                raise ValueError(f'options file contains key {k} which is not in cellpose model')
    merged_options = {**default_options, **options_yml}

        
    cellpose_output = cellpose_model.eval(
        image, channels=channels, channel_axis=0,
        **merged_options
    )
    cellpose_output = cellpose_output[0]

    if clear_border is True:

        cellpose_output = [skimage.segmentation.clear_border(im) for im in cellpose_output]
        cellpose_output = [skimage.segmentation.relabel_sequential(im)[0] for im in cellpose_output]
    
    # save output
    for im, im_m, p in zip(cellpose_output, image_measure, image_path):
        
        props=None
        if len(properties) > 0:
            props = compute_props(
                    label_image=im,
                    intensity_image=im_m,
                    output_path=output_path,
                    image_name=p,
                    properties=properties,
                    channel_names=channel_measure_names
                    )

        if output_path is not None:
            output_path = Path(output_path)
            save_path = output_path.joinpath(p.stem+'_mask.tif')
            skimage.io.imsave(save_path, im, check_contrast=False)

    return cellpose_output, props


def compute_props(
    label_image, intensity_image, output_path=None,
    image_name=None, properties=None, channel_names=None):
    """Compute properties of segmented image.
    
    Parameters
    ----------
    label_image : array
        image with labeled cells
    intensity_image : array
        image with intensity values
    output_path : str or Path
        path to output folder
    image_name : str or Path
        either path to image or image name
    properties = list of str, default None
        list of types of properties to compute. Any of 'intensity', 'perimeter', 'shape', 'position', 'moments'
    channel_names: list of str, default None
        names of channel(s) in which to measure intensity
    """
    
    if (image_name is not None) and (output_path is not None):
        image_name = Path(image_name)
        output_path = Path(output_path).joinpath('tables')
        if not output_path.exists():
            output_path.mkdir(parents=True)
    
    if properties is None:
        properties = []
        
    if intensity_image is None:
        if "intensity" in properties:
            warnings.warn("Computing intensity features but no intensity image provided. Result will be zero.")
        intensity_image = np.zeros(label_image.shape)[:,:,np.newaxis]
        
    props = regionprops_table(
        image=intensity_image[:,:,-1], labels=label_image,
        size='size' in properties,
        perimeter='perimeter' in properties,
        shape='shape' in properties,
        position='position' in properties,
        moments='moments' in properties,
        intensity=False,
        )

    if 'intensity' in properties:
        intensity_measure = sk_regionprops_table(
            label_image=label_image, intensity_image=intensity_image,
            properties=['max_intensity', 'mean_intensity', 'min_intensity'])
        intensity_measure = pd.DataFrame(intensity_measure)
        if channel_names is not None:
            for ind, c in enumerate(channel_names):
                intensity_measure.rename(
                    columns={
                        f'mean_intensity-{ind}': f'mean_intensity-{c}',
                        f'min_intensity-{ind}': f'min_intensity-{c}',
                        f'max_intensity-{ind}': f'max_intensity-{c}'}, inplace=True)
        props = pd.concat([props, intensity_measure], axis=1)

    if output_path is not None:
        props.to_csv(output_path.joinpath(image_name.stem+'_props.csv'), index=False)

    return props


def load_props(output_path, image_name):
    """Load properties for an analyzed image.
    
    Parameters
    ----------
    output_path : str or Path
        path to output folder
    image_name : str or Path
        either path to image or image name
    
    Returns
    -------
    props : pandas dataframe
        properties of segmented cells
    """

    # get file name
    image_name = Path(image_name)
    output_path = Path(output_path).joinpath('tables')

    # load properties
    props_path = Path(output_path).joinpath(image_name.stem+'_props.csv')
    props=None
    if props_path.exists():
        props = pd.read_csv(props_path)

    return props

def load_allprops(output_path):
    """Load all properties files for a given folder.
    
    Parameters
    ----------
    output_path : str or Path
        path to output folder

    Returns
    -------
    all_props : pandas dataframe
        properties of segmented cells in all images
    
    """

    # get file name
    output_path = Path(output_path).joinpath('tables')
    if not output_path.exists():
        return None
    table_names = list(output_path.glob('*_props.csv'))
    
    all_props = []
    for p in table_names:
        props = pd.read_csv(p)
        props['name'] = p.stem
        all_props.append(props)
    all_props = pd.concat(all_props)

    all_props.to_csv(output_path.joinpath('summary.csv'), index=False)

    return all_props

