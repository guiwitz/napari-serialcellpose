from pathlib import Path
import skimage.io
import skimage.segmentation
from skimage.measure import regionprops_table
import pandas as pd

def run_cellpose(image_path, cellpose_model, output_path, scaling_factor=1,
                 diameter=None, flow_threshold=0.4, cellprob_threshold=0.0, clear_border=True):
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
    clear_border : bool
        remove cells touching border

    Returns
    -------
    cellpose_output : list of arrays
        list of segmented images
    """


    if not isinstance(image_path, list):
        image_path = [image_path]
    image = [skimage.io.imread(x) for x in image_path]
    
    for i in range(len(image)):
        if image[i].ndim == 3:
            image_gray = skimage.color.rgb2gray(image[i])
            image_gray = skimage.util.img_as_ubyte(image_gray)
            image[i] = image_gray
        if scaling_factor != 1:
            image[i] = image[i][::scaling_factor, ::scaling_factor]

    output_path = Path(output_path)
    
    # run cellpose
    cellpose_output = cellpose_model.eval(image, channels = [[0,0]], diameter=diameter, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold)
    cellpose_output = cellpose_output[0]

    if clear_border is True:

        cellpose_output = [skimage.segmentation.clear_border(im) for im in cellpose_output]
        cellpose_output = [skimage.segmentation.relabel_sequential(im)[0] for im in cellpose_output]
    
    # save output
    for im, p in zip(cellpose_output, image_path):
        save_path = output_path.joinpath(p.stem+'_mask.tif')
        skimage.io.imsave(save_path, im, check_contrast=False)

        compute_props(
            label_image=im,
            output_path=output_path,
            image_name=p
            )

    return cellpose_output


def compute_props(label_image, output_path, image_name):
    """Compute properties of segmented image.
    
    Parameters
    ----------
    label_image : array
        image with labeled cells
    output_path : str or Path
        path to output folder
    image_name : str or Path
        either path to image or image name

    """
    
    image_name = Path(image_name)
    output_path = Path(output_path).joinpath('tables')
    if not output_path.exists():
        output_path.mkdir(parents=True)

    props = regionprops_table(label_image,
                          properties=('label', 'area', 'eccentricity', 'solidity', 'feret_diameter_max'))
        
    props['eccentricity'] = props['eccentricity'].round(2)
    props['feret_diameter_max'] = props['feret_diameter_max'].round(2)
    props['solidity'] = props['solidity'].round(2)

    props = pd.DataFrame(props)
    props.to_csv(output_path.joinpath(image_name.stem+'_props.csv'), index=False)


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
    table_names = list(output_path.glob('*_props.csv'))
    
    all_props = []
    for p in table_names:
        props = pd.read_csv(p)
        props['name'] = p.stem
        all_props.append(props)
    all_props = pd.concat(all_props)

    return all_props

