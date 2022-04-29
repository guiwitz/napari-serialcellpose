from pathlib import Path
import skimage.io
import skimage.segmentation
from skimage.measure import regionprops_table
import pandas as pd

def run_cellpose(image_path, cellpose_model, output_path, scaling_factor=1,
                 diameter=None, clear_border=True):
    """Run cellpose on image"""


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
    cellpose_output = cellpose_model.eval(image, channels = [[0,0]], diameter=diameter)
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
    """Compute cell properties"""
    
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
    """Load properties from file"""

    # get file name
    image_name = Path(image_name)
    output_path = Path(output_path).joinpath('tables')

    # load properties
    props_path = Path(output_path).joinpath(image_name.stem+'_props.csv')
    if props_path.exists():
        props = pd.read_csv(props_path)

    return props

def load_allprops(output_path):
    """Load properties from file"""

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

