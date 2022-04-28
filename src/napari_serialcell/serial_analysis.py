from pathlib import Path
import skimage.io

def run_cellpose(image_path, cellpose_model, output_path, scaling_factor=1,
                 diameter=None):
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
    cellpose_output = cellpose_model.eval(image, channels = [[0,0]])
    cellpose_output = cellpose_output[0]
    
    # save output
    for im, p in zip(cellpose_output, image_path):
        save_path = output_path.joinpath(p.stem+'_mask.tif')
        skimage.io.imsave(save_path, im, check_contrast=False)

    return cellpose_output


