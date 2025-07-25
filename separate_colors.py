#!/usr/bin/env python
#
# Moritz Blumer | 2025-07-24
#
# Plot different color ranges separately below the original image.


## FILE INFO
__author__ = 'Moritz Blumer, 2025'
__email__ = 'lmb215@cam.ac.uk'



## CONFIG BLOCK

# set channel ranges
HSV_RANGES = {
    'red': [
        ((0,   30, 56), (15,  255, 255)),
        ((145, 30, 56), (180, 255, 255))
    ],
    'yellow': [
        ((16,  30, 80), (40,  255, 255)),
    ],
    'blue': [
        ((60,  1, 1), (140, 255, 255)),
    ],
    'shadows': [
        ((0,  0, 0),   (180, 255, 5)),
    ],
}

# defaults for red/yellow/blue panels
BRIGHTNESS = 1.3            # increase brightness by x
CONTRAST = 1.1              # increase contrast by x
SATURATION = 1.1            # increase saturation by x
DARKS_THRESHOLD = 80        # select "dark regions" (value > x)
BRIGHTENING_FACTOR = 1.5    # brighten those "dark regions" by x

# defaults for "shades" panels
S_BRIGHTNESS = 0.9          # increase brightness by x
S_CONTRAST = 1.5            # increase contrast by x



## SETUP

#  import packages
import argparse
import cv2
import numpy as np



## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    global input_path, output_path

    parser = argparse.ArgumentParser(description="Plot different color" \
        " ranges separately below the original image.")

    # add arguments
    parser.add_argument('input_path', type=str, help='input PNG or JPG')
    parser.add_argument('output_path', type=str, help='input PNG or JPG')

    # parse
    args = parser.parse_args()

    # reassign variable names
    input_path, output_path = args.input_path, args.output_path



## FUNCTIONS


def read_image(image_path):
    '''
    Read image including alpha channel (or adding fully opaque alpha channel 
    if missing)
    '''

    image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    if image.shape[2] == 4:
        bgr = image[:, :, :3]
        alpha = image[:, :, 3]
    else:
        bgr = image
        alpha = np.full(bgr.shape[:2], 255, dtype=np.uint8)

    return(bgr, alpha)


def adjust_bgr(bgr, contrast=1.1, brightness=1.3):
    """
    Increase contrast and/or brightness in a BGR image
    """

    # convert to float
    bgr_float = bgr.astype(np.float32)

    # adjust contrast
    bgr_float = (bgr_float - 127.5) * contrast + 127.5

    # adjust brightness
    bgr_float = bgr_float * brightness

    # clip
    bgr_adj = np.clip(bgr_float, 0, 255).astype(np.uint8)

    return bgr_adj


def hsv_thresholding_mask(hsv, hsv_ranges):
    '''
    Create mask based on inpout HSV ranges
    '''

    # create empty mask
    mask = np.zeros(hsv.shape[:2], dtype=np.uint8)

    # apply thresholding
    for lower, upper in hsv_ranges:
        mask |= cv2.inRange(hsv, np.array(lower), np.array(upper))

    return mask


def modify_shades(hsv, darks_threshold, brightening_factor):
    '''
    Modify HSV image by adjusting shades (value < darks_threshold) by
    multiplying value channel by brightening_factor
    '''

    # create copy
    hsv_adj = hsv.copy()

    # extract value channel
    values = hsv_adj[:, :, 2]

    # extract dark areas where value < darks_threshold
    low_value_mask = values < darks_threshold

    # adjust value and clip
    values[low_value_mask] = np.clip(
        values[low_value_mask] * brightening_factor, 0, 255
    ).astype(np.uint8)

    # update values
    hsv_adj[:, :, 2] = values

    return hsv_adj


def modify_saturation(hsv, saturation_factor=1.1):
    """
    Adjust the saturation channel of an hsv
    """

    # create copy
    hsv_adj = hsv.copy()

    # extract saturation channel
    saturation = hsv_adj[:, :, 1].astype(np.float32)

    # adjust saturation
    saturation *= saturation_factor

    # clip
    saturation = np.clip(saturation, 0, 255)

    # update values
    hsv_adj[:, :, 1] = saturation.astype(np.uint8)

    return hsv_adj


def resize_img(img, width_new):
    '''
    Resize an image
    '''

    height, width = img.shape[:2]
    height_new = int(height * width_new / width)
    img_resized = cv2.resize(
        img,
        (width_new, height_new),
        interpolation=cv2.INTER_AREA
    )
    return img_resized


def compile_composite_image(img_a, img_b1, img_b2, img_b3, img_b4):
    '''
    Compile a compound image of the original in full size and 4 modified
    versions in a single row below
    '''

    # resize bottom row images
    bottom_imgs = [
        resize_img(img, int(img_a.shape[:2][1] / 4)) \
            for img in [img_b1, img_b2, img_b3, img_b4]
    ]

    # stack bottom row images horizontally
    bottom_row = bottom_imgs[0]
    for i in range(1, 4):
        bottom_row = np.hstack((bottom_row, bottom_imgs[i]))

    bottom_row = resize_img(
        bottom_row,
        img_a.shape[1],
    )

    composite_img = np.vstack((img_a, bottom_row))

    return composite_img



## MAIN

def main():

    # parse arguments
    cli()

    bgr_adj_dct = {}

    # read image
    bgr, alpha = read_image(input_path)

    # processred/yellow/blue ranges
    for i in ['red', 'yellow', 'blue']:
        bgr_adj = adjust_bgr(bgr, CONTRAST, BRIGHTNESS)
        hsv = cv2.cvtColor(bgr_adj, cv2.COLOR_BGR2HSV)
        mask = hsv_thresholding_mask(hsv, HSV_RANGES[i])
        hsv_adj = modify_shades(hsv, DARKS_THRESHOLD, BRIGHTENING_FACTOR)
        hsv_adj = modify_saturation(hsv_adj, SATURATION)
        bgr_adj_dct[i] = cv2.cvtColor(hsv_adj, cv2.COLOR_HSV2BGR)
        bgr_adj_dct[i][mask == 0] = [255, 255, 255]

    # process shadows range separately
    bgr_adj = adjust_bgr(bgr, S_CONTRAST, S_BRIGHTNESS)
    hsv = cv2.cvtColor(bgr_adj, cv2.COLOR_BGR2HSV)
    mask = hsv_thresholding_mask(hsv, HSV_RANGES['shadows'])
    bgr_adj_dct['shadows'] = cv2.cvtColor(hsv_adj, cv2.COLOR_HSV2BGR)
    bgr_adj_dct['shadows'][mask == 0] = [255, 255, 255]

    # create composite image
    composite_img = compile_composite_image(
        bgr,
        bgr_adj_dct['red'],
        bgr_adj_dct['yellow'],
        bgr_adj_dct['blue'],
        bgr_adj_dct['shadows'],
    )

    # write image
    cv2.imwrite(output_path, composite_img)


if __name__ == '__main__':
    main()

