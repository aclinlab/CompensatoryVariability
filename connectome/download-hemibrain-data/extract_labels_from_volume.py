from collections.abc import Mapping, Iterable, Iterator

import numpy as np
import pandas as pd


def extract_labels_from_volume(points_df, volume, box_zyx=None, vol_scale=0, label_names=None):
    """
    Given a list of point coordinates and a label volume, assign a
    label to each point based on its position in the volume.

    Note:
        Works in-place.  (Adds columns to the input DataFrame.)
    
    Extracting values from an array in numpy is simple.
    In the simplest case, this is equivalent to:
    
        coords = points_df[['z', 'y', 'x']].values.transpose()
        points_df['label'] = volume[(*coords,)]

    But this function supports extra features:
    
    - Points outside the volume extents are handled gracefully (they remain unlabeled).
    - The volume can be offset from the origin (doesn't start at (0,0,0)).
    - The volume can be provided in downscaled form, in which case the
      given points will be downscaled before sampling is performed.
    - Both label values (ints) and label names are output, if the label names were specified.
    
    Args:
        points_df:
            DataFrame with at least columns ['x', 'y', 'z'].
            The points in this DataFrame should be provided at SCALE-0,
            regardless of vol_scale.
            This function appends two additional columns to the DataFrame, IN-PLACE.
        
        volume:
            3D ndarray of label voxels
        
        box_zyx:
            The (min,max) coordinates in which the volume resides in the point coordinate space.
            It is assumed that this box is provided at the same scale as vol_scale,
            (i.e. it is not necessarily given using scale-0 coordiantes).
        
        vol_scale:
            Specifies the scale at which volume (and box_zyx) were provided.
            The coordinates in points_df will be downscaled accordingly.
            (For example, if vol_scale=3, then it is assumed that the provided volume 
            is downsampled by 8x in all dimensions.)

        label_names:
            Optional.  Specifies how label IDs map to label names.
            If provided, a new column 'label_name' will be appended to
            points_df in addition to the 'label' column.

            Must be either:
            - a mapping of `{ label_id: name }` (or `{ name : label_id }`),
              indicating each label ID in the output image, or
            - a list label names in which case the mapping is determined automatically
              by enumerating the labels in the given order (starting at 1).
    
    Returns:
        None.  Results are appended to the points_df as new column(s).
    """
    if box_zyx is None:
        box_zyx = np.array(([0]*volume.ndim, volume.shape))

    assert ((box_zyx[1] - box_zyx[0]) == volume.shape).all() 

    assert points_df.index.duplicated().sum() == 0, \
        "This function doesn't work if the input DataFrame's index has duplicate values."

    downsampled_coords_zyx = (points_df[['z', 'y', 'x']] // (2**vol_scale)).astype(np.int32)

    # Drop everything outside the combined_box
    min_z, min_y, min_x = box_zyx[0] #@UnusedVariable
    max_z, max_y, max_x = box_zyx[1] #@UnusedVariable
    dc = downsampled_coords_zyx
    downsampled_coords_zyx = dc.loc[   (dc['z'] >= min_z) & (dc['z'] < max_z)
                                     & (dc['y'] >= min_y) & (dc['y'] < max_y)
                                     & (dc['x'] >= min_x) & (dc['x'] < max_x) ]
    del dc

    downsampled_coords_zyx -= box_zyx[0]

    points_df.drop(columns=['label', 'label_name'], errors='ignore', inplace=True)
    points_df['label'] = volume.dtype.type(0)
    points_df.loc[downsampled_coords_zyx.index, 'label'] = volume[tuple(downsampled_coords_zyx.values.transpose())]

    if label_names is not None:
        if isinstance(label_names, Mapping):
            # We need a mapping of label_ids -> names.
            # If the user provided the reverse mapping,
            # then flip it.
            (k,v) = next(iter(label_names.items()))
            if isinstance(k, str):
                # Reverse the mapping
                label_names = { v:k for k,v in label_names.items() }
        else:
            label_names = dict(enumerate(label_names, start=1))
        
        name_set = ['<unspecified>', *label_names.values()]
        default_names = ['<unspecified>']*len(points_df)
        # FIXME: More than half of the runtime of this function is spent on this line!
        #        Is there some way to speed this up?
        points_df['label_name'] = pd.Categorical( default_names,
                                                  categories=name_set,
                                                  ordered=False )
        for label, name in label_names.items():
            rows = points_df['label'] == label
            points_df.loc[rows, 'label_name'] = name
