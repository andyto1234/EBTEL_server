import numpy as np

def bresenham(x1, y1, x2, y2):
    """
    Returns an array of all pixel coordinates which the line defined by `x1, y1` and
    `x2, y2` crosses. Uses Bresenham's line algorithm to enumerate the pixels along
    a line. This was adapted from ginga

    Parameters
    ----------
    x1, y1, x2, y2 :`int`

    References
    ----------
    | https://github.com/ejeschke/ginga/blob/master/ginga/BaseImage.py#L387
    | http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    | https://ginga.readthedocs.org/en/latest/
    """
    for x in [x1, y1, x2, y2]:
        if type(x) not in (int, np.int64):
            raise TypeError('All pixel coordinates must be of type int')
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx - dy
    res = []
    x, y = x1, y1
    while True:
        res.append((x, y))
        if (x == x2) and (y == y2):
            break
        e2 = 2 * err
        if e2 > -dy:
            err = err - dy
            x += sx
        if e2 <  dx:
            err = err + dx
            y += sy
    return np.array(res)


def get_intersecting_pixels(coord, image_wcs):
    # Find pixels between each loop segment
    px, py = image_wcs.world_to_pixel(coord)
    px = np.round(px).astype(int)
    py = np.round(py).astype(int)
    loop_pix = []
    for i in range(px.shape[0]-1):
        b = bresenham(px[i], py[i], px[i+1], py[i+1])
        # Pop the last one, unless this is the final entry because the first point
        # of the next section will be the same
        if i < px.shape[0]-2:
            b = b[:-1]
        loop_pix.append(b)
    return np.vstack(loop_pix)


def filter_pix(loop_coord, ref_map):
    pix = get_intersecting_pixels(loop_coord, ref_map.wcs)
    pix_x = pix[:,0][(pix[:,0]>=0) & (pix[:,0]<len(ref_map.data[0,:])) & (pix[:,1]>=0) & (pix[:,1]<len(ref_map.data[:,0]))]
    pix_y = pix[:,1][(pix[:,0]>=0) & (pix[:,0]<len(ref_map.data[0,:])) & (pix[:,1]>=0) & (pix[:,1]<len(ref_map.data[:,0]))]

    return pix_x, pix_y

