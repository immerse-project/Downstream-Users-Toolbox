import numpy as np
import datetime as dt

def define_centroid_cross_section(cross_section_pts, fr_obj, res=50, search_res=100):

    # First turn the line into a set of points
    line_pts = []
    for i in np.arange(0, len(cross_section_pts) -1):
        this_sect = [cross_section_pts[i], cross_section_pts[i+1]]
        this_dirn = np.asarray([this_sect[1][0] - this_sect[0][0], this_sect[1][1] - this_sect[0][1]])
        this_sect_len = np.sqrt(this_dirn[0]**2 + this_dirn[1]**2)
        this_dirn = this_dirn/this_sect_len    
        tot_pts = np.floor(this_sect_len/res) + 1
        this_line_sect = [this_sect[0] + j*res*this_dirn for j in np.arange(0,tot_pts)]
        line_pts.append(np.asarray(this_line_sect))
    line_pts = np.vstack(line_pts)

    # Exclude any outside the domain
    #in_domain = fr_obj.in_domain(line_pts[:,0], line_pts[:,1], cartesian=True)
    #end_pts = line_pts[[np.min(np.where(in_domain)) - 1, np.max(np.where(in_domain)) +1], :]
    #line_pts = line_pts[in_domain,:]
    
    # Reduce the number of centroids to search in
    mm_x = [np.min(line_pts[:,0]) - search_res, np.max(line_pts[:,0]) + search_res]
    mm_y = [np.min(line_pts[:,1]) - search_res, np.max(line_pts[:,1]) + search_res]

    node_nos = np.arange(0, len(fr_obj.variables['xc']))
    node_x = fr_obj.variables['xc']
    node_y = fr_obj.variables['yc']

    keep_pts = np.logical_and(np.logical_and(node_x >= mm_x[0], node_x <= mm_x[1]), np.logical_and(node_y >= mm_y[0], node_y <= mm_y[1]))
    
    node_x = node_x[keep_pts]
    node_y = node_y[keep_pts]
    node_nos = node_nos[keep_pts]

    # Find the nearest to each of the line subsample points
    close_node = []

    for this_pt in line_pts:
        dists = np.sqrt((this_pt[0] - node_x)**2 + (this_pt[1] - node_y)**2)
        close_node.append(node_nos[np.argmin(dists)])

    close_node_clean = [close_node[0]]

    for this_node in close_node[1:]:
        if this_node != close_node_clean[-1]:
            close_node_clean.append(this_node)

    if len(close_node_clean) == len(np.unique(close_node_clean)):
        return close_node_clean

def unstructured_grid_depths(h, zeta, sigma, nan_invalid=False):
    """
    Calculate the depth time series for cells in an unstructured grid.

    Parameters
    ----------
    h : np.ndarray
        Water depth
    zeta : np.ndarray
        Surface elevation time series
    sigma : np.ndarray
        Sigma vertical distribution, range 0-1 (`siglev' or `siglay' in FVCOM)
    nan_invalid : bool, optional
        Set values shallower than the mean sea level (`h') to NaN. Defaults to not doing that.

    Returns
    -------
    depths : np.ndarray
        Time series of model depths.

    """

    if nan_invalid:
        invalid = -zeta > h
        zeta[invalid] = np.nan

    abs_water_depth = zeta + h
    # Add zeta again so the range is surface elevation (`zeta') to mean water depth rather (`h') than zero to water
    # depth (`h' + `zeta') which is much more useful for plotting.
    depths = abs_water_depth[:, np.newaxis, :] * sigma[np.newaxis, ...] + zeta[:, np.newaxis, :]

    return depths

def cross_sect_trapeze(cross_centres, zeta_centre, grid_info):
    siglev_deps = unstructured_grid_depths(grid_info['h_center'][cross_centres], zeta_centre, grid_info['siglev_center'][:, cross_centres])
    siglay_deps = unstructured_grid_depths(grid_info['h_center'][cross_centres], zeta_centre, grid_info['siglay_center'][:, cross_centres])

    cross_x = grid_info['xc'][cross_centres]
    cross_y = grid_info['yc'][cross_centres]
    cross_z = siglev_deps

    mid_x = cross_x[0:-1] + (cross_x[1:] - cross_x[0:-1])/2
    mid_y = cross_y[0:-1] + (cross_y[1:] - cross_y[0:-1])/2

    end1_x = mid_x[0] + (cross_x[0] - mid_x[0])*2
    end1_y = mid_y[0] + (cross_y[0] - mid_y[0])*2

    end2_x = mid_x[-1] + (cross_x[-1] - mid_x[-1])*2
    end2_y = mid_y[-1] + (cross_y[-1] - mid_y[-1])*2

    mid_z = siglay_deps[:,:,0:-1] + (siglay_deps[:,:,1:] - siglay_deps[:,:,0:-1])/2

    cross_z_1 = cross_z[:,0:-1,0] + (cross_z[:,1:,0] - cross_z[:,0:-1,0])/2
    cross_z_2 = cross_z[:,0:-1,-1] + (cross_z[:,1:,-1] - cross_z[:,0:-1,-1])/2

    end1_z = mid_z[:,:,0] + (cross_z_1 - mid_z[:,:,0])*2 
    end2_z = mid_z[:,:,-1] + (cross_z_2 - mid_z[:,:,-1])*2

    mid_x = np.hstack([end1_x, mid_x, end2_x])
    mid_y = np.hstack([end1_y, mid_y, end2_y])
    mid_z = np.concatenate([end1_z[:,:,np.newaxis].T, mid_z.T, end2_z[:,:,np.newaxis].T]).T
    # Extend points outside grid

    end1_x = cross_x[0] - (cross_x[1] - cross_x[0])
    end1_y = cross_y[0] - (cross_y[1] - cross_y[0])

    end2_x = cross_x[-1] + (cross_x[-1] - cross_x[-2]) 
    end2_y = cross_y[-1] + (cross_y[-1] - cross_y[-2])

    cross_x = np.hstack([end1_x, cross_x, end2_x])
    cross_y = np.hstack([end1_y, cross_y, end2_y])

    end1_z = cross_z[:,:,0]
    end2_z = cross_z[:,:,-1]
    
    cross_z = np.concatenate([end1_z[:,:,np.newaxis].T, cross_z.T, end2_z[:,:,np.newaxis].T]).T

    end1_z = mid_z[:,0,:] + (mid_z[:,0,:] - mid_z[:,1,:]) 
    end2_z = mid_z[:,-1,:] + (mid_z[:,-2,:] - mid_z[:,-1,:])

    mid_z = np.concatenate([end1_z[:,np.newaxis,:], mid_z, end2_z[:,np.newaxis,:]], axis=1)

    # Get control trapezoid corners        

    t_step_tot = cross_z.shape[0]
    lr_ind = cross_z.shape[2] -1
    ud_ind = mid_z.shape[1] -1

    all_t = []
    for t_step in np.arange(0,t_step_tot):
        cols = []
        for i in np.arange(0, lr_ind):
            row = []
            for j in np.arange(0, ud_ind):
                line_a = [[mid_x[i], mid_y[i], mid_z[t_step,j,i]], [mid_x[i], mid_y[i], mid_z[t_step,j+1,i]]]
                line_b = [[cross_x[i], cross_y[i], cross_z[t_step, j,i]], [cross_x[i+1], cross_y[i+1], cross_z[t_step, j,i+1]]]
                pp = bisect_lines_pt(np.asarray(line_a), np.asarray(line_b))
                row.append(pp)
            cols.append(np.asarray(row))
        all_t.append(np.asarray(cols))

    # Collect trapezoids

    all_sqrs = []
    all_norms = []
    all_areas = []

    for cols in all_t:
        sqrs = []
        norms = []
        areas = []

        for i in np.arange(0, cols.shape[0] - 1):
            this_row = []
            this_row_norms = []
            this_row_areas = []
            for j in np.arange(0, cols.shape[1] -1):
                this_sqr = [cols[i,j,:], cols[i+1,j,:], cols[i+1,j+1,:], cols[i, j+1, :]]            
                this_row.append(this_sqr)
               
                this_norm = [cols[i,j,1] - cols[i+1,j,1], cols[i+1,j,0] - cols[i,j,0]] 
                this_row_norms.append(this_norm/np.linalg.norm(this_norm))
                
                this_row_areas.append(get_quad_area(this_sqr))

            sqrs.append(np.asarray(this_row))
            norms.append(np.asarray(this_row_norms))
            areas.append(np.asarray(this_row_areas))

        all_norms.append(np.asarray(norms))
        all_sqrs.append(np.asarray(sqrs))
        all_areas.append(np.asarray(areas))

    return np.asarray(all_sqrs), np.asarray(all_norms), np.asarray(all_areas), [mid_x, mid_y, mid_z]

def get_area_heron(s1, s2, s3):
    """
    Calculate the area of a triangle/set of triangles based on side length (Herons formula). Could tidy by combining
    with get_area.

    Parameters
    ----------
    s1, s2, s3 : tuple, list (float, float)
        Side lengths of the three sides of a triangle. Can be 1D arrays of lengths or lists of lengths.

    Returns
    -------
    area : tuple, np.ndarray
        Area of the triangle(s). Units of v0, v1 and v2.

    """

    s1 = np.asarray(s1)
    s2 = np.asarray(s2)
    s3 = np.asarray(s3)

    p = 0.5 * (s1 + s2 + s3)

    area = np.sqrt(p * (p - s1) * (p - s2) * (p - s3))

    return abs(area)

def get_quad_area(points):
    s1_1 = line_3d_len(points[0], points[2])
    s2_1 = line_3d_len(points[0], points[1])
    s3_1 = line_3d_len(points[1], points[2])    

    s1_2 = s1_1
    s2_2 = line_3d_len(points[2], points[3])
    s3_2 = line_3d_len(points[3], points[0])    
    

    a1 = get_area_heron(s1_1, s2_1, s3_1)
    a2 = get_area_heron(s1_2, s2_2, s3_2) 

    return a1+a2

def dist_along_cross(mid_x, mid_y, origin):
    dists = [np.sqrt((origin[0] - mid_x[0])**2 + (origin[1] - mid_y[0])**2)]
    for i in np.arange(1,len(mid_x)):
        this_dist = np.sqrt((mid_x[i] - mid_x[i-1])**2 + (mid_y[i] - mid_y[i-1])**2)
        dists.append(dists[-1] + this_dist)
    return np.asarray(dists)

def line_3d_len(p1, p2):
    return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2)

def bisect_lines_pt(lineA, lineB):
    A = lineA[0]
    B = lineB[0]
    dA = (lineA[1] - lineA[0])/np.linalg.norm(lineA[1] - lineA[0])
    dB = (lineB[1] - lineB[0])/np.linalg.norm(lineB[1] - lineB[0])

    u = A - B
    b = np.dot(dA, dB)
    d = np.dot(dA, u)
    e = np.dot(dB, u)
    t_intersect = (b * e - d) / (1 - b * b)
    P = A + t_intersect * dA
    return P

from matplotlib.tri import LinearTriInterpolator
from matplotlib.tri.triangulation import Triangulation

def node_to_centre(field, filereader):
    """
    TODO: docstring

    """
    tt = Triangulation(filereader.variables['x'][:], filereader.variables['y'][:], np.asarray(filereader.variables['nv'][:]-1, dtype=int).T)
    interped_out = []
    if len(field.shape) == 1:
        field = field[np.newaxis, :]

    for this_t in field:
        ct = LinearTriInterpolator(tt, this_t)
        interped_out.append(ct(filereader.variables['xc'][:], filereader.variables['yc'][:]))

    return np.asarray(interped_out)

def get_plot_sqr(sqrs):
    all_out = []
    for this_row in np.arange(0, sqrs.shape[0]):
        row_out = []
        for this_col in np.arange(0, sqrs.shape[1]):
            this_sqr = np.asarray(sqrs[this_row, this_col, :,:])
            plot_x = np.hstack([this_sqr[:,0], this_sqr[0,0]])
            plot_y = np.hstack([this_sqr[:,2], this_sqr[0,2]])
            row_out.append([plot_x, plot_y])

        all_out.append(row_out)

    return all_out

