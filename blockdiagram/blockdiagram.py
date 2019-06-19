import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mayavi import mlab
import time, sys
import scipy
from PIL import Image, ImageDraw

##################
# MAIN FUNCTIONS #
##################

def create_block_diagram(strat,dx,ve,xoffset,yoffset,scale,ci,strat_switch,contour_switch,bottom,topo_min,topo_max):
    """function for creating a 3D block diagram in Mayavi
    strat - input array with stratigraphic surfaces
    dx - size of gridcells in the horizontal direction in 'strat'
    ve - vertical exaggeration
    offset - offset in the y-direction relative to 0
    scale - scaling factor
    ci - contour interval
    strat_switch - 1 if you want to plot stratigraphy on the sides; 0 otherwise
    contour_switch - 1 if you want to plot contours on the top surface; 0 otherwise
    bottom - elevation value for the bottom of the block"""

    r,c,ts = np.shape(strat)
    z = scale*strat[:,:,ts-1]
    if strat_switch == 1:
        z1 = strat[:,:,0]
    else:
        z1 = strat[:,:,-1]

    X1 = scale*(xoffset + np.linspace(0,r-1,r)*dx)
    Y1 = scale*(yoffset + np.linspace(0,c-1,c)*dx)
    mlab.surf(X1,Y1,z,warp_scale=ve,colormap='gist_earth',vmin=scale*topo_min,vmax=scale*topo_max) #, line_width=5.0, representation='wireframe')
    if contour_switch == 1:
        contours = list(np.arange(vmin,vmax,ci*scale)) # list of contour values
        mlab.contour_surf(X1,Y1,z,contours=contours,warp_scale=ve,color=(0,0,0),line_width=1.0)

    gray = (0.6,0.6,0.6) # color for plotting sides
    
    # updip side:
    vertices, triangles = create_section(z1[:,0],dx,bottom) 
    x = scale*(xoffset + vertices[:,0])
    y = scale*(yoffset + np.zeros(np.shape(vertices[:,0])))
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=gray)
    
    # downdip side:
    vertices, triangles = create_section(z1[:,-1],dx,bottom) 
    x = scale*(xoffset + vertices[:,0])
    y = scale*(yoffset + (c-1)*dx*np.ones(np.shape(vertices[:,0])))
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=gray)

    # left edge (looking downdip):
    vertices, triangles = create_section(z1[0,:],dx,bottom) 
    x = scale*(xoffset + np.zeros(np.shape(vertices[:,0])))
    y = scale*(yoffset + vertices[:,0])
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=gray)
    
    # right edge (looking downdip):
    vertices, triangles = create_section(z1[-1,:],dx,bottom) 
    x = scale*(xoffset + (r-1)*dx*np.ones(np.shape(vertices[:,0])))
    y = scale*(yoffset + vertices[:,0])
    z = scale*ve*vertices[:,1]
    mlab.triangular_mesh(x,y,z,triangles,color=gray)
    
    # bottom face of block:
    vertices = dx*np.array([[0,0],[r-1,0],[r-1,c-1],[0,c-1]])
    triangles = [[0,1,3],[1,3,2]]
    x = scale*(xoffset + vertices[:,0])
    y = scale*(yoffset + vertices[:,1])
    z = scale*bottom*np.ones(np.shape(vertices[:,0]))
    mlab.triangular_mesh(x,y,ve*z,triangles,color=gray)

def add_stratigraphy_to_block_diagram(strat,facies,h,thalweg_z,dx,ve,xoffset,yoffset,scale,layers_switch,color_mode,colors,line_thickness,export):
    """function for adding stratigraphy to the sides of a block diagram
    colors layers by relative age
    strat - input array with stratigraphic surfaces
    facies - 1D array of facies codes for layers
    h - channel depth (height of point bar)
    thalweg_z - array of thalweg elevations for each layer
    dx - size of gridcells in the horizontal direction in 'strat'
    ve - vertical exaggeration
    offset - offset in the y-direction relative to 0
    scale - scaling factor
    layers_switch - if equals 1, stratigraphic boundaries will be plotted on the sides as black lines
    color_mode - determines what kind of plot is created; can be 'property', 'time', or 'facies'
    colors - colors scheme for facies (list of RGB values)
    line_thickness - tube radius for plotting layers on the sides
    export - if equals 1, the display can be saved as a VRML file for use in other programs (e.g., 3D printing)""" 
    r,c,ts=np.shape(strat)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=ts-1)
    cmap = matplotlib.cm.get_cmap('viridis')
    for layer_n in range(ts-1): # main loop
        update_progress(layer_n/(ts-1))
        vmin = scale*thalweg_z[layer_n] # minimum elevation (for colormap)
        vmax = vmin + scale*h # maximum elevation (for colormap)

        top = strat[:,0,layer_n+1]  # updip side
        base = strat[:,0,layer_n]
        if layers_switch == 1:
            X1 = scale*(xoffset + dx*np.arange(0,r))
            Y1 = scale*(yoffset + np.zeros(np.shape(base)))
            Z1 = ve*scale*base
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if np.max(top-base)>0:
            Points,Inds = triangulate_layers(top,base,dx)
            for i in range(len(Points)):
                vertices = Points[i]
                triangles, scalars = create_triangles(vertices)
                X1 = scale*(xoffset + vertices[:,0])
                Y1 = scale*(yoffset + dx*0*np.ones(np.shape(vertices[:,0])))
                Z1 = scale*vertices[:,1]
                plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scale*scalars,cmap,norm,vmin,vmax,export)

        top = strat[:,-1,layer_n+1]  # downdip side
        base = strat[:,-1,layer_n]
        if layers_switch == 1:
            X1 = scale*(xoffset + dx*np.arange(0,r))
            Y1 = scale*(yoffset + dx*(c-1)*np.ones(np.shape(base)))
            Z1 = ve*scale*base
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if np.max(top-base)>0:
            Points,Inds = triangulate_layers(top,base,dx)
            for i in range(len(Points)):
                vertices = Points[i]
                triangles, scalars = create_triangles(vertices)
                X1 = scale*(xoffset + vertices[:,0])
                Y1 = scale*(yoffset + dx*(c-1)*np.ones(np.shape(vertices[:,0])))
                Z1 = scale*vertices[:,1]
                plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scale*scalars,cmap,norm,vmin,vmax,export)

        top = strat[0,:,layer_n+1]  # left edge (looking downdip)
        base = strat[0,:,layer_n]
        if layers_switch == 1:
            X1 = scale*(xoffset + np.zeros(np.shape(base)))
            Y1 = scale*(yoffset + dx*np.arange(0,c))
            Z1 = ve*scale*base
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if np.max(top-base)>0:
            Points,Inds = triangulate_layers(top,base,dx)
            for i in range(len(Points)):
                vertices = Points[i]
                triangles, scalars = create_triangles(vertices)
                X1 = scale*(xoffset + dx*0*np.ones(np.shape(vertices[:,0])))
                Y1 = scale*(yoffset + vertices[:,0])
                Z1 = scale*vertices[:,1]
                plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scale*scalars,cmap,norm,vmin,vmax,export)

        top = strat[-1,:,layer_n+1] # right edge (looking downdip)
        base = strat[-1,:,layer_n]
        if layers_switch == 1:
            X1 = scale*(xoffset + dx*(r-1)*np.ones(np.shape(base)))
            Y1 = scale*(yoffset + dx*np.arange(0,c))
            Z1 = ve*scale*base
            mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
        if np.max(top-base)>0:
            Points,Inds = triangulate_layers(top,base,dx)
            for i in range(len(Points)):
                vertices = Points[i]
                triangles, scalars = create_triangles(vertices)
                X1 = scale*(xoffset + dx*(r-1)*np.ones(np.shape(vertices[:,0])))
                Y1 = scale*(yoffset + vertices[:,0])
                Z1 = scale*vertices[:,1]
                plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scale*scalars,cmap,norm,vmin,vmax,export)

def create_exploded_view(strat,facies,topo,h,nx,ny,gap,dx,ve,scale,strat_switch,layers_switch,contour_switch,color_mode,colors,line_thickness,bottom,export):
    """function for creating an exploded-view block diagram
    inputs:
    strat - stack of stratigraphic surfaces
    facies - 1D array of facies codes for layers
    topo - stack of topographic surfaces
    nx - number of blocks in x direction
    ny - number of blocks in y direction
    gap - gap between blocks (number of gridcells)
    dx - gridcell size
    ve - vertical exaggeration
    scale - scaling factor (for whole model)
    strat_switch - if equals 1, the stratigraphy will be plotted on the sides of the blocks
    layers_switch - if equals 1, the stratigraphic surfaces will be plotted on the sides (adds a lot of triangles - not good for 3D printing)
    contour_swicth - if equals 1, contours will be plotted on the top surface
    color_mode - determines what kind of plot is created; can be 'property', 'time', or 'facies'
    colors - colors scheme for facies (list of RGB values)
    line_thickness - - tube radius for plotting layers on the sides
    bottom - elevation value for the bottom of the block
    export - if equals 1, the display can be saved as a VRML file for use in other programs (e.g., 3D printing)"""
    r,c,ts=np.shape(strat)
    thalweg_z = []
    for layer_n in range(ts-1):
        t = layer_n - np.mod(layer_n,3)
        thalweg_z.append(np.min(topo[:,:,int(t+t/3)]))
    topo_min = np.min(strat[:,:,-1])
    topo_max = np.max(strat[:,:,-1])
    count = 0
    for i in range(nx):
        for j in range(ny):
            x1 = i*int(c/nx)
            x2 = (i+1)*int(c/nx)
            y1 = j*int(r/ny)
            y2 = (j+1)*int(r/ny)
            xoffset = (y1+j*gap)*dx
            yoffset = (x1+i*gap)*dx
            create_block_diagram(strat[y1:y2,x1:x2,:],dx,ve,xoffset,yoffset,scale,5.0,strat_switch,contour_switch,bottom,topo_min,topo_max)
            add_stratigraphy_to_block_diagram(strat[y1:y2,x1:x2,:],facies,h,thalweg_z,dx,ve,xoffset,yoffset,scale,layers_switch,color_mode,colors,line_thickness,export)
            count = count+1
            print("block "+str(count)+" done, out of "+str(nx*ny)+" blocks")

def create_fence_diagram(strat,facies,topo,h,nx,ny,gap,dx,ve,scale,layers_switch,color_mode,colors,line_thickness,bottom,export):
    """function for creating a fence diagram
    inputs:
    strat - stack of stratigraphic surfaces
    facies - 1D array of facies codes for layers
    topo - stack of topographic surfaces
    nx - number of strike sections
    ny - number of dip sections
    gap - gap between blocks (number of gridcells)
    dx - gridcell size
    ve - vertical exaggeration
    scale - scaling factor (for whole model)
    layers_switch - if equals 1, the stratigraphic surfaces will be plotted on the sides (adds a lot of triangles - not good for 3D printing)
    color_mode - determines what kind of plot is created; can be 'property', 'time', or 'facies'
    colors - colors scheme for facies (list of RGB values)
    line_thickness - - tube radius for plotting layers on the sides
    bottom - elevation value for the bottom of the block
    export - if equals 1, the display can be saved as a VRML file for use in other programs (e.g., 3D printing)"""
    r,c,ts=np.shape(strat)
    gray = (0.6,0.6,0.6)
    thalweg_z = []
    for layer_n in range(ts-1):
        t = layer_n - np.mod(layer_n,3)
        thalweg_z.append(np.min(topo[:,:,int(t+t/3)]))
    topo_min = np.min(strat[:,:,-1])
    topo_max = np.max(strat[:,:,-1])
    cmap = matplotlib.cm.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=ts-1)
    for nsec in range(1,nx+1): # strike sections
        x1 = nsec*int(c/(nx+1))
        vertices, triangles = create_section(strat[:,x1,0],dx,bottom) 
        x = scale*(vertices[:,0])
        y = scale*(x1*dx+np.zeros(np.shape(vertices[:,0])))
        z = scale*ve*vertices[:,1]
        mlab.triangular_mesh(x,y,z,triangles,color=gray)
        for layer_n in range(ts-1): # main loop
            update_progress(layer_n/(ts-1))
            vmin = scale*thalweg_z[layer_n] # minimum elevation (for colormap)
            vmax = vmin + scale*h # maximum elevation (for colormap)
            top = strat[:,x1,layer_n+1]  
            base = strat[:,x1,layer_n]
            if layers_switch == 1:
                X1 = scale*(dx*np.arange(0,r))
                Y1 = scale*(x1*dx+np.zeros(np.shape(base)))
                Z1 = ve*scale*base
                mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
            if np.max(top-base)>0:
                Points,Inds = triangulate_layers(top,base,dx)
                for i in range(len(Points)):
                    vertices = Points[i]
                    triangles, scalars = create_triangles(vertices)
                    X1 = scale*(vertices[:,0])
                    Y1 = scale*(x1*dx+dx*0*np.ones(np.shape(vertices[:,0])))
                    Z1 = scale*vertices[:,1]
                    plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scale*scalars,cmap,norm,vmin,vmax,export)
        print('done with section '+str(nsec)+' of '+str(nx)+' strike sections')
    for nsec in range(1,ny+1): # dip sections
        y1 = nsec*int(r/(ny+1))
        vertices, triangles = create_section(strat[y1,:,0],dx,bottom) 
        x = scale*(y1*dx+np.zeros(np.shape(vertices[:,0])))
        y = scale*(vertices[:,0])
        z = scale*ve*vertices[:,1]
        mlab.triangular_mesh(x,y,z,triangles,color=gray)
        for layer_n in range(ts-1): # main loop
            update_progress(layer_n/(ts-1))
            vmin = scale*thalweg_z[layer_n] # minimum elevation (for colormap)
            vmax = vmin + scale*h # maximum elevation (for colormap)
            top = strat[y1,:,layer_n+1]  
            base = strat[y1,:,layer_n]
            if layers_switch == 1:
                X1 = scale*(y1*dx+np.zeros(np.shape(base)))
                Y1 = scale*(dx*np.arange(0,c))
                Z1 = ve*scale*base
                mlab.plot3d(X1,Y1,Z1,color=(0,0,0),tube_radius=line_thickness)
            if np.max(top-base)>0:
                Points,Inds = triangulate_layers(top,base,dx)
                for i in range(len(Points)):
                    vertices = Points[i]
                    triangles, scalars = create_triangles(vertices)
                    X1 = scale*(y1*dx + dx*0*np.ones(np.shape(vertices[:,0])))
                    Y1 = scale*(vertices[:,0])
                    Z1 = scale*vertices[:,1]
                    plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scale*scalars,cmap,norm,vmin,vmax,export)
        print('done with section '+str(nsec)+' of '+str(ny)+' dip sections')

########################
# ADDITIONAL FUNCTIONS #
########################

def update_progress(progress):
    """progress bar from https://stackoverflow.com/questions/3160699/python-progress-bar
    update_progress() : Displays or updates a console progress bar
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%"""
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

def triangulate_layers(top,base,dx):
    """function for creating vertices of polygons that describe one layer"""
    x = dx * np.arange(0,len(top))
    ind1 = np.argwhere(top-base>0).flatten()
    ind2 = np.argwhere(np.diff(ind1)>1)
    ind2 = np.vstack((np.array([[-1]]),ind2))
    ind2 = np.vstack((ind2,np.array([[len(top)]])))
    Points = [] # list for points to be triangulated
    Inds = []
    for i in range(len(ind2)-1):
        ind3 = ind1[int(ind2[i])+1:int(ind2[i+1])+1]
        if (ind3[0] != 0) & (ind3[-1] != len(top)-1):
            ind3 = np.hstack((ind3[0]-1,ind3))
            ind3 = np.hstack((ind3,ind3[-1]+1)) 
            top1 = top[ind3][:-1]
            base1 = base[ind3][1:]
            x1 = np.concatenate((x[ind3][:-1], x[ind3][::-1][:-1]))
            inds = np.concatenate((ind3[:-1], ind3[::-1][:-1]))
        if (ind3[0] == 0) & (ind3[-1] != len(top)-1):
            ind3 = np.hstack((ind3,ind3[-1]+1))
            top1 = top[ind3][:-1]
            base1 = base[ind3]
            x1 = np.concatenate((x[ind3][:-1], x[ind3][::-1]))
            inds = np.concatenate((ind3[:-1], ind3[::-1]))
        if (ind3[0] != 0) & (ind3[-1] == len(top)-1):
            ind3 = np.hstack((ind3[0]-1,ind3))
            top1 = top[ind3]
            base1 = base[ind3][1:]
            x1 = np.concatenate((x[ind3], x[ind3][::-1][:-1]))
            inds = np.concatenate((ind3, ind3[::-1][:-1]))
        if (ind3[0] == 0) & (ind3[-1] == len(top)-1):
            top1 = top[ind3]
            base1 = base[ind3]
            x1 = np.concatenate((x[ind3], x[ind3][::-1]))
            inds = np.concatenate((ind3, ind3[::-1]))
        npoints = len(top1)+len(base1)
        y = np.hstack((top1,base1[::-1]))
        vertices = np.vstack((x1,y)).T
        Points.append(vertices)
        Inds.append(inds)
    return Points,Inds

def create_triangles(vertices):
    """function for creating list of triangles from vertices
    inputs:
    vertices - 2xn array with coordinates of polygon
    returns:
    triangles - indices of the 'vertices' array that from triangles (for triangular mesh)
    scalars - 'fake' elevation values for each vertex of the polygon, used for coloring (relies on the base of the polygon)"""
    n = len(vertices[:,0])
    Z1 = vertices[:,1]
    triangles = []
    if (np.mod(n,2)==0) & (vertices[int((n-1)/2),0] != vertices[int((n-1)/2+1),0]): # if polygon is in the interior of the block
        triangles.append([0,1,n-1])
        for i in range(1,int(n/2-1)):
            triangles.append([i,i+1,n-i])
            triangles.append([i+1,n-i,n-i-1])
        triangles.append([int(n/2-1),int(n/2),int(n/2+1)])
        scalars = np.hstack((Z1[0],Z1[int(n/2):][::-1],Z1[int(n/2)+1:]))
    if (np.mod(n,2)==0) & (vertices[int((n-1)/2),0] == vertices[int((n-1)/2+1),0]): # if polygon touches both sides of the block
        for i in range(0,int(n/2-1)):
            triangles.append([i,i+1,n-i-1])
            triangles.append([i+1,n-i-1,n-i-2])
        scalars = np.hstack((Z1[int(n/2):][::-1],Z1[int(n/2):]))
    if np.mod(n,2)!=0: # if polygon has one segment on the side of the block
        if vertices[int((n-1)/2),0] == vertices[int((n-1)/2+1),0]: # if polygon touches the right side of the block
            triangles.append([0,1,n-1])
            for i in range(1,int((n-1)/2)):
                triangles.append([i,i+1,n-i])
                triangles.append([i+1,n-i,n-i-1])
            scalars = np.hstack((Z1[0],Z1[int((n+1)/2):][::-1],Z1[int((n+1)/2):]))
        else:
            for i in range(0,int((n-1)/2)-1): # if polygon touches the left side of the block
                triangles.append([i,i+1,n-i-1])
                triangles.append([i+1,n-i-1,n-i-2])
            triangles.append([int((n-1)/2-1),int((n-1)/2),int((n-1)/2+1)])
            scalars = np.hstack((Z1[int((n+1)/2)-1:][::-1],Z1[int((n+1)/2):]))
    return triangles, scalars

def create_section(profile,dx,bottom):
    """function for creating a cross section from a top surface
    inputs:
    profile - elevation data for top surface
    dx - gridcell size
    bottom - elevation value for the bottom of the block
    returns:
    vertices - coordinates of vertices
    triangles - indices of the 'vertices' array that from triangles (for triangular mesh)
    """
    x1 = dx*np.linspace(0,len(profile)-1,len(profile))
    x = np.hstack((x1,x1[::-1]))
    y = np.hstack((profile,bottom*np.ones(np.shape(x1))))
    vertices = np.vstack((x,y)).T
    n = len(x)
    triangles = []
    for i in range(0,int((n-1)/2)):
        triangles.append([i,i+1,n-i-1])
        triangles.append([i+1,n-i-1,n-i-2])
    return vertices, triangles

def plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scalars,cmap,norm,vmin,vmax,export):
    """function for plotting layers on one side of a block
    inputs:
    layer_n - layer number
    facies - 1D array of facies codes for layers
    color_mode - determines what kind of plot is created; can be 'property', 'time', or 'facies'
    colors - list of RGB values used if color_mode is 'facies'
    X1,Y1,Z1 - coordinates of mesh vertices
    ve - vertical exaggeration
    triangles - indices of triangles used in mesh
    vertices - coordinates of the vertices
    scalars - scalars used for coloring the mesh in 'property' mode (= z-value of the base of current layer)
    cmap - colormap used for layers in 'time' mode
    norm - color normalization function used in 'time' mode
    export - if equals 1, the display can be saved as a VRML file for use in other programs (e.g., 3D printing)
    """
    if color_mode == 'time':
        mlab.triangular_mesh(X1,Y1,ve*Z1,triangles,color=cmap(norm(layer_n))[:3])
    if color_mode == 'property':
        if facies[layer_n] == 1:
            if export == 1:
                vmin = ve*vmin
                vmax = ve*vmax
                mesh = mlab.triangular_mesh(X1,Y1,ve*Z1,triangles,scalars=ve*scalars,colormap='YlOrBr',vmin=vmin,vmax=vmax)
                cmapf = matplotlib.cm.get_cmap('YlOrBr',256)
                normf = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
                z_range = np.linspace(np.min(ve*Z1),np.max(ve*Z1),256)
                mesh.module_manager.scalar_lut_manager.lut.table = (np.array(cmapf(normf(z_range)))*255).astype('uint8')
            else:
                mesh = mlab.triangular_mesh(X1,Y1,ve*Z1,triangles,scalars=scalars,colormap='YlOrBr',vmin=vmin,vmax=vmax)
        else:
            mlab.triangular_mesh(X1,Y1,ve*Z1,triangles,color=tuple(colors[int(facies[layer_n])]))
    if color_mode == 'facies':
        mlab.triangular_mesh(X1,Y1,ve*Z1,triangles,color=tuple(colors[int(facies[layer_n])]))

def create_random_section_2_points(strat,facies,thalweg_z,h,scale,ve,color_mode,colors,x1,x2,y1,y2,s1,dx,bottom,export):
    r, c, ts = np.shape(strat)
    dist = dx*((x2-x1)**2 + (y2-y1)**2)**0.5
    s2 = s1*dx+dist
    num = int(dist/float(dx))
    cmap = matplotlib.cm.get_cmap('viridis')
    Xrand, Yrand, Srand = np.linspace(x1,x2,num), np.linspace(y1,y2,num), np.linspace(s1*dx,s2,num)
    base = scipy.ndimage.map_coordinates(strat[:,:,0], np.vstack((Yrand,Xrand)))
    vertices, triangles = create_section(base,dx,bottom) 
    gray = (0.6,0.6,0.6) # color for plotting basal part of panel
    mlab.triangular_mesh(scale*np.hstack((dx*Xrand,dx*Xrand[::-1])),scale*np.hstack((dx*Yrand,dx*Yrand[::-1])),scale*ve*vertices[:,1],triangles,color=gray)
    for layer_n in range(0,ts-1):
        update_progress(layer_n/(ts-1))
        vmin = thalweg_z[layer_n] # minimum elevation (for colormap)
        vmax = vmin + h # maximum elevation (for colormap)
        top = scipy.ndimage.map_coordinates(strat[:,:,layer_n+1], np.vstack((Yrand,Xrand)))
        base = scipy.ndimage.map_coordinates(strat[:,:,layer_n], np.vstack((Yrand,Xrand)))
        if np.max(top-base)>1e-6:
            Points, Inds = triangulate_layers(top,base,dx)
            for i in range(len(Points)):
                vertices = Points[i]
                inds = Inds[i]
                triangles, scalars = create_triangles(vertices)
                X1 = scale*dx*Xrand[inds]
                Y1 = scale*dx*Yrand[inds]
                Z1 = scale*vertices[:,1]
                plot_layers_on_one_side(layer_n,facies,color_mode,colors,X1,Y1,Z1,ve,triangles,vertices,scalars,cmap,norm,vmin,vmax,export)
        
def create_random_section_n_points(strat,facies,topo,h,scale,ve,color_mode,colors,x1,x2,y1,y2,dx,bottom,export):
    r, c, ts = np.shape(strat)
    thalweg_z = []
    for layer_n in range(ts-1):
        t = layer_n - np.mod(layer_n,3)
        thalweg_z.append(np.min(topo[:,:,int(t+t/3)]))
    if len(x1)==1:
        create_random_section_2_points(strat,facies,thalweg_z,h,scale,ve,color_mode,colors,x1,x2,y1,y2,0,dx,bottom,export)
    else:
        count = 0
        dx1,dy1,ds1,s1 = compute_derivatives(x1,y1)
        for i in range(len(x1)):
            create_random_section_2_points(strat,facies,thalweg_z,h,scale,ve,color_mode,colors,x1[i],x2[i],y1[i],y2[i],s1[i],dx,bottom,export)
            count = count+1
            print("panel "+str(count)+" done, out of "+str(len(x1))+" panels")

def create_random_cookie(strat,facies,topo,h,scale,ve,color_mode,colors,x1,x2,y1,y2,dx,bottom,export):
    r, c, ts = np.shape(strat)
    thalweg_z = []
    for layer_n in range(ts-1):
        t = layer_n - np.mod(layer_n,3)
        thalweg_z.append(np.min(topo[:,:,int(t+t/3)]))
    count = 0
    dx1,dy1,ds1,s1 = compute_derivatives(x1,y1)
    for i in range(len(x1)):
        create_random_section_2_points(strat,facies,thalweg_z,h,scale,ve,color_mode,colors,x1[i],x2[i],y1[i],y2[i],s1[i],dx,bottom,export)
        count = count+1
        print("panel "+str(count)+" done, out of "+str(len(x1)+1)+" panels")
    create_random_section_2_points(strat,facies,thalweg_z,h,scale,ve,color_mode,colors,x2[-1],x1[0],y2[-1],y1[0],s1[-1]+np.sqrt((x1[0]-x2[-1])**2+(y1[0]-y2[-1])**2),dx,bottom,export)
    polygon = []
    for i in range(len(x1)):
        polygon.append((x1[i]+0.5,y1[i]+0.5))
    polygon.append((x2[-1]+0.5,y2[-1]+0.5))
    img = Image.fromarray(strat[:,:,-1])
    ImageDraw.Draw(img).polygon(polygon, outline=0, fill=1)
    mask = np.array(img)
    mask[mask!=1] = np.nan
    mask[mask==1] = strat[:,:,-1][mask==1]
    r,c = np.shape(strat[:,:,-1])
    Y1 = scale*(np.linspace(0,r-1,r)*dx)
    X1 = scale*(np.linspace(0,c-1,c)*dx)
    topo_min = np.min(strat[:,:,-1])
    topo_max = np.max(strat[:,:,-1])
    mlab.surf(X1,Y1,scale*mask.T,warp_scale=ve,colormap='gist_earth',vmin=scale*topo_min,vmax=scale*topo_max)
        
def compute_derivatives(x,y):
    dx = np.diff(x) # first derivatives
    dy = np.diff(y)   
    ds = np.sqrt(dx**2+dy**2)
    s = np.hstack((0,np.cumsum(ds)))
    return dx, dy, ds, s

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()

def select_random_section(strat):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.imshow(strat[:,:,-1],cmap='viridis')
    plt.tight_layout()
    ax.set_title('click to build line segments')
    line, = ax.plot([], [])  # empty line
    linebuilder = LineBuilder(line)
    xcoords = linebuilder.xs
    ycoords = linebuilder.ys
    return xcoords, ycoords