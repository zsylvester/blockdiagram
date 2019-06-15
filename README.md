
<img src="https://github.com/zsylvester/blockdiagram/blob/master/blockdiagram_logo.png" width="400">

## Description

'blockdiagram' is a Python module for creating block diagrams and other three-dimensional displays from stratigraphic models. It is designed to work with [meanderpy](https://github.com/zsylvester/meanderpy), but it should work with any model that consists of a stack of stratigraphic surfaces.

## Requirements

numpy  
matplotlib  
mayavi  
scipy  
PIL

## Usage

The main function in 'blockdiagram' is 'create_exploded_view'. It can either be used to generate a normal block diagram or an exploded-view block diagram, in which the model is split into several smaller blocks so that more stratigraphic detail is visible. Here is a typical set of input parameters:

```python
mlab.figure(bgcolor=(1,1,1)) 
# parameters
ve = 15.0 # vertical exaggeration
scale = 0.1 # scaling of diagram (important for 3D printing)
strat_switch = 1 # equals 1 if you want stratigraphy displayed on the sides
layers_switch = 0 # equals 1 if you want stratigraphic boundaries displayed on the sides
contour_switch = 0 # equals 1 if you want contours displayed on the top surface
dx = 10.0 # cell size for display
bottom = np.min(chb_3d.strat) - 1.5 # elevation of bottom side of diagram
color_mode = 'property' # determines how the stratigraphy will be colored; can be 'property', 'facies', or 'time'
colors = [[0.5,0.25,0],[0.9,0.9,0],[0.5,0.25,0]] # colors for 'facies' display
line_thickness = 1.0 # thickness of lines if 'layers_switch' is 1
gap = 20 # distance between exploded blocks (if any; in number of gridcells)
h = 5.0 # channel depth (m)
nx = 1 # number of blocks in x direction
ny = 1 # number of blocks in y direction

bd.create_exploded_view(chb_3d.strat,chb_3d.facies,chb_3d.topo,h,nx,ny,gap,dx,ve,scale,strat_switch,
                        layers_switch,contour_switch,color_mode,colors,line_thickness,bottom)
```
If the command above is run with nx=1 and ny=1 (the number of blocks in the x and y directions), a simple block diagram is displayed:
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_1.png" width="800">

Changing nx and ny to 3 results in something like this:
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_2.png" width="800">

Both of the models above have been colored using the 'property' setting for the 'color_mode' parameter, so that the change from yellow to brown in the point bars reflects the change in grain size (and porosity/permeability). This setting can also be set to 'facies' (when each facies, e.g., point bar vs. overbank, gets its own color) or to 'time', when layers are colored according to their relative age - see example below.
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_5.png" width="800">

Another functionality is to create a 'random' section from the model. In order to do that, the location of the section has to be selected on a map of the top surface, using the 'select_random_section' function:

```python
xcoords, ycoords = bd.select_random_section(chb_3d.strat) # define x and y coordinates for random section
mlab.figure(bgcolor=(1,1,1))
color_mode = 'property'
bd.create_random_section_n_points(chb_3d.strat,chb_3d.facies,chb_3d.topo,h,scale,ve,color_mode,colors,
                                   xcoords[:-1],xcoords[1:],ycoords[:-1],ycoords[1:],dx,bottom)
```
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_6.png" width="300">
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_7.png" width="800">

Finally, you can also cut a "cookie" from the model, using the 'create_random_cookie' option:

```python
xcoords, ycoords = bd.select_random_section(chb_3d.strat) # define x and y coordinates for random section
mlab.figure(bgcolor=(1,1,1))
bd.create_random_cookie(chb_3d.strat,chb_3d.facies,chb_3d.topo,h,scale,ve,color_mode,colors,xcoords[:-1],xcoords[1:],
                        ycoords[:-1],ycoords[1:],dx,bottom)
```
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_4.png" width="300">
<img src="https://github.com/zsylvester/blockdiagram/blob/master/fluvial_model_example_5.png" width="800">

## License

'blockdiagram' is licensed under the Apache License 2.0

Copyright 2019 Zoltan Sylvester
