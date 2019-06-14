
<img src="https://github.com/zsylvester/blockdiagram/blob/master/blockdiagram_logo.png" width="400">

## Description

'blockdiagram' is a Python module for creating block diagrams and other three-dimensional displays from stratigraphic models. It is designed to work well with [meanderpy](https://github.com/zsylvester/meanderpy), but it should work with any model that consists of a stack of stratigraphic surfaces.

## Usage

```python
ve = 15.0
scale = 0.1
strat_switch = 1
layers_switch = 0
contour_switch = 0
dx = 10.0
bottom = np.min(chb_3d.strat) - 10
colors = [[0.5,0.25,0],[0.9,0.9,0],[0.5,0.25,0]]
line_thickness = 1.0
gap = 20
color_mode = 'property'
h = 5.0

bd.create_exploded_view(chb_3d.strat,chb_3d.facies,chb_3d.topo,h,3,3,gap,dx,ve,scale,strat_switch,layers_switch,contour_switch,color_mode,colors,line_thickness,bottom)
```

## License

'blockdiagram' is licensed under the Apache License 2.0
