# InsoleTilter
Created insole and mold models for specified tilt angle.

## Requirements
These are the packages that this application has been developed with.
No testing with other versions has been done, but it may work.
- Python - 3.7.1
- PyQt - 5.12.3
- PyQtGraph - 0.10.0
- Trimesh 2.38.17
- Numpy - 1.16.3
- Splipy - 1.3.1
- PyOpenGL - 3.1.1a1

## Usage
To launch the application type `python src/main.py` in the directory you cloned the repo. It will load a basic sole model. Use the mouse to navigate around the model.

- You can load a new mesh by using the `Load mesh` button. Make sure your mesh is a surface and not a volume.
- Use the buttons to move the spline's control points. Click `Recompute Spline` to recompute and display the new spline. Click `Cut Top` to cut the top part along the spline. Then click `Cut Front` to cut the front part between the two edges of the spline.
- Use the spinner to choose the portion of the mesh where the rotation linearly decreases.
- Use the radio buttons to choose between a normal (Full) rotation and a stretched rotation.
- Use the slider to determine the rotation angle.
- When you're satisfied with your parameters, click the `Compute Insole and Mold models` checkbox to compute the models. Changing parameters with this option checked will drastically diminish the performances of the software, so it is recommended to uncheck it before making modifications.
- Use the export buttons to save the stl files you need.
- Use the dropdown menu to choose the shader you want for your display.
