# Target Volume Predictive Model for Adaptive Radiation Therapy

Given an input of a time series of dimensional points P, the PredictiveModel function performs linear quadratic estimation (Kalman filtering) on each point, one dimension at a time, through all samples. Assumes input arguments maintains row correspondence for each point.

The output produced is:

- Python: Set of points in a ```.roi``` file defining a Kalman filtered volume
- MATLAB implementation: raw_data cell array defining a Kalman filtered volume

Mathematically:

![forumulas](https://raw.githubusercontent.com/yazanobeidi/art-predictive-model/master/images/formulas.jpg)

## Usage

**MATLAB**:

Runtime example can be executed with ```runPredictiveModel``` script from command.

**PYTHON**:

- Define your series of volumes similar to ```plan.roi```.
- Place your ```plan.roi``` and ```predModel.py``` in the same directory. 
- Run ```predModel.py```. 
- See ```output.txt``` or ```newartplan.roi``` for output.

## Clinical application

This algorithm was designed as part of adaptive radiation therapy for clinical
target volume determination. Here are some images of it in use:

![vol](https://raw.githubusercontent.com/yazanobeidi/art-predictive-model/master/images/predictive_model_1.jpg)


![vol2](https://raw.githubusercontent.com/yazanobeidi/art-predictive-model/master/images/predictive_model_2.jpg)


![vol3](https://raw.githubusercontent.com/yazanobeidi/art-predictive-model/master/images/predictive_model_3.jpg)


## Disclaimer

The Python code is quite un-Pythonic. I may provide an updated version some day however it does what it needs to do and I rather spend my time doing new things so probably not. So take this project for the idea rather than its current implementation. If you feel like refactoring some code feel free to create a pull request!

## License

Copyright 2015, 2016 Yazan Obeidi.

You may use any of this for educational purposes.
