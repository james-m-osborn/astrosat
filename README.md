# astrosat

Astrosat is a Python package which calculates which satellites can be seen by a given observer in a given field of view at a given observation time and observation duration.
This includes the geometry of the satellite and observer but also estimates the expected apparent brightness of the satellite to aid astronomers in assessing the impact on their observations.

## Installation
### PyPI
Astrosat is available on the Python Package Index - [PyPI](https://pypi.org/project/astrosat/) and can be installed using:

`pip install astrosat`

### Installation From Source
Alternatively astrosat can be installed from source by cloning this project and running:

`python setup.py install`

## Usage
Astrosat reads in an observation from a configuration file, then predicts any intersections between the observation and
satellites using online tools.
The easiest way to get started is to use the example files provided with this project.
The example configuration file is `params.yaml` and an example script for generating a plot of the FOV of an observation
with the satellite trails overlaid is given in `tel_fov.py`.
