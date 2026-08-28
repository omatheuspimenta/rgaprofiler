import os
import re
from setuptools import setup, find_packages

PATH_ROOT = os.path.dirname(__file__)


with open(os.path.join(PATH_ROOT, "README.md")) as f:
    readme = f.read()

with open(os.path.join(PATH_ROOT, "requirements.txt")) as f:
    requirements = f.read()




setup(
    name="signalp6",
    version='6.0+h',
    description="SignalP 6.0 signal peptide prediction tool",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="healthtech.dtu.dk",
    author="Felix Teufel",
    author_email="felix.teufel@gmail.com",
    packages=['signalp', 'signalp.conversion_utils'],
    package_data = {'signalp': ['model_weights/*', 'model_weights/sequential_models_signalp6*']},
    python_requires=">=3.6, <4",
    install_requires=requirements,
    scripts=['signalp/conversion_utils/signalp6_convert_models'],
    entry_points = {"console_scripts":['signalp6=signalp:predict']},

)