import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blockdiagram",
    version="0.1.0",
    author="Zoltan Sylvester",
    author_email="zoltan.sylvester@beg.utexas.edu",
    description="Module for creating block diagrams and other three-dimensional displays from stratigraphic models",
    keywords = 'geomorphology, stratigraphy, 3D visualization',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zsylvester/blockdiagram",
    packages=['blockdiagram'],
    scripts=['/Users/zoltan/Dropbox/Git/blockdiagram/blockdiagram/blockdiagram.py'],
    install_requires=['numpy','matplotlib','mayavi','scipy','pillow'],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
