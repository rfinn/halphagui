from setuptools import setup, find_packages

setup(
    name="hatools",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "astropy",
        "photutils",
        "statmorph",
        "matplotlib",
        "scipy",
    ],
    python_requires=">=3.9",
    include_package_data=True,
    package_data={"hatools": ["filter_traces/*.fits"]},  #
)
