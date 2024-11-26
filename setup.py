from setuptools import find_packages, setup

setup(
    name="docktprep",
    version="0.0.1.dev0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "docktprep=docktprep.main:main",
        ],
    },
)
