from setuptools import setup

setup(
    # Your setup arguments
    python_requires='>=3.6',  # Your supported Python ranges
    name = "molcompview",
    version = "1.0.0",
    description = "MolCompass Viewer",
    author = "Sergey Sosnin <serg.sosnin@gmail.com>",
    include_package_data=False,
    install_requires=[
        'molcomplib',
        'rdkit',
    ],
    packages=["molcompview"],
    entry_points={
        'console_scripts': [
            'mcv = molcompview.main:main',
        ],
    },
    license = "MIT",
)
