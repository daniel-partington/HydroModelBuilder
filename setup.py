try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'HydroModelBuilder',
    'author': 'Daniel Partington, Takuya Iwanaga, Michael Asher',
    'url': 'https://github.com/daniel-partington/HydroModelBuilder',
    'download_url': 'https://github.com/daniel-partington/HydroModelBuilder/archive/master.zip',
    'author_email': 'dpartington1982@gmail.com',
    'version': '0.1',
    'setup_requires':['numpy', 'gdal'],  # known issue in numpy - have to have it in setup_requires as well
    'install_requires':['nose', 'flopy', 'numpy', 'pandas', 'gdal'],
    'packages': ['HydroModelBuilder'],
    'scripts': [],
    'name': 'hydro_model_builder'
    }

setup(**config)	