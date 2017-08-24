try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup
    import re
    import os

    def find_packages(path='.'):
        ret = []
        for root, dirs, files in os.walk(path):
            if '__init__.py' in files:
                ret.append(re.sub('^[^A-z0-9_]+', '', root.replace('/', '.')))
        return ret
    # End find_packages()
# End try

config = {
    'description': 'HydroModelBuilder',
    'author': 'Daniel Partington, Takuya Iwanaga, Michael Asher',
    'url': 'https://github.com/daniel-partington/HydroModelBuilder',
    'download_url': 'https://github.com/daniel-partington/HydroModelBuilder/archive/master.zip',
    'author_email': 'dpartington1982@gmail.com',
    'version': '0.1',
    'setup_requires': ['numpy', 'gdal'],  # known issue in numpy - have to have it in setup_requires as well
    'install_requires': ['nose', 'flopy', 'numpy', 'pandas', 'gdal'],
    'packages': find_packages(),
    'scripts': [],
    'name': 'hydro_model_builder'
}

setup(**config)
