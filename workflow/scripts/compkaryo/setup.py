from setuptools import setup, find_packages

setup(name='compkaryo',
     version='0.1',
     description='computational karyotyping for An. coluzzii and gambiae',
     url='https://github.com/rrlove/compkaryo',
     author='Becca Love',
     author_email='rachelrebeccalove@gmail.com',
     license='GNU General Public License v3.0',
     packages=['compkaryo'],
     install_requires=['numpy','scikit-allel'],
     zip_safe=False,
     include_package_data=True
     )
