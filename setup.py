from setuptools import setup, find_packages

setup(
    name='DynamicHTVS_Screener',
    version='2.1',
    packages=find_packages(),
    package_data={
        'DynamicHTVS_lib.VMD': ['getContacts.tcl'],
    },
)
