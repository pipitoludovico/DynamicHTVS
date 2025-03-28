from setuptools import setup, find_packages

setup(
    name='DynamicHTVS_Screener',
    version='3.4',
    packages=find_packages(),
    package_data={
        'DynamicHTVS_lib.VMD': ['getContacts.tcl'],
        'DynamicHTVS_lib.parameters' : ["*.*"]
    },
)
