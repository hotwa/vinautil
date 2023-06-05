from setuptools import setup # auto find packages please import find_packages
from pathlib import Path
import yaml

with open('meta.yaml', 'r') as file:
    doc = yaml.load(file, Loader=yaml.FullLoader)
    version = doc['package']['version']

here = Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name='vinautil',
    version=version,
    description=(
    	'a tool for autodock vina'
    	 ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/hotwa/vinautil.git',
    author='lingyu zeng',
    author_email='pylyzeng@gmail.com',
    license='MIT License',
    classifiers=[
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords=["tools", "autodock vina"],
    package_dir={"": "vinautil"},  # Optional
    packages=['vinautil', 'vinautil.vutils', 'vinautil.pymolutils'],
    python_requires='>=3.9',
    install_requires=[],
    dependency_links=[
        "https://github.com/openbabel/openbabel"
    ],
    platforms=["linux-64"],
    entry_points={
        'console_scripts': [
            'scardock = vinautil.scardock:SCARdock',  # Replace 'my_module' with the actual module that contains SCARdock function
            'scardocktest = vinautil.scardock:SCARdocktest'
        ],
    },
    package_data={
        'vinautil': ['test/*'],
    }
)
