from setuptools import setup # auto find packages please import find_packages
from pathlib import Path

here = Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name='vinautil',
    version='1.0.0',
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
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords=["tools", "autodock vina"],
    package_dir={"": "vinautil"},  # Optional
    packages=['vinautil', 'vinautil.utils', 'vinautil.pymol'],
    python_requires='>=3.8',
    install_requires=[],
    dependency_links=[
        "https://github.com/openbabel/openbabel"
    ],
    platforms=["win-64", "osx-64", "linux-64", "osx-arm64", "any"],
)
