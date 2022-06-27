from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

VERSION = '0.0.1'
DESCRIPTION = 'Provide kernels of Alchemical Integral Transform'

# Setting up
setup(
    name="pyalchemy",
    version=VERSION,
    author="Simon León Krug",
    author_email="<simonleon.krug@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    url="https://github.com/SimonLeonKrug/pyalchemy"
    packages=find_packages(),
    install_requires=[],
    keywords=['python', 'alchemy', 'APDFT', 'AIT', 'transform'],
    classifiers=[
        "Intended Audience :: Computational Chemists/Materials Physicists",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Inependent",
    ]
)
