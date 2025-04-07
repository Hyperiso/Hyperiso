from setuptools import setup, find_packages

setup(
    name="Phyperiso",
    version="0.0.1",
    packages=find_packages(where="core"),
    package_dir={"": "core"},
    install_requires=[
    "pandas",
    "numpy"],
    python_requires=">=3.9"
)
