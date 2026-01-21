import setuptools


setuptools.setup(
    name="pcctool",
    version="1.0.0",
    author="kangLab yijiyang",
    author_email="15623109189@163.com",
    description="A chain tool for confident position lift between hg38 and yao",
    url="",
    install_requires=[],
    python_requires='>=3.10',
    packages=setuptools.find_packages(),
    entry_points={'console_scripts': ['pcctool = yaolinkhg38chain.pcc:main'], },
)
