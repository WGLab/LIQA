import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="liqa",
    version="1.1.20",
    author="Yu Hu",
    author_email="huyu999999@gmail.com",
    description="A statistical tool to quantify isoform-specific expression using long-read RNA-seq",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WGLab/LIQA",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["pysam", "lifelines"],
    #python_requires='>=3.6',
    #packages=["liqa_src","liqa_bin"],
    #package_dir={"liqa_src":"liqa_src","liqa_bin":"liqa_src/liqa_bin"},
    packages=["liqa_src"],
    package_dir={"liqa_src":"liqa_src"},
    scripts=["liqa_src/group_process.pl","liqa_src/testDAS.R", "liqa_src/PreProcess.pl", "liqa_src/PreProcess_gtf.pl"],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "liqa=liqa_src.liqa:main",
        ]
    },
)
