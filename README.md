<h3 align="center">
  <img src="https://raw.githubusercontent.com/bcb-unl/run_dbcan/master/docs/source/_static/img/run_dbcan_v5_logo.png" width="400" alt="dbCAN-logo"/><br/>

  run_dbcan - Standalone Tool of <a href="http://bcb.unl.edu/dbCAN2/">dbCAN3</a>

</h3>

<p align="center">
  <a href="https://github.com/bcb-unl/run_dbcan/stargazers"><img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/bcb-unl/run_dbcan?style=for-the-badge&logo=starship&labelColor=363a4f&color=b7bdf8"></a>
  <a href="https://pypi.org/p/dbcan/"><img alt="PyPI - Version" src="https://img.shields.io/pypi/v/dbcan?style=for-the-badge&logo=pypi&color=74c7ec"></a>
  <a href="https://anaconda.org/bioconda/dbcan"><img alt="Conda Downloads" src="https://img.shields.io/conda/dn/bioconda/dbcan?style=for-the-badge&logo=anaconda&labelColor=363a4f&color=a6da95"></a>
  <a href="https://run-dbcan.readthedocs.io/en/latest/"><img alt="Read the Docs" src="https://img.shields.io/readthedocs/dbcan?style=for-the-badge&logo=Read%20the%20Docs&labelColor=363a4f&color=cba6f7"></a>
  <a href="https://github.com/bcb-unl/run_dbcan/issues"><img alt="GitHub Issues or Pull Requests" src="https://img.shields.io/github/issues/bcb-unl/run_dbcan?style=for-the-badge&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNTYgMjU2Ij4KPHBhdGggZD0iTTIxNiwzMlYxOTJhOCw4LDAsMCwxLTgsOEg3MmExNiwxNiwwLDAsMC0xNiwxNkgxOTJhOCw4LDAsMCwxLDAsMTZINDhhOCw4LDAsMCwxLTgtOFY1NkEzMiwzMiwwLDAsMSw3MiwyNEgyMDhBOCw4LDAsMCwxLDIxNiwzMloiIHN0eWxlPSJmaWxsOiAjQ0FEM0Y1OyIvPgo8L3N2Zz4%3D&labelColor=363a4f&color=f5a97f"></a>
  <br/>
  <a href="#"><img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/dbcan?style=for-the-badge&logo=python&labelColor=363a4f&color=99d1db"></a>
  <a href="https://github.com/bcb-unl/run_dbcan/releases/latest"><img alt="GitHub Release" src="https://img.shields.io/github/v/release/bcb-unl/run_dbcan?style=for-the-badge&logo=github&labelColor=363a4f&color=89dceb"></a>
  <a href="https://github.com/bcb-unl/run_dbcan/blob/master/LICENSE"><img alt="GitHub License" src="https://img.shields.io/github/license/bcb-unl/run_dbcan?style=for-the-badge&labelColor=363a4f&color=eba0ac"></a>
  <a href="https://github.com/bcb-unl/run_dbcan/actions/workflows/build_dbcan_docker.yml"><img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/bcb-unl/run_dbcan/build_dbcan_docker.yml?branch=master&style=for-the-badge&logo=github&labelColor=363a4f&color=f2cdcd"></a>
  <a href="https://github.com/bcb-unl/run_dbcan/actions/workflows/test_dbcan.yml"><img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/bcb-unl/run_dbcan/test_dbcan.yml?branch=master&style=for-the-badge&logo=github&label=Test&labelColor=363a4f&color=a6d189"></a>
</p>

## Update
10/20/2025:
1. **SignalP6.0 Topology Annotation**: Added support for SignalP6.0 signal peptide prediction. Use `--run_signalp` flag in `CAZyme_annotation` command to enable topology annotation. Results are automatically added to the overview.tsv file.
2. **Global Logging System**: Implemented comprehensive logging system with `--log-level`, `--log-file`, and `--verbose` options for better debugging and monitoring.
3. **Database Download Command**: Added new `database` command for easy database downloading. Supports both HTTP and AWS S3 sources (use `--aws_s3` flag for faster downloads). Use `--cgc/--no-cgc` to control CGC-related database downloads.
4. **Code Structure Improvements**: Continued refactoring with object-oriented programming, improved modularity, and centralized configuration management.

5/12/2025:
`dev-dbcan` branch is used to test new functions and fix issues. After testing, this branch will be merged into the main branch and update docker/conda/pypi. If you want to use those beta functions, please replace the code folder (dbcan) with your current package.

3/16/2025:
1. Rewrite the structure of <a href="https://github.com/linnabrown/run_dbcan">run_dbcan 4.0</a> (suggested by Haidong), using object-oriented programming (OOP) to improve maintainability and readability.
2. Added new function: cgc_circle, which can visualize CGC in genome.

**Future plans** Add prediction of food consumption through CAZyme. If you have new suggestions, please contact Dr. Yanbin Yin (yyin@unl.edu), Xinpeng Zhang (xzhang55@huskers.unl.edu), and Dr. Haidong Yi (hyi@stjude.org).

## Introduction

**Notice**

This is the updated version of <a href="https://github.com/linnabrown/run_dbcan">run_dbcan 4.0</a>. Many changes have been made and described in https://run-dbcan.readthedocs.io/en/latest/. From now on, this repo is the official run_dbcan site, and the site at <a href="https://github.com/linnabrown/run_dbcan">run_dbcan 4.0</a> will be no longer maintained.

**run_dbcan** is the standalone version of the [dbCAN3](http://bcb.unl.edu/dbCAN2/) annotation tool for automated CAZyme annotation. This tool, known as `run_dbcan`, incorporates pyHMMER (replacing HMMER for better performance), Diamond, and dbCAN_sub for annotating CAZyme families, and integrates CAZyme Gene Clusters (CGCs) and substrate predictions.

## Main Commands

The tool provides the following main commands:

- `database` - Download dbCAN databases (supports HTTP and AWS S3)
- `CAZyme_annotation` - Annotate CAZymes using Diamond, pyHMMER, and dbCAN-sub
- `gff_process` - Generate GFF files for CGC identification
- `cgc_finder` - Identify CAZyme Gene Clusters (CGCs)
- `substrate_prediction` - Predict substrate specificities of CGCs
- `cgc_circle_plot` - Generate circular plots for CGCs
- `easy_CGC` - Complete CGC analysis pipeline (annotation + GFF processing + CGC identification)
- `easy_substrate` - Complete CGC analysis with substrate prediction
- `Pfam_null_cgc` - Annotate null genes in CGCs using Pfam

All commands support global logging options: `--log-level`, `--log-file`, and `--verbose`.

For usage discussions, visit our [issue tracker](https://github.com/bcb-unl/run_dbcan/issues). To learn more, read the [dbcan doc]. If you're interested in contributing, whether through issues or pull requests, please review our contribution guide.

## Reference

Please cite the following `dbCAN` publications if you use `run_dbcan` in your research:

> **dbCAN3: automated carbohydrate-active enzyme and substrate annotation**
>
> Jinfang Zheng, Qiwei Ge, Yuchen Yan, Xinpeng Zhang, Le Huang, Yanbin Yin,
>
> Nucleic Acids Research, 2023;, gkad328, doi: [10.1093/nar/gkad328](https://doi.org/10.1093/nar/gkad328).

> **dbCAN2: a meta server for automated carbohydrate-active enzyme annotation**
>
> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin
>
> Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, doi: [10.1093/nar/gky418](https://doi.org/10.1093/nar/gky418).

> **dbCAN-seq: a database of carbohydrate-active enzyme (CAZyme) sequence and annotation**
>
> Le Huang, Han Zhang, Peizhi Wu, Sarah Entwistle, Xueqiong Li, Tanner Yohe, Haidong Yi, Zhenglu Yang, Yanbin Yin
>
> Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D516–D521, doi: [10.1093/nar/gkx894\*](https://doi.org/10.1093/nar/gkx894*).

[dbcan doc]: https://run-dbcan.readthedocs.io/en/latest/
