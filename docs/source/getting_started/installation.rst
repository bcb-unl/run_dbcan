Installation
=============

We support multiple methods for installing run_dbCAN to accommodate different user preferences and requirements.

Installation Options
---------------------

.. admonition:: Available Methods

   - **Conda** (recommended)
   - **Pip**
   - **Docker**

Conda Installation (Recommended)
---------------------------------

We recommend using Conda for installation as it provides the most stable environment and simplifies dependency management.

Before installation, we recommend verifying your environment configuration:

.. code-block:: shell

   # Check conda version
   conda --version

   # Update all packages to the latest versions
   conda update --all

   # Verify Python and pip paths (please also double-check it after the installation)
   which python
   which pip

To install run_dbCAN using Conda, follow these steps:

1. Download the environment file from our GitHub repository:

   `run_dbCAN Conda Environment Files <https://github.com/bcb-unl/run_dbcan/tree/master/envs>`_

2. Create and activate the environment:

   .. code-block:: shell

      conda env create -f environment.yml
      conda activate run_dbcan_env

   ``envs/environment.yml`` creates an environment named ``run_dbcan_env``; use the name given in the
   ``name:`` line of the file you downloaded.

.. hint::
    If you encounter any issues during installation, please refer to the `Troubleshooting` section.


Pip Installation
------------------

For users who prefer pip, the package is available on PyPI. However, please note that you will need to install **Diamond** separately as it is not available through PyPI.

.. code-block:: shell

   pip install dbcan

.. note::
   Before using the pip installation, ensure that **Diamond** is properly installed in your environment.

Docker Installation
-------------------

We also provide a Docker image for users who prefer containerized environments. You can pull the image from the GitHub Container Registry (tags: ``latest`` or a release version such as ``5.2.9``):

.. code-block:: shell

   docker pull ghcr.io/bcb-unl/run_dbcan:latest
   docker run --rm -it ghcr.io/bcb-unl/run_dbcan:latest run_dbcan --help

Troubleshooting
-------------------

If you encounter any issues during installation, please refer to the following troubleshooting tips:

1. Check the conda version
2. Update all packages to the latest versions
3. Verify the python and pip paths, please make sure the path is activated with `current env` not others.
4. Check the installation logs for any error messages
5. Check :doc:`../FAQ/troubleshooting` for known limitations and common errors.
6. If you still meet any issues, please feel free to contact us either via `GitHub issues <https://github.com/bcb-unl/run_dbcan/issues>`_ or by email.



.. warning::
   If you have multiple Python/pip installations, ensure that you're using the correct versions from your conda environment. This is **especially important** when installing additional dependencies.


Why We Use PyPI, bioconda, and docker
--------------------------------------

1. We've uploaded run_dbCAN to PyPI to simplify the installation process. This eliminates the need to clone the entire repository from GitHub. Users only need to download the environment files, which are available in the GitHub repository.
2. We also set the `automatic deployment` via `Github Workflow` to PyPI/docker, which means that the latest version will always be available for installation.
3. Bioconda also provides the `autobump` to upload the latest version to the bioconda channel.
