Installation
============

run_dbCAN installs the commands ``run_dbcan``, ``dbcan_plot``, ``dbcan_utils`` and ``dbcan_asmfree``.
Pick one of the three methods below; conda is recommended because it also installs DIAMOND.

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Method
     - Command
   * - **Conda** (recommended)
     - ``conda env create -f environment.yml`` with the file from
       `envs/ <https://github.com/bcb-unl/run_dbcan/tree/master/envs>`_ (bioconda package: ``dbcan``)
   * - Pip
     - ``pip install dbcan``; install `DIAMOND <https://github.com/bbuchfink/diamond>`_ yourself (not on PyPI)
   * - Docker
     - ``docker pull ghcr.io/bcb-unl/run_dbcan:latest`` (or a release tag such as ``5.2.9``)

Conda
-----

.. code-block:: shell

   # download envs/environment.yml from the GitHub repository, then
   conda env create -f environment.yml
   conda activate run_dbcan_env      # the name in the "name:" line of the file
   run_dbcan --help

Pip
---

.. code-block:: shell

   pip install dbcan
   diamond version        # DIAMOND must be on PATH

Docker
------

.. code-block:: shell

   docker pull ghcr.io/bcb-unl/run_dbcan:latest
   docker run --rm -it ghcr.io/bcb-unl/run_dbcan:latest run_dbcan --help

Mount your data and database folders with ``-v`` (for example ``-v $PWD:/data -w /data``).

If something goes wrong
-----------------------

* Check that ``which python`` and ``which pip`` point into the active conda environment, not to a system Python.
  This matters most when you install extra packages.
* Update conda (``conda update -n base conda``) and read the error in the installation log.
* See :doc:`../FAQ/troubleshooting`, or open a `GitHub issue <https://github.com/bcb-unl/run_dbcan/issues>`_.

.. note::

   New releases are published automatically to PyPI and the GitHub Container Registry (GitHub workflow), and
   to bioconda by its autobump bot, so all three methods give the latest version.

Next: :doc:`Download the databases <../user_guide/prepare_the_database>`.
