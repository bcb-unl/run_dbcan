.. _nextflow-hpc:

Nextflow Pipeline: HPC Configuration
====================================

This document describes the special configurations needed when running the Nextflow pipeline on HPC clusters, particularly UNL HCC (Holland Computing Center). It covers key settings for Nextflow and nf-core pipelines to ensure stable and efficient execution on Slurm-based systems.

Introduction
------------

Running Nextflow pipelines on HPC clusters requires specific configurations to ensure that:

- Jobs are properly submitted through the job scheduler (Slurm)
- Resources are requested appropriately
- Container environments work correctly
- I/O performance is optimized
- The scheduler is not overloaded

This guide provides best practices and configuration examples based on UNL HCC and nf-core recommendations.

UNL HCC Characteristics
------------------------

Important Notes
~~~~~~~~~~~~~~~

1. **Login Node Restrictions**
   - Login nodes should only be used for lightweight tasks (editing files, submitting jobs, monitoring, etc.)
   - **Do not run Nextflow pipelines on login nodes**, as they may be terminated
   - Jobs must be submitted through Slurm

2. **Resource Requests**
   - Must explicitly specify runtime (``--time``), memory (``--mem`` or ``--mem-per-cpu``), and CPU count
   - Default resources may be very limited, causing job failures
   - Maximum runtime is **7 days**

3. **Storage Locations**
   - Output should be placed in ``/work/group/user`` or ``$WORK`` directory
   - Avoid storing large amounts of data in ``$HOME``
   - Using ``/scratch`` space can improve I/O performance

4. **QoS (Quality of Service)**
   - Using ``--qos=ac_<group>`` can improve scheduling priority
   - Requires research group to have acknowledgement credit

Key Nextflow HPC Configurations
-------------------------------

Executor Configuration
~~~~~~~~~~~~~~~~~~~~~~

**Must use Slurm executor**, otherwise tasks will run on login nodes:

.. code-block:: groovy

   process {
       executor = 'slurm'
       queue = 'batch'  // Adjust according to your cluster
   }

Resource Limits
~~~~~~~~~~~~~~~

Set reasonable resource limits to prevent requesting excessive resources:

.. code-block:: groovy

   process {
       resourceLimits = [
           cpus: 32,
           memory: 128.GB,
           time: 7.d
       ]
   }

Submission Rate Limiting
~~~~~~~~~~~~~~~~~~~~~~~~

Prevent Nextflow from submitting too many tasks, which can overload the scheduler:

.. code-block:: groovy

   executor {
       name = 'slurm'
       queueSize = 50              // Limit number of simultaneously submitted tasks
       submitRateLimit = '5 sec'    // Submit one task every 5 seconds
       pollInterval = '1 min'       // Polling interval
   }

Work Directory and Temporary Directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: groovy

   // Work directory - use work storage
   workDir = '/work/$USER/nf-work'

   // Temporary directory - use scratch space
   process {
       scratch = '/scratch/$USER/$SLURM_JOB_ID'
   }

Container Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

For Singularity/Apptainer:

.. code-block:: groovy

   singularity {
       enabled = true
       autoMounts = true
       cacheDir = '/work/$USER/.singularity_cache'
       // Bind necessary directories if needed
       // runOptions = '--bind /work:/work --bind /scratch:/scratch'
   }

Environment Variables
~~~~~~~~~~~~~~~~~~~~~

Set in Slurm script:

.. code-block:: bash

   export NXF_WORK=/work/$USER/nf-work
   export NXF_SINGULARITY_CACHEDIR=/work/$USER/.singularity_cache
   export SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
   export TMPDIR=/scratch/$USER/nf_test_tmp
   export NXF_OPTS="-Xms1g -Xmx4g"  # Java memory settings

nf-core Special Settings
-------------------------

Profile System
~~~~~~~~~~~~~~

nf-core pipelines use a profile system to manage configurations for different environments:

.. code-block:: bash

   # Use multiple profiles
   nextflow run pipeline.nf -profile test,singularity,slurm

Common profiles:

- ``test`` - Use test data
- ``singularity`` - Use Singularity containers
- ``slurm`` - Use Slurm executor
- ``conda`` - Use Conda environment (not recommended on HPC)

Custom Configuration Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can create custom configuration files (e.g., ``conf/slurm.config``) and reference them:

.. code-block:: bash

   nextflow run main.nf -c conf/slurm.config -profile test,singularity

Parameter Validation
~~~~~~~~~~~~~~~~~~~~

nf-core pipelines automatically validate parameters:

.. code-block:: groovy

   params {
       validate_params = true  // Enabled by default
   }

Reports and Timeline
~~~~~~~~~~~~~~~~~~~~

nf-core pipelines automatically generate:

- Execution timeline
- Execution report
- Execution trace
- Pipeline DAG

These files are located in the ``pipeline_info/`` directory.

Common Issues and Solutions
---------------------------

Issue 1: Tasks Running on Login Nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: Slurm executor not configured

**Solution**: Set ``executor = 'slurm'`` in configuration file

Issue 2: Jobs Killed by System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: Requested resources exceed limits or job runtime too long

**Solution**: 

- Check ``resourceLimits`` settings
- Ensure individual task time does not exceed 7 days
- Consider splitting long-running tasks

Issue 3: Scheduler Overload
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: Nextflow submitting tasks too quickly

**Solution**: Set ``submitRateLimit`` and ``queueSize``

Issue 4: Poor I/O Performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: Using network storage as work directory

**Solution**: 

- Use ``/scratch`` as temporary directory
- Set ``scratch`` option
- Use ``stageInMode`` and ``stageOutMode``

Issue 5: Container Permission Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: Singularity mount issues

**Solution**: 

- Set ``autoMounts = true``
- Use ``runOptions`` to bind necessary directories
- Check file permissions

Issue 6: Insufficient Memory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause**: Java heap memory setting too small

**Solution**: Set ``NXF_OPTS`` environment variable:

.. code-block:: bash

   export NXF_OPTS="-Xms1g -Xmx4g"

Testing Workflow
----------------

1. **Small-scale Test**

   .. code-block:: bash

      sbatch run_nextflow_test.slurm

2. **Check Logs**

   .. code-block:: bash

      tail -f logs/nf-test-<jobid>.out
      tail -f logs/nf-test-<jobid>.err

3. **Verify Job Submission**

   .. code-block:: bash

      squeue -u $USER

4. **Check Output**

   .. code-block:: bash

      ls -lh test_results/

Example Slurm Script
--------------------

A complete example Slurm submission script (``run_nextflow_test.slurm``) is provided in the pipeline repository. Key components include:

- Resource requests (CPU, memory, time)
- Module loading (Nextflow, Singularity)
- Environment variable setup
- Pipeline execution command

Example Configuration File
---------------------------

A complete Slurm executor configuration file (``conf/slurm.config``) is provided in the pipeline repository. This file includes:

- Slurm executor settings
- Resource limits
- Submission rate limiting
- Singularity container configuration
- Work and temporary directory settings

Reference Resources
-------------------

- `Nextflow HPC Documentation <https://www.nextflow.io/docs/latest/executor.html#slurm-executor>`_
- `nf-core Config Repository <https://github.com/nf-core/configs>`_
- `UNL HCC Documentation <https://hcc.unl.edu/docs/>`_
- `nf-core HPC Best Practices <https://nf-co.re/docs/usage/hpc>`_
