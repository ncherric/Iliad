.. _projectinfo/faq:

==========================
Frequently Asked Questions
==========================

.. contents::

How do I provide pre-existing file, such as a sorted BAM that I already have and want to perform variant calling on it?
-----------------------------------------------------------------------------------------------------------------------

Use the appropriate Snakemake flags. In this case, if you move the sorted BAM file to ``./Iliad/results/sortedBam/`` directory, then you must add ``--ignore-incomplete`` as long as you are certain the file is complete, there is an associated index file tagging along with it, and you are aware of the risks of running that flag.

In variant calling stage, I receive the error: ``[E:: faidx_adjust_position] The sequence "7" was not found.``?
---------------------------------------------------------------------------------------------------------------

Every once in a while, the conditional rules may fail. We believe this is a timing issue, but the quick fix is to let the jobs that are already running complete. Snakemake handles errors but shutting down after all currently running jobs complete. Start the run again but add the ``--rerun-incomplete`` flag.


I just setup conda, mamba, or micromamba package manager, I run a dry-run of workflow and no errors. But when I execute the workflow I get an error that reads: ``/usr/bin/bash: line 1: conda: command not found``
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Likely there are some missing exports in your .bashrc file (A hidden file located in your home directory). Try adding the following exports to your .bashrc as mentioned in this `link. <https://github.com/snakemake/snakemake/issues/1120#issuecomment-950735467>`_

.. code-block:: console

    export -f conda
    export -f __conda_activate
    export -f __conda_reactivate
    export -f __add_sys_prefix_to_path
    export -f __conda_hashr
    export -f __conda_exe