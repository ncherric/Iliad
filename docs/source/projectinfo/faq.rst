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
