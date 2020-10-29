.. _install:

Install
=======

In order to install hprfe2, you need to clone its git repository
(here assuming it is installed in the ``~/apps/`` directory)::

   >>> cd ~/apps
   >>> git clone https://github.com/marandra/hprfe2.git

For convenience, make sure ``~/apps/hprfe2/hprfe2/`` is in the
``$PATH`` environment variable.

..
   tip::
     There is a simple tab-completion bash script located in the ``utils``
     directory that you can source for enabling tab-completion (in bash):
  
     ``echo "source $HOME/apps/hprfe2/utils/hprfe2-completions.bash" >> ~HOME/.bashrc``

HPRFE2 dependencies
-------------------

HPRFE2 requires a working installation of KratosMultiphysics compiled and 
MultiscaleROMApplication.

KratosMultiphysics
~~~~~~~~~~~~~~~~~~

For KratosMultiphysics, refer to its
`installation page <https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md>`_.

MultiscaleROMApplication
~~~~~~~~~~~~~~~~~~~~~~~~

Clone the git repository::

    $ cd ~/apps/Kratos/applications
    $ git clone https://github.com/marandra/MultiscaleROMApplication.git

Add the following line to ...::

    To be completed...
    ...
    ...

and compile Kratos.

Check the installation (there should not be error messages)::
    $ python -m KratosMultiphysics
    $ python -m KratosMultiphysics.StructuralMechanicsApplication
    $ python -m KratosMultiphysics.MultiscaleROMApplication

Python modules
~~~~~~~~~~~~~~

Required python modules are in ``requirements.txt``, in the project directory.
Check that they are intalled, or install them if missing::
    $python -m pip install -r ~/apps/hprfe2/requirements.txt --user

In case you prefer to use a virtual environment, to keep things tidy::
    $python -m venv venv-hprfe2
    $source venv-hprfe2/bin/activate
    $python -m pip install -r ~/apps/hprfe2/requirements.txt install

Usage
=====

Usage.

Configuration
-------------

The default configuration file of each material is ``configuration.json``.
It must be present and located at the root directory of the material.

An initial configuration file with default parameters can be generated with::
$ hprfe2 config

This step generates a ``config.json`` similar to::
  
    {
      "config_data": {
        "cases_test_dataset": [ 5 ],
        "rve_data_points": [ 200, 400 ],
        "rve_data_points_range_list": [[100, 1600, 20], [1600, 2600, 100]],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30, 40, 50, 60],
      }
    }
    

Sampling
--------
Sampling usage.


Bases generation
----------------
Bases generation usage.

