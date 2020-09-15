Quick start
===========

In this section we will follow a step-by-step process to review the several stages involved in the generation of a sample material dataset.
For this, we will create a project for the material ``COMPOSITE_01`` and sample a unit cell of this material using 9 trajectories.

To follow along, make sure the :ref:`installation <install>` is complete before proceding, no need to configure it just yet.

Let's get started.

Sampling
--------

Create root directory for our material::

  >>> mkdir COMPOSITE_01
  >>> cd COMPOSITE_01

The next step generates an initial configuration file.
Also, it creates a base directory structure and copies template case files
from a specified location. For this tutorial, we will use the template files
bundled with the installation files in the ``sample`` directory of this project::

  >>> pwd
  COMPOSITE_01
  >>> python hpr.py init .../hprfe2_project/sample/template_case
  Written configuration file configuration.json.
  Created sampling directory sampling
  Template files copied to sampling directory
  >>> ls
  configuration.json sampling
  >>> ls sampling
  MainKratos.py   model.mdpa              ProjectParameters_quiet.json
  materials.json  ProjectParameters.json  _training_strain_set.dat

The script also created a ``sampling`` directory and populate it with our
Kratos case, which includes ``MainKratos.py``, ``model.mdpa`` (unit cell
discretization), ``materials.json`` (COMPOSITE_01 constituve model and
material parameters), ``ProjectParameters.json`` (case configuration
for Kratos), and ``_training_strain_set.dat`` with the list of strains for the
sampling process.

The configuration file contains the at least following needed parameters:

.. code-block:: json
  :emphasize-lines: 4,8

    {
      "config_data": {
        "cases_test_dataset": [ 5 ],
        "rve_data_points": [ 200, 400 ],
        "rve_data_points_range_list": [[100, 1600, 20], [1600, 2600, 100]],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30, 40, 50, 60],
        "strain_svd_cutoff": 0.1,
      }
    }


In the following step, we generate the sampling directories::

  >>> python hpr.py deploy
  case_0 [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  case_1 [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
  case_2 [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
  case_3 [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
  case_4 [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
  case_5 [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  case_6 [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
  case_7 [1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
  case_8 [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]

This process create the file structure for each trajectory,
and populates them with the Kratos case.
The only difference between them is the strain tensor value in their respective
``ProjectParameters.json`` files.

At this point, we should have the following file structure (here showing only ``case_0``, as it is the same for the other directories)::

  COMPOSITE_01
  ├── configuration.json
  └── training
      ├── case_0
      │   ├── MainKratos.py
      │   ├── materials.json
      │   ├── model.mdpa
      │   └── ProjectParameters.json
      ├── case_0
      ├── case_1
      ├── case_2
      ├── case_3
      ├── case_4
      ├── case_5
      ├── case_6
      ├── case_7
      ├── case_8
      ├── MainKratos.py
      ├── materials.json
      ├── model.mdpa
      ├── ProjectParameters.json
      └── _training_strain_set.dat

(only the files relevant to this turorial are shown, there are more auxiliar
files in this directories for more complex use.)

We must now run every case.
In this tutorial, we just enter to each directory and run Kratos::

  >>> pwd
  COMPOSITE_01
  >>> cd sampling
  >>> cd case_0
  >>> python3 MainKratos.py
  >>> cd ..
  >>> cd case_1
  >>> python3 MainKratos.py
  >>> cd ..
  ...

but in real-life cases we should have our own script for managing the jobs (more on this later).

Basis generation
----------------
.. note::
    Add configuration option here.

Aca una muestra de codigo::

  >>> python hola.py
  Hola.
  >>> ls
  lsout

Datasets generation
-------------------


