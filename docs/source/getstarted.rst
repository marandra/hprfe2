Quick start
===========

In this section we will follow a step-by-step process to review the several
stages involved in the generation of a material dataset.
For this, we will create a project for the material ``COMPOSITE_01`` and
sample a unit cell of this material using 9 trajectories.

To follow along, make sure the :ref:`installation <install>` is complete,
and (recommended) that the location of the ``hprfe2`` executable is in the
``$PATH`` environment variable before proceding.

Let's get started.

Set up
-------

.. todo:: 
   Repeat all tutorial from a nested directory


Create root directory for our material:

.. code-block:: console

   $ mkdir COMPOSITE_01
   $ cd COMPOSITE_01

The next step generates an initial configuration file.
Also, it creates a base directory structure and copies template case files from
a specified location.
For this tutorial, we will use a test case bundled in the project files in the 
``utils`` directory:

.. code-block:: console

  $ hprfe2 init /path/to/hprfe2_project/utils/template_case
  Created directory sampling
  Created directory bases
  Created directory datasets
  Written configuration file config.json.
  Template files copied to sampling directory
  $ ls
  bases  config.json  datasets  hprfe2.log  sampling
  $ ls sampling
  MainKratos.py   model.mdpa              ProjectParameters_quiet.json
  materials.json  ProjectParameters.json  strain_set.dat

The script creates the required directories, and
populate the ``sampling`` directory with our Kratos case, which includes
``MainKratos.py``, ``model.mdpa`` (unit cell discretization), 
``materials.json`` (COMPOSITE_01 constituve model and material parameters),
``ProjectParameters.json`` (case configuration for Kratos),
and ``strain_set.dat`` with the list of strains for the sampling process.

The configuration file contains the defautl values for the most frequently
set parameters:

.. code-block:: json

    {
      "config_data": {
        "rve_data_points": [150, 200],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30],
        "reconstruction_pairs": [[20, 150], [30, 200]]
      }
    }

No need to modify it just yet, we will continue the tutorial with the default
values.

Sampling
--------

In the following step, we generate the sampling directories:

.. code-block:: console

  $ hprfe2 deploy
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

At this point, we should have the following file structure (here showing only
``case_0``, as it is the same for the other directories):

.. code-block:: console

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
      └── strain_set.dat

(only the files relevant to this turorial are shown, there are more auxiliar
files in this directories for more complex use.)

We must now run every case.
In this tutorial, we just enter to each directory and run Kratos:

.. code-block:: console

  $ pwd
  COMPOSITE_01
  $ cd sampling
  $ cd case_0
  $ python3 MainKratos.py
  ... (Kratos output omitted) ...
  $ cd ..
  $ cd case_1
  $ python3 MainKratos.py
  ... (Kratos output omitted) ...
  $ cd ..
  ...

but in real-life cases we should have our own script for managing the jobs 
(more on this later).

Bases generation
----------------

Now we generate the modal bases and the ROC integration points reduced sets,
as well as other auxiliar files:

.. code-block:: console

  $ hprfe2 generate
  ... (dense output omitted) ...

This step will write output files in the ``bases`` directory:

.. code-block:: console

  $ ls bases
  bases_ENERGY_FREE_284m.npy
  bases_R_VALUE_31m.npy
  bases_STRAIN_FLUCTUANT_255m.npy
  roc_150ip
  roc_200ip
  roc_ROMip
  sv_ENERGY_FREE_elastic.dat
  sv_ENERGY_FREE_inelastic.dat
  sv_R_VALUE_elastic.dat
  sv_R_VALUE_inelastic.dat
  sv_STRAIN_FLUCTUANT_elastic.dat
  sv_STRAIN_FLUCTUANT_inelastic.dat

``bases_STRAIN_FLUCTUANT_255m.npy``,
``bases_ENERGY_FREE_284m.npy``,
``bases_R_VALUE_31m.npy``
are modal bases for strain, energy and r-value, respectively, with 
``sv_*.dat`` files being the modes' corresponding singular values.
``roc_150ip``, ``roc_200ip`` are the reduced sets of 150 and 200
integration points, and ``roc_ROMip`` being the complete set for ROM analysis
(as required in the configuration file).

Datasets generation
-------------------

In this step, we will create the datasets for all the combinations of number
of points (150, 200, ROM) and modes (20, 30) required in the configuration.

.. code-block:: console

  $ hprfe2 pack
  ...
  Generating datasets/rve_20m_150ip.json
  Generating datasets/rve_30m_150ip.json
  Generating datasets/rve_20m_200ip.json
  Generating datasets/rve_30m_200ip.json
  Generating datasets/rve_20m_ROMip.json
  Generating datasets/rve_30m_ROMip.json

The output files are written in the ``datasets`` directory:

.. code-block:: console

  $ ls datasets
  rve_20m_150ip.json
  rve_20m_200ip.json
  rve_20m_ROMip.json
  rve_30m_150ip.json
  rve_30m_200ip.json
  rve_30m_ROMip.json

Reconstruction data
------------------------------

The following step is generate the correlation matrices required for later
reconstruction of RVE fields, for the combination of modes and points specified
in the configuration:

.. code-block:: json
  :emphasize-lines: 6

    {
      "config_data": {
        "rve_data_points": [150, 200],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30],
        "reconstruction_pairs": [[20, 150], [30, 200]]
      }
    }

.. code-block:: console

  $ hprfe2 resources
  ... (dense output omitted) ...

Besides the previous files in the ``bases`` directory, we now find the
correlation matrices 

.. code-block:: console

  $ ls datasets
  ...
  correlation_r_value_20m_150ip.npy
  correlation_r_value_20m_200ip.npy
  correlation_r_value_30m_200ip.npy
  correlation_strain_20m.npy       
  correlation_strain_30m.npy       
  ...
