Quick start
===========

In this section we will follow a step-by-step process to review the several
stages involved in the generation of a material dataset.
For this, we will create a project for a material ``COMPOSITE_01`` and
sample a unit cell of this material using 9 trajectories.

To follow along, make sure the :ref:`installation <install>` is complete,
and (recommended) that the location of the ``hprfe2`` executable is in the
``$PATH`` environment variable before proceding.

Let's get started.

Set up
-------

Create root directory for our material:

.. code-block:: console

   $ mkdir COMPOSITE_01
   $ cd COMPOSITE_01

The next step generates an initial configuration file.
Also, it creates a base directory structure and copies template case files from
a specified location.
For this tutorial, we will use a test case bundled in the project files in the 
``utils`` directory, assuming a project directory ``~/apps/hprfe2``:

.. code-block:: console

  $ hprfe2 init ~/apps/hprfe2/utils/sampling_case ~/apps/hprfe2/utils/validation_case
  Created directory sampling
  Created directory bases
  Created directory datasets
  Created directory validation
  Written configuration file config.json.
  Template files copied to sampling directory
  Template files copied to validation directory
  $ ls
  bases  config.json  datasets  hprfe2.log  sampling  validation
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
        "rve_data_points": [100, 200],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30],
        "reconstruction_pairs": [[20, 100], [30, 200]]
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
  └── sampling
      ├── case_0
      │   ├── MainKratos.py
      │   ├── materials.json
      │   ├── model.mdpa
      │   ├── ProjectParameters.json
      │   └── ProjectParameters_quiet.json
      ├── case_0
      ├── case_1
      ├── ...
      ├── case_7
      ├── case_8
      ├── MainKratos.py
      ├── materials.json
      ├── model.mdpa
      ├── ProjectParameters.json
      ├── ProjectParameters_quiet.json
      └── strain_set.dat

(only the files relevant to this turorial are shown, there are also auxiliar
files in these directories for more complex use cases.)

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
  roc_100ip
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
``roc_100ip``, ``roc_200ip`` are the reduced sets of 100 and 200
integration points, and ``roc_ROMip`` being the complete set for ROM analysis
(as required in the configuration file).

Datasets generation
-------------------

In this step, we will create the datasets for all the combinations of number
of points (100, 200, ROM) and modes (20, 30) required in the configuration.

.. code-block:: console

  $ hprfe2 pack
  ...
  Generating datasets/rve_20m_100ip.json
  Generating datasets/rve_30m_100ip.json
  Generating datasets/rve_20m_200ip.json
  Generating datasets/rve_30m_200ip.json
  Generating datasets/rve_20m_ROMip.json
  Generating datasets/rve_30m_ROMip.json

The output files are written in the ``datasets`` directory:

.. code-block:: console

  $ ls datasets
  rve_20m_100ip.json
  rve_20m_200ip.json
  rve_20m_ROMip.json
  rve_30m_100ip.json
  rve_30m_200ip.json
  rve_30m_ROMip.json

Datasets for reconstruction
---------------------------

The following step is generate the correlation matrices required for later
reconstruction of RVE fields, for the combination of modes and points specified
in the configuration:

.. code-block:: json
  :emphasize-lines: 6

    {
      "config_data": {
        "rve_data_points": [100, 200],
        "rve_data_points_rom": true,
        "rve_data_modes": [20, 30],
        "reconstruction_pairs": [[20, 100], [30, 200]]
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
  correlation_r_value_20m_100ip.npy
  correlation_r_value_30m_200ip.npy
  correlation_strain_20m.npy       
  correlation_strain_30m.npy       
  ...

1-ip multiscale simulation
--------------------------

The ``validate`` directory is created and populated at the initialization step.
The directory contains template files for the simulation of a macrostructure
consisting in a 1-ip tetrahedron, using the reduced datasets computed previously.
It reproduces the strain state of the validation cases selected in the configuration file.

The ``validate`` module creates a file structure for testing each dataset generated.

.. code-block:: console

  $ hprfe2 validate
  case_0 _20m_100ip
  case_0 _20m_200ip
  case_0 _20m_ROMip
  case_0 _30m_100ip
  case_0 _30m_200ip
  case_0 _30m_ROMip

The file structure will be:

.. code-block:: console

  validation
  ├── case_0
  │   ├── _20m_100ip
  │   │   ├── macro_materials.json
  │   │   ├── macro_model.mdpa
  │   │   ├── MainKratos.py
  │   │   ├── ProjectParameters.json
  │   │   └── ProjectParameters_quiet.json
  │   ├── _20m_200ip
  │   ├── _20m_ROMip
  │   ├── _30m_100ip
  │   ├── _30m_200ip
  │   └── _30m_ROMip
  ├── macro_materials.json
  ├── macro_model.mdpa
  ├── MainKratos.py
  ├── ProjectParameters.json
  ├── tmp_case_0_20m_100ip.bash
  ├── tmp_case_0_20m_200ip.bash
  ├── tmp_case_0_20m_ROMip.bash
  ├── tmp_case_0_30m_100ip.bash
  ├── tmp_case_0_30m_200ip.bash
  └── tmp_case_0_30m_ROMip.bash

  
The final step in at this stage is to run the desired script(s).
In this case we will test the smaller (and faster) case:

.. code-block:: console

  $ cd validation
  $ bash tmp_case_0_20m_100ip.bash
  ... (dense output ommited) ...

All the generated files will be inside the corresponding directory.
For this case:

.. code-block:: console

  $ ll case_0/_20m_100ip/
  homogenized_stress.dat
  macro_materials.json
  macro_model.mdpa
  MainKratos.py
  Multiscale_0.post.msh
  Multiscale_0.post.res
  outMainKratos
  outMainKratos_quiet
  ProjectParameters.json
  ProjectParameters_quiet.json
  rve_runtime_data_el1_ip0.json
  time.dat
  time_quiet.dat
  vtk_output

being ``homogenized_stress.dat`` homogenized strain-stress data,
``Multiscale_0.post.*`` GiD visualization files,
``outMainKratos`` Kratos output,
``vtk_ouput`` Paraview visualization files,
``time.dat`` time information,
and ``*_quiet`` are the equivalent files but run without writing output,
useful for speedup measurements.

Reconstruction of fields
------------------------

.. warning::
    The module for fields' reconstruction is not yet integrated in HPRFE2.
    The following are temporary directions using transition scripts.

The reconstruction of the RVE fields requires datasets for each validation case.
The script will look for the required files in the ``resources`` directory,
located in the ``case_0`` directory (in this tutorial).
Resources dataset are: model and material RVE data, reduced dataset,
strain bases, strain and damage variable correlation matrices.

Resources must by collected manually for now.

.. code-block:: console

  $ mkdir case_0/resources
  $ cp ../sampling/model.mdpa  case_0/resources
  $ cp ../sampling/materials.json case_0/resources
  $ cp ../datasets/rve_20m_100ip.json case_0/resources
  $ cp ../bases/bases_STRAIN_FLUCTUANT_255m.npy case_0/resources
  $ cp ../bases/correlation_strain20m.npy case_0/resources
  $ cp ../bases/correlation_r_value_20m_100ip.npy case_0/resources

We now ``cd`` into our case's directory and run the reconstruction script,
which if located in the ``offline_scrips`` directory of the ``MultiscaleROMApplication`` project.
The script takes the project's root path, the runtime generated dataset,
and the resources file location.

.. code-block:: console

  $ python ~/offline_scripts/reconstruct_rve_variables.py ~/COMPOSITE_01 rve_runtime_data_el1_ip0.json ../resources/

The scripts writes the ``rve_reconstructed.h5`` and ``rve_reconstructed.xdmf`` files for visualization in Paraview.

