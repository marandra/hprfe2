.. _install:

Install
=======

HPRFE2 needs to be installed inside a previous Kratos installation.

.. todo::
        Make sure to add definition of required environment variables.

.. todo::
        Add reference to Kratos installation.

.. note::
        Assuming Kratos already installed

Basic steps::

   >>> cd /path/to/Kratos/applications
   >>> git checkout MultiscaleROMApplication
   >>> cd MultiscaleROMApplication

Using venv
----------

It is more convenient to create a virtual environment::
    >>> python -m pyvenv venv
    >>> source venv/bin/activate.sh
    >>> pip install -r requirements.txt

 

Usage
=====

Usage.

Configuration
-------------

The default configuration file of each material is ``configuration.json``.
It must be present and located at the root directory of the material.

An initial configuration file with default parameters can be generated with 
>>> python offline_common.py --init

This step generates a ``configuration.json`` similar to::
  
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
    

Sampling
--------
Sampling usage.


Bases generation
----------------
Bases generation usage.

