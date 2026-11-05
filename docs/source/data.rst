.. _Data:

Data
====

The :class:`ExpansionPlanningData<gtep.gtep_data.ExpansionPlanningData>` class stores data for the
IDAES-GTEP model. It includes a function to load data
structured using the Prescient data loader.

.. code-block::

   def load_prescient(self, data_path, options_dict=None):

      self.data_type = "prescient"
      options_dict = {
          "data_path": data_path,
          "input_format": "rts-gmlc",
          "start_date": "01-01-2020",
          "num_days": 365,
          "sced_horizon": 1,
          "sced_frequency_minutes": 60,
          "ruc_horizon": 36,
      }

      prescient_options = PrescientConfig()
      prescient_options.set_value(options_dict)

All the data is stored within the :class:`ExpansionPlanningModel<gtep.gtep_model.ExpansionPlanningModel>` class and
used for the solution of the model. Table 1 shows a detailed
description of the relevant parameters included.

.. table:: Table 1: Data needed per component in the
           :class:`ExpansionPlanningModel<gtep.gtep_model.ExpansionPlanningModel>` class
   :widths: 25 10 10 30
	    
   ============================ ========= ============= ============================================
   Grid Components              Type      Units         Description                                 
   ============================ ========= ============= ============================================
   `loads`                      Parameter               Demand at each bus
   `lossRate`                   Parameter               Loss rate for each transmission line
   `extensionMultiplier`        Parameter               Cost of life extension multiplier
   `minOperatingReserve`        Parameter               Minimum operating reserve
   `minSpinningReserve`         Parameter               Minimum spinning reserve
   `peakLoad`                   Parameter MW            Maximum level of power demand
   `reserveMargin`              Parameter MW            Unused power determined as a fraction of
                                                        `peakLoad`
   `weights`                    Parameter dimensionless 
   `investmentFactor`           Parameter dimensionless
   `deficitPenalty`             Parameter MW            Generation deficits
   ============================ ========= ============= ============================================

.. table:: 
   :widths: 25 10 10 30

   ============================ ========= ============= ============================================
   Generator Components
   ============================ ========= ============= ============================================
   `lifetimes`                  Parameter               Lifetime of each generator              
   `startupCost`                Parameter               Startup cost for each generator
   `capitalMultiplier`          Parameter               Multiplier for new generator investments
   `startFuel`                  Parameter               Fuel required to be consumed for startup
                                                        process
   `fuelCost`                   Parameter USD           Cost per unit of fuel at each generator
   `emissionsFactor`            Parameter dimensionless :math:`CO_{2}` emission factor for each
                                                        generator
   `generatorInvestmentCost`    Parameter               Cost of investment in each new generator
   `maxSpinningReserve`         Parameter               Maximum spinning reserve available for each
                                                        generator
   `maxQuickstartReserve`       Parameter               Maximum quickstart reserve available for
                                                        each generator
   `rampUpRates`                Parameter               Ramp up rates for each generator
   `rampDownRates`              Parameter               Ramp down rates for each generator
   `gensAtRegion`               Parameter               Matching for each generator to its
                                                        respective region
   `fixedOperatingCost`         Parameter USD           Operating costs for each generator
   **Thermal Specific**
   `thermalCapacity`            Parameter               Maximum output of each thermal generator
   `thermalMin`                 Parameter               Minimum output of each thermal generator
   `spinningReserveFraction`    Parameter               Maximum fraction of maximum thermal
                                                        generation output as spinning reserve       
   `quickstartReserveFraction`  Parameter               Maximum fraction of maximum thermal
                                                        generation output as quickstart reserve     
                                                        transmission line
   **Renewable Specific**
   `renewableCapacity`          Parameter MW            Maximum capacity of each renewable generator
   `renewableCapacityValue`     Parameter dimensionless Fraction of `renewableCapacity` that can be
                                                        counted towards planning reserve requirement
   `renewableQuota`             Parameter MW
   `curtailmentCost`            Parameter USD/MW        Cost of curtailed renewable energy
   `loadShedCost`               Parameter USD           Cost of load shedding
   ============================ ========= ============= ============================================

.. table:: 
   :widths: 25 10 10 30

   ============================= ========= ============= ============================================
   Transmission Line Components 
   ============================= ========= ============= ============================================
   `transmissionCapacity`        Parameter               Long term thermal capacity of each
                                                         transmission line
   `distance`                    Parameter               Distance between terminal buses for each
                                                         transmission line
   ============================= ========= ============= ============================================

.. table:: 
   :widths: 25 10 10 30

   ============================ ========= ============= ============================================
   Branch Components             
   ============================ ========= ============= ============================================
   `branchInvestmentCost`       Parameter USD           Cost of investment in each new branch
   `branchInvestmentCost`       Parameter USD
   `branchCapitalMultiplier`    Parameter
   `branchExtensionMultiplier`  Parameter
   ============================ ========= ============= ============================================


Optional hydro inputs
---------------------

Two optional files give hydro units a minimum output and an inflow
(energy) budget. They are read by
:meth:`load_hydro_timeseries<gtep.gtep_data.ExpansionPlanningData.load_hydro_timeseries>`,
which ``load_prescient`` calls automatically, and are used only when the
model options ``advanced_hydro`` and ``hydro_budget`` are both ``True``
(see the advanced hydropower section of the modeling page). The files
are still read and checked when ``hydro_budget`` is ``False``, but the
model ignores them.

.. table:: Table 2: Optional hydro input files
   :widths: 25 20 45

   ============================== ==================== =================================================================
   File                           Units                Description
   ============================== ==================== =================================================================
   ``DAY_AHEAD_hydro_min.csv``    MW                   Minimum output of each listed hydro unit in each period.
                                                       Stored as ``hydro_p_min``.
   ``DAY_AHEAD_hydro_budget.csv`` MW (average)         Inflow to each listed unit's reservoir, as average MW of
                                                       generation-equivalent energy in each period (for hourly data,
                                                       the same number as MWh per period). Stored as ``hydro_budget``.
   ============================== ==================== =================================================================

Both files use the layout of ``DAY_AHEAD_renewables.csv``: the columns
``Year``, ``Month``, ``Day`` and ``Period`` (counting from 1 within each
day), followed by one column per hydro generator ``GEN UID``:

.. code-block:: text

   Year,Month,Day,Period,1_HYDRO
   2020,1,1,1,4.2
   2020,1,1,2,3.8

- Each file may list any subset of the ``HYDRO`` units. A column that is
  not a ``HYDRO`` generator in ``gen.csv`` raises an error, as do blank
  or negative values and dates that do not cover the simulation.
- Periods have the length given by the ``DAY_AHEAD`` ``Period_Resolution``
  in ``simulation_objects.csv`` (3600 s, hourly, in the bundled cases),
  and ``Period`` counts from 1 within each day. Timestamps are built the
  same way as EGRET's RTS-GMLC parser, which Prescient uses for the other
  ``DAY_AHEAD`` files. Values are then aligned to the model's time steps
  (averaged when the files are finer, held when they are coarser).
- A unit without a column uses its ``p_max`` series instead, for the
  minimum, the inflow, or both. With neither file a unit runs at
  ``p_max``, the same as with ``hydro_budget`` set to ``False``.
- The files are listed, with their required columns, in
  ``gtep/data/required_data_files.csv``.

.. note::

   The existing hydro ``p_max`` and ``p_min`` series are filled in by
   Prescient from ``REAL_TIME_renewables.csv`` (the value at the start of
   each model time step), not from ``DAY_AHEAD_renewables.csv``.
   ``timeseries_pointers.csv`` in the bundled cases points both
   ``PMin MW`` and ``PMax MW`` at the same column, which is why hydro
   output was fixed before these inputs were added.


.. currentmodule:: gtep.gtep_data


.. automodule:: gtep.gtep_data
    :members: