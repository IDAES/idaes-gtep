.. _Data:

Data
====

The `ExpansionPlanningData()` class storages standard data for the
IDAES-GTEP model. Within, it includes a function to load data
structure using the Prescient data loader.

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

All the data is stored within the `ExpansionPlanningModel()` class and
used for the solution of the model. Table 1 shows a detailed
description of the relevant parameters included.

.. table:: Table 1: Data needed per component in the
           `ExpansionPlanningModel()` class
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


.. currentmodule:: gtep.gtep_data

.. automodule:: gtep
    :members:

.. automodule:: gtep.gtep_data
    :members: 
    
