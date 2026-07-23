# Data {#Data}

[![image](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gsclaros/idaes-gtep.git/HEAD?urlpath=%2Fdoc%2Ftree%2Fdocs%2Fsource%2Fdata.md)

The
`ExpansionPlanningData<gtep.gtep_data.ExpansionPlanningData>`{.interpreted-text
role="class"} class stores data for the IDAES-GTEP model. It includes a
function to load data structured using the Prescient data loader.

``` 
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
```

All the data is stored within the
`ExpansionPlanningModel<gtep.gtep_model.ExpansionPlanningModel>`{.interpreted-text
role="class"} class and used for the solution of the model. Table 1
shows a detailed description of the relevant parameters included.

+-------------------------+---------+---------+------------------------------+
| Generator Components    |         |         |                              |
+=========================+=========+=========+==============================+
| [lifetimes]{.title-ref} | Pa      |         | Lifetime of each generator   |
| [s                      | rameter |         | Startup cost for each        |
| tartupCost]{.title-ref} | Pa      |         | generator Multiplier for new |
| [capital                | rameter |         | generator investments Fuel   |
| Multiplier]{.title-ref} | Pa      |         | required to be consumed for  |
| [startFuel]{.title-ref} | rameter |         | startup process              |
|                         | Pa      |         |                              |
|                         | rameter |         |                              |
+-------------------------+---------+---------+------------------------------+
| [fuelCost]{.title-ref}  | Pa      | USD     | Cost per unit of fuel at     |
|                         | rameter |         | each generator               |
+-------------------------+---------+---------+------------------------------+
| [emiss                  | Pa      | dimens  | $CO_{2}$ emission factor for |
| ionsFactor]{.title-ref} | rameter | ionless | each generator Cost of       |
|                         |         |         | investment in each new       |
| [generatorInve          | Pa      |         | generator Maximum spinning   |
| stmentCost]{.title-ref} | rameter |         | reserve available for each   |
| [maxSpinn               | Pa      |         | generator Maximum quickstart |
| ingReserve]{.title-ref} | rameter |         | reserve available for each   |
|                         |         |         | generator Ramp up rates for  |
| [maxQuickst             | Pa      |         | each generator Ramp down     |
| artReserve]{.title-ref} | rameter |         | rates for each generator     |
|                         |         |         | Matching for each generator  |
| [r                      | Pa      |         | to its respective region     |
| ampUpRates]{.title-ref} | rameter |         |                              |
| [ram                    | Pa      |         |                              |
| pDownRates]{.title-ref} | rameter |         |                              |
| [ge                     | Pa      |         |                              |
| nsAtRegion]{.title-ref} | rameter |         |                              |
+-------------------------+---------+---------+------------------------------+
| [fixedOpe               | Pa      | USD     | Operating costs for each     |
| ratingCost]{.title-ref} | rameter |         | generator                    |
| **Thermal Specific**    |         |         |                              |
| [therm                  | Pa      |         | Maximum output of each       |
| alCapacity]{.title-ref} | rameter |         | thermal generator Minimum    |
| [                       | Pa      |         | output of each thermal       |
| thermalMin]{.title-ref} | rameter |         | generator Maximum fraction   |
| [spinningReser          | Pa      |         | of maximum thermal           |
| veFraction]{.title-ref} | rameter |         | generation output as         |
|                         |         |         | spinning reserve Maximum     |
| [quickstartReser        | Pa      |         | fraction of maximum thermal  |
| veFraction]{.title-ref} | rameter |         | generation output as         |
|                         |         |         | quickstart reserve           |
| **Renewable Specific**  |         |         | transmission line            |
+-------------------------+---------+---------+------------------------------+
| [renewab                | Pa      | MW      | Maximum capacity of each     |
| leCapacity]{.title-ref} | rameter |         | renewable generator          |
+-------------------------+---------+---------+------------------------------+
| [renewableCap           | Pa      | dimens  | Fraction of                  |
| acityValue]{.title-ref} | rameter | ionless | [re                          |
|                         |         |         | newableCapacity]{.title-ref} |
| [rene                   | Pa      | MW      | that can be counted towards  |
| wableQuota]{.title-ref} | rameter |         | planning reserve requirement |
+-------------------------+---------+---------+------------------------------+
| [curta                  | Pa      | USD/MW  | Cost of curtailed renewable  |
| ilmentCost]{.title-ref} | rameter |         | energy                       |
+-------------------------+---------+---------+------------------------------+
| [lo                     | Pa      | USD     | Cost of load shedding        |
| adShedCost]{.title-ref} | rameter |         |                              |
+-------------------------+---------+---------+------------------------------+

+-------------------------+---------+---------+------------------------------+
| Transmission Line       |         |         |                              |
| Components              |         |         |                              |
+=========================+=========+=========+==============================+
| [transmissi             | Pa      |         | Long term thermal capacity   |
| onCapacity]{.title-ref} | rameter |         | of each transmission line    |
|                         |         |         | Distance between terminal    |
| [distance]{.title-ref}  | Pa      |         | buses for each transmission  |
|                         | rameter |         | line                         |
+-------------------------+---------+---------+------------------------------+

+-------------------------+---------+---------+------------------------------+
| Branch Components       |         |         |                              |
+=========================+=========+=========+==============================+
| [branchInve             | Pa      | USD USD | Cost of investment in each   |
| stmentCost]{.title-ref} | rameter |         | new branch                   |
| [branchInve             | Pa      |         |                              |
| stmentCost]{.title-ref} | rameter |         |                              |
| [branchCapital          | Pa      |         |                              |
| Multiplier]{.title-ref} | rameter |         |                              |
| [branchExtension        | Pa      |         |                              |
| Multiplier]{.title-ref} | rameter |         |                              |
+-------------------------+---------+---------+------------------------------+

::: currentmodule
gtep.gtep_data
:::

::: {.automodule members=""}
gtep.gtep_data
:::
