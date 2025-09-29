IDAES-GTEP Tutorial Notebook
============================

Presented & last updated 9/19/24
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This notebook is intended as an introductory tutorial to using the
IDAES-GTEP tool. It walks through loading a small test case (PJM 5-bus)
and solves expansion planning models with a few different assumptions on
the network. It demonstrates some basic result visualizations on the
investment options and grid operations.

.. code:: ipython3

    from gtep.gtep_model import ExpansionPlanningModel
    from gtep.gtep_data import ExpansionPlanningData
    from gtep.gtep_solution import ExpansionPlanningSolution
    from pyomo.core import TransformationFactory
    from pyomo.contrib.appsi.solvers.highs import Highs


.. parsed-literal::

    WARNING: DEPRECATED: pyomo.core.expr.current is deprecated.  Please import
    expression symbols from pyomo.core.expr  (deprecated in 6.6.2) (called from
    <frozen importlib._bootstrap>:241)
    Interactive Python mode detected; using default matplotlib backend for plotting.


Loads default set of representative days – #TODO allow non defaults by
Tuesday

.. code:: ipython3

    data_path = "./data/5bus"
    data_object = ExpansionPlanningData()
    data_object.load_prescient(data_path)

Builds expansion planning object but not specific model yet – #TODO note
issues that can occur with num_reps too large. Also, make config
overwrite these periods for the distinct times.

.. code:: ipython3

    mod_object = ExpansionPlanningModel(
        stages=1,
        data=data_object.md,
        num_reps=1,
        len_reps=1,
        num_commit=24,
        num_dispatch=4,
    )

.. code:: ipython3

    mod_object.create_model()

.. code:: ipython3

    TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
    TransformationFactory("gdp.bigm").apply_to(mod_object.model)

.. code:: ipython3

    opt = Highs()
    mod_object.results = opt.solve(mod_object.model)

#TODO – demonstrate capabilities to save & load solution info

.. code:: ipython3

    sol_object = ExpansionPlanningSolution()
    sol_object.load_from_model(mod_object)
    sol_object.dump_json("./gtep_solution.json")
    sol_object.import_data_object(data_object)
    sol_object.plot_levels(save_dir="./plots/")
