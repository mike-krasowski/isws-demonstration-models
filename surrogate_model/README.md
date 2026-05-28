# isws-demonstration-models
# --> surrogate_model

This set of scripts attempts to have an artificial neural network (ANN) mimic a MODFLOW model simulation. Running the **surrogate.py** file will generate a MODFLOW model of a "toy" model and an ANN that attempts to replicate that same "toy" model. The script generates plots for both the MODFLOW model (both results and verification plots) and the ANN. Nothing should need to be changed or updated in this file to accomplish this work.

The MODFLOW "toy" model is a 25x25 model with 3 layers and 5 stress periods, with the first stress period steady-state and the rest being transient. The "toy" model also includes a well that withdraws a decreasing rate of over model run (i.e., 100 L^3/t in SP1, 90 L^3/t in SP2, 80 L^3/t in SP3, 70 L^3/t in SP4, and 60 L^3/t in SP5). The "toy" model also contains 9 RIV cells that are clustered to replicate as a lake. Lake water levels are 5.5 L (unit "Length") for each modeled stress period. A constant recharge of 0.005 L^3/t are applied to all active cells in the model's top layer. More details on the "toy" model can be found in the `run_to_model()` function of the **toy_model.py** script.

In **surrogate.py**, model data is correlated into a large dataframe for ANN training. Each row in the data frame includes information for each model cell during a given modeled stress period. A row in the dataframe includes:

1. stress period information,
2. cell coordinates,
3. a cell's distance from the well,
4. well pumpage,
5. a cell's distance from surface water cells (i.e., RIV cells),
6. the stage of the nearest SW cell, and
7. the head pressure at each cell for ANN testing.

The ANN is then trained on a random selection of 80% of the aggregated model data (for all 5 stress periods), specifically relating how parameters 1-6 (independent variables) correlate with item 7 (dependent variable). The ANN then evaluates its predictive capabilities against the remaining 20% of the aggregated data.
