# isws-demonstration-models
# --> layered_groundwater_model_modflow

To run the desired scenarios and produce the desired animations, run **main.py**. Common user inputs are found in this 
script between the comment blocks labeled "common user inputs". These include the *max_workers* variable for setting how
many models/animations you would like to run in parallel. This should be equal to or less than the number of CPU cores 
your computer has. The next variable *scenarios* should be the file names of your scenarios schedules. These should be 
found in the repo at */inputs/scenarios/* and should be csvs formatted similarly to the ones already present.

To modify the models' behavior, edit the function *run_the_models()* found in **modules.py**. To modify the animation, 
edit the function *make_the_animation()* found in **modules.py**. The parameter being plotted in the animation can be 
changed via the *parameter* input to the *make_the_animation()* function (can be set to "head" "concentration" or 
"temperature").

### scenarios
scenario 1 - **s1.csv** - is a long-running flow animation meant to be an "idle" flow animation of the model

scenario 2 - **s2.csv** - is a flow and transport model meant to simulate the video of dye injection at well 2 and dye
uptake at well 6

### notes

 - The scenario csvs currently only work for flow and transport phenomena. The energy model inputs are currently hard coded and do not read from the csvs.
 - I (mike-krasowski) have intermittently gotten errors with *sim = flopy.mf6.MFSimulation()* associated with a "structure" or "block". This has twice been solved by what I can only describe as "unintentional iterative troubleshooting." For instance, I added a pumping well for one stress period to **s1.csv** and the error went away. I reverted the change and the error did not return.
   - Another time this was solved by restarting my IDE. That seems to be a reliable way to fix this error. this makes me think it's an issue with file/folder accessing?
 