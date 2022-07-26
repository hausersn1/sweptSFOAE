Stimulus frequency otoacoustic emissions are responses from the cochlea to pure tone stimuli. Since the stimulus is the same frequency as the emission, we have to use some tricks to separate the stimulus and the OAE. 

`makeSFOAEstim.m` generates a sweep of a probe tone and a suppressor, as well as sets some other stimulus parameters. This saves a `stim` struct so which is read by `runSFOAEswept`

`runSFOAEswept` generates the stimuli and plays them with the ER10X. The `stim` struct contains the data and the stimulus parameters. 

`analyzeSFOAEswept` will analyze a loaded stimuli structure that contains the data and plot the amplitude and phase. 