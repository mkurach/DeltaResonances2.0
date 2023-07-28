# DeltaResonances
Analysis I am working on in HADES Collaboration, GSI Darmstadt.
I am focusing on Delta0 and Delta++ reconstruction from their pion-proton and pion+proton decay channels.

**makeHistTermBasic**/**compHistTerm**
- produce 4 basic histograms (from article)
- input: Jędrzejowe files from *lustre*
- output: root files with clear histograms in  *outputBasic*

**compTermBasic**/**compTerm**
- draw pictures with plots for 4 basic histograms
- input: root file from *outputBasic*
- output: png in *outputBasic*

**makeHistTermPt**
- produce histograms of pt for different mass ranges
- input: Jędrzejowe files from *lustre*
- output: root files with clear histograms *outputPt*

**compPtMFit**
- fit function to the plots, draw T_eff vs. M plot (delta only or additinal particles)
- input: root file from *outputPt*
- output: png in *outputPt* and pdf in *ladne*

**compPtYFit**
- fit funtion to the pt for differen y ranges, draw T_eff vs. y plot
- input: histograms from *outputBasic*
- output: png in *outputPt* and pdf in *ladne*

**makeParticlesPt**
- produce pt histograms for 3 stable particles
- input: Jędrzejowe files from *lustre*
- output: root files in *outputPart*

**makeParticlesMt**
- produce mt histograms for 3 stable particles
- input: Jędrzejowe files from *lustre*
- output: root files in *outputPart*

**compParticlesFitPt**
- fit funtion to pt histograms for stable particles
- input: root files from *outputPart*
- output: root files with fitted function in *outputPart*

**compParticlesFitMt**
- fit funtion to mt histograms for stable particles
- input: root files from *outputPart*
- output: root files with fitted function in *outputPart*

**makeHistEvMix**
- perform event mixing on deltas
- input: Jędrzejowe files from *lustre*
- output: root files in *outputEvMix*

**compEvMix**
- compare event mixing and true Monte Carlo delta
- input: root files from *outputEvMix* and from *outputBasic*
- output: png in *outputEvMix*

**compHistSig**
- create signal, background and subtracion from event mixing
- input: root files from *outputEvMix*
- output: root file and comparison plots in *outputSig*

**diffFactor**
- compare different division factors in event mixing method
- input: signal histograms from *outputEvMix*
- output: root files in *outputSig*

**compMFit**
- perform fit to mass histograms
- input: histograms from *outputBasic*
- output: png in *outputM* and pdf in *ladne*

**compMtMFit**
- perform fit for mt for different m ranges
- input: histograms in *outputMt*
- output: root files and png in *outputMt*

**compMFitNew**
- perform fit for mass for Manley and Monitz parametrisation
- input:  root files from *outputBasic*
- output: root files in *outputM*