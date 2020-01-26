% 493 volumes

mglRetinotopy('displayName=fMRIProjFlex','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks',3,'doEyeCalib=0');

%% also 493 volumes (for fMRIProj32)

warning('Did you check the # of ignored triggers?');
mglRetinotopy('displayName=fMRIproj32','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks',3,'doEyeCalib=0','easyFixTask',-0.5,'barWidth',1,'elementVelocity',2,'elementSize',0.75,'barSweepExtent',14);


%% fMRIProjFlex 176 frames and 1.5 s (mux3) sequence

warning('Did you change the # of ignored stimuli?');
mglRetinotopy('displayName=fMRIproj32','bars=1','fixedRandom=1','stimulusPeriod=26','stepsPerCycle',20,'blanks',3,'doEyeCalib=0','easyFixTask',-0.5,'barWidth',1,'elementVelocity',2,'elementSize',0.75,'barSweepExtent',14);
