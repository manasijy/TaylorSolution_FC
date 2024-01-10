This repository contains three sub-repositories: TaylorFC, TaylorG and TaylorLPcystal
TaylorFC calculates texture evolution using taylor full contraint model with full dislocation glide. It solves taylor ambiguity using minimum plastic spin criteria. There is a funciton to determine the slip system solutions with repective shears.
It gives a clear perspective for a learner about how different crystal plasticity equations are implemented to determine plastic spin and updating grain rotations.
The TaylorG is a general script to calculate texture for a general deformation mode involving any of all of full, partial slip and/or twin. Same method is followed here.
The Taylor single crystal script is to caluclate the rotation of individual orientation. This could be done using the previous scripts too but to reduce the computation time for a single orientation this script was developed.
