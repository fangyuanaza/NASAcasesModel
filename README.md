# NASAcasesModel

The folder contains two explict models used in DAY2 02 Fang challenge.pdf in 
2022 Symposium on Turbulence Modeling: Roadblocks, and the Potential for Machine Learning
(https://turbmodels.larc.nasa.gov/turb-prs2022.html)

Two models are presented. Generally, both model 1 and 2 improve subsonic jet and wall-bounded hump
while remain similar performance as baseline k-Omega SST for plate, channel and airfoil cases.
The difference lies in model 1 improve the simulation of plate when comparing to theorical data 
while not contribute to hump too much (improve 31%). For model 2, the plate results lie within the
uncertainty range of experimental data while can contribute to hump much more (improve 69%).

The expression can be found in the article ()

Based on our simulation result, we recently found the simulation result of only keep the leading term
Vij of Rij (no aij(no correction to Reynolds stress (you can set nonlinearStress_=2*k_*0*T1))) and can 
achieve similar results as the whole expression.

For a posteriori test on NACA 0012 cases with different angles of attack, the two models give similar 
results as baseline k-Omega SST and they both simulate stall phenomenon at 19 degree.

Email: yuanfang1@student.unimelb.edu.au
