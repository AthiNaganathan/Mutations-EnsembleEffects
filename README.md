# Mutations-EnsembleEffects
Includes data on protein ensemble properties (coupling free energies from the WSME model) that exhibit an exponential distance dependence on introducing mutations

Download the coupling free energies data for PDZ (https://tinyurl.com/5xxnvrnc), CheY (https://tinyurl.com/2s49dczy), and CypA (https://tinyurl.com/2p8jxv3r). The 'mat' files include the differential coupling matrices and incdices which are used to generate the ouputs using the MATLAB script below.

Execute 'Coupling_AlaMut.m' - line 3 can be modified to choose from among the three proteins studied (Pdz, CheY and CypA). Line 4 specifies the residue number whose effect on the ensemble (via alanine substitution) can be observed. 

Files_Variables.docx provides information on the various parameters.
