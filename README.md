ToML
==================================================================================

All the codes in the folders: "BIM_solver", "DBIM_solver", "CSI_solver", "MR_CSI_solver" are written by Lei Guo, from School of ITEE, UQ

Email: l.guo3@uq.edu.au

Instruction for using the code named "process_data" in all the four solvers' folder
-----------------------------------------------------------------------------------

You have four variables: "eps_data", "sigma_data", "freq_low", "freq_high", "show_flag".

You should generate "eps_data" and "sigma_data" from the solvers (they are the calcualated permittivity and conductivity from the tomography solver).

The size of "eps_data" and "sigma_data" is 116 x 100 x 84....
116 is the number of pixles along x-direction, and 100 is the number of pixels along y-direction...
84 is the number of frequency samples....
You have the tomography images under 84 different frequency samples, from 701MHz to 1.199 GHz, with the frequency step of 6MHz

"freq_low" is the lowest frequency bin that you want to use for the post-processing, and "freq_high" is the highest frequency bin.

For example, when you choose freq_low = 1, freq_high = 84, then you selected the tomography images under the frequencies from 701MHz to 1.199GHz to do the post-processing.