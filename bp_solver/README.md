# Read Me

* Call the func_BP_CST to get the eps_BP (reconstructed permittivity) and the 	 (reconstructed conductivity)
* func_BP_CST has three parameters. The parameter 'freq_num' is the index of frequency sample
* Totally 13 frequency samples are supplied in func_BP_CST, they are 700MHz, 750MHz, 800MHz, 850MHz, 900MHz, 950MHz,
* 1000MHz, 1050MHz, 1100MHz, 1150MHz, 1200MHz, 1250MHz, and 1300MHz
* freq_num = 1 means the first frequency sample is selected, which is 700MHz
* freq_num = 2 means the second frequency sample is selected, which is 750MHz
* freq_num = 3 means the third frequency sample is selected, which is 800MHz, and so on...
* target_filename is the name of the S-parameter files. For example, target_filename = '101309_1.s16p' means we build the permittiivity
* and conductivity for the file "101309_1.s16p"
* To get the results from all the brain models, you need to write your own loop to visit all the S-parameter files 
* The parameter resolution specifies the resolution of the reconstructed images. It has to be a string variable
* Only two options are supplied for the parameter of resolution, they are '2mm' and '4mm'
* Therefore, use resolution = '2mm', or use resolution = '4mm'
* The folder "Parameters" includes all the variables that you need to use in func_BP_CST
* Please include the folder "Parameters" in your current working path before using this code
* I suggest to use the frequency samples of 700MHz ~ 950MHz. 
* I suggest do not go to the frequency samples higher than 950 MHz because based on the BP results, they look not very good

* Codes are written by Dr. Lei Guo from the University of Queensland
* Email: l.guo3@uq.edu.au

-------------------------------------------



Alex's comments: 
* According to the **Sparams** that I supplied, there's no such fraction frequencies as the ones used in the function.
* Why is it that empty system is used, despite the fact that it is not practically measured.
* Lei usually uses 4mm and then upsample the output image to increase the resolution by 2 (making it look like 2mm).
* How to generate results with the same meshing specs as done in DeepHead?
* What are the values of permittivity outside the domain of interest?
* What is the geometry of the results? I.e. where is antenna one located in the outcome image when shown as is (coordinates 0, 0 at the top)

func_BP_CST_with_empty
This is the same as func_BP_CST except that it uses the empty system in addition to two cals. Notice that in reality, at least at the time of writing this, we don't have empty domain measurement so this is not practical. We want the tuning process in reality be similar to that done in the experiment. This function with this name is not used. If you want to make it used remove _with_empty suffix and avoid name collission with the function alredy with that name.



Lei's comments (07/01/2021):
* A new function named 'func_BP_CST_dynamic.m' is added. This is an interface function of the BP solver.
  You still need to write your own code to visit all data from winner simulation and pass them into this interface function 
* You need five input parameters for 'func_BP_CST_dynamic.m'. 
  * The first one is the frequency (you can input any frequency such 856.35 MHz now). 
  * The second one is the file name of target.
  * The third and fourth ones are the file names of the calibration phantoms (they either from CST or measurement)
  * and the last one is the resolution (although you can specify any value of resolution, 
    but I don't suggest to go beyond 2 because it will take huge time to calculate the Green functions)
* A new function named 'MoM_GreenFunc.m' is added. You don't need to worry about this function because it is only called within 'func_BP_CST_dynamic.m' 
* The code named 'test_newBP_solver.m' is only used for testing 'func_BP_CST_dynamic.m'. You can delete or hide it if you want
* A new function named 'build_Exp_Cal_model.m' is added. You don't need to worry about this function because it is only called in 'MoM_GreenFunc.m'
