These files were written using MATLAB2023b and EasySpin 6.0.6, use other version at your own peril. (EasySpin is trivially easy to download just google it)
The main file is TT_esr_fluctuating_J.m, which calls fluctuating_J_populations.m which in turn calls SLE_code.m, so you'll need to download all three to do anything. 
Currently the code loads a dataset called ESR_1.5us.mat, which imports a variable called d; remove this if you just want to simulate.
