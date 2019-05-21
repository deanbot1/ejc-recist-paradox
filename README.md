# ejc-recist-paradox

# Disclaimer
THIS CODE IS AS-IS without any representation of validation beyond reproducing the graphs of the EJC manuscript. <<paste in license>>
# License

# About the code
This is code supporting the publication: 
Directional inconsistency between Response Evaluation Criteria in Solid Tumors (RECIST) time to progression and response speed and depth, EJC 109 (2019) 
https://doi.org/10.1016/j.ejca.2018.11.008
Please cite and acknowledge if you use this code in original or modified form. We would also welcome collaborations.

# How to use this code:
1. Launch MATLAB and change directory to the code/ subdirectory of this project
2. Never change directories within MATLAB. All paths are relative to the code directory, so our code won't work if you change working directory

# To generate the results in the article, from the matlab prompt:

>> run % runs the parameter estimation, looping thru all indications

>> post_process % generates article figures: 6,7,3a,3b,4c,5a (see titles and/or figure names)

>> analyze_paradox % generates article figures: 1,4a,4b,5,and graphical abstract left figure

>> make_results_table  % generates article figure 2 and table 1 

# Questions 
please contact Dean Bottino: dean.bottino@takeda.com, or LinkedIn

"The derived data sets in the data/ subdirectory are based on research using information obtained from www.projectdatasphere.org, which is maintained by Project Data Sphere, LLC. 
Neither Project Data Sphere, LLC nor the owner(s) of any information from the web site have contributed to, approved or are in any way responsible for the contents of this code and data set.”
 
