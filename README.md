# ejc-recist-paradox

# License & Disclaimer
MIT License

Copyright (c) 2019 Millennium Pharmaceuticals, Inc

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"The derived data sets in the data/ subdirectory are based on research using information obtained from www.projectdatasphere.org, which is maintained by Project Data Sphere, LLC. 
Neither Project Data Sphere, LLC nor the owner(s) of any information from the web site have contributed to, approved or are in any way responsible for the contents of this code and data set.”  Per the project data sphere user agreement, only registered users of projectdatasphere.org may access the datasets in this repository. 

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


 
