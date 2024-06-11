# Class project
## Titile : Counting single cell-encapsulated alginate microparticles from microscopic images
-----
### Introduction
#### Background

##### 1. Single Cell Analysis

* The interesting technique the can give the precised results without cell-to-cell variation and cell-to-cell interaction. It can be used with since bacteria cells to animal cells. This single cell analysis can be used for exploring more information of minor population. Therefore, it will be useful in various aspects such as small molecule detection, whole genomr sequencing and antibiotic resistance detection.

<img width="994" alt="Screenshot 2024-04-08 at 9 24 16 PM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/e3262ca5-c59e-4b5b-a236-dcd03d4cbe58">

##### 2. Alginate Microparticles

* Encapsulation technique is necessary thing for single cell studying. This technique can protect cells and  avoid the contamination. Alginate is the one of the best polymer used for encapsulation of cells because of its biocompartibility and low cost. Alginate will be gelatinized by calcium ion (called egg-box model) and normally formed as spherical shape.

<img width="994" alt="Screenshot 2024-04-08 at 9 20 44 PM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/1bc59ce0-9061-41a8-a728-938ea328bdce">


* For using the alginate microparticles in the aspect of single cell analysis, size and the number of microparticles is neccessary point and should be considered

-----

#### Question to be solved by class project
* As mention above, morphology, size, and also number of alginate microparticles is important factors for using in single cell analysis and also other works. Normally, morphology, size, and number of alginate microparticles can be calculated by ImageJ from microscopic capture. However, it can analyze only one picture per process and the gained result willl be analyzed again for making "histogram" or tendency graph by other applications.

* Therefore, it will be better if some program can be used to analyze alginate microparticles in the aspect of size, number, and circularity of whole sample pictures as once and display the results as histogram of size of microparticle or tendency graph. In this study, R programm will be used for improving the analyzing alginate microparticles process as workflow below :
  ##### Workflow

<img width="994" alt="Screenshot 2024-04-08 at 10 09 11 PM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/9c81b681-324d-4c1c-a8d0-09ab67a070ea">

-----
# BeSIDE
#### BeSIDE - Bead Size Distribution Estimator
BeSIDE (Bead Size Distribution Estimator) is the program created for estimate the size of microparticles by using many microscopic picture for analyzing before visualize in th form of various histograms and raw data.csv.  

<img width="892" alt="Screenshot 2024-06-11 at 9 33 14 AM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/1ad52d73-f8ea-4531-af07-570abee20aaf">

## Running procedure 
#### 1. Loading "BeSIDE.R" that is provided in datafolder
#### 2. Setting folder for pictures that you want to analyze. 
Each picture has to rename as Img1, Img2, Img3,... before analyzing step. Further, inside this folder, creating one folder named "result" to prepare the place for gaining analyzed results. Nite that this version can be applied with microscopic picture with 10X or magnitude equal as 100 only.
 
<img width="936" alt="Screenshot 2024-06-11 at 10 01 18 AM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/534cf197-be79-4d7c-bce2-3eda0791507b">


 #### 3. Run "BeSIDE.R" for opening application for analyzing results
 * 3.1 Select color that you want to design the histogram graph. For assuring that your selected color is set, click "SET" button. If the meassage "Your update is complete mean your selected color is already set. If you want to reset color as default, you can just click "RESET" button.
   
 * 3.2 For selecting folder "CLICK RUN BUTTON ONE TIME" for choosing the folder that has your microscopic pictures inside. After choosing, system will run autonomically. Please waiting until text  "Analysis is done" happened. At this step all analyzed results will be collected in folder named "result" that you created at the beginning.

<img width="980" alt="Screenshot 2024-06-11 at 10 08 47 AM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/e776f6a7-c899-4321-9edc-727b4dad8ae5">

 * 3.3 If you want to change color of design histogram graph, you can select color again and click "UPDATE" button for updating the new set color. The new design will be saved again and saved in the "result" folder instead of previous version autonomically.
     
<img width="1018" alt="Screenshot 2024-06-11 at 10 22 00 AM" src="https://github.com/igchoi/IBT618-SystemsBiotechnology/assets/164750221/724b8610-d187-4c01-93bb-7be922829e3e">

-----
#### Tools
* Tool named "pliman" will be used and adapted for this study.
#### Links
1. https://github.com/TiagoOlivoto/pliman












