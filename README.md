# Undergrad-Honours-Research-Project
Welcome! This repository contains my undergrad honours thesis work, and the work I've been completing afterwards to be able to write a paper. 

The skills and techniques demonstrated in this project are: 
-

-----------------------------------------------------------------------------------------------|

In summary, my project seeks to analyze wound healing responses in bone cells (osteoblasts). When osteoblasts encounter physical damage, ATP is released into the extracellular space. This ATP acts upon receptors on osteoblasts, which triggers calcium from the extracelluar space, and from internal cell stores, to flood inside the osteoblasts, leading to a wound healing response. This calcium response can be measured and visualized experimentally using calcium dyes, leading to the production of Calcium (concentration) vs. time(s) response curves. It turns out, there is a wide heterogeniety of calcium responses that are seen, and it depends on factors such as extracelluar [ATP] concentration, and the calcium receptors on the cell. My projects seeks to biophysically model the calcium reponses using an ordinary differential equation (ODE) model formalism. Specifically, my work saught to explore what produces the fast oscialltios vs the slow oscillations seen on the calcium curves. Further, a clustering analysis was implemented to highlight the heterogeneity of calcium responses, and to see if certain responsee curves occured more frequently. This clustering analysis was coded by one of our collaborators: Nicholas Mikolajeicz, and I implemented a few changes. 

My final paper for the thesis can be read in the file: "Phgy_Project_Final_Paper_AnthonyQ_WithEdits.pdf". The MATLAB code which contains and runs the ODE model can be found in the folder: "Matlab_Model_Github". The R code for the clustering analysis can be found in the file: "calciumSignalling_signatureClustering_060420_Anthony'sEdits.Rmd". The results of the clustering analysis, when adjsuted to various levels of curve smoothing (i.e. 'span'), can be found in the folders: "calciumSignatureClusters_span=3%", "calciumSignatureClusters_span=5%", and "calciumSignatureClusters_span=10%". 



