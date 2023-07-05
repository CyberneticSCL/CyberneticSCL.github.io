Setup and run instructions:

1) Run "Utility_Battery_Programmer.exe" to initialize the setup wizard and proceed
by pressing "Next". The wizard installs the package and Matlab's Runtime Compiler.

2) After the installation is complete, make sure your user account has the premission
to read from and write to all the files in "/Utility_Battery_Programmer/application/Output"
and "/Utility_Battery_Programmer/application/load_data".

3) To run the package, run "/Utility_Battery_Programmer/application/Utility_Battery_Programmer.exe".
A dialoge window should open. Answer the questions by typing a value or name in the 
box below each question. 

4) If your answer to Question 1 is NOT "3", then leave the box below Question "2" 
empty. If your answer to Question 1 is "3", then type down the name of your ".csv"
file which contains the hourly retail load data (e.g., retail_load_15_to_17-July.csv).
The datafile must have a single column containing the hourly retail load values
in MW for any number of days. Example files can be found in the Folder: 
"/Optimal_Battery_Programmerapplication/load_data".

5) In Question 10, the data file containing the hourly market prices must have the
same format as described in Item 4.

6) Please make sure that the data files for Questions 2 and 9 are located in the 
path: "/Utility_Battery_Programmer/application/load_data".

7) At this time, the algorithm may exhibit instability or take longer to converge
for the cases where:
battery capacity < 100 MWh, charging rate < 30 MW, maximum arbitrage < 50 MW. Indeed
for some of these cases the problem can be infeasible.

