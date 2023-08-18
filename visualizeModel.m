% init instance with your parameters
network_model = initModel(100,30,1,0.5,1,200000);

% visualization of the generated instance
gplot(network_model.adjMat,network_model.locMat,'-o');

% Our instance generator contains randomization, i.e., it generates different instances with the same parameters.

% You can use different parameters or run the algorithm repeatedly to generate multiple different instances. 

% If you want to save instances, please design your own function to save single instance information to a file like txt.