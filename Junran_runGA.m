function [x,fval,exitflag,output,population,score] = Junran_runGA(nvars,...
                                                            lb,...
                                                            ub,...
                                                            PopInitRange_Data,...
                                                            PopulationSize_Data,...
                                                            Generations_Data,...
                                                            InitialPopulation_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopInitRange', PopInitRange_Data);
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'Generations', Generations_Data);
options = gaoptimset(options,'InitialPopulation', InitialPopulation_Data);
options = gaoptimset(options,'Display', 'diagnose');
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv });
options = gaoptimset(options,'Vectorized', 'off');
options = gaoptimset(options,'UseParallel',1);
options = gaoptimset(options,'TolCon',1e-3);
[x,fval,exitflag,output,population,score] = ...
ga(@Junran_spmet,nvars,[],[],[],[],lb,ub,@Junran_constaints,[],options);
% ga(@Junran_OCV_OBJ,nvars,[],[],[],[],lb,ub,[],[],options);


