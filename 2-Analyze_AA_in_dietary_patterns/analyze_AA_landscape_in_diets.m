%----- Analyze landscape of amino acids in human dietary patterns -------%
% Input datasets:
%   1. ../data/AA_variables.mat: variables resulting from the previous step
%   (1-AA_Variability_Foods)
%   2. list_animal_dairy_fat.csv: list of animal and dairy products in the
%   food category "Fats"
%   3. list_animal_fat.csv: list of animal products (dairy products not
%   included) in the food category "Fats"
%   4. USDA_Goals.mat: MATLAB workspace file consisting of variables
%   describing the ranges of nutrient intakes recommended by the USDA
%   2015-2020 dietary guidelines
%
%  Output dataset:
%   ../data/AA_variables.mat
%------------------------------------------------------------------------%

addpath(genpath('../functions'));
load ../data/AA_variables.mat;
define_diets; %Generate the mathematical forms of 10 human dietary patterns

% Compute ranges of absolute AA intake in dietary patterns
get_abs_AA_intake_ranges;

% Sample relative AA abundance in dietary patterns
sample_AA_composition_in_diets;

% Save variables to the workspace
save ../data/AA_variables.mat