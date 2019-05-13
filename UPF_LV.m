clc;clear;close all;
load('BilinModel_LV_subjects90pct_phases90pct_svd05it_proc05it\BilinModel_LV_subjects90pct_phases90pct_svd05it_proc05it.mat');
TR = triangulation(model.faces,model.mean);
g = trisurf(TR);axis image;
g.EdgeColor = 'w';
% g.CData = model.pointdata.epiendo;
g.CData = model.celldata.part;
Y = model.W*model.B;
