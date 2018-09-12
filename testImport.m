clc;close all;clear all; 

% filename to import
filename = 'S03-Trial-Walk-1-JointCenters.csv';

% specify the names of markers to export from the vicon file, with one
% marker name per ROW, not column
markerNames = {'LAJC';
               'LHJC';
               'LKJC';
               'RAJC';
               'RHJC';
               'RKJC'};

markers = importViconMarkers(...
     'path2file',fullfile(pwd,filename),...
     'markerNames',markerNames);
 
% pairs of marker names to draw lines between 
markerPairs = {'LAJC','LKJC';
               'LHJC','LKJC';
               'LHJC','RHJC';
               'RAJC','RKJC';
               'RHJC','RKJC';               
                };

% color of lines drawn between marker names            
markerSetColour = {'b','b','g','r','r'};            
            
% plot and visualise in an interactive graphical user interface
[hFig ] = animViconMarkersV2('markerData',markers,...
                           'markerSet',markerPairs,...
                           'markerSetColour',markerSetColour);