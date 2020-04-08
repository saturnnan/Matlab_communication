% Skull information needs to be added.
clear
close all;clc
%% create the computational grid(kgrid)
Nx = 200;           % number of grid points in the x (row) direction
Ny = 200;           % number of grid points in the y (column) direction
Nz = 320;
dx = 0.2e-3;        % grid point spacing in the x direction [m]
dy = dx;              % grid point spacing in the y direction [m]
dz = dx;              % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%% Data load 
% load('Skull_9_CT_ShearW_1119.mat');
% https://www.itis.ethz.ch/virtual-population/tissue-properties/database/
TissuePar.Water.Type = 'Water';
TissuePar.Water.Dens = 994;
TissuePar.Water.HeatCap = 4178;
TissuePar.Water.ThermalCond = 0.60;
TissuePar.Water.SoundSpeed = 1482.3;
AcousticPar.Freq = 0.8e6;                                                          % [Hz]
medium.sound_speed_compression = ...
    TissuePar.Water.SoundSpeed * ones(Nx, Ny, Nz);               % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny, Nz);              % [m/s]
medium.density = TissuePar.Water.Dens * ones(Nx, Ny, Nz); % [kg/m^3]
% % skull medium (According to the skull coordinate matrix)
% hounsfield2density                             % see k-Wave MATLAB toolbox
% medium.sound_speed_compression
% medium.sound_speed_shear  % [m/s]
% medium.density
%% PML design
kgrid_PML.Alpha_X = 10;
kgrid_PML.Alpha_Y = 10;
kgrid_PML.Alpha_Z = 10;
kgrid_PML.Size_X   =  5;
kgrid_PML.Size_Y   =  5;
kgrid_PML.Size_Z   =  5;
%% Emission source and sensor
ArrayPar.EleR             =   4.5e-3;     % Element r [m] 
ArrayPar.R_Sensor         =   15e-3;      % Sensor r [m] 
ArrayPar.position_ArrayR  = - 25e-3;      % Element position in z [m]
ArrayPar.position_Sensor  =   25e-3;      % Sensor  position in z [m]
[kgrid.t_array, kgrid.dt] = makeTime(kgrid, medium.sound_speed_compression);
% % Extend reception time of sensor ensuring a stable signal
kgrid.t_array = 0:kgrid.dt: 1.8 * size(kgrid.t_array,2)* kgrid.dt;
%% define source element and a time varying sinusoidal source
% Custom elements or makeDisc(k-Wave MATLAB toolbox)
source.u_mask = disk_element(kgrid.Nx, kgrid.Ny, kgrid.Nz,...
    ArrayPar.EleR, ArrayPar.position_ArrayR, kgrid.dx);
Ele_U = 2.08e-08;     % Vibration speed
source_mag_u = Ele_U; % 800kHz D9mm
source.uz = source_mag_u*sin(2*pi*AcousticPar.Freq*kgrid.t_array);
% Custom sensor
[sensor.mask,SensorCoo] = square_element(kgrid,ArrayPar);
sensor.record = {'p'};
% GPU simulation
input_args = {'DisplayMask', source.u_mask,...
    'PMLAlpha',[kgrid_PML.Alpha_X,kgrid_PML.Alpha_Y,kgrid_PML.Alpha_Z],...
    'PMLSize', [kgrid_PML.Size_X, kgrid_PML.Size_Y, kgrid_PML.Size_Z],...
    'PlotSim',false,'DataCast', 'gpuArray-single'};
% ShearWave GPU
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});
