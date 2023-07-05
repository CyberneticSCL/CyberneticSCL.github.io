% Installation Script for PIETOOLS
% This script adds SOSTOOLS, MULTIPOLY and PIETOOLS to the path overriding
% previous installations. 
% The script also installs latest stable version of SeDuMi overriding
% previous installations.
% The script does not delete any previous versions of sostools and sedumi, 
% but just changes path priority.
% Date - 6/9/2020

clc;
disp('----------------------------------------------');
disp('Installation of PIETOOLS.');
disp('----------------------------------------------');
disp(' ');
fprintf(['Choose the installation directory.\n',...
      'A folder "PIETOOLS" is going to be created in the specified location.\n',...
      'If not specified current directory will be chosen as the installation directory/\n']);
  

% get the installation folder
default_dir = pwd;
c = uigetdir(pwd);
if isequal(c,0)
    fprintf(['No directory selected.\n',... 
        'Installing in the current directory "%s"?\n'],default_dir);
    c = default_dir;
end
 
% create a new directory in that folder
d = [c,filesep,'PIETOOLS'];
if isequal(exist(d,'dir'),7)
    error('The installation directory "%s" already exists.\nPlease, remove or rename the folder or change the installation path.',d);
end
disp('Creating the directory "PIETOOLS".');
out = mkdir(d);
if ~out
    error(['An error appear when trying to create the folder "%s".\n',...
          'Please, check administrative access to modify the destination directory.'],c); 
end

% enter that directory
cd(d);


% download SeDuMi
disp(' ');
disp('Downloading the SeDuMi Master from Github.');
[f, c] = urlwrite('https://github.com/sqlp/sedumi/archive/master.zip', 'sedumi_master.zip');
rehash;

% unzip SeDuMi
disp("Extracting SeDuMi files");
unzip('sedumi_master.zip',[d,filesep,'SeDuMi']);
delete sedumi_master.zip;

if isequal(c,0)
    error('Could not download SeDuMi from the internet. Check internet access or website status.');
end

% Downloading SOSTOOLS, MULTIPOlY and PIETOOLS  %%%Change Path to PIETOOLS zip%%%

disp(' ');
disp('Downloading the PIETOOLS from control.asu.edu.');
[f, c] = urlwrite('http://control.asu.edu/Software/PIETOOLS.zip', 'PIETOOLS_master.zip');
rehash;

% unzip PIETOOLS
disp("Extracting PIETOOLS files");
unzip('PIETOOLS_master.zip',[d,filesep,'PIETOOLS']);
delete PIETOOLS_master.zip;

if isequal(c,0)
    error('Could not download PIETOOLS from the internet. Check internet access or website status.');
end

% get back to the original directory
cd(default_dir);

% add path to PIETOOLS+SeDuMi
disp(' ');
disp('Adding PIETOOLS path to Matlab.');
addpath(d);

% save path for future
disp(' ');
disp('Saving path for future sessions.');
status = savepath;

if status
    fprintf('Could not save the path to a default location,\nplease provide a location where you want to save the path.');
    cn = uigetdir(pwd);
    if isequal(cn,0)
        disp(' ');
        fprintf('No directory specified, saving the path to the current directory "%s".\n\n',default_dir);
        cn = default_dir;
    end
    sn = savepath([cn,filesep,'pathdef.m']);
    if sn
        error(['Could not save the path automatically.\n',...
            'Please, open the "Set Path" button in the Matlab menu and save the path manually to some location.']);
    end
end

disp(' ');
disp('Installation finished.');