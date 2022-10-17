% A simple example code that saves out FFST Meyer filter bank coefficients 
% to a file.
% 
% Parameters are specifically set for 512x512 input image, but can be
% adjusted according to user needs.
%
% For more details on how the filterbank is generated, check the
% documentation.

%Parameters for filterbank generation
l = [512 512];
numOfScales = floor(0.5 * log2(max(l)));
shearlet_spect = 'meyerShearletSpect';
shearlet_arg = [];

%Get filterbank
Psi = scalesShearsAndSpectra ( l, numOfScales, shearlet_spect, shearlet_arg );

%Drop all coefficients smaller than 1e-6
Psi(abs(Psi) < 1e-6) = 0;
n = sum(sum(sum(abs(Psi) > 0))); %total nr of non-zero coefficients

%Write it out to file
fileID = fopen('FFST_512x512_table.txt','w');
fprintf(fileID,'%d\n',n);
for i=1:l(1)
    for j=1:l(2)
        for k=1:size(Psi,3)
            if (abs(Psi(i,j,k)) > 0)
                fprintf(fileID,'%d %d %d %f\n',i, j, k, Psi(i,j,k));
            end
        end
    end
end
fclose(fileID);