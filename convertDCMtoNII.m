function [  ] = convertDCMtoNII( filenames , outputFolder )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

matlabbatch{1}.spm.util.import.dicom.data = filenames;

matlabbatch{1}.spm.util.import.dicom.root = 'flat';
matlabbatch{1}.spm.util.import.dicom.outdir = {outputFolder};
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

spm_jobman('run',matlabbatch);

end

