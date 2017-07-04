function save_spuriousC_files(rS, St, wc, eq, filename)
% save_spuriousC_files(rS, St, wc, eq, filename)
%
% saves four files in the local directory: 
%
%   filename.rS, filename.St, filename.wc, filename.eq
%
% which can be read by spuriousC.
%
%   DNAdesign package (Winfree lab, Caltech)

fid=fopen([filename '.rS'],'w');
fprintf(fid,'%s',rS);
fclose(fid);

fid=fopen([filename '.St'],'w');
fprintf(fid,'%s',St);
fclose(fid);

fid=fopen([filename '.wc'],'w');
fprintf(fid,'%d ',wc);
fclose(fid);

fid=fopen([filename '.eq'],'w');
fprintf(fid,'%d ',eq);
fclose(fid);

