function K = convert13ks(fpath)
flist = dir([fpath,filesep,'*.csv']);
N = 60;
frame = 1:N;
K = zeros(13,N);
for k = 1:length(flist)
    mat = readmatrix([fpath,filesep,flist(k).name]);
    K(k,:) = interp1(mat(:,1),mat(:,2),frame,'spline');
end
end
