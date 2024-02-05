function slowVideo(file, FR)
Vobj = VideoReader(file);
n = Vobj.NumFrames;

v = VideoWriter([file(1:end-4) '_slow.avi']);
v.FrameRate = FR;

open(v);
%%
for i= 1:n
    matrix_image = read(Vobj, i);
    writeVideo(v,matrix_image); 
end
