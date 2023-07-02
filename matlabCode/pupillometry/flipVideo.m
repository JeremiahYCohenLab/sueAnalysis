function flipVideo(file, dim)
Vobj = VideoReader(file);
n = Vobj.NumFrames;

v = VideoWriter([file(1:end-4) '_flip.avi']);
v.FrameRate = Vobj.FrameRate;

open(v);
%%
for i= 1:n
    matrix_image = read(Vobj, i);
    if dim == 1
        matrix_image_rev = flipud(matrix_image);
    else
        if dim == 2
            matrix_image_rev = fliplr(matrix_image);
        end
    end
 
  writeVideo(v,matrix_image_rev);  
 
end
%%
close(v)

