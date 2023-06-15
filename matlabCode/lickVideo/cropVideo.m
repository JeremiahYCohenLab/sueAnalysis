% videoName = 'F:\lickSampleVidoes\testCamera\testEdmound\19-25-16.451';
% videoName = 'F:\lickSampleVidoes\testCamera\test2\18-51-24.980';
% videoName = 'F:\lickSampleVidoes\testCamera\test3\13-46-43.923';
videoName = 'F:\lickSampleVidoes\testCamera\test4\20-13-31.296';
vid1=VideoReader([videoName '.avi']);
n=vid1.NumberOfFrames;
writerObj1 = VideoWriter([videoName '-crop.avi']);
writerObj1.FrameRate = 10;
open(writerObj1);
%%
for i= 76:130
    % 7615:7665
    % 7506:7545
    % 5525:5583
  im=read(vid1,i);
  imc=imcrop(im,[180 230 330 150]);% The dimention of the new video 
  % Edmund 50 150 550 250
  % 2 100 250 330 150
  % 3 150 200 374 170
 
  writeVideo(writerObj1,imc);  
 
end
%%
close(writerObj1)
%%
a = 1/10;
b = 2*a/(1-a^2);
