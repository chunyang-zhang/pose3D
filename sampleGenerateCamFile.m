%%Specify the num of routings
nFrames = 2086;
y = 1.0;


%notice the axis

sz = 0.2;
cR = 5;
cw = 640;
ch = 480;
fx = 700;
fy = 700;
cx = cw/2;
cy = ch/2;

%%only z and y , yaw moves. Demo circle here

yawP = 0:2*pi/(nFrames-1):2*pi;

fileId = fopen('./models/18.cam','w');
fMat = zeros(nFrames,16);
fprintf(fileId,':camera\n');
fprintf(fileId,'%f %f %f %f %f %f %f\n',cw,ch,fx,fy,cx,cy,sz);
for i=1:1:2086
    crotMat = eul2rotm([0,-yawP(i)-pi-pi/2,0],'ZYX');
    cvec = [ cR*cos(yawP(i)) y cR*sin(yawP(i))];
    callMat = [crotMat cvec';0 0 0 1];
    callVec = reshape(callMat,[16,1]);
    fprintf(fileId,'%f ',callVec);
    fprintf(fileId,'\n');
end
fclose(fileId);