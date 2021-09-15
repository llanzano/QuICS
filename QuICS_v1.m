%%   User-fiendly code based on the following article:
%$£   EVALUATION OF STED SUPER-RESOLUTION IMAGE QUALITY BY IMAGE CORRELATION SPECTROSCOPY (QuICS)
%$£  
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£ 
%$£                                                                     %$£
%$£                         Luca Lanzanò & Diaspro-Lab                  %$£
%$£      University of Catania - Department of Physics and Astronomy    %$£
%$£         Istituto Italiano di Tecnologia - Nanoscopy Department      %$£
%$£                      User-Friendly Version (26-08-2021)             %$£
%$£                                                                     %$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£

% Main input file is a single .TIF image 
%(or a .TIF stack of 2 independent replicates of the same image, if available)

% optional input file is a .TIF mask image that define a ROI of analysis

% output parameters R, B, N
% R = average FWHM of the 'particles' in the image, extracted by ACF
% B = brightness defined as G0*<I> 
% N = the relative intensity variance due to the noise


function ResultsICS=QuICS_v1(filename );

% default values
px=0.020; %px size in um
mode='c';
maxlagd=12;

%open files ... 
[filename,pathname, filterindex] = uigetfile({'*.tif'},'Select the .TIF image ');
filenamefull = [pathname, filename];   
A=simple_ICCS_readfiletif(filenamefull);
A(isnan(A))=0;
T=size(A,3);
X=size(A,1);
Y=size(A,2);

%load mask
Mask=ones(X,Y);
maskch=1;
[filenamemask,pathnamemask, filterindex] = uigetfile({'*.tif'},'Select Mask or Cancel');
if filenamemask ~= 0
filenamefull3 = [pathnamemask, filenamemask];   
Mask=simple_ICCS_readfiletif(filenamefull3);
Thrmask=0;
Mask=simpleICCS_Threshold(Mask,Thrmask);
% Mask(Mask>0)=1;
end

xch1=A(:,:,1);
if T==2
A1=xch1;
A2=A(:,:,2);
else
% downsampling from checkerboard
X1=floor(X/2);
Y1=floor(Y/2);
for i=1:X1
    for j=1:Y1
A1d(i,j)=mean( [ A(2*i-1,2*j-1), A(2*i,2*j) ]) ;
A2d(i,j)=mean( [ A(2*i,2*j-1), A(2*i-1,2*j) ]) ;  
Mask1(i,j)=mean( Mask(2*i-1:2*i,2*j-1:2*j) ,'all' ) ; 
    end
end

% oversample
[A1]=simple_oversample_image(A1d, 1.9999);
[A2]=simple_oversample_image(A2d, 1.9999);
[Mask_do]=simple_oversample_image(Mask1, 1.9999);
end   
    

%padding (adding pixels with average value outside the ROI)
Extra=0;
[x1pad,Np,Mask,~,Iav]=simple_PadForICS_fromMask(xch1,Extra,Mask);
% [x1pad_e,Np_e,Mask_even,~,Iav_e]=simple_PadForICS_fromMask(xch1,Extra,Mask_even);
% [x1pad_o,Np_o,Mask_odd,~,Iav_o]=simple_PadForICS_fromMask(xch1,Extra,Mask_odd);
if T==2
    Mask_c=Mask;
else
    Mask_c=Mask_do;
end

[x1pad_e,Np_e,Mask_c,~,Iav_e]=simple_PadForICS_fromMask(A1,Extra,Mask_c);
[x1pad_o,Np_o,Mask_c,~,Iav_o]=simple_PadForICS_fromMask(A2,Extra,Mask_c);

% x1=cat(3,x1pad_e,x1pad_o);
% x1=cat(3,x1pad_e(:,:,1),x1pad_o(:,:,2));

maxlag=min( [ ceil(X/4), maxlagd ] );
[m,n]=size(x1pad_e);

% all function are corrected for padding
%calculate 'noise-free' corr function (checkerboard, downsampled)
CCF=simple_ICCS_CCFmean(x1pad_e,x1pad_o)*m*n./(m*n-Np_e);

%calculate standard auto corr function
[m,n]=size(x1pad);
ACF1=simple_ICCS_CCFmean(x1pad,x1pad)*m*n./(m*n-Np);  % (m*n - Np) is the # pixels of the ROI=area ROI

AreaRatio1=( m*n-Np ) / ( m*n ) ; 

% AreaRatio2=( m*n-Np ) / ( m*n ) ;

maxlag=min(maxlag,length(ACF1));

% %set gamma factor
% gamma=0.5;

%set parameters for fit
lag00=1; %first spatial lag to fit
if mode=='c'
lag00=0;
end
lagfit0=maxlag-lag00; %points to fit
figcheck2=1;
figure
subplot(2,2,3);
imagesc(xch1);
axis image
str1 = { 'R=' , 'B=' , 'N='};
dim1 = [.60 0.15 .7 .76 ];
an1=annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
lag0=[lag00 ];
lagfit=[lagfit0 ];
while figcheck2==1
% AreaPSFRatio1=m*n/(pi*2*((FWHM1nm/1000/px)*0.25*sqrt(2/log(2)))^2);
%fit of the 'noise-free' ACF 
subplot(2,2,1)
switch mode
    case 'c'
% (cross of odd and even ACFs)
[param3, fit12, chisq12]=simple_ICCS_Fit_ICS_1Dsingle(lag0(1):lag0(1)+lagfit(1)-1,CCF(lag0(1)+1:lag0(1)+lagfit(1)),4,0,'CCF');
    case 'a'
% or (auto with selected range)
[param3, fit12, chisq12]=simple_ICCS_Fit_ICS_1Dsingle(lag0(1):lag0(1)+lagfit(1)-1,ACF1(lag0(1)+1:lag0(1)+lagfit(1)),4,0,'ACF');
end
xplot=0:lag0(1)+lagfit(1)-1;
plot(xplot,ACF1(xplot+1),'s',xplot,CCF(xplot+1),'*', xplot, (param3(1)+param3(2).*exp(-((xplot-0).*(xplot-0)./(param3(3)*param3(3)) ))) ,'-r' );
legend('ACF','CCF', 'noise-free(fit)')
% xplot=0:lag0(1)+lagfit(1)-1;
% plot(xplot,ACF1(1:lag0(1)+lagfit(1)),'o', )
axis auto


%calculation of parameters:

FWHM=param3(3)*1.18*px*1000;
Br=param3(2)*Iav(1);
rnv=( ACF1(1)-param3(1)-param3(2) )/ param3(2);

str1 = {['R=',num2str(FWHM,3),' nm'], ['B=',num2str(Br,2)] , ['N=', num2str(rnv,2) ]};
set(an1,'String',str1);


subplot(2,2,4)
imagesc(Mask);
axis image


prompt2 = {'Min lag:','Points for fit:','Pixel size (um)','Fit auto (a) or cross (c)'}; 
dlg_title2 = 'Set paramaters for fit'; 
num_lines = 1;
def2 = {num2str(lag0),num2str(lagfit),num2str(px),mode}; 
answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck2=~isempty(answer2); 
if figcheck2==1
lag0=str2num(answer2{1});
lagfit=str2num(answer2{2});
px=str2num(answer2{3});
mode=answer2{4};
end
end


Results(:,1)=xplot;
Results(:,2)=ACF1(xplot+1);
Results(:,3)=CCF(xplot+1); fit12;
Results(:,4)=(param3(1)+param3(2).*exp(-((xplot-0).*(xplot-0)./(param3(3)*param3(3)) )))  ;

R=FWHM;
B=Br;
N=rnv;

ResultsICS= [R B N];

filenameout=filenamefull(1:end-4) ;

answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1
%save data in Matlab
save([filenameout,'_ICS.mat'] );
% writing files txt
dlmwrite([filenameout,'.txt'],Results,'delimiter',';','precision',4);
% dlmwrite([filenameout,'_Norm.txt'],ResultsNorm,'delimiter',';','precision',4);
end


end



%required funtions

function [Aover]=simple_oversample_image(A, factor)

[sx,sy] = size(A);
xq = (0:1/factor:sx)';
yq = (0:1/factor:sy)';
F = griddedInterpolant(double(A));
Aover = uint16(F({xq,yq}));

end

function A=simple_ICCS_readfiletif(fname)
info = imfinfo(fname);
nslice = numel(info);
A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end

end

function y=simpleICCS_smooth_simple(M,sm,n)
y=M;
if sm>0
filt = (1/(8+1/sm))*[1 1 1; 1 1/sm 1; 1 1 1]; % sm factor <=1 
    for i=1:n
    y = filter2(filt,y);
    end
end
    
end

function [y,Np, varargout]=simple_PadForICS_fromMask(x1,Extra,Mask)
[m,n,p]=size(x1);
Mask=double(Mask);
%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
end
%% adding average on zeros
for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 || isnan(Mask(i,j))
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function [y,Np, varargout]=simple_PadForICS_sm(x1,Extra,Thr, sm)
[m,n,p]=size(x1);
for kk=1:p
x1s(:,:,kk)=simpleICCS_smooth_simple(x1(:,:,kk),sm,2);
end

%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
ys=y;
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
    ys(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1s(:,:,k);
end

%% adding average on zeros
Mask=ys(:,:,1);
Mask2=ys(:,:,2);
Mask(Mask<=Thr(1)& Mask2<=Thr(2))=0;
Mask(Mask2>Thr(2)| Mask>Thr(1) )=1;

for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function  Output=simple_ICCS_CCFmean(x1,x2)

NumberOfAngles=180;
[X,Y]=size(x1);
%ACF=conv2(x1,x2,'same');
F1=fft2(x1);
F2=fft2(x2);
ACF= F1.*conj(F2);
G=((sum(sum(x1)))*(sum(sum(x2)))/X/Y);
ACF= ifft2(ACF);
ACF= fftshift(ACF)./G-1;

[R, C]=size(ACF);
if iseven(R)
r0=R/2+1;
else
r0=(R+1)/2;
end
if iseven(C)
c0=C/2+1;
% Radius=min(R/2,C/2);
else
c0=(C+1)/2;
% Radius=min((R-1)/2,(C-1)/2);
end
Radius=min(r0-1,c0-1);

if NumberOfAngles==1
    Output=ACF(r0,c0:end);
else
ACF1=flipud(ACF(1:r0-1,c0:end));
ACF2=ACF(r0:end,c0:end);
ProfMat=zeros(NumberOfAngles*2,Radius);

for j=1:2
    if j==1
        y=ACF1';
    else
        y=ACF2;
    end
    
% CALCULATION OF ROTATIONAL MEAN
% Definition of angles
t=(pi/NumberOfAngles/2:pi/NumberOfAngles/2:pi/2);
   
% Matrix
y=y(1:Radius,1:Radius);
% Cycle between the 2nd and 2nd to last angles
[~, y1y]=size(y);

for i=1:NumberOfAngles
   rt=ceil(cos(t(i))*(1:Radius));
   ct=ceil(sin(t(i))*(1:Radius));
   profile=y((rt-1).*y1y+ct);

   if j==1
   ProfMat(NumberOfAngles+i,:)=profile;
   else
%    ProfMat(i,:)=profile;
   ProfMat(i,:)=[profile(2:end),profile(end)];  % excluding the central ACF point (0,0)
   end   
end

end


Output=[double(ACF(r0,c0)) sum(ProfMat)./(2*NumberOfAngles)];
% 
% OrientedProfiles=min_fw_Profile;
% OrientedProfiles(2,:)=max_fw_Profile;
% Angles=[min_th,max_th];

end


end

function bool=iseven(x)

if mod(x,2) == 0
bool=1;
else
bool=0;
end
end

function [param, fval, chisq]=simple_ICCS_Fit_ICS_1Dsingle(x,y,w0,Display, title1)
%fixed 

my=min(y);
My=max(y);
fun = @(Param) sum( (( (Param(1)+Param(2).*exp(-((x-0).*(x-0)./(Param(3)*Param(3)) ))) -y ).^2)./(abs(y))  ) ;
[param, chisqpar]=fminsearch( fun,[my My w0]);
param(3)=abs(param(3));
fval=(param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) )));
chisq=sum( (( (param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) ))) -y ).^2)./((param(2)^2)  )) ;

if Display==1
%     figure;
    plot(x,y,'o')
    hold on
    plot(x, fval, '-k')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat(title1, '  w=',num2str(param(3),2),'   G0=',num2str(param(2),2) ));
end
end

function B=simpleICCS_Threshold(A,thr)
  
if length(thr)==1
B=A;
B(B<=thr)=0;
B(B>thr)=1;
else
B=A;
B(B>thr(2))=0;
B(B<=thr(1))=0;
B(B>0)=1;
end

end


